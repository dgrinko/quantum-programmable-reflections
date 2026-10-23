# Numerically stable Clebsch–Gordan coefficients for SUNRepresentations.jl v0.3.
#
# SUNRepresentations.jl (≤ v0.4) fixes the highest-weight vectors of s1 ⊗ s2 → s3 and then
# obtains all other CGCs by applying lowering operators weight by weight. That recursion
# amplifies rounding errors exponentially in the length of the weight chain: for SU(2) the
# CGCs of spin j ⊗ j are off by 1e-9 at j = 10, by O(1) at j = 20 and by many orders of
# magnitude beyond, which makes the qubit transition matrix wrong from n ≈ 70 qubits on
# (and the d ≥ 3 ones visibly inaccurate already for moderate n).
#
# `include` this file right after `using SUNRepresentations`. It replaces the lowering step
# by one that imposes the lowering *and* raising equations in the least-squares sense, which
# is stable (SU(2) CGCs agree with the Racah formula to ~1e-14 for every j ≤ 40), and it keeps
# the CGCs in memory only, so that CGCs cached on disk by the unstable version are never used.
#
# The lowering routine below is adapted from SUNRepresentations.jl
# (https://github.com/QuantumKitHub/SUNRepresentations.jl), which is distributed under the
# MIT License, Copyright (c) 2020 maartenvd.
using LinearAlgebra
using SUNRepresentations

@eval SUNRepresentations begin
    function lower_weight_CGC!(CGC, s1::I, s2::I, s3::I) where {I <: SUNIrrep{N}} where {N}
        N123 = size(CGC, 4)
        T = eltype(CGC)

        Jm_list1 = annihilation(s1)
        Jm_list2 = annihilation(s2)
        Jm_list3 = annihilation(s3)

        map1 = weightmap(basis(s1))
        map2 = weightmap(basis(s2))
        map3 = weightmap(basis(s3))
        weights1 = map(weight, basis(s1))
        weights2 = map(weight, basis(s2))

        # reverse lexographic order: so all relevant parents should come earlier
        # and should thus have been solved
        w3list = sort(collect(keys(map3)); rev = true)

        # basis |m₁, m₂⟩ of the subspace of s1 ⊗ s2 with the same weight as w3 in s3
        wshift = div(sum(s1.I) + sum(s2.I) - sum(s3.I), N)
        productbases = Dict{NTuple{N, Int}, Vector{CartesianIndex{2}}}()
        productbasis(w3) = get!(productbases, w3) do
            m1m2list = CartesianIndex{2}[]
            for (w1, m1list) in map1
                for m2 in get(map2, w3 .- w1 .+ wshift, _emptyindexlist), m1 in m1list
                    push!(m1m2list, CartesianIndex(m1, m2))
                end
            end
            return m1m2list
        end
        lower(w, l) = Base.setindex(Base.setindex(w, w[l] - 1, l), w[l + 1] + 1, l + 1)
        Lrows = Int[]
        Lvals = T[]

        # Solving only the lowering equations J⁻₃ |m₃⟩ = (J⁻₁ ⊗ 𝕀 + 𝕀 ⊗ J⁻₂) |m₁, m₂⟩ weight by
        # weight is unstable: rounding errors in the components along other irreps s3′ ⊂ s1 ⊗ s2
        # get amplified at every step (for SU(2) by c_{J′}(M) / c_J(M) > 1 when J′ > J), and the
        # CGCs of large irreps become garbage. We therefore also impose the raising equations
        # (J⁺₁ ⊗ 𝕀 + 𝕀 ⊗ J⁺₂) X = Y J⁺₃ (with J⁺ = J⁻ᵀ) and solve both in the least-squares
        # sense, which damps instead of amplifies these errors. Restricted to the weight space
        # of w3, with Lₗ and Kₗ the blocks of J⁻ₗ on s1 ⊗ s2 and on s3 coming from the parent
        # weight w3 + αₗ, and Yₗ the (already computed) CGC block of that parent, this reads
        #   min_X ∑ₗ ‖X Kₗ - Lₗ Yₗ‖² + ‖Lₗᵀ X - Yₗ Kₗᵀ‖²
        # with normal equations X A + B X = R, A = ∑ₗ Kₗ Kₗᵀ, B = ∑ₗ Lₗ Lₗᵀ, R = 2 ∑ₗ Lₗ Yₗ Kₗᵀ.
        # A is positive definite (all states of s3 except its highest weight are reached by
        # lowering) and B is positive semidefinite, so the solution is unique.
        for w3 in view(w3list, 2:length(w3list))
            m3list = map3[w3]
            m1m2list = productbasis(w3)
            m1m2index = Dict(Im1m2 => i for (i, Im1m2) in enumerate(m1m2list))

            A = zeros(T, length(m3list), length(m3list))
            B = zeros(T, length(m1m2list), length(m1m2list))
            R = zeros(T, length(m1m2list), length(m3list), N123)
            for (l, (J⁻₁, J⁻₂, J⁻₃)) in enumerate(zip(Jm_list1, Jm_list2, Jm_list3))
                w3′ = Base.setindex(w3, w3[l] + 1, l)
                w3′ = Base.setindex(w3′, w3[l + 1] - 1, l + 1)
                m3′list = get(map3, w3′, _emptyindexlist)
                isempty(m3′list) && continue
                m1m2′list = productbasis(w3′)

                K = T[J⁻₃[m3, m3′] for m3 in m3list, m3′ in m3′list]
                mul!(A, K, K', true, true)
                YK = [T[CGC[Im1m2′, m3′, α] for Im1m2′ in m1m2′list, m3′ in m3′list] * K' for α in 1:N123]

                # Lₗ is very sparse: accumulate B and R column by column
                for (i′, Im1m2′) in enumerate(m1m2′list)
                    m1′, m2′ = Tuple(Im1m2′)
                    empty!(Lrows)
                    empty!(Lvals)
                    # apply J⁻₁
                    for m1 in get(map1, lower(weights1[m1′], l), _emptyindexlist)
                        v = J⁻₁[m1, m1′]
                        iszero(v) && continue
                        push!(Lrows, m1m2index[CartesianIndex(m1, m2′)])
                        push!(Lvals, v)
                    end
                    # apply J⁻₂
                    for m2 in get(map2, lower(weights2[m2′], l), _emptyindexlist)
                        v = J⁻₂[m2, m2′]
                        iszero(v) && continue
                        push!(Lrows, m1m2index[CartesianIndex(m1′, m2)])
                        push!(Lvals, v)
                    end
                    for (i1, v1) in zip(Lrows, Lvals)
                        for (i2, v2) in zip(Lrows, Lvals)
                            B[i2, i1] += v1 * v2
                        end
                        for α in 1:N123, j in axes(R, 2)
                            R[i1, j, α] += 2 * v1 * YK[α][i′, j]
                        end
                    end
                end
            end

            # solve X A + B X = R: in the eigenbasis A = U Diagonal(μ) Uᵀ of the small matrix A,
            # column k of X U solves a linear system with the positive definite matrix B + μₖ 𝕀
            μ, U = eigen!(Symmetric(A))
            RU = [view(R, :, :, α) * U for α in 1:N123]
            XU = similar(R)
            for (k, μₖ) in enumerate(μ)
                F = cholesky!(Symmetric(B + UniformScaling(μₖ)))
                for α in 1:N123
                    XU[:, k, α] = F \ RU[α][:, k]
                end
            end
            for α in 1:N123
                X = view(XU, :, :, α) * U'

                # fill in CGC
                # loop over X in column major order, CGC is hashmap anyways
                @inbounds for (j, m3) in enumerate(m3list), (i, Im1m2) in enumerate(m1m2list)
                    CGC[Im1m2, m3, α] = X[i, j]
                end
            end
        end
        return CGC
    end

    # keep CGCs in RAM only: never read from (or write to) the on-disk cache
    @noinline function _get_CGC(::Type{T}, @nospecialize(key)) where {T}
        d::SparseArray{T, 4} = get!(() -> _CGC(T, key...), CGC_CACHE, key)
        return d
    end
end

# drop anything computed before this file was included
empty!(SUNRepresentations.CGC_CACHE)
