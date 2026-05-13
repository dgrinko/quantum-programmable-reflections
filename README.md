# Numerics for *Quantum Programmable Reflections*

This repository contains the numerical code accompanying the paper

> E. Schoute, D. Grinko, Y. Subaşı, T. Volkoff, *Quantum Programmable Reflections*.

The code reproduces the numerical evidence used in the section
**"Optimization of the Holevo information"** and the associated conjectures.
Concretely, for a uniform ensemble of reflections $e^{i\pi|\psi\rangle\langle\psi|}$
applied to an $n$-copy entangled probe state, it computes the Holevo information
$S\left(\widetilde{\mathcal{T}(\rho)}\right)$ in the mixed Schur basis and compares
it to the information-theoretic upper bound $2\log_2 D(n,d)$, where
$D(n,d)=\binom{n+d-1}{d-1}^2$.

## Contents

| Path | Language | What it does |
|------|----------|--------------|
| `lower_bound_optimization.ipynb` | Julia | Builds the block-diagonal reflected ensemble in the mixed Schur basis, assembles the transition matrix $K$ from $U(d)$ Clebsch–Gordan coefficients, and **maximizes** the Holevo information over the probe-state symmetry weights $\{q_\lambda\}$ with `Ipopt`/`JuMP`. Includes the qubit special case and the comparison with the Yang–Renner–Chiribella (YRC) probe state. |
| `symm_subspace_holevo_information_d_gtreq_3.ipynb` | Julia | Computes the Holevo information for the *suboptimal* probe state supported only on the symmetric irrep $\theta_n$ (i.e. $q_\lambda=\delta_{\lambda,\theta_n}$). In this case the reduced Clebsch–Gordan coefficients have the closed combinatorial form of Vilenkin & Klimyk, Vol. 3, §18.2.6 Eq. (6), so the calculation reaches much larger $n$ and $d$ ($d=3,\dots,20$, $n$ up to ~40). |
| `holevo_information_data_d_3_20_n_1_41.pkl` | — | Pre-computed output of the symmetric-subspace notebook: an $18\times 40$ array of Holevo-information values for $d=3,\dots,20$ and $n=1,\dots,40$. |
| `log_log_holevo_information_plot.ipynb` | Python | Loads the `.pkl` data and produces the log–log plot of the gap $1-r$ (with $r$ the ratio of the Holevo information to the upper bound) versus $n$. This is the figure 4 in the paper; the gap scales like $c(d)\,n^{-\alpha}$ with $\alpha\approx 0.3$ roughly independent of $d$. |
| `Project.toml`, `Manifest.toml` | — | Pinned Julia environment for the two Julia notebooks above. |
| `requirements.txt` | — | Python dependencies for `log_log_holevo_information_plot.ipynb`. |

## Requirements

### Julia (`lower_bound_optimization.ipynb`, `symm_subspace_holevo_information_d_gtreq_3.ipynb`)

The environment is pinned in `Project.toml` / `Manifest.toml`.  Each Julia notebook
activates it automatically in its first cell (it searches the current directory and its
parents for `Project.toml`), so you only need to install the dependencies once — from this
directory:

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Key packages: [`SUNRepresentations.jl`](https://github.com/QuantumKitHub/SUNRepresentations.jl)
(GT patterns and $SU(d)$ Clebsch–Gordan coefficients), `JuMP` + `Ipopt` (entropy
optimization), `Combinatorics`, `Memoize`, `SpecialFunctions`, `Plots`.

You also need a Julia Jupyter kernel: `import Pkg; Pkg.add("IJulia")`.

### Python (`log_log_holevo_information_plot.ipynb`)

```bash
pip install -r requirements.txt
```

## Running

```bash
jupyter lab        # then open any of the notebooks and run the cells top to bottom
```

The Julia notebooks are CPU-bound: the cost is dominated by evaluating $SU(d)$
Clebsch–Gordan coefficients, which is why the general optimization
(`lower_bound_optimization.ipynb`) is limited to small $(n,d)$ while the
symmetric-subspace specialization reaches the ranges stored in the `.pkl` file.

## Notation map (code ↔ paper)

- `n` — number of program copies; `d` — dimension of the reflected state.
- `λirreps(n,d)` — partitions $\lambda\vdash_d n$ labelling the $SU(d)$ irreps that appear in the probe state.
- `mixed_irreps(n,d)` / `ν_irrep(k,d)` — the irreps $\nu_k=(k,\vec 0_{d-2},-k)$ of the partially transposed permutation algebra.
- `LinSystem(n,d)` — returns `(K, target)` where `K` maps the weights $q_\lambda$ to the squared norms $\lVert v_\nu\rVert^2$, and `target` is the distribution $d_{\nu_k}/D(n,d)$ that would saturate the upper bound.
- `maximize_entropy(n,d)` — solves $\max_{\{q_\lambda\}}\; -\sum_\nu d_\nu \lVert v_\nu\rVert^2 \log_2\lVert v_\nu\rVert^2$.
- `D(n,d) = binomial(n+d-1,d-1)^2` — the quantity $D(n,d)$; the upper bound on the Holevo information is $2\log_2 D(n,d)$.
- `symmetric_subspace_entropy(n,d)` — the Holevo information for the probe state $q_\lambda=\delta_{\lambda,\theta_n}$ (built from the Vilenkin–Klimyk reduced CG coefficients, `ScalarFactor`).
