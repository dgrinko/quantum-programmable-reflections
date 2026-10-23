# Numerics for *Quantum Programmable Reflections*

This repository contains the numerical code accompanying the paper

> E. Schoute, D. Grinko, Y. Subaşı, T. Volkoff, *Quantum Programmable Reflections*.

The code reproduces the numerical evidence used in Section 4.3,
**"Optimization of the Holevo information"**, and the associated conjectures.
Concretely, for a uniform ensemble of reflections $e^{i\pi|\psi\rangle\langle\psi|}$
applied to an $n$-copy entangled probe state, it computes the Holevo information
$S\left(\widetilde{\mathcal{T}(\rho)}\right)$ in the mixed Schur basis and compares
it to the upper bound $\log_2 D(n,d) = 2\log_2\binom{n+d-1}{d-1}$, where
$D(n,d)=\binom{n+d-1}{d-1}^2$ bounds the rank of the reflected ensemble.

## Contents

| Path | Language | What it does |
|------|----------|--------------|
| `lower_bound_optimization.ipynb` | Julia | Builds the block-diagonal reflected ensemble in the mixed Schur basis, assembles the transition matrix $K$ from $SU(d)$ Clebsch–Gordan coefficients, and **maximizes** the Holevo information over the probe-state symmetry weights $\{q_\lambda\}$ with `Ipopt`/`JuMP` (a convex problem). Includes the qubit case: a check of Conjecture 4.6 up to $n = 80$ and the comparison with the Yang–Renner–Chiribella (YRC) probe state. |
| `symm_subspace_holevo_information_d_gtreq_3.ipynb` | Julia | Computes the Holevo information for the *suboptimal* probe state supported only on the symmetric irrep $\theta_n$ (i.e. $q_\lambda=\delta_{\lambda,\theta_n}$). In this case the reduced Clebsch–Gordan coefficients have the closed combinatorial form of Vilenkin & Klimyk, Vol. 3, §18.2.6 Eq. (6), so the calculation reaches much larger $n$ and $d$ ($d=3,\dots,20$, $n$ up to 40). The closed form is evaluated in 256-bit arithmetic because its alternating sums cancel strongly. |
| `holevo_information_data_d_3_20_n_1_41.pkl` | — | Pre-computed output of the symmetric-subspace notebook: an $18\times 40$ array of Holevo-information values for $d=3,\dots,20$ and $n=1,\dots,40$. |
| `log_log_holevo_information_plot.ipynb` | Python | Loads the `.pkl` data and produces the log–log plot of the gap $1-r$ (with $r$ the ratio of the Holevo information to the upper bound) versus $n$. This is Figure 4 of the paper, drawn at its printed one-column size (`plots/ratfig.pdf`); the gap scales like $c(d)\,n^{-\alpha}$ with fitted exponents $\alpha\approx 0.2$–$0.3$. |
| `stable_cgc.jl` | Julia | Numerically stable replacement for the Clebsch–Gordan routine of `SUNRepresentations.jl` (see below). Included by `lower_bound_optimization.ipynb`. |
| `plot_style.jl` | Julia | Shared plot style for the two Julia notebooks: colours, the blue ramp for curves indexed by $d$, direct labels at curve ends and log-axis ticks. |
| `plots/` | PDF | The figures produced by the notebooks; see [Figures](#figures) below. |
| `Project.toml`, `Manifest.toml` | — | Pinned Julia environment for the two Julia notebooks above. |
| `requirements.txt` | — | Python dependencies for `log_log_holevo_information_plot.ipynb`. |

## Figures

The notebooks write their figures to `plots/` (each run overwrites its own files). Throughout,
$r = S / (\text{upper bound})$ is the Holevo information $S$ of the reflected ensemble divided by
its upper bound: $\log_2\binom{n+2}{2}$ for $d = 2$ and $\log_2 D(n,d) = 2\log_2\binom{n+d-1}{d-1}$
for $d \ge 3$. Curves for different $d$ are blue, light for small $d$ and dark for large $d$, and are
labelled at their ends; gray lines are the upper bound.

| File | What it shows |
|------|---------------|
| `holevo_ratio_optimized.pdf` | The Holevo information maximized over the probe weights $\{q_\lambda\}$, as the ratio $r$ to the upper bound, against the number of copies $n$, for $d = 2, \dots, 11$ ($n$ up to `max_copies[d]`, which is limited by the cost of the CG coefficients). For $d = 2$ the ratio is exactly 1 for every $n$: the flat spectrum of Conjecture 4.6 is reached, as in Corollary 4.7. For $d \ge 3$ it stays below 1, increases with $n$ after a small dip at small $n$ for $d = 3, 4$, and decreases with $d$: at the largest $n$ computed it is 0.99 for $d = 3$ and 0.49 for $d = 11$. |
| `holevo_gain_optimized_vs_symmetric.pdf` | How much optimizing over $\{q_\lambda\}$ gains over the probe supported only on the symmetric irrep ($q_\lambda = \delta_{\lambda,\theta_n}$): the difference of the two Holevo informations in units of the upper bound. Only $d = 2, 3, 4$ appear, because for $d \ge 5$ (in the range computed) the optimum *is* the symmetric irrep. The gain is at most about 0.03 of the bound. |
| `qubit_probe_weights.pdf` | The qubit probe weights $q_j$ against the spin $j$ for $n = 80$: in blue the solution of the linear system of Conjecture 4.6, which gives a flat spectrum and hence the maximal Holevo information $\log_2\binom{n+2}{2} = 11.697$ bits; in orange the Yang–Renner–Chiribella weights $q_j \propto (2j+1)^2$, which give 11.607 bits. The flat-spectrum weights are concentrated much more strongly on the largest spins. |
| `symmetric_probe_holevo_vs_bound.pdf` | The Holevo information $S$ of the symmetric-irrep probe (from the closed-form reduced CG coefficients of Vilenkin & Klimyk) and the upper bound $\log_2 D(n,d)$, in bits, against $n = 1, \dots, 22$ for $d = 10$ (the notebook's defaults `copies = 1:22`, `dims = 10:10`). Both are measured from zero, so the vertical distance shows how far the probe is from saturating the bound. |
| `symmetric_probe_gap.pdf` | The same data as the gap $1 - r$ on log–log axes. Both axes span the same number of decades, so the slope is not exaggerated. A power law $c\,n^{-\alpha}$ would be a straight line; for $d = 10$ the gap falls from 0.50 at $n = 1$ to 0.28 at $n = 22$, and the curve bends downward (its local slope steepens from $-0.08$ to $-0.31$). |
| `ratfig.pdf` | Figure 4 of the paper: the gap $1 - r$ of the symmetric-irrep probe for $d = 4, \dots, 20$ and $n = 1, \dots, 40$, read from `holevo_information_data_d_3_20_n_1_41.pkl`, on log–log axes. The colour bar gives $d$, and the extremes $d = 4$ and $d = 20$ are labelled. At $n = 40$ the gap ranges from about 0.04 ($d = 4$) to 0.47 ($d = 20$). $d = 3$ is omitted because $r = 1$ exactly at $n = 1$. The notebook also fits $c(d)\,n^{-\alpha}$ for $n \ge 11$, giving $\alpha \approx 0.2$–$0.3$. Drawn at its printed one-column size. |

## Clebsch–Gordan coefficients

`SUNRepresentations.jl` (≤ v0.4) computes the CGCs of $\lambda\otimes\mu\to\nu$ by fixing the
highest-weight vectors and then applying lowering operators weight by weight. That recursion
amplifies rounding errors exponentially in the length of the weight chain: the $SU(2)$ CGCs of
spin $j\otimes j$ are off by $O(1)$ for $j\approx 20$, the qubit transition matrix breaks down from
$n\approx 70$ (the solution of the linear system of Conjecture 4.6 stops being a probability
vector), and for $d=3$, $n=16$ the CGCs are already off by $\sim 10^{-6}$.
`stable_cgc.jl` replaces the lowering step by one that imposes the lowering *and* the raising
equations in the least-squares sense, which is stable: the $SU(2)$ CGCs agree with the Racah
formula to $\sim 10^{-14}$ for all $j\le 40$, and the $SU(3)$–$SU(5)$ CGCs agree with a 256-bit
computation to $\sim 10^{-14}$. It also keeps the CGCs in memory only, so that CGCs stored in
`SUNRepresentations`' on-disk cache by the unstable version are never used.

## Requirements

### Julia (`lower_bound_optimization.ipynb`, `symm_subspace_holevo_information_d_gtreq_3.ipynb`)

The environment is pinned in `Project.toml` / `Manifest.toml`. Each Julia notebook
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
optimization), `Combinatorics`, `Memoize`, `Plots`.

You also need a Julia Jupyter kernel: `import Pkg; Pkg.add("IJulia")`.

### Python (`log_log_holevo_information_plot.ipynb`)

```bash
pip install -r requirements.txt
```

## Running

```bash
jupyter lab        # then open any of the notebooks and run the cells top to bottom
```

The Julia notebooks are CPU-bound. In `lower_bound_optimization.ipynb` the cost is dominated by
evaluating $SU(d)$ Clebsch–Gordan coefficients, which is why the general optimization is limited
to small $(n,d)$, while the symmetric-subspace specialization reaches the ranges stored in the
`.pkl` file.

## Notation map (code ↔ paper)

- `n` — number of program copies; `d` — dimension of the reflected state.
- `λirreps(n,d)` — partitions $\lambda\vdash_d n$ labelling the $U(d)$ irreps that appear in the probe state; for $d = 2$ the probe weights are labelled by the spin $j$ instead (`qubit_spins(n)`).
- `ν_irrep(k,d)` — the blocks $\nu_k=(k,\vec 0_{d-2},-k)$ of the reflected ensemble; for $d = 2$ the blocks are the spins $J\equiv n \pmod 2$ (`qubit_blocks(n)`).
- `LinSystem(n,d)` — returns `(K, target)` where `K` maps the weights $q_\lambda$ to the block weights $p_\nu = d_\nu\lVert v_\nu\rVert^2$, and `target` is the flat distribution that would saturate the upper bound ($d_{\nu_k}/D(n,d)$ for $d\ge 3$, $(2J+1)/\binom{n+2}{2}$ for $d=2$).
- `maximize_entropy(n,d)` — solves $\max_{\{q_\lambda\}}\; -\sum_\nu d_\nu \lVert v_\nu\rVert^2 \log_2\lVert v_\nu\rVert^2$.
- `D(n,d) = binomial(n+d-1,d-1)^2` — the quantity $D(n,d)$ of Eq. (195); `holevo_bound(n,d)` is the upper bound $\log_2 D(n,d)$ for $d\ge 3$ and $\log_2\binom{n+2}{2}$ for $d=2$.
- `symmetric_subspace_entropy(n,d)` — the Holevo information for the probe state $q_\lambda=\delta_{\lambda,\theta_n}$ (built from the Vilenkin–Klimyk reduced CG coefficients, `ScalarFactor`).
