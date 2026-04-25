# WENO-5 advection — MMS convergence (Layer B)

Manufactured solution: `f(x) = sin(2πx + 1)` on `[0, 1]` with periodic boundary
conditions. Cell averages are integrated analytically from the antiderivative
so the residual is purely the WENO5 left-biased reconstruction error at the
right cell edge `q_{i+1/2}^L`.

The phase shift places both critical points of `f` away from every dyadic cell
face at `n ∈ {32, 64, 128, 256}`. Without it, the WENO5-JS nonlinear weights
fail to recover the linear weights `ω_k → d_k` at `f'(x_{i+1/2}) = 0`, producing
a spurious accuracy dip typically reported as "3rd-order at critical points"
(Henrick, Aslam & Powers, JCP 2005). The phase shift sidesteps this regime so
the rule's advertised 5th-order flux reconstruction is exercised cleanly.

| `n`  | `dx`         | L∞ error     | observed order |
| ---: | :----------- | :----------- | :------------- |
|  32  | 0.03125000   | 3.4137e-05   | —              |
|  64  | 0.01562500   | 1.0656e-06   | 5.002          |
| 128  | 0.00781250   | 3.3254e-08   | 5.002          |
| 256  | 0.00390625   | 1.0377e-09   | 5.002          |

Theoretical asymptotic order: **5.0** (Jiang & Shu 1996, smooth regions).
Acceptance threshold: **min(observed order) ≥ 4.7**.

The 4.7 floor leaves headroom for the small accuracy hit from the
`ε = 1e-6` regularisation of the nonlinear weights; empirically the orders
sit essentially at 5.0 for this solution.

These numbers are reproduced by the
`weno5_advection MMS convergence` test item in
`test/test_weno5_advection_rule.jl`.
