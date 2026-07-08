# PPM reconstruction — MMS convergence (Layer B)

Manufactured solution: `f(x) = sin(2πx) + 0.3·cos(4πx)` on `[0, 1]` with periodic
boundary conditions.

Cell averages are computed from the analytic antiderivative (no quadrature
error), so the residual is purely the unlimited PPM reconstruction error
(CW84 §1, eqs. 1.6–1.10). The parabola is sampled at 5 sub-points per cell
and the L∞ error is taken over all samples and all cells.

| `n`  | `dx`        | L∞ error    | observed order |
| ---: | :---------- | :---------- | :------------- |
|  16  | 0.0625000   | 2.9052e-03  | —              |
|  32  | 0.0312500   | 2.2013e-04  | 3.722          |
|  64  | 0.0156250   | 2.1838e-05  | 3.334          |
| 128  | 0.0078125   | 2.6265e-06  | 3.056          |

Theoretical asymptotic order: **3.0** (CW84 §1).
Acceptance threshold: **min(observed order) ≥ 2.7**.

The observed orders converge from above toward 3 — the leading-order error
constant is dominant on the coarsest grid and the interpolation is locally
4th-order at edge values, so the per-cell L∞ peaks slightly faster than the
asymptote at low resolutions before settling onto the 3rd-order trajectory.

These numbers are reproduced by the `ppm_reconstruction MMS convergence`
test item in `test/test_ppm_reconstruction_rule.jl`.
