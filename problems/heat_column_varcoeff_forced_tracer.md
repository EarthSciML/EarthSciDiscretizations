---
esm: "2.3.0"
license: "MIT"
name: "heat_column_varcoeff_forced_tracer"
references:
  - citation: "Roache, P. J. (2002). Code verification by the method of manufactured solutions. J. Fluids Eng. 124(1), 4-10."
    doi: "10.1115/1.1436090"
tags:
  - "esd:problem"
  - "grid:column_nonuniform_1d_tracer"
  - "equation:heat"
  - "mms"
---

```esm:description
MMS driver for grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer with a GENUINELY VARIABLE interface coefficient and UNEQUAL, time-dependent per-tracer face fluxes -- the tracer-axis twin of problems/heat_column_varcoeff_forced, built so that the index slips the constant-K, equal-flux problems/heat_column_varcoeff_flux_bc_tracer cannot see are O(1) errors here. u_t = d/dz(K du/dz) + S on the unit column [0, 1] for TWO tracers, u shaped [lev, tracer] with NTRACER = 2 bound at the rule's import edge, on the grid's shared quadratic-stretch MMS column, with the shared interface coefficient K(z) = 1 + z (K[kn] = 1 + ze[kn]). Tracer s carries its own manufactured solution u_s(z,t) = exp(-t) ( p_s z + q_s sin(pi z) + r_s cos(pi z) ) with p = (1, -2), q = (1, 0.5), r = (0, 1) (spelled p_s = 4 - 3s, q_s = 1.5 - s/2, r_s = s - 1; model-local templates mms_p / mms_q / mms_r / mms_exact). Then K u_z = (1+z) exp(-t) ( p_s + pi q_s cos(pi z) - pi r_s sin(pi z) ) (mms_flux), so the prescribed face fluxes are the [tracer] observeds kdudz_bot[s] = mms_flux(0, s, t) = exp(-t) (p_s + pi q_s) and kdudz_top[s] = mms_flux(1, s, t) = 2 exp(-t) (p_s - pi q_s): about (4.14, -0.43) exp(-t) at the surface and (-4.28, -7.14) exp(-t) at the top -- different at the two faces of each tracer and different between the tracers. The forcing is S_s = u_t - d/dz(K u_z) = -u_s - exp(-t) ( p_s + q_s (pi cos(pi z) - pi^2 (1+z) sin(pi z)) - r_s (pi sin(pi z) + pi^2 (1+z) cos(pi z)) ) (mms_forcing), evaluated at the runtime layer centres zc. Detected slips: reading K[n-1] for K[n] (K varies by dz across an interface, an O(dz) relative change of every flux, an O(1) change of the divergence); swapping the two face fluxes (bottom and top differ in sign and size); reading another tracer's face flux or field (the tracers' solutions and fluxes differ). Nothing is discretely exact and the column totals are not conserved (the faces carry net flux and S is nonzero), so the t = 0.1 relative L2 error over the [lev, tracer] field is the full second-order spatial error of the scheme. The mesh follows the column grid's edge-array contract (ze = mms_column_ze of grids/column_nonuniform_1d/mesh/mms_column.esm; zc and dz grid-derived). The ic applies the closed per-layer centre template mms_column_zc (grids/column_nonuniform_1d/mesh/mms_column_centre.esm) inside its aggregate over [lev, tracer]; the 6.6.5 test references spell the layer centre inline from literals and NLEV, because the inline-test oracle does not resolve template references.
```

::model[Heat]

:var[$K$]{description="Interface diffusivity K = 1 + z on the interfaces, K[kn] = 1 + ze[kn], shared by both tracers." shape='["lev_nodes"]' units="m^2/s"}

:var[$S$]{description="Manufactured forcing mms_forcing(zc, s, t), per tracer." shape='["lev", "tracer"]' units="K/s"}

:var[$\mathit{dz}$]{description="Layer thicknesses, grid-derived by column_cell_thickness." shape='["lev"]' units="m"}

:var[$\mathit{kdudz\_bot}$]{description="Prescribed K du/dz on the surface interface, per tracer: mms_flux(0, s, t) = exp(-t) (p_s + pi q_s)." shape='["tracer"]' units="K*m/s"}

:var[$\mathit{kdudz\_top}$]{description="Prescribed K du/dz on the top interface, per tracer: mms_flux(1, s, t) = 2 exp(-t) (p_s - pi q_s)." shape='["tracer"]' units="K*m/s"}

:var[$u$]{description="Two tracers over the lev (layer) axis." shape='["lev", "tracer"]' units="K"}

:var[$\mathit{zc}$]{description="Layer centres, grid-derived by column_cell_center." shape='["lev"]' units="m"}

:var[$\mathit{ze}$]{description="Interface coordinates of the quadratic-stretch column, the shared mms_column_ze." shape='["lev_nodes"]' units="m"}

$$
u(t_{0}) = \left[\sum\operatorname{mms\_exact}(\operatorname{mms\_column\_zc}(i = i), s, 0)\right]_{i \in \mathit{lev}, s \in \mathit{tracer}}
$$

$$
\frac{\mathrm{d} u}{\mathrm{d} t} = \frac{\partial}{\partial\mathit{lev}} \left(K \frac{\partial u}{\partial\mathit{lev}}\right) + S
$$

$$
K = 1 + \mathit{ze}
$$

$$
S = \left[\sum\limits\operatorname{mms\_forcing}(\mathit{zc}_{i}, s, t)\right]_{i \in \mathit{lev}, s \in \mathit{tracer}}
$$

$$
\mathit{kdudz\_bot} = \left[\sum\operatorname{mms\_flux}(0, s, t)\right]_{s \in \mathit{tracer}}
$$

$$
\mathit{kdudz\_top} = \left[\sum\operatorname{mms\_flux}(1, s, t)\right]_{s \in \mathit{tracer}}
$$

$$
\mathit{ze} = \operatorname{mms\_column\_ze}
$$

$$
\mathit{zc} = \operatorname{column\_cell\_center}(\mathit{ze} = \mathit{ze})
$$

$$
\mathit{dz} = \operatorname{column\_cell\_thickness}(\mathit{ze} = \mathit{ze})
$$

:::template{description="Manufactured solution of tracer s: u_s(z,t) = exp(-t) (p_s z + q_s sin(pi z) + r_s cos(pi z))."}
$$
\operatorname{mms\_exact}(z, s, t) \coloneqq \exp(-t) \left(\operatorname{mms\_p}(s) z + \operatorname{mms\_q}(s) \sin(3.141592653589793 z) + \operatorname{mms\_r}(s) \cos(3.141592653589793 z)\right)
$$
:::

:::template{description="Diffusive flux of tracer s, K u_z = (1+z) exp(-t) (p_s + pi q_s cos(pi z) - pi r_s sin(pi z)), positive upward."}
$$
\operatorname{mms\_flux}(z, s, t) \coloneqq \left(1 + z\right) \exp(-t) \left(\operatorname{mms\_p}(s) + 3.141592653589793 \operatorname{mms\_q}(s) \cos(3.141592653589793 z) - 3.141592653589793 \operatorname{mms\_r}(s) \sin(3.141592653589793 z)\right)
$$
:::

:::template{description="Manufactured forcing of tracer s: S_s = u_t - d/dz((1+z) u_z) = -u_s - exp(-t) ( p_s + q_s (pi cos(pi z) - pi^2 (1+z) sin(pi z)) - r_s (pi sin(pi z) + pi^2 (1+z) cos(pi z)) )."}
$$
\operatorname{mms\_forcing}(z, s, t) \coloneqq -\operatorname{mms\_exact}(z, s, t) - \exp(-t) \left(\operatorname{mms\_p}(s) + \operatorname{mms\_q}(s) \left(3.141592653589793 \cos(3.141592653589793 z) - 9.869604401089358 \left(1 + z\right) \sin(3.141592653589793 z)\right) - \operatorname{mms\_r}(s) \left(3.141592653589793 \sin(3.141592653589793 z) + 9.869604401089358 \left(1 + z\right) \cos(3.141592653589793 z)\right)\right)
$$
:::

:::template{description="Ramp coefficient of tracer s: p_s = 4 - 3s, (1, -2)."}
$$
\operatorname{mms\_p}(s) \coloneqq 4 - 3 s
$$
:::

:::template{description="Sine coefficient of tracer s: q_s = 1.5 - s/2, (1, 0.5)."}
$$
\operatorname{mms\_q}(s) \coloneqq 1.5 - 0.5 s
$$
:::

:::template{description="Cosine coefficient of tracer s: r_s = s - 1, (0, 1)."}
$$
\operatorname{mms\_r}(s) \coloneqq s - 1
$$
:::

```esm:expression_template_imports
- bindings:
    NTRACER: 2
  ref: "../grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer.md"
- ref: "../grids/column_nonuniform_1d/mesh/mms_column.esm"
- ref: "../grids/column_nonuniform_1d/mesh/mms_column_centre.esm"
```

```esm:tests
- assertions:
    - expected: 0
      reduce: "L2_error"
      reference: $\left[\sum1 \left(\left(4 - 3 s\right) \left(0.5 \left(\left(0.6 \frac{i - 1}{\mathit{NLEV}} + 0.4 \frac{i - 1}{\mathit{NLEV}} \frac{i - 1}{\mathit{NLEV}}\right) + \left(0.6 \frac{i}{\mathit{NLEV}} + 0.4 \frac{i}{\mathit{NLEV}} \frac{i}{\mathit{NLEV}}\right)\right)\right) + \left(1.5 - 0.5 s\right) \sin(3.141592653589793 \left(0.5 \left(\left(0.6 \frac{i - 1}{\mathit{NLEV}} + 0.4 \frac{i - 1}{\mathit{NLEV}} \frac{i - 1}{\mathit{NLEV}}\right) + \left(0.6 \frac{i}{\mathit{NLEV}} + 0.4 \frac{i}{\mathit{NLEV}} \frac{i}{\mathit{NLEV}}\right)\right)\right)) + \left(s - 1\right) \cos(3.141592653589793 \left(0.5 \left(\left(0.6 \frac{i - 1}{\mathit{NLEV}} + 0.4 \frac{i - 1}{\mathit{NLEV}} \frac{i - 1}{\mathit{NLEV}}\right) + \left(0.6 \frac{i}{\mathit{NLEV}} + 0.4 \frac{i}{\mathit{NLEV}} \frac{i}{\mathit{NLEV}}\right)\right)\right))\right)\right]_{i \in \mathit{lev}, s \in \mathit{tracer}}$
      time: 0
      tolerance:
        abs: 1.0e-9
      variable: "u"
    - expected: 0
      reduce: "L2_error"
      reference: $\left[\sum0.9048374180359595 \left(\left(4 - 3 s\right) \left(0.5 \left(\left(0.6 \frac{i - 1}{\mathit{NLEV}} + 0.4 \frac{i - 1}{\mathit{NLEV}} \frac{i - 1}{\mathit{NLEV}}\right) + \left(0.6 \frac{i}{\mathit{NLEV}} + 0.4 \frac{i}{\mathit{NLEV}} \frac{i}{\mathit{NLEV}}\right)\right)\right) + \left(1.5 - 0.5 s\right) \sin(3.141592653589793 \left(0.5 \left(\left(0.6 \frac{i - 1}{\mathit{NLEV}} + 0.4 \frac{i - 1}{\mathit{NLEV}} \frac{i - 1}{\mathit{NLEV}}\right) + \left(0.6 \frac{i}{\mathit{NLEV}} + 0.4 \frac{i}{\mathit{NLEV}} \frac{i}{\mathit{NLEV}}\right)\right)\right)) + \left(s - 1\right) \cos(3.141592653589793 \left(0.5 \left(\left(0.6 \frac{i - 1}{\mathit{NLEV}} + 0.4 \frac{i - 1}{\mathit{NLEV}} \frac{i - 1}{\mathit{NLEV}}\right) + \left(0.6 \frac{i}{\mathit{NLEV}} + 0.4 \frac{i}{\mathit{NLEV}} \frac{i}{\mathit{NLEV}}\right)\right)\right))\right)\right]_{i \in \mathit{lev}, s \in \mathit{tracer}}$
      time: 0.1
      tolerance:
        abs: 0.00025
      variable: "u"
  description: "Two manufactured solutions u_s(z,t) = exp(-t) (p_s z + q_s sin(pi z) + r_s cos(pi z)), p = (1, -2), q = (1, 0.5), r = (0, 1), with K = 1 + z on the interfaces, the matching per-tracer forcing and the per-tracer prescribed face fluxes K u_z at z = 0 and z = 1 (unequal at the two faces and between tracers), on the quadratic-stretch column. The t=0 assertion pins the ic wiring; the t=0.1 assertion bounds the relative L2 error over the [lev, tracer] field at the default NLEV = 64 (second-order; per-resolution values pinned in tests/conformance/convergence/heat_column_varcoeff_forced_tracer/golden/errors.json)."
  id: "mms_column_forced_tracer"
  time_span:
    end: 0.1
    start: 0
```
