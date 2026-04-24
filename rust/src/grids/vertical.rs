//! Vertical grid family (1D column) — accessor runtime.
//!
//! Conforms to the cross-binding contract in `docs/GRIDS_API.md` §2.5, §3.3, §7.
//!
//! Per the 2026-04-20 scope correction (mayor), the `.esm` lowering is a small
//! declarative config (family + coordinate kind + interface levels + optional
//! hybrid coefficients), NOT a serialized geometry blob. Cell centers and
//! widths are derived on demand from the interface `levels` via pure
//! arithmetic.
//!
//! Mirrors the Python / Julia / TypeScript sibling bindings
//! (`python/src/earthsci_toolkit/grids/vertical.py`,
//! `src/grids/vertical.jl`, `typescript/src/grids/vertical.ts`). Indices are
//! 0-based in Rust (Julia is 1-based); semantics and accessor outputs are
//! otherwise identical.
//!
//! Supported coordinate kinds:
//!
//! | Kind                 | Value domain                 | Required options      |
//! |----------------------|------------------------------|-----------------------|
//! | `Sigma`              | [0, 1]; 1 = surface, 0 = top | `nz` or `levels`      |
//! | `Eta`                | hybrid sigma-pressure (NCAR) | `ak`, `bk`            |
//! | `Z`                  | geometric altitude (m)       | `levels` (increasing) |
//! | `Theta`              | potential temperature (K)    | `levels` (increasing) |
//! | `HybridSigmaTheta`   | blended sigma→theta          | `nz` or `levels`      |
//! | `ZStar`              | generalized height           | `levels` (increasing) |

use serde_json::{json, Value};

use crate::{Dtype, Grid, GridError, Result};

const API_VERSION: &str = "1.0.0";

/// Vertical coordinate kind. See module docs for semantics.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VerticalCoordinate {
    Sigma,
    Eta,
    Z,
    Theta,
    HybridSigmaTheta,
    ZStar,
}

impl VerticalCoordinate {
    /// Canonical wire-form name (snake_case) used in `.esm` output and schema.
    pub fn name(self) -> &'static str {
        match self {
            VerticalCoordinate::Sigma => "sigma",
            VerticalCoordinate::Eta => "eta",
            VerticalCoordinate::Z => "z",
            VerticalCoordinate::Theta => "theta",
            VerticalCoordinate::HybridSigmaTheta => "hybrid_sigma_theta",
            VerticalCoordinate::ZStar => "z_star",
        }
    }

    /// Parse the canonical wire-form name.
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "sigma" => VerticalCoordinate::Sigma,
            "eta" => VerticalCoordinate::Eta,
            "z" => VerticalCoordinate::Z,
            "theta" => VerticalCoordinate::Theta,
            "hybrid_sigma_theta" => VerticalCoordinate::HybridSigmaTheta,
            "z_star" => VerticalCoordinate::ZStar,
            _ => return None,
        })
    }
}

/// Scalar metric fields on a vertical grid. See [`VerticalGrid::metric_eval`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VerticalMetricName {
    /// Layer thickness (native units, always positive).
    Dz,
    /// Cell-centre coordinate value (native units).
    Z,
    /// Sigma at cell centre. Valid only for sigma-like coordinates
    /// (`Sigma`, `HybridSigmaTheta`, `Eta`).
    Sigma,
    /// Reference pressure at cell centre (`p = ak + bk * p0` averaged across
    /// the layer's two interfaces). Requires hybrid coefficients.
    Pressure,
    /// Hybrid `ak` averaged across the layer's two interfaces.
    Ak,
    /// Hybrid `bk` averaged across the layer's two interfaces.
    Bk,
}

impl VerticalMetricName {
    /// Parse the canonical wire-form name.
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "dz" => VerticalMetricName::Dz,
            "z" => VerticalMetricName::Z,
            "sigma" => VerticalMetricName::Sigma,
            "pressure" => VerticalMetricName::Pressure,
            "ak" => VerticalMetricName::Ak,
            "bk" => VerticalMetricName::Bk,
            _ => return None,
        })
    }
}

/// Axis-aligned vertical neighbours of a layer. Top and bottom layers drop
/// the out-of-range side (`None`).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub struct VerticalNeighbors {
    /// Layer `k-1` (closer to surface), if any.
    pub down: Option<usize>,
    /// Layer `k+1` (closer to top), if any.
    pub up: Option<usize>,
}

/// 1D vertical column.
///
/// Stores only the declarative parameters plus the pre-tabulated interface
/// `levels` (length `nz + 1`); `centers` and `widths` are cached from `levels`
/// at build time so accessor calls are `O(1)`. Per the §6 schema, no
/// per-cell derived arrays appear in the `.esm` wire form.
#[derive(Clone, Debug)]
pub struct VerticalGrid {
    coordinate: VerticalCoordinate,
    levels: Vec<f64>,
    centers: Vec<f64>,
    widths: Vec<f64>,
    ak: Vec<f64>,
    bk: Vec<f64>,
    p0: f64,
    transition: Option<f64>,
    ghosts: usize,
    dtype: Dtype,
}

impl VerticalGrid {
    /// Vertical coordinate kind.
    #[inline]
    pub fn coordinate(&self) -> VerticalCoordinate {
        self.coordinate
    }

    /// Interface values (length `nz + 1`).
    #[inline]
    pub fn levels(&self) -> &[f64] {
        &self.levels
    }

    /// Cell-centre values (length `nz`).
    #[inline]
    pub fn centers(&self) -> &[f64] {
        &self.centers
    }

    /// Layer thicknesses (length `nz`, always positive).
    #[inline]
    pub fn widths(&self) -> &[f64] {
        &self.widths
    }

    /// Hybrid `ak` at interfaces (empty when the coordinate has no hybrid
    /// coefficients).
    #[inline]
    pub fn ak(&self) -> &[f64] {
        &self.ak
    }

    /// Hybrid `bk` at interfaces (empty when the coordinate has no hybrid
    /// coefficients).
    #[inline]
    pub fn bk(&self) -> &[f64] {
        &self.bk
    }

    /// Reference surface pressure (Pa).
    #[inline]
    pub fn p0(&self) -> f64 {
        self.p0
    }

    /// `hybrid_sigma_theta`-only transition sigma in `(0, 1)`, if set.
    #[inline]
    pub fn transition(&self) -> Option<f64> {
        self.transition
    }

    /// Halo layer width.
    #[inline]
    pub fn ghosts(&self) -> usize {
        self.ghosts
    }

    /// Number of layers (== `n_cells`).
    #[inline]
    pub fn nz(&self) -> usize {
        self.centers.len()
    }

    /// Arity — always `1` for the vertical family.
    #[inline]
    pub fn ndim(&self) -> usize {
        1
    }

    /// Cell count (`nz`).
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nz()
    }

    /// Vertex count (`nz + 1`).
    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.nz() + 1
    }

    /// Edge count (`nz`).
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.nz()
    }

    /// Topology string (always `"column"` for this family per §7).
    #[inline]
    pub fn topology(&self) -> &'static str {
        "column"
    }

    /// Scalar cell-centre value at layer `k`.
    pub fn cell_center(&self, k: usize) -> Result<f64> {
        self.check_layer(k)?;
        Ok(self.centers[k])
    }

    /// Scalar layer thickness at layer `k`.
    pub fn cell_width(&self, k: usize) -> Result<f64> {
        self.check_layer(k)?;
        Ok(self.widths[k])
    }

    /// Axis-aligned vertical neighbours of layer `k`. Bottom layer has no
    /// `down` side; top layer has no `up` side.
    pub fn neighbors(&self, k: usize) -> Result<VerticalNeighbors> {
        self.check_layer(k)?;
        let nz = self.nz();
        let down = if k > 0 { Some(k - 1) } else { None };
        let up = if k + 1 < nz { Some(k + 1) } else { None };
        Ok(VerticalNeighbors { down, up })
    }

    /// Evaluate a scalar metric field at layer `k`.
    pub fn metric_eval(&self, name: VerticalMetricName, k: usize) -> Result<f64> {
        self.check_layer(k)?;
        match name {
            VerticalMetricName::Dz => Ok(self.widths[k]),
            VerticalMetricName::Z => Ok(self.centers[k]),
            VerticalMetricName::Sigma => match self.coordinate {
                VerticalCoordinate::Sigma
                | VerticalCoordinate::HybridSigmaTheta
                | VerticalCoordinate::Eta => Ok(self.centers[k]),
                other => Err(GridError::InvalidOption(
                    "metric_name",
                    format!("'sigma' undefined for coordinate '{}'", other.name()),
                )),
            },
            VerticalMetricName::Pressure => {
                if self.ak.is_empty() || self.bk.is_empty() {
                    return Err(GridError::InvalidOption(
                        "metric_name",
                        format!(
                            "'pressure' requires hybrid ak/bk (coordinate '{}' has none)",
                            self.coordinate.name()
                        ),
                    ));
                }
                let p_lo = self.ak[k] + self.bk[k] * self.p0;
                let p_hi = self.ak[k + 1] + self.bk[k + 1] * self.p0;
                Ok(0.5 * (p_lo + p_hi))
            }
            VerticalMetricName::Ak => {
                if self.ak.is_empty() {
                    return Err(GridError::InvalidOption(
                        "metric_name",
                        "'ak' unavailable (no hybrid coefficients)".to_string(),
                    ));
                }
                Ok(0.5 * (self.ak[k] + self.ak[k + 1]))
            }
            VerticalMetricName::Bk => {
                if self.bk.is_empty() {
                    return Err(GridError::InvalidOption(
                        "metric_name",
                        "'bk' unavailable (no hybrid coefficients)".to_string(),
                    ));
                }
                Ok(0.5 * (self.bk[k] + self.bk[k + 1]))
            }
        }
    }

    /// Evaluate a scalar metric by its canonical wire-form name.
    pub fn metric_eval_by_name(&self, name: &str, k: usize) -> Result<f64> {
        let m = VerticalMetricName::from_name(name).ok_or_else(|| {
            GridError::InvalidOption("metric_name", format!("unknown metric: {name:?}"))
        })?;
        self.metric_eval(m, k)
    }

    /// Provenance block per `docs/GRIDS_API.md` §6.4. Matches the Python /
    /// Julia / TypeScript sibling bindings' key set so cross-binding byte
    /// equality is achievable by the conformance harness after stripping
    /// `binding` and `binding_version`.
    pub fn provenance(&self) -> Value {
        json!({
            "binding": "rust",
            "binding_version": env!("CARGO_PKG_VERSION"),
            "family": "vertical",
            "version": API_VERSION,
            "coordinate": self.coordinate.name(),
            "dtype": dtype_str(self.dtype),
        })
    }

    // internals -----------------------------------------------------------

    fn check_layer(&self, k: usize) -> Result<()> {
        let nz = self.nz();
        if k >= nz {
            return Err(GridError::InvalidOption(
                "k",
                format!("layer index {k} out of range [0, {nz})"),
            ));
        }
        Ok(())
    }
}

impl Grid for VerticalGrid {
    fn family(&self) -> &'static str {
        "vertical"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// Declarative `.esm` lowering per the 2026-04-20 scope correction.
    ///
    /// Emits `family` + `topology` + `dtype` + `ndim` + `ghosts` + cell /
    /// vertex / edge counts + an `options` bag (coordinate kind, `nz`, full
    /// `levels` array, and hybrid `ak`/`bk` + `p0` when applicable) +
    /// provenance + schema version. Derived arrays (centres, widths) do not
    /// appear in the wire form.
    fn to_esm(&self) -> Value {
        let mut options = serde_json::Map::new();
        options.insert(
            "coordinate".to_string(),
            Value::String(self.coordinate.name().to_string()),
        );
        options.insert("nz".to_string(), Value::from(self.nz()));
        options.insert("levels".to_string(), json!(self.levels));
        if !self.ak.is_empty() {
            options.insert("ak".to_string(), json!(self.ak));
        }
        if !self.bk.is_empty() {
            options.insert("bk".to_string(), json!(self.bk));
        }
        if matches!(
            self.coordinate,
            VerticalCoordinate::Eta | VerticalCoordinate::HybridSigmaTheta
        ) || !self.ak.is_empty()
            || !self.bk.is_empty()
        {
            options.insert("p0".to_string(), json!(self.p0));
        }
        if let Some(t) = self.transition {
            options.insert("transition".to_string(), json!(t));
        }

        json!({
            "family": "vertical",
            "topology": self.topology(),
            "dtype": dtype_str(self.dtype),
            "ndim": 1,
            "ghosts": self.ghosts,
            "n_cells": self.n_cells(),
            "n_vertices": self.n_vertices(),
            "n_edges": self.n_edges(),
            "options": Value::Object(options),
            "provenance": self.provenance(),
            "schema_version": API_VERSION,
        })
    }
}

fn dtype_str(d: Dtype) -> &'static str {
    match d {
        Dtype::F64 => "float64",
        Dtype::F32 => "float32",
    }
}

// --- builder ---------------------------------------------------------------

/// Create a new [`VerticalBuilder`] with default options.
///
/// ```
/// use earthsci_grids::{vertical, Grid};
/// let g = vertical::builder()
///     .coordinate(vertical::VerticalCoordinate::Sigma)
///     .nz(16)
///     .build()
///     .unwrap();
/// assert_eq!(g.family(), "vertical");
/// assert_eq!(g.n_cells(), 16);
/// ```
pub fn builder() -> VerticalBuilder {
    VerticalBuilder::default()
}

/// Builder for [`VerticalGrid`]. Required options are validated at
/// [`VerticalBuilder::build`] time so partial builders can be reused in tests.
#[must_use]
#[derive(Clone, Debug, Default)]
pub struct VerticalBuilder {
    coordinate: Option<VerticalCoordinate>,
    nz: Option<usize>,
    levels: Option<Vec<f64>>,
    ak: Option<Vec<f64>>,
    bk: Option<Vec<f64>>,
    p0: Option<f64>,
    transition: Option<f64>,
    dtype: Option<Dtype>,
    ghosts: Option<usize>,
}

impl VerticalBuilder {
    /// Vertical coordinate kind. Required.
    pub fn coordinate(mut self, coordinate: VerticalCoordinate) -> Self {
        self.coordinate = Some(coordinate);
        self
    }

    /// Number of layers. Required for uniform `Sigma` / `HybridSigmaTheta`
    /// when `levels` is omitted. Must be `>= 1`.
    pub fn nz(mut self, nz: usize) -> Self {
        self.nz = Some(nz);
        self
    }

    /// Explicit interface values (length `nz + 1`).
    ///
    /// Required for `Z`, `Theta`, `ZStar`. Optional for `Sigma` /
    /// `HybridSigmaTheta` (overrides `nz` if both given). Forbidden for `Eta`
    /// (use `ak`/`bk` instead).
    pub fn levels(mut self, levels: Vec<f64>) -> Self {
        self.levels = Some(levels);
        self
    }

    /// Hybrid A coefficients at interfaces (length `nz + 1`). Required for
    /// `Eta`; optional for `HybridSigmaTheta`.
    pub fn ak(mut self, ak: Vec<f64>) -> Self {
        self.ak = Some(ak);
        self
    }

    /// Hybrid B coefficients at interfaces (length `nz + 1`). Required for
    /// `Eta`; optional for `HybridSigmaTheta`.
    pub fn bk(mut self, bk: Vec<f64>) -> Self {
        self.bk = Some(bk);
        self
    }

    /// Reference surface pressure (Pa). Defaults to `1.0e5`. Must be positive
    /// and finite.
    pub fn p0(mut self, p0: f64) -> Self {
        self.p0 = Some(p0);
        self
    }

    /// `hybrid_sigma_theta`-only: sigma value in `(0, 1)` at which the
    /// coordinate transitions from pure-sigma to pure-theta.
    pub fn transition(mut self, transition: f64) -> Self {
        self.transition = Some(transition);
        self
    }

    /// Element precision. Defaults to [`Dtype::F64`].
    pub fn dtype(mut self, dtype: Dtype) -> Self {
        self.dtype = Some(dtype);
        self
    }

    /// Halo layer width. Defaults to 0.
    pub fn ghosts(mut self, ghosts: usize) -> Self {
        self.ghosts = Some(ghosts);
        self
    }

    /// Validate options and construct the grid.
    pub fn build(self) -> Result<VerticalGrid> {
        let coordinate = self
            .coordinate
            .ok_or(GridError::MissingOption("coordinate"))?;
        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);
        let p0 = self.p0.unwrap_or(1.0e5);
        if !(p0.is_finite() && p0 > 0.0) {
            return Err(GridError::InvalidOption(
                "p0",
                format!("must be positive and finite, got {p0}"),
            ));
        }
        if let Some(nz) = self.nz {
            if nz < 1 {
                return Err(GridError::InvalidOption(
                    "nz",
                    format!("must be >= 1, got {nz}"),
                ));
            }
        }

        let levels: Vec<f64>;
        let ak: Vec<f64>;
        let bk: Vec<f64>;
        let mut transition: Option<f64> = None;

        match coordinate {
            VerticalCoordinate::Sigma => {
                if let Some(lv) = self.levels {
                    let lv = coerce_levels(lv, "levels", MonoDir::Decreasing, Some((0.0, 1.0)))?;
                    check_nz_consistency(self.nz, lv.len(), "sigma")?;
                    levels = lv;
                } else {
                    let nz = self
                        .nz
                        .ok_or_else(|| missing("sigma requires `nz` or `levels`"))?;
                    levels = uniform_sigma_levels(nz);
                }
                ak = Vec::new();
                bk = Vec::new();
            }
            VerticalCoordinate::Z | VerticalCoordinate::Theta | VerticalCoordinate::ZStar => {
                let lv_opt = self.levels.ok_or_else(|| {
                    missing(match coordinate {
                        VerticalCoordinate::Z => "'z' requires explicit `levels`",
                        VerticalCoordinate::Theta => "'theta' requires explicit `levels`",
                        _ => "'z_star' requires explicit `levels`",
                    })
                })?;
                let lv = coerce_levels(lv_opt, "levels", MonoDir::Increasing, None)?;
                check_nz_consistency(self.nz, lv.len(), coordinate.name())?;
                levels = lv;
                ak = Vec::new();
                bk = Vec::new();
            }
            VerticalCoordinate::Eta => {
                let ak_in = self
                    .ak
                    .ok_or_else(|| missing("'eta' requires `ak` (length nz+1)"))?;
                let bk_in = self
                    .bk
                    .ok_or_else(|| missing("'eta' requires `bk` (length nz+1)"))?;
                if ak_in.len() != bk_in.len() {
                    return Err(GridError::InvalidOption(
                        "eta",
                        format!(
                            "ak/bk must have equal length; got {} vs {}",
                            ak_in.len(),
                            bk_in.len()
                        ),
                    ));
                }
                let nz_eff = self.nz.unwrap_or(ak_in.len().saturating_sub(1));
                if nz_eff < 1 {
                    return Err(GridError::InvalidOption(
                        "nz",
                        format!("must be >= 1, got {nz_eff}"),
                    ));
                }
                let ak_arr = coerce_hybrid(ak_in, nz_eff + 1, "ak")?;
                let bk_arr = coerce_hybrid(bk_in, nz_eff + 1, "bk")?;
                let mut sigma = Vec::with_capacity(ak_arr.len());
                for i in 0..ak_arr.len() {
                    sigma.push(ak_arr[i] / p0 + bk_arr[i]);
                }
                for i in 0..(sigma.len() - 1) {
                    if sigma[i + 1] >= sigma[i] {
                        return Err(GridError::InvalidOption(
                            "eta",
                            "synthesized sigma (ak/p0 + bk) must be strictly decreasing"
                                .to_string(),
                        ));
                    }
                }
                if self.levels.is_some() {
                    return Err(GridError::InvalidOption(
                        "levels",
                        "'eta' derives `levels` from ak/bk; pass only `ak` and `bk`".to_string(),
                    ));
                }
                levels = sigma;
                ak = ak_arr;
                bk = bk_arr;
            }
            VerticalCoordinate::HybridSigmaTheta => {
                let lv = if let Some(lv) = self.levels {
                    let lv = coerce_levels(lv, "levels", MonoDir::Decreasing, Some((0.0, 1.0)))?;
                    check_nz_consistency(self.nz, lv.len(), "hybrid_sigma_theta")?;
                    lv
                } else {
                    let nz = self
                        .nz
                        .ok_or_else(|| missing("'hybrid_sigma_theta' requires `nz` or `levels`"))?;
                    uniform_sigma_levels(nz)
                };
                if let Some(t) = self.transition {
                    if !(t.is_finite() && t > 0.0 && t < 1.0) {
                        return Err(GridError::InvalidOption(
                            "transition",
                            format!("hybrid_sigma_theta `transition` must be in (0, 1); got {t}"),
                        ));
                    }
                    transition = Some(t);
                }
                ak = match self.ak {
                    Some(a) => coerce_hybrid(a, lv.len(), "ak")?,
                    None => Vec::new(),
                };
                bk = match self.bk {
                    Some(b) => coerce_hybrid(b, lv.len(), "bk")?,
                    None => Vec::new(),
                };
                levels = lv;
            }
        }

        let (centers, widths) = centers_and_widths(&levels);

        Ok(VerticalGrid {
            coordinate,
            levels,
            centers,
            widths,
            ak,
            bk,
            p0,
            transition,
            ghosts,
            dtype,
        })
    }
}

// --- helpers ---------------------------------------------------------------

#[derive(Copy, Clone)]
enum MonoDir {
    Increasing,
    Decreasing,
}

fn missing(msg: &str) -> GridError {
    GridError::InvalidOption("coordinate", msg.to_string())
}

fn coerce_levels(
    values: Vec<f64>,
    label: &'static str,
    dir: MonoDir,
    domain: Option<(f64, f64)>,
) -> Result<Vec<f64>> {
    if values.len() < 2 {
        return Err(GridError::InvalidOption(
            label,
            format!("must have >= 2 entries; got {}", values.len()),
        ));
    }
    for (i, &v) in values.iter().enumerate() {
        if !v.is_finite() {
            return Err(GridError::InvalidOption(
                label,
                format!("entry {i} is not finite ({v})"),
            ));
        }
    }
    if let Some((lo, hi)) = domain {
        for (i, &v) in values.iter().enumerate() {
            if v < lo || v > hi {
                return Err(GridError::InvalidOption(
                    label,
                    format!("entries must lie in [{lo}, {hi}]; index {i} = {v}"),
                ));
            }
        }
    }
    for i in 0..(values.len() - 1) {
        let d = values[i + 1] - values[i];
        let ok = match dir {
            MonoDir::Increasing => d > 0.0,
            MonoDir::Decreasing => d < 0.0,
        };
        if !ok {
            let want = match dir {
                MonoDir::Increasing => "increasing",
                MonoDir::Decreasing => "decreasing",
            };
            return Err(GridError::InvalidOption(
                label,
                format!("must be strictly {want} (violation at index {i})"),
            ));
        }
    }
    Ok(values)
}

fn coerce_hybrid(values: Vec<f64>, expected_len: usize, label: &'static str) -> Result<Vec<f64>> {
    if values.len() != expected_len {
        return Err(GridError::InvalidOption(
            label,
            format!(
                "must have length nz+1 = {expected_len}; got {}",
                values.len()
            ),
        ));
    }
    for (i, &v) in values.iter().enumerate() {
        if !v.is_finite() {
            return Err(GridError::InvalidOption(
                label,
                format!("entry {i} is not finite ({v})"),
            ));
        }
    }
    Ok(values)
}

fn check_nz_consistency(nz: Option<usize>, levels_len: usize, coord: &str) -> Result<()> {
    if let Some(nz) = nz {
        if levels_len != nz + 1 {
            return Err(GridError::InvalidOption(
                "nz",
                format!("{coord}: nz={nz} inconsistent with levels length {levels_len}"),
            ));
        }
    }
    Ok(())
}

fn uniform_sigma_levels(nz: usize) -> Vec<f64> {
    (0..=nz).map(|k| 1.0 - (k as f64) / (nz as f64)).collect()
}

fn centers_and_widths(levels: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = levels.len() - 1;
    let mut centers = Vec::with_capacity(n);
    let mut widths = Vec::with_capacity(n);
    for k in 0..n {
        centers.push(0.5 * (levels[k] + levels[k + 1]));
        widths.push((levels[k + 1] - levels[k]).abs());
    }
    (centers, widths)
}

// --- tests -----------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn sigma_n16() -> VerticalGrid {
        builder()
            .coordinate(VerticalCoordinate::Sigma)
            .nz(16)
            .build()
            .unwrap()
    }

    fn z_l4() -> VerticalGrid {
        builder()
            .coordinate(VerticalCoordinate::Z)
            .levels(vec![0.0, 100.0, 300.0, 1000.0, 2000.0])
            .build()
            .unwrap()
    }

    fn eta_l4() -> VerticalGrid {
        // ak + bk*p0 produces strictly decreasing sigma = ak/p0 + bk.
        builder()
            .coordinate(VerticalCoordinate::Eta)
            .ak(vec![0.0, 2000.0, 5000.0, 10000.0, 0.0])
            .bk(vec![1.0, 0.7, 0.3, 0.05, 0.0])
            .build()
            .unwrap()
    }

    // --- builder validation ---

    #[test]
    fn builder_requires_coordinate() {
        let err = builder().nz(4).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("coordinate")));
    }

    #[test]
    fn sigma_requires_nz_or_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("coordinate", _)));
    }

    #[test]
    fn z_requires_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Z)
            .nz(4)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("coordinate", _)));
    }

    #[test]
    fn theta_requires_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Theta)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("coordinate", _)));
    }

    #[test]
    fn z_star_requires_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::ZStar)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("coordinate", _)));
    }

    #[test]
    fn eta_requires_ak_and_bk() {
        let err = builder()
            .coordinate(VerticalCoordinate::Eta)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("coordinate", _)));
        let err = builder()
            .coordinate(VerticalCoordinate::Eta)
            .ak(vec![0.0, 1.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("coordinate", _)));
    }

    #[test]
    fn eta_rejects_unequal_ak_bk_lengths() {
        let err = builder()
            .coordinate(VerticalCoordinate::Eta)
            .ak(vec![0.0, 1.0, 2.0])
            .bk(vec![1.0, 0.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("eta", _)));
    }

    #[test]
    fn eta_rejects_nondecreasing_synth_sigma() {
        // Constant sigma = ak/p0 + bk violates strictness.
        let err = builder()
            .coordinate(VerticalCoordinate::Eta)
            .ak(vec![0.0, 0.0, 0.0])
            .bk(vec![1.0, 0.5, 0.5])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("eta", _)));
    }

    #[test]
    fn eta_rejects_explicit_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Eta)
            .ak(vec![0.0, 1000.0, 0.0])
            .bk(vec![1.0, 0.5, 0.0])
            .levels(vec![1.0, 0.5, 0.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("levels", _)));
    }

    #[test]
    fn sigma_rejects_out_of_domain_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .levels(vec![1.2, 0.5, 0.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("levels", _)));
    }

    #[test]
    fn sigma_rejects_nondecreasing_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .levels(vec![1.0, 0.5, 0.5, 0.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("levels", _)));
    }

    #[test]
    fn z_rejects_nonincreasing_levels() {
        let err = builder()
            .coordinate(VerticalCoordinate::Z)
            .levels(vec![0.0, 100.0, 50.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("levels", _)));
    }

    #[test]
    fn levels_must_have_two_entries() {
        let err = builder()
            .coordinate(VerticalCoordinate::Z)
            .levels(vec![0.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("levels", _)));
    }

    #[test]
    fn levels_must_be_finite() {
        let err = builder()
            .coordinate(VerticalCoordinate::Z)
            .levels(vec![0.0, f64::NAN, 100.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("levels", _)));
    }

    #[test]
    fn nz_inconsistent_with_levels_is_rejected() {
        let err = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .nz(2)
            .levels(vec![1.0, 0.5, 0.25, 0.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nz", _)));
    }

    #[test]
    fn p0_must_be_positive() {
        let err = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .nz(4)
            .p0(-1.0)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("p0", _)));
    }

    #[test]
    fn transition_must_be_in_unit_open_interval() {
        let err = builder()
            .coordinate(VerticalCoordinate::HybridSigmaTheta)
            .nz(4)
            .transition(1.0)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("transition", _)));
    }

    #[test]
    fn nz_must_be_positive() {
        let err = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .nz(0)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nz", _)));
    }

    // --- uniform sigma ---

    #[test]
    fn uniform_sigma_is_surface_to_top() {
        let g = sigma_n16();
        assert_eq!(g.n_cells(), 16);
        assert_eq!(g.n_vertices(), 17);
        assert_eq!(g.n_edges(), 16);
        assert_eq!(g.ndim(), 1);
        assert_eq!(g.topology(), "column");
        // levels[0] = 1 (surface), levels[nz] = 0 (top).
        assert_relative_eq!(g.levels()[0], 1.0, epsilon = 1e-15);
        assert_relative_eq!(g.levels()[16], 0.0, epsilon = 1e-15);
        // centers midway between interfaces.
        assert_relative_eq!(g.centers()[0], 1.0 - 0.5 / 16.0, epsilon = 1e-15);
        // widths are uniform.
        for w in g.widths() {
            assert_relative_eq!(*w, 1.0 / 16.0, epsilon = 1e-15);
        }
    }

    // --- z / theta / z_star ---

    #[test]
    fn z_levels_stored_verbatim() {
        let g = z_l4();
        assert_eq!(g.levels(), &[0.0, 100.0, 300.0, 1000.0, 2000.0]);
        assert_relative_eq!(g.centers()[0], 50.0, epsilon = 1e-15);
        assert_relative_eq!(g.widths()[0], 100.0, epsilon = 1e-15);
        assert_relative_eq!(g.widths()[3], 1000.0, epsilon = 1e-15);
    }

    #[test]
    fn z_star_accepts_increasing_levels() {
        let g = builder()
            .coordinate(VerticalCoordinate::ZStar)
            .levels(vec![0.0, 10.0, 25.0])
            .build()
            .unwrap();
        assert_eq!(g.coordinate(), VerticalCoordinate::ZStar);
        assert_eq!(g.n_cells(), 2);
    }

    #[test]
    fn theta_accepts_increasing_levels() {
        let g = builder()
            .coordinate(VerticalCoordinate::Theta)
            .levels(vec![290.0, 310.0, 400.0])
            .build()
            .unwrap();
        assert_relative_eq!(g.widths()[0], 20.0, epsilon = 1e-15);
        assert_relative_eq!(g.widths()[1], 90.0, epsilon = 1e-15);
    }

    // --- eta ---

    #[test]
    fn eta_synthesizes_sigma_from_ak_bk() {
        let g = eta_l4();
        assert_eq!(g.coordinate(), VerticalCoordinate::Eta);
        assert_eq!(g.levels().len(), 5);
        // First interface: ak=0, bk=1 -> sigma=1. Last: ak=0, bk=0 -> sigma=0.
        assert_relative_eq!(g.levels()[0], 1.0, epsilon = 1e-15);
        assert_relative_eq!(g.levels()[4], 0.0, epsilon = 1e-15);
        // Synthesized sigma must be strictly decreasing.
        for i in 0..(g.levels().len() - 1) {
            assert!(g.levels()[i + 1] < g.levels()[i]);
        }
    }

    // --- hybrid_sigma_theta ---

    #[test]
    fn hybrid_sigma_theta_defaults_to_uniform_sigma() {
        let g = builder()
            .coordinate(VerticalCoordinate::HybridSigmaTheta)
            .nz(8)
            .transition(0.3)
            .build()
            .unwrap();
        assert_eq!(g.coordinate(), VerticalCoordinate::HybridSigmaTheta);
        assert_eq!(g.n_cells(), 8);
        assert_eq!(g.transition(), Some(0.3));
    }

    // --- accessors ---

    #[test]
    fn neighbors_interior_has_both_sides() {
        let g = sigma_n16();
        let n = g.neighbors(5).unwrap();
        assert_eq!(n.down, Some(4));
        assert_eq!(n.up, Some(6));
    }

    #[test]
    fn neighbors_bottom_has_no_down() {
        let g = sigma_n16();
        let n = g.neighbors(0).unwrap();
        assert_eq!(n.down, None);
        assert_eq!(n.up, Some(1));
    }

    #[test]
    fn neighbors_top_has_no_up() {
        let g = sigma_n16();
        let n = g.neighbors(15).unwrap();
        assert_eq!(n.down, Some(14));
        assert_eq!(n.up, None);
    }

    #[test]
    fn neighbors_rejects_out_of_range() {
        let g = sigma_n16();
        assert!(g.neighbors(16).is_err());
    }

    #[test]
    fn cell_center_and_width_match_slices() {
        let g = z_l4();
        for k in 0..g.n_cells() {
            assert_relative_eq!(g.cell_center(k).unwrap(), g.centers()[k], epsilon = 1e-15);
            assert_relative_eq!(g.cell_width(k).unwrap(), g.widths()[k], epsilon = 1e-15);
        }
        assert!(g.cell_center(g.n_cells()).is_err());
    }

    // --- metrics ---

    #[test]
    fn metric_dz_equals_width() {
        let g = z_l4();
        let dz = g.metric_eval(VerticalMetricName::Dz, 2).unwrap();
        assert_relative_eq!(dz, 700.0, epsilon = 1e-15);
    }

    #[test]
    fn metric_z_equals_center() {
        let g = z_l4();
        let z = g.metric_eval(VerticalMetricName::Z, 0).unwrap();
        assert_relative_eq!(z, 50.0, epsilon = 1e-15);
    }

    #[test]
    fn metric_sigma_rejected_for_z() {
        let g = z_l4();
        let err = g.metric_eval(VerticalMetricName::Sigma, 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn metric_sigma_ok_for_sigma_like() {
        let g = sigma_n16();
        let s = g.metric_eval(VerticalMetricName::Sigma, 0).unwrap();
        assert_relative_eq!(s, g.centers()[0], epsilon = 1e-15);
        let g = eta_l4();
        let s = g.metric_eval(VerticalMetricName::Sigma, 0).unwrap();
        assert_relative_eq!(s, g.centers()[0], epsilon = 1e-15);
    }

    #[test]
    fn metric_pressure_requires_hybrid() {
        let g = sigma_n16();
        let err = g.metric_eval(VerticalMetricName::Pressure, 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn metric_pressure_averages_interfaces() {
        let g = eta_l4();
        let p = g.metric_eval(VerticalMetricName::Pressure, 0).unwrap();
        // Interface 0: ak=0, bk=1 -> p_lo = 0 + 1*p0.
        // Interface 1: ak=2000, bk=0.7 -> p_hi = 2000 + 0.7*p0.
        let p0 = g.p0();
        let expected = 0.5 * ((0.0 + 1.0 * p0) + (2000.0 + 0.7 * p0));
        assert_relative_eq!(p, expected, epsilon = 1e-12);
    }

    #[test]
    fn metric_ak_bk_require_hybrid() {
        let g = sigma_n16();
        assert!(g.metric_eval(VerticalMetricName::Ak, 0).is_err());
        assert!(g.metric_eval(VerticalMetricName::Bk, 0).is_err());
    }

    #[test]
    fn metric_ak_bk_average_interfaces() {
        let g = eta_l4();
        let ak = g.metric_eval(VerticalMetricName::Ak, 1).unwrap();
        let bk = g.metric_eval(VerticalMetricName::Bk, 1).unwrap();
        assert_relative_eq!(ak, 0.5 * (2000.0 + 5000.0), epsilon = 1e-12);
        assert_relative_eq!(bk, 0.5 * (0.7 + 0.3), epsilon = 1e-12);
    }

    #[test]
    fn metric_eval_by_name_parses_wire_form() {
        let g = z_l4();
        let dz = g.metric_eval_by_name("dz", 0).unwrap();
        assert_relative_eq!(dz, 100.0, epsilon = 1e-15);
    }

    #[test]
    fn metric_eval_by_name_rejects_unknown() {
        let g = z_l4();
        let err = g.metric_eval_by_name("foo", 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn metric_eval_rejects_out_of_range() {
        let g = z_l4();
        assert!(g.metric_eval(VerticalMetricName::Dz, g.n_cells()).is_err());
    }

    // --- to_esm ---

    #[test]
    fn to_esm_sigma_matches_golden_shape() {
        let g = sigma_n16();
        let doc = g.to_esm();
        assert_eq!(doc["family"], "vertical");
        assert_eq!(doc["topology"], "column");
        assert_eq!(doc["dtype"], "float64");
        assert_eq!(doc["ndim"], 1);
        assert_eq!(doc["ghosts"], 0);
        assert_eq!(doc["n_cells"], 16);
        assert_eq!(doc["n_vertices"], 17);
        assert_eq!(doc["n_edges"], 16);
        assert_eq!(doc["schema_version"], "1.0.0");
        let opts = &doc["options"];
        assert_eq!(opts["coordinate"], "sigma");
        assert_eq!(opts["nz"], 16);
        assert!(opts["levels"].is_array());
        assert_eq!(opts["levels"].as_array().unwrap().len(), 17);
        // No hybrid coefficients + pure sigma => no p0 key.
        assert!(opts.get("p0").is_none());
        assert!(opts.get("ak").is_none());
        assert!(opts.get("bk").is_none());

        // Byte-level checks mirror the committed golden fixture
        // `discretizations/grids/vertical/sigma_uniform_n16.esm`. The
        // conformance harness canonicalizes sort-key JSON + strips the
        // `binding*` provenance fields, so semantic equality here suffices.
        let levels = opts["levels"].as_array().unwrap();
        assert_eq!(levels[0], 1.0);
        assert_eq!(levels[16], 0.0);
        assert_eq!(levels[8], 0.5);

        // Derived arrays must not leak into the wire form.
        let wire = serde_json::to_string(&doc).unwrap();
        assert!(!wire.contains("\"centers\""));
        assert!(!wire.contains("\"widths\""));
        assert!(wire.len() < 2_000, "wire too large: {} bytes", wire.len());
    }

    #[test]
    fn to_esm_eta_includes_ak_bk_p0() {
        let g = eta_l4();
        let doc = g.to_esm();
        let opts = &doc["options"];
        assert_eq!(opts["coordinate"], "eta");
        assert!(opts["ak"].is_array());
        assert!(opts["bk"].is_array());
        assert!(opts["p0"].is_number());
    }

    #[test]
    fn to_esm_z_omits_hybrid_and_p0() {
        let g = z_l4();
        let doc = g.to_esm();
        let opts = &doc["options"];
        assert_eq!(opts["coordinate"], "z");
        assert_eq!(opts["nz"], 4);
        assert!(opts.get("ak").is_none());
        assert!(opts.get("bk").is_none());
        assert!(opts.get("p0").is_none());
    }

    #[test]
    fn to_esm_hybrid_sigma_theta_carries_transition_when_set() {
        let g = builder()
            .coordinate(VerticalCoordinate::HybridSigmaTheta)
            .nz(4)
            .transition(0.25)
            .build()
            .unwrap();
        let doc = g.to_esm();
        let opts = &doc["options"];
        assert_eq!(opts["transition"], 0.25);
        // p0 is still carried because hybrid_sigma_theta is sigma-pressure hybrid.
        assert!(opts["p0"].is_number());
    }

    #[test]
    fn to_esm_roundtrips_via_json() {
        let g = sigma_n16();
        let doc = g.to_esm();
        let text = serde_json::to_string(&doc).unwrap();
        let reparsed: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(reparsed["family"], "vertical");
        assert_eq!(reparsed["options"]["nz"], 16);
    }

    #[test]
    fn provenance_identifies_rust_binding() {
        let g = sigma_n16();
        let doc = g.to_esm();
        let prov = &doc["provenance"];
        assert_eq!(prov["binding"], "rust");
        assert_eq!(prov["binding_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(prov["family"], "vertical");
        assert_eq!(prov["version"], "1.0.0");
        assert_eq!(prov["coordinate"], "sigma");
    }

    #[test]
    fn dtype_f32_propagates_to_esm() {
        let g = builder()
            .coordinate(VerticalCoordinate::Sigma)
            .nz(4)
            .dtype(Dtype::F32)
            .build()
            .unwrap();
        assert_eq!(g.dtype(), Dtype::F32);
        assert_eq!(g.to_esm()["dtype"], "float32");
    }

    // --- coordinate name round-trip ---

    #[test]
    fn coordinate_name_round_trips() {
        for c in [
            VerticalCoordinate::Sigma,
            VerticalCoordinate::Eta,
            VerticalCoordinate::Z,
            VerticalCoordinate::Theta,
            VerticalCoordinate::HybridSigmaTheta,
            VerticalCoordinate::ZStar,
        ] {
            assert_eq!(VerticalCoordinate::from_name(c.name()), Some(c));
        }
        assert_eq!(VerticalCoordinate::from_name("bogus"), None);
    }

    #[test]
    fn metric_name_round_trips() {
        for m in [
            VerticalMetricName::Dz,
            VerticalMetricName::Z,
            VerticalMetricName::Sigma,
            VerticalMetricName::Pressure,
            VerticalMetricName::Ak,
            VerticalMetricName::Bk,
        ] {
            let name = match m {
                VerticalMetricName::Dz => "dz",
                VerticalMetricName::Z => "z",
                VerticalMetricName::Sigma => "sigma",
                VerticalMetricName::Pressure => "pressure",
                VerticalMetricName::Ak => "ak",
                VerticalMetricName::Bk => "bk",
            };
            assert_eq!(VerticalMetricName::from_name(name), Some(m));
        }
        assert_eq!(VerticalMetricName::from_name("bogus"), None);
    }
}
