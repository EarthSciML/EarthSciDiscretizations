---
title: "Reprojection"
description: "Coordinate-transform template fragments (longitude-latitude, Lambert conformal, …) for use inside regridding and coupling expressions."
---

Coordinate-transform template fragments (longitude-latitude, Lambert conformal, …) for use inside regridding and coupling expressions.

These are forward/inverse **scalar** templates over the evaluable core — but they also apply **in-model over coordinate arrays**: invoke them inside a `faq` node whose bindings mix indexed reads, literals, and parameter references, and every `apply_expression_template` inlines (esm-spec §9.6.2 Option A) into a closed scalar AST at each cell. The `tests/conformance/reprojection/lcc_grid_roundtrip` case pins exactly this faq-mapped lowering: byte-identical expanded AST across all five bindings, and a working round-trip simulation on Julia + Python + Rust.

### `lambert_conformal`

<p class="rule-meta"><code>crs:lambert_conformal</code></p>

Source: [`reprojection/lambert_conformal.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/reprojection/lambert_conformal.esm)

Spherical Lambert Conformal Conic (PROJ +proj=lcc +R=...), as a template library (esm-spec 9.7) following the reprojection signature convention established by reprojection/longlat.esm: four public scalar templates lambert_conformal_forward_x / _forward_y (inputs lon, lat in DEGREES plus the projection constants; outputs easting/northing x, y in METRES) and lambert_conformal_inverse_lon / _inverse_lat (inputs x, y in metres plus the same constants; outputs degrees), plus five stem-prefixed helper templates (_t, _n, _F, _rho, _theta) that factor Snyder's intermediate quantities -- helper invocations inline at registration time by pure substitution (9.7.3), so every public template lowers to a closed scalar AST over the evaluable core (tan, atan, atan2, log, ^, sqrt, sign, sin, cos; all esm-spec 4.2 elementary functions). The projection constants {lat_1, lat_2 (standard parallels, deg), lat_0 (origin latitude, deg), lon_0 (central meridian, deg), R (sphere radius, m)} are TEMPLATE PARAMETERS, not metaparameters: metaparameters are load-time integers (9.7.6) while these are real-valued, and template parameters bind arbitrary expressions at each apply site -- numeric literals for a fixed CRS, or parameter references so one document drives several parameter sets (the WRF CONUS and NEI2016 sets of the conformance manifest flow through this ONE library with different bindings). THE MATH (Snyder 1987, PP 1395, pp. 104-110, sphere; the same closed form as archive/reprojection/lambert_conformal.esm, re-expressed in the 0.8.0 library form): with phi = lat*pi/180 etc., cone constant n = ln(cos phi_1 / cos phi_2) / ln(tan(pi/4 + phi_2/2) / tan(pi/4 + phi_1/2)) (eq. 15-3), F = cos phi_1 * tan(pi/4 + phi_1/2)^n / n (15-2), rho = R*F / tan(pi/4 + phi/2)^n (15-1, and rho_0 = rho at phi_0), theta = n*(lambda - lambda_0) (14-4); FORWARD x = rho*sin(theta) (14-1), y = rho_0 - rho*cos(theta) (14-2); INVERSE rho = sign(n)*sqrt(x^2 + (rho_0 - y)^2) (14-10), theta = atan2(x, rho_0 - y) (14-11, n > 0 branch: for n < 0 Snyder reverses the signs of x, y, and rho_0 - y first -- both documented parameter sets have n > 0, approx. 0.7156 for WRF 30/60 and 0.6305 for NEI2016 33/45), lat = (2*atan((R*F/rho)^(1/n)) - pi/2)*180/pi (15-5), lon = lon_0 + (theta/n)*180/pi (14-9). inverse(forward(...)) is the algebraic identity wherever the forward map is defined (singular only at the far pole phi = -90 deg for n > 0). Radian-degree conversions are inlined as the double-precision literals pi/180 = 0.017453292519943295, 180/pi = 57.29577951308232, pi/4 = 0.7853981633974483, pi/360 = 0.008726646259971648, pi/2 = 1.5707963267948966. Independent analytic reference points (both parameter sets, forward + inverse + round-trip) live in tests/conformance/reprojection/lambert_conformal/golden/points.json.

<details><summary>Templates (9)</summary>

**`lambert_conformal_t`**

Helper: the conformal tangent t(lat) = tan(pi/4 + phi/2) with phi = lat*pi/180 (so the half-angle is lat*pi/360). The building block of Snyder eqs. 15-1..15-3.

Parameters: `lat`

<div class="math-block">tan(0.7853981633974483 + lat * 8.726646259971648e-3)</div>

**`lambert_conformal_n`**

Helper: cone constant n = ln(cos phi_1 / cos phi_2) / ln(t(lat_2)/t(lat_1)) (Snyder eq. 15-3). Depends only on the standard parallels: approx. 0.715567 for WRF (30/60), 0.630478 for NEI2016 (33/45). Requires lat_1 != lat_2 (for equal parallels Snyder's limit is n = sin phi_1, eq. 15-4 -- not expressed here).

Parameters: `lat_1`, `lat_2`

<div class="math-block">log(cos(lat_1 * 0.017453292519943295) / cos(lat_2 * 0.017453292519943295)) / log(lambert_conformal_t&lt;lat=lat_2&gt; / lambert_conformal_t&lt;lat=lat_1&gt;)</div>

**`lambert_conformal_F`**

Helper: F = cos phi_1 * t(lat_1)^n / n (Snyder eq. 15-2), the constant tying the cone to the unit sphere.

Parameters: `lat_1`, `lat_2`

<div class="math-block">cos(lat_1 * 0.017453292519943295) * lambert_conformal_t&lt;lat=lat_1&gt;^lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt; / lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt;</div>

**`lambert_conformal_rho`**

Helper: rho = R*F / t(lat)^n (Snyder eq. 15-1), the cone radius of the parallel at lat (metres). rho_0 is this template invoked with lat bound to lat_0. Singular at the far pole (t -> 0 for n > 0 as lat -> -90).

Parameters: `lat`, `lat_1`, `lat_2`, `R`

<div class="math-block">R * lambert_conformal_F&lt;lat_1=lat_1, lat_2=lat_2&gt; / lambert_conformal_t&lt;lat=lat&gt;^lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt;</div>

**`lambert_conformal_theta`**

Helper: theta = n * (lon - lon_0) * pi/180 (Snyder eq. 14-4), the polar angle on the developed cone (radians).

Parameters: `lon`, `lon_0`, `lat_1`, `lat_2`

<div class="math-block">lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt; * (lon - lon_0) * 0.017453292519943295</div>

**`lambert_conformal_forward_x`**

FORWARD, easting (metres): x = rho(lat) * sin(theta) (Snyder eq. 14-1). At lon = lon_0, theta = 0 so x = 0 (the central meridian is the y-axis). lat_0 is unused (x is independent of the origin latitude; uniform reprojection signature).

Parameters: `lon`, `lat`, `lat_1`, `lat_2`, `lat_0`, `lon_0`, `R`

<div class="math-block">lambert_conformal_rho&lt;R=R, lat=lat, lat_1=lat_1, lat_2=lat_2&gt; * sin(lambert_conformal_theta&lt;lat_1=lat_1, lat_2=lat_2, lon=lon, lon_0=lon_0&gt;)</div>

**`lambert_conformal_forward_y`**

FORWARD, northing (metres): y = rho_0 - rho(lat) * cos(theta) (Snyder eq. 14-2), with rho_0 = rho at the origin latitude (the _rho helper invoked with lat bound to lat_0). At (lon_0, lat_0): theta = 0 and rho = rho_0, so y = 0 -- the projection origin maps to the false origin (0, 0).

Parameters: `lon`, `lat`, `lat_1`, `lat_2`, `lat_0`, `lon_0`, `R`

<div class="math-block">lambert_conformal_rho&lt;R=R, lat=lat_0, lat_1=lat_1, lat_2=lat_2&gt; - lambert_conformal_rho&lt;R=R, lat=lat, lat_1=lat_1, lat_2=lat_2&gt; * cos(lambert_conformal_theta&lt;lat_1=lat_1, lat_2=lat_2, lon=lon, lon_0=lon_0&gt;)</div>

**`lambert_conformal_inverse_lon`**

INVERSE, longitude (degrees): lon = lon_0 + (theta/n)*180/pi with theta = atan2(x, rho_0 - y) (Snyder eqs. 14-11, 14-9; n > 0 branch of 14-11 -- for n < 0 Snyder reverses the signs of x, y, rho_0 - y first, a case no documented parameter set exercises). atan2(rho*sin theta, rho*cos theta) = theta for rho > 0, so this inverts the forward polar angle exactly.

Parameters: `x`, `y`, `lat_1`, `lat_2`, `lat_0`, `lon_0`, `R`

<div class="math-block">lon_0 + atan2(x, lambert_conformal_rho&lt;R=R, lat=lat_0, lat_1=lat_1, lat_2=lat_2&gt; - y) / lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt; * 57.29577951308232</div>

**`lambert_conformal_inverse_lat`**

INVERSE, latitude (degrees): lat = (2*atan((R*F/rho)^(1/n)) - pi/2)*180/pi (Snyder eq. 15-5), with rho = sign(n)*sqrt(x^2 + (rho_0 - y)^2) (14-10, taking the sign of n) and rho_0 = rho(lat_0). Since R*F/rho = t(lat)^n on the forward image, this recovers the input latitude: the round-trip identity. lon_0 is unused (uniform reprojection signature).

Parameters: `x`, `y`, `lat_1`, `lat_2`, `lat_0`, `lon_0`, `R`

<div class="math-block">(2 * atan((R * lambert_conformal_F&lt;lat_1=lat_1, lat_2=lat_2&gt; / (sign(lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt;) * sqrt(x^2 + (lambert_conformal_rho&lt;R=R, lat=lat_0, lat_1=lat_1, lat_2=lat_2&gt; - y)^2)))^(1 / lambert_conformal_n&lt;lat_1=lat_1, lat_2=lat_2&gt;)) - 1.5707963267948966) * 57.29577951308232</div>

</details>

**References**

- Snyder, J. P. (1987). Map Projections -- A Working Manual. U.S. Geological Survey Professional Paper 1395, pp. 104-110: Lambert Conformal Conic, sphere. Forward eqs. 14-1, 14-2, 14-4, 15-1, 15-2, 15-3; inverse eqs. 14-9, 14-10, 14-11, 15-5; numerical check example pp. 295-296.

### `longlat`

<p class="rule-meta"><code>crs:longlat</code></p>

Source: [`reprojection/longlat.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/reprojection/longlat.esm)

Geographic lon-lat 'projection' (PROJ +proj=longlat; identity / plate carree in degrees), as a template library (esm-spec 9.7). This is the CONVENTION-ESTABLISHING degenerate base case for the reprojection/ family: every esd:reproject library declares four match-less scalar templates -- <crs>_forward_x and <crs>_forward_y taking the geographic point (lon, lat, in degrees) plus the projection constants and returning the projected coordinate (x resp. y), and <crs>_inverse_lon and <crs>_inverse_lat taking the projected point (x, y) plus the SAME constants and returning the geographic coordinate. A coordinate transform is a POINT-WISE SCALAR FORMULA, so each template body is a scalar expression; a consumer maps it over a coordinate array with an ordinary aggregate (bind lon/lat to indexed reads), or evaluates it at scalars. Projection constants are TEMPLATE PARAMETERS, not metaparameters: metaparameters are load-time integers (9.7.6) while lon_0 etc. are real-valued, and template parameters bind arbitrary expressions (numeric literals included) at each apply_expression_template site. Every template takes the full signature even where an argument is unused (longlat's y is independent of lon and lon_0) -- uniform signatures are what let a consumer swap CRS libraries without rewiring bindings; esm-spec 9.6.1 explicitly allows a parameter to occur zero times in the body. The math (archive/reprojection/longlat.esm, re-expressed in the 0.8.0 library form): FORWARD x = lon - lon_0, y = lat; INVERSE lon = x + lon_0, lat = y. inverse(forward(lon, lat)) = (lon - lon_0) + lon_0 = lon and lat = lat, so the round-trip is the identity by construction; with lon_0 = 0 the forward map degenerates to the pure geographic identity, and a non-zero lon_0 keeps the parameter load-bearing. Units convention: lon/lat and x/y are all DEGREES for this CRS (a projected CRS such as lambert_conformal returns metres -- the units are part of each library's documented contract, not of the shared signature).

<details><summary>Templates (4)</summary>

**`longlat_forward_x`**

FORWARD, x component (degrees): x = lon - lon_0. lat is unused (uniform reprojection signature).

Parameters: `lon`, `lat`, `lon_0`

<div class="math-block">lon - lon_0</div>

**`longlat_forward_y`**

FORWARD, y component (degrees): y = lat. lon and lon_0 are unused (uniform reprojection signature).

Parameters: `lon`, `lat`, `lon_0`

<div class="math-block">lat</div>

**`longlat_inverse_lon`**

INVERSE, longitude (degrees): lon = x + lon_0. Composes with longlat_forward_x to the identity: (lon - lon_0) + lon_0 = lon.

Parameters: `x`, `y`, `lon_0`

<div class="math-block">x + lon_0</div>

**`longlat_inverse_lat`**

INVERSE, latitude (degrees): lat = y. x and lon_0 are unused (uniform reprojection signature).

Parameters: `x`, `y`, `lon_0`

<div class="math-block">y</div>

</details>

**References**

- Snyder, J. P. (1987). Map Projections -- A Working Manual. U.S. Geological Survey Professional Paper 1395, pp. 90-91 (Equirectangular projection, eqs. 12-1..12-6 with phi_1 = 0: the degenerate plate carree case, here kept in degrees).

