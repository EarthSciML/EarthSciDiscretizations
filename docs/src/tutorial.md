# Tutorial: Advection on the Sphere

This tutorial demonstrates a simple advection simulation on the cubed-sphere grid using EarthSciDiscretizations.jl. We transport a cosine bell initial condition under solid-body rotation and visualize the results.

## Step 1: Create Grid and Initial Condition

We create a C16 cubed-sphere grid and project a cosine bell centered at (lon=0, lat=0) onto it.

```@example tutorial
using EarthSciDiscretizations
using EarthSciDiscretizations: evaluate_arrayop

Nc = 16
R = 1.0
grid = CubedSphereGrid(Nc; R=R)

# Cosine bell initial condition
function cosine_bell(lon, lat; lon0=0.0, lat0=0.0, r0=1.0)
    # Great-circle distance
    d = acos(clamp(sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(lon - lon0), -1.0, 1.0))
    return d < r0 ? 0.5 * (1.0 + cos(pi * d / r0)) : 0.0
end

q0 = project_initial_condition((lon, lat) -> cosine_bell(lon, lat; r0=1.0), grid)
println("Initial condition range: [$(minimum(q0)), $(maximum(q0))]")
println("Grid: $(Nc)x$(Nc) per panel, 6 panels")
```

## Step 2: Set Up Solid-Body Rotation

We define a solid-body rotation velocity field. The zonal wind is $u = \omega \cos(\text{lat})$ and the meridional wind is $v = 0$. We approximate the Courant number in computational coordinates.

```@example tutorial
# Angular velocity for one full rotation in T_rot time units
T_rot = 2pi  # one rotation period
omega = 2pi / T_rot

# Courant numbers in computational coordinates
# For solid-body rotation, the xi-courant is approximately omega * cos(lat) * dt / dxi
dt = 0.01  # time step
dxi = grid.dξ

courant_xi = zeros(6, Nc, Nc)
courant_eta = zeros(6, Nc, Nc)
for p in 1:6, i in 1:Nc, j in 1:Nc
    courant_xi[p, i, j] = omega * cos(grid.lat[p, i, j]) * dt / dxi
end

println("Max Courant number: $(maximum(abs.(courant_xi)))")
println("Time step: $(dt)")
println("Rotation period: $(T_rot)")
```

## Step 3: Forward Euler Time Integration

We use a simple forward Euler scheme to advance the tracer field. The `transport_2d` operator returns the tendency (time derivative), which we evaluate and apply.

Note: the transport operator works on the interior cells [6, Nc-2, Nc-2], so we update only the interior and keep boundary cells fixed.

```@example tutorial
q = copy(q0)
nsteps = 100

for step in 1:nsteps
    ao = transport_2d(q, courant_xi, courant_eta, grid)
    tendency = evaluate_arrayop(ao)

    # Update interior cells: index (p, i, j) in tendency maps to
    # physical cell (p, i+1, j+1)
    for p in 1:6, i in 1:size(tendency, 2), j in 1:size(tendency, 3)
        q[p, i + 1, j + 1] += dt * tendency[p, i, j]
    end
end

println("After $(nsteps) steps:")
println("  q range: [$(minimum(q)), $(maximum(q))]")
println("  Mass (sum q*area): $(sum(q .* grid.area))")
println("  Initial mass:      $(sum(q0 .* grid.area))")
```

## Step 4: Visualization

We plot the initial and final tracer fields on a Robinson projection using CairoMakie and GeoMakie.

```@example tutorial
using CairoMakie
using GeoMakie

fig = Figure(size=(900, 400))

for (col, (data, title)) in enumerate(zip([q0, q], ["Initial condition", "After $(nsteps) steps"]))
    ax = GeoAxis(fig[1, col];
        dest="+proj=robin",
        title=title,
    )

    for p in 1:6
        lons = rad2deg.(grid.lon[p, :, :])
        lats = rad2deg.(grid.lat[p, :, :])
        vals = data[p, :, :]
        surface!(ax, lons[:], lats[:], zeros(length(lons[:]));
            color=vals[:], colormap=:viridis,
            colorrange=(0.0, 1.0), shading=NoShading)
    end
end

fig
```

## Step 5: Animation

We create an animation of the advection process by running the time loop and recording frames.

```@example tutorial
q_anim = copy(q0)
nframes = 50
steps_per_frame = 2

fig_anim = Figure(size=(600, 350))
ax_anim = GeoAxis(fig_anim[1, 1];
    dest="+proj=robin",
    title="Advection on the Sphere",
)

# Initial plot - collect all panel data
all_lons = Float64[]
all_lats = Float64[]
all_vals = Observable(Float64[])

for p in 1:6
    append!(all_lons, rad2deg.(vec(grid.lon[p, :, :])))
    append!(all_lats, rad2deg.(vec(grid.lat[p, :, :])))
end
all_vals[] = vcat([vec(q_anim[p, :, :]) for p in 1:6]...)

surface!(ax_anim, all_lons, all_lats, zeros(length(all_lons));
    color=all_vals, colormap=:viridis,
    colorrange=(0.0, 1.0), shading=NoShading)

record(fig_anim, "advection.gif", 1:nframes; framerate=10) do frame
    for _ in 1:steps_per_frame
        ao = transport_2d(q_anim, courant_xi, courant_eta, grid)
        tendency = evaluate_arrayop(ao)
        for p in 1:6, i in 1:size(tendency, 2), j in 1:size(tendency, 3)
            q_anim[p, i + 1, j + 1] += dt * tendency[p, i, j]
        end
    end
    all_vals[] = vcat([vec(q_anim[p, :, :]) for p in 1:6]...)
end

nothing
```

![Advection animation](advection.gif)
