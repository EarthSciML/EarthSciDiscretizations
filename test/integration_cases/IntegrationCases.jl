module IntegrationCases

using JSON

export run_case

# Dispatcher contract: `run_case` returns a `(outcome::Symbol, message::String)`
# tuple. Outcome ∈ (:pass, :skip, :fail). A plain tuple keeps the walker free
# of cross-module binding lookups (Julia 1.12 world-age tightening makes those
# fragile for files loaded via `Base.include` at runtime).

include("cartesian_full_pipeline.jl")
include("cubed_sphere_advection.jl")

"""
    run_case(case_spec::AbstractDict, base_dir::AbstractString) -> Tuple{Symbol,String}

Dispatch on `case_spec["kind"]` to a case runner. The manifest path is resolved
relative to `base_dir` (the directory holding the `cases.json` that listed the
case). Returns `(outcome, message)` where outcome ∈ (:pass, :skip, :fail).
"""
function run_case(case_spec, base_dir::AbstractString)
    name = String(get(case_spec, "name", "<unnamed>"))
    kind = String(get(case_spec, "kind", ""))
    manifest_rel = String(get(case_spec, "manifest", ""))
    if isempty(manifest_rel)
        return (:fail, "$name: case missing 'manifest' field")
    end
    manifest_path = isabspath(manifest_rel) ? manifest_rel : joinpath(base_dir, manifest_rel)
    if !isfile(manifest_path)
        return (:fail, "$name: manifest not found at $manifest_rel")
    end
    manifest = JSON.parse(read(manifest_path, String))

    if kind == "stub"
        reason = String(get(manifest, "skip_reason", "stub case (no runner wired yet)"))
        return (:skip, "$name: $reason")
    elseif kind == "cartesian_full_pipeline"
        return run_cartesian_full_pipeline(name, manifest)
    elseif kind == "cubed_sphere_advection"
        return run_cubed_sphere_advection(name, manifest)
    else
        return (:fail, "$name: unknown case kind '$kind'")
    end
end

end # module
