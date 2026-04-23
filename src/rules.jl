"""
    RuleFile

Reference to a discretization rule JSON file on disk.

Fields:
- `family::Symbol` — scheme family (e.g. `:finite_volume`)
- `name::String`   — rule name (filename without `.json`)
- `path::String`   — absolute path to the JSON file

Parsing and validation against the EarthSciSerialization §7 schema will be
delegated to the ESS rule engine once it lands.
"""
struct RuleFile
    family::Symbol
    name::String
    path::String
end

"""
    load_rules(path::AbstractString) -> Vector{RuleFile}

Discover discretization rule files under `path` (typically the repo's
`discretizations/` directory). The first-level subdirectories are treated as
scheme families; each `*.json` file inside is a rule.

This is a parser delegator stub. It returns `RuleFile` references only; the
file contents are not parsed here. Once the EarthSciSerialization rule
engine lands, this function will delegate to ESS for schema validation and
materialization.
"""
function load_rules(path::AbstractString)
    isdir(path) || throw(ArgumentError("load_rules: not a directory: $path"))
    rules = RuleFile[]
    for family_entry in sort(readdir(path))
        family_dir = joinpath(path, family_entry)
        isdir(family_dir) || continue
        for file in sort(readdir(family_dir))
            endswith(file, ".json") || continue
            push!(
                rules, RuleFile(
                    Symbol(family_entry),
                    String(chop(file; tail = length(".json"))),
                    abspath(joinpath(family_dir, file)),
                )
            )
        end
    end
    return rules
end
