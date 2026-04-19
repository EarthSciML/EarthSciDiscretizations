using Test
using TestItems

@testitem "load_rules discovers rule files by family" begin
    using EarthSciDiscretizations: load_rules, RuleFile

    mktempdir() do root
        mkdir(joinpath(root, "finite_volume"))
        mkdir(joinpath(root, "finite_difference"))
        write(joinpath(root, "finite_volume", "ppm.json"), "{}")
        write(joinpath(root, "finite_volume", "muscl.json"), "{}")
        write(joinpath(root, "finite_difference", "central_4th.json"), "{}")
        write(joinpath(root, "finite_volume", "notes.txt"), "ignored")

        rules = load_rules(root)

        @test length(rules) == 3
        @test all(r -> r isa RuleFile, rules)
        families = sort(unique([r.family for r in rules]))
        @test families == [:finite_difference, :finite_volume]
        fv_names = sort([r.name for r in rules if r.family == :finite_volume])
        @test fv_names == ["muscl", "ppm"]
        @test all(r -> endswith(r.path, ".json"), rules)
        @test all(r -> isabspath(r.path), rules)
    end
end

@testitem "load_rules errors on missing directory" begin
    using EarthSciDiscretizations: load_rules

    @test_throws ArgumentError load_rules(joinpath(tempdir(), "esd_nonexistent_xyz"))
end

@testitem "load_rules returns empty for empty catalog" begin
    using EarthSciDiscretizations: load_rules

    mktempdir() do root
        @test isempty(load_rules(root))
    end
end
