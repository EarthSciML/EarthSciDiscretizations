@testsnippet VertRemapSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Vertical remap placeholder" setup=[VertRemapSetup] tags=[:vertical] begin
    q = randn(10)
    dp_old = ones(10)
    dp_new = ones(10)
    @test_throws ErrorException vertical_remap_tendency(q, dp_old, dp_new, 10)
end
