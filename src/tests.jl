using Test
include("./search.jl")

@testset "eucDist" begin
    @test euc_dist([0,0],[0,0]) ==  0
end