using Test
include("./search.jl")

@testset "swapNodes" begin
    vars = read_input("input/21_4_1.vrp")
    sol =  read_test_vrp("test.vrp",vars)
    solcp = deepcopy(sol)
    @test isequal(sol , solcp)
end