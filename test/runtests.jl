using Test
using NumericalMethods: derivative

@testset "derivatives" begin
    @test derivative((x) -> 2*x, 2) ≈ 2
    @test derivative((x) -> 2*x, 2, "forward") ≈ 2
    @test derivative((x) -> 2*x, 2, "backward") ≈ 2
    @test derivative((x) -> x^2, 2) ≈ 4
end
