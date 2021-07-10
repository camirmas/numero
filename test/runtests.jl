using Test
using NumericalMethods

@testset "derivatives" begin
    @test derivative((x) -> 2*x, 2) ≈ 2
    @test derivative((x) -> 2*x, 2, "forward") ≈ 2
    @test derivative((x) -> 2*x, 2, "backward") ≈ 2
    @test derivative((x) -> x^2, 2) ≈ 4
end

@testset "root finding" begin
    @test bisect((x) -> x^2-1, 0, 2, .5, 10, 0) ≈ 1
    @test false_position((x) -> x^2-1, 0, 2, 0, 100, .5) ≈ 1
end
