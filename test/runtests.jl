using Test
using NumericalMethods

@testset "derivatives" begin
    @test derivative((x) -> 2*x, 2) ≈ 2
    @test derivative((x) -> 2*x, 2, "forward") ≈ 2
    @test derivative((x) -> 2*x, 2, "backward") ≈ 2
    @test derivative((x) -> x^2, 2) ≈ 4
end

@testset "bracketing methods" begin
    @test bisect((x) -> x^2-1, 0, 2, .5, 10, 0) ≈ 1
    @test false_position((x) -> x^2-1, 0, 2, 0, 100, .5) ≈ 1
end

@testset "open methods" begin
    @test fixed_point((x) -> 1, 0, 0, 10) == 1 # f(x) = x^2-1
    @test fixed_point((x) -> exp(-x), 0, .5, 10) ≈ .5648793473910495 # f(x) = e^(-x)-x
    @test newton_raphson((x -> exp(-x)-x), 0, .5, 10) ≈ .5671431598525681
end

@testset "Gauss elimination" begin
    m = [3 -.1 -.2
		 .1 7 -.3
	     .3 -.2 10]
	b = [7.85, -19.3, 71.4]
    @test naive_gauss(m, b) ≈ [3.0, -2.5, 7.0]
end

@testset "One-Dimensional Unconstrained Optimization" begin
    f(x) = 2sin(x)-x^2/10
	@test gold(0, 4, 8, .01, f) ≈ 1.4427190999915878
end
