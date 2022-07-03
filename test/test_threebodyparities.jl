using ThreeBodyDecay
using Test

Ps = ThreeBodyParities('+', '-', '+'; P0='+')

@testset "Three Body Parities structure" begin
    # 
    @test Ps.P1 == Ps[1] == '+'
    @test Ps.P2 == Ps[2] == '-'
    @test Ps.P3 == Ps[3] == '+'
    @test Ps.P0 == Ps[4] == '+'
    # 
    @test_throws BoundsError Ps[5]
    @test_throws ErrorException ThreeBodyParities('+', '+', '+')
end

@testset "operations and interate" begin
    @test sum(i == '+' for i in Ps) == 3
end
