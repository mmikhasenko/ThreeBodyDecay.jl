using ThreeBodyDecay
using Test

Ps = ThreeBodyParities('+', '-', '+'; P0='+')

@testset "Three Body Parities structure" begin
    # 
    @test Ps.P1 == Ps[1] == '+'
    @test Ps.P2 == Ps[2] == '-'
    @test Ps.P3 == Ps[3] == '+'
    @test Ps.P0 == Ps[0] == Ps[4] == '+'
    # 
    @test_throws ErrorException Ps[5]
    @test_throws ErrorException ThreeBodyParities(0, 1, 0; two_js0=1)
    @test_throws ErrorException ThreeBodyParities('+', '+', '+')
end

@testset "operations and interate" begin
    @test sum(i == '+' for i in Ps) == 3
end
