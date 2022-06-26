using ThreeBodyDecay
using Test

two_js = ThreeBodySpins(1,0,2; two_h0=3)

@testset "Three Body Masses structure" begin
    # 
    @test two_js.two_h1 == two_js[1] == 1
    @test two_js.two_h2 == two_js[2] == 0
    @test two_js.two_h3 == two_js[3] == 2
    @test two_js.two_h0 == two_js[4] == 3
    # 
    @test_throws ErrorException ThreeBodySpins(0,1,1; two_js0=1)
    @test_throws ErrorException ThreeBodySpins(0,1,1)
    @test_throws BoundsError two_js[5]
end

@testset "operations and interate" begin
    @test iseven(sum(two_js))
    @test div(sum(Tuple(two_js) .|> x2), 2) == sum(two_js)
end