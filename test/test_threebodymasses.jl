using ThreeBodyDecay
using Test

ms = ThreeBodyMasses(1.0, 2.0, 3.0; m0=120.0)

@testset "Three Body Masses structure" begin
    # 
    @test ms.m1 == ms[1] == 1
    @test ms.m2 == ms[2] == 2
    @test ms.m3 == ms[3] == 3.0
    @test ms.m0 == ms[4] == 120
    # 
    @test_throws BoundsError ms[5]
    @test_throws ErrorException ThreeBodyMasses(0, 1, 1; m0=1)
    @test_throws UndefKeywordError ThreeBodyMasses(0, 1, 1)
end

@testset "operations and interate" begin
    ms² = ms^2
    @test sum(ms²) == sum(abs2, ms)
    # 
    @test ms²[1] == ms[1]^2
    @test ms²[2] == ms[2]^2
    @test ms²[3] == ms[3]^2
    @test ms²[4] == ms[4]^2
end

# @testset "Limits for ThreeBodyMasses" begin
let
    m1 = 0.938
    m2 = 0.49367
    m3 = 0.13957
    m0 = 2.46867
    ms = ThreeBodyMasses(m1, m2, m3; m0=m0)
    #
    @test lims1(ms)[1] == (m2 + m3)^2
    @test lims2(ms)[1] == (m3 + m1)^2
    @test lims3(ms)[1] == (m1 + m2)^2
    #
    @test lims1(ms)[2] == (m0 - m1)^2
    @test lims2(ms)[2] == (m0 - m2)^2
    @test lims3(ms)[2] == (m0 - m3)^2
    #
    @test lims1(ms) == lims(1, ms)
    @test lims2(ms) == lims(2, ms)
    @test lims3(ms) == lims(3, ms)
end
