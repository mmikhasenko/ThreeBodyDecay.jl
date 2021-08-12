using ThreeBodyDecay
using Test

ms = ThreeBodyMasses(1,2,3.0; m0=120)

@testset "Three Body Masses structure" begin
    # 
    @test ms.m1 == ms[1] == 1
    @test ms.m2 == ms[2] == 2
    @test ms.m3 == ms[3] == 3.0
    @test ms.m0 == ms[0] == ms[4] == 120
    # 
    @test_throws ErrorException ms[5]
    @test_throws ErrorException ThreeBodyMasses(0,1,1; m0=1)
    @test_throws ErrorException ThreeBodyMasses(0,1,1)
end

@testset "operations and interate" begin
    @test sum(ms^2) == sum(abs2, ms)
end

@testset "Limits for ThreeBodyMasses" begin
    m1 = 0.938; m2 = 0.49367; m3 = 0.13957; m0 = 2.46867
    ms = ThreeBodyMasses(m1=m1,m2=m2,m3=m3,m0=m0)
    #
    @test lims1(ms)[1] == (m2+m3)^2
    @test lims2(ms)[1] == (m3+m1)^2
    @test lims3(ms)[1] == (m1+m2)^2
    #
    @test lims1(ms)[2] == (m0-m1)^2
    @test lims2(ms)[2] == (m0-m2)^2
    @test lims3(ms)[2] == (m0-m3)^2
    #
    @test lims1(ms) == lims(1, ms)
    @test lims2(ms) == lims(2, ms)
    @test lims3(ms) == lims(3, ms)
    # 
    σs = randomPoint(ms)
    @test ThreeBodyDecay.inrange(σs[1], lims1(ms))
    @test ThreeBodyDecay.inrange(σs[2], lims2(ms))
    @test ThreeBodyDecay.inrange(σs[3], lims3(ms))
    @test inphrange(σs, ms)
end
