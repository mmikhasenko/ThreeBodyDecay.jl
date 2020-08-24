using Test
using ThreeBodyDecay

@testset "Definition of ThreeBodySystem" begin
    m1 = 0.938; m2 = 0.49367; m3 = 0.13957; m0 = 2.46867
    tbs = ThreeBodySystem(m1,m2,m3,m0=m0)
    #
    @test lims1(tbs.ms)[1] == (m2+m3)^2
    @test lims2(tbs.ms)[1] == (m3+m1)^2
    @test lims3(tbs.ms)[1] == (m1+m2)^2
    #
    @test lims1(tbs.ms)[2] == (m0-m1)^2
    @test lims2(tbs.ms)[2] == (m0-m2)^2
    @test lims3(tbs.ms)[2] == (m0-m3)^2
    #
end
