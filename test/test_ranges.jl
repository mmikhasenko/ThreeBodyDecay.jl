using Test
using ThreeBodyDecay

@testset "Definition of ThreeBodySystem" begin
    m1 = 0.938; m2 = 0.49367; m3 = 0.13957; m0 = 2.46867
    tbs = ThreeBodySystem(m0,m1,m2,m3)
    #
    @test tbs.mthsq[1] == (m2+m3)^2
    @test tbs.mthsq[2] == (m3+m1)^2
    @test tbs.mthsq[3] == (m1+m2)^2
    #
    @test tbs.sthsq[1] == (m0-m1)^2
    @test tbs.sthsq[2] == (m0-m2)^2
    @test tbs.sthsq[3] == (m0-m3)^2
    #
end
