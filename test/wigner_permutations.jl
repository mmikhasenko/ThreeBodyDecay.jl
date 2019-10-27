using ThreeBodyDecay
using Test

@testset "Wigner angle permutations" begin

    mp = 0.938; mK = 0.49367; mpi  = 0.13957; mXic = 2.46867
    tbs = ThreeBodySystem(mXic,mp,mK,mpi)
    dpp = randomPoint(tbs)

    # (23) cosζ31_for1 = cosζ23_for1
    @test cosζk1_for1(2,[dpp.σ123[1],dpp.σ123[2],dpp.σ123[3]],[tbs.msq[1],tbs.msq[2],tbs.msq[3],s(tbs)]) ==
          cosζk1_for1(3,[dpp.σ123[1],dpp.σ123[3],dpp.σ123[2]],[tbs.msq[1],tbs.msq[3],tbs.msq[2],s(tbs)])
    # (31) cosζ21_for2 = cosζ32_for2
    @test cosζk2_for2(3,[dpp.σ123[1],dpp.σ123[2],dpp.σ123[3]],[tbs.msq[1],tbs.msq[2],tbs.msq[3],s(tbs)]) ==
          cosζk2_for2(1,[dpp.σ123[3],dpp.σ123[2],dpp.σ123[1]],[tbs.msq[3],tbs.msq[2],tbs.msq[1],s(tbs)])
    # (12) cosζ32_for3 = cosζ13_for3
    @test cosζk3_for3(1,[dpp.σ123[1],dpp.σ123[2],dpp.σ123[3]],[tbs.msq[1],tbs.msq[2],tbs.msq[3],s(tbs)]) ==
          cosζk3_for3(2,[dpp.σ123[2],dpp.σ123[1],dpp.σ123[3]],[tbs.msq[2],tbs.msq[1],tbs.msq[3],s(tbs)])

    # (123) cosζ31_for1 = cosζ21_for2
    @test cosζk2_for2(1,[dpp.σ123[1],dpp.σ123[2],dpp.σ123[3]],[tbs.msq[1],tbs.msq[2],tbs.msq[3],s(tbs)]) ==
          cosζk1_for1(3,[dpp.σ123[2],dpp.σ123[3],dpp.σ123[1]],[tbs.msq[2],tbs.msq[3],tbs.msq[1],s(tbs)])
    # (123)^2 cosζ31_for1 = cosζ23_for3
    @test cosζk3_for3(2,[dpp.σ123[1],dpp.σ123[2],dpp.σ123[3]],[tbs.msq[1],tbs.msq[2],tbs.msq[3],s(tbs)]) ==
          cosζk1_for1(3,[dpp.σ123[3],dpp.σ123[1],dpp.σ123[2]],[tbs.msq[3],tbs.msq[1],tbs.msq[2],s(tbs)])
    #
end
