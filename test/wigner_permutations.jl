using ThreeBodyDecay
using Parameters
using Test

@testset "Wigner angle permutations" begin

    mp = 0.938; mK = 0.49367; mpi  = 0.13957; mXic = 2.46867
    tbs = ThreeBodySystem(mXic,mp,mK,mpi)
    σs = randomPoint(tbs.ms)

    @unpack m1,m2,m3,m0 = tbs.ms
    # (23) cosζ31_for1 = cosζ23_for1
    @test cosζk1_for1(2,[σs.σ1,σs.σ2,σs.σ3],[m1^2,m2^2,m3^2,m0^2]) ==
          cosζk1_for1(3,[σs.σ1,σs.σ3,σs.σ2],[m1^2,m3^2,m2^2,m0^2])
    # (31) cosζ21_for2 = cosζ32_for2
    @test cosζk2_for2(3,[σs.σ1,σs.σ2,σs.σ3],[m1^2,m2^2,m3^2,m0^2]) ==
          cosζk2_for2(1,[σs.σ3,σs.σ2,σs.σ1],[m3^2,m2^2,m1^2,m0^2])
    # (12) cosζ32_for3 = cosζ13_for3
    @test cosζk3_for3(1,[σs.σ1,σs.σ2,σs.σ3],[m1^2,m2^2,m3^2,m0^2]) ==
          cosζk3_for3(2,[σs.σ2,σs.σ1,σs.σ3],[m2^2,m1^2,m3^2,m0^2])

    # (123) cosζ31_for1 = cosζ21_for2
    @test cosζk2_for2(1,[σs.σ1,σs.σ2,σs.σ3],[m1^2,m2^2,m3^2,m0^2]) ==
          cosζk1_for1(3,[σs.σ2,σs.σ3,σs.σ1],[m2^2,m3^2,m1^2,m0^2])
    # (123)^2 cosζ31_for1 = cosζ23_for3
    @test cosζk3_for3(2,[σs.σ1,σs.σ2,σs.σ3],[m1^2,m2^2,m3^2,m0^2]) ==
          cosζk1_for1(3,[σs.σ3,σs[1],σs[2]],[m3^2,m1^2,m2^2,m0^2])
    #
end
