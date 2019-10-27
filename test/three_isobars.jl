using ThreeBodyDecay
using Test
using StaticArrays

@testset "Three isobars in general decay" begin

    mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
    tbs = ThreeBodySystem(mLb,mJψ,mp,mK; two_JP = (1,'+'), two_jps = ([2,1,0],['-','+','-']))
    dpp = randomPoint(tbs)

    two_Λ = 1; two_λs = (0,1,0)
    @test Zsτ(1, 1,1, two_Λ,two_λs, dpp, tbs) != 0.0 # Lambda
    @test Zsτ(2, 2,0, two_Λ,two_λs, dpp, tbs) != 0.0 # JψK
    @test Zsτ(3, 3,1, two_Λ,two_λs, dpp, tbs) != 0.0 # Jψp
end
