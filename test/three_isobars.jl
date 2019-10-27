using ThreeBodyDecay
using Test
# using StaticArrays

@testset "Three isobars in general decay" begin
    mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
    tbs = ThreeBodySystem(mLb,mJψ,mp,mK; two_jps = ([2,1,0,1],['-','+','-','+']))
    dpp = randomPoint(tbs); dpp.two_λs[1]

    # Lambda
    two_s1 = 3; two_τ1 = rand(-two_s1:2:two_s1) # hat angle 11 is equal to 0
    @test (Zksτ(1, two_s1,two_τ1, dpp, tbs) != 0.0) || (two_Λ(dpp)!=two_τ1-dpp.two_λs[1])
    # JψK
    two_s2 = 2; two_τ2 = rand(-two_s2:2:two_s2)
    @test Zksτ(2, two_s2,two_τ2, dpp, tbs) != 0.0
    # Jψp
    two_s3 = 1; two_τ3 = rand(-two_s3:2:two_s3)
    @test Zksτ(3, two_s3,two_τ3, dpp, tbs) != 0.0
end
