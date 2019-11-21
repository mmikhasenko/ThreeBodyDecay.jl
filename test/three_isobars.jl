using ThreeBodyDecay
using Test
# using StaticArrays

@testset "Three isobars in general decay" begin
    mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
    tbs = ThreeBodySystem(mJψ,mp,mK,mLb; two_jps = ([2,1,0,1],['-','+','-','+']))
    dpp = randomPoint(tbs);

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


# mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
# tbs = ThreeBodySystem(mJψ,mp,mK,mLb; two_jps = ([2,1,0,1],['-','+','-','+']))
# dpp = randomPoint(tbs);
# dpp.two_λs[1]
# two_s1 = 3; two_τ1 = rand(-two_s1:2:two_s1) # hat angle 11 is equal to 0
# # coupling scheme
# CScheme1 = coupling_scheme23((two_s1,'+'),tbs)
# # couplings
# cs1 = rand(length(CS1))
# # helicity couplings
# hτ = sum(c*clebsch_for_chaink(1, (two_s1, two_τ1), CouSch, dpp.two_λs, tbs.two_js) for (c,CS) in zip(coupls1,CScheme1))
#
# # Lambda
# @test (Zksτ(1, two_s1,two_τ1, dpp, tbs) != 0.0) || (two_Λ(dpp)!=two_τ1-dpp.two_λs[1])
#     # JψK
#     two_s2 = 2; two_τ2 = rand(-two_s2:2:two_s2)
#     @test Zksτ(2, two_s2,two_τ2, dpp, tbs) != 0.0
#     # Jψp
#     two_s3 = 1; two_τ3 = rand(-two_s3:2:two_s3)
#     @test Zksτ(3, two_s3,two_τ3, dpp, tbs) != 0.0
# end
