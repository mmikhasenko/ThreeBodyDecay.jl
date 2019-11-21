using ThreeBodyDecay
using Test

let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    tbs = ThreeBodySystem(m1,m2,m3,m0)
    #
    dpp = randomPoint(tbs)
    #
    cosζ13_for1_ = cosζ13_for1(dpp.σs,tbs.msq)
    cosζ21_for1_ = cosζ21_for1(dpp.σs,tbs.msq)
    cosζ23_for1_ = cosζ23_for1(dpp.σs,tbs.msq)
    #
    ζ13_for1 = acos(cosζ13_for1_)
    ζ21_for1 = acos(cosζ21_for1_)
    ζ23_for1 = acos(cosζ23_for1_)
    # test of angles
    @test ζ23_for1 ≈ ζ21_for1+ζ13_for1
    #
    # test of phi angle
    prod_13 = [sum(wignerd(1,ν,λ,cosζ23_for1_) *
        wignerd(1,λ,0,cosζ21_for1_) * phase(2,1,2*λ,0)
        for λ=-1:1) for ν=-1:1]
    just_13 = [wignerd(1,ν,0,cosζ13_for1_) for ν=-1:1]
    @test sum(prod_13 .≈ just_13) == 3
end
