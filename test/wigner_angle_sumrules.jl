using ThreeBodyDecay
using ThreeBodyDecay.PartialWaveFunctions
using Test

let
    ms = ThreeBodyMasses(m1=0.938, m2=0.49367, m3=0.13957, m0=2.46867)
    #
    σs = randomPoint(ms)
    #
    cosζ13_for1_ = cosζ(wr(1, 3, 1), σs, ms^2)
    cosζ21_for1_ = cosζ(wr(2, 1, 1), σs, ms^2)
    cosζ23_for1_ = cosζ(wr(2, 3, 1), σs, ms^2)
    #
    ζ13_for1 = acos(cosζ13_for1_)
    ζ21_for1 = acos(cosζ21_for1_)
    ζ23_for1 = acos(cosζ23_for1_)
    # test of angles
    @test ζ23_for1 ≈ ζ21_for1 + ζ13_for1
    #
    # sum rules for the combination of d functions
    prod_13 = [sum(
        wignerd(1, ν, λ, cosζ23_for1_) *
        wignerd(1, λ, 0, cosζ21_for1_) * minusone()^λ
        for λ = -1:1) for ν = -1:1]
    just_13 = [wignerd(1, ν, 0, cosζ13_for1_) for ν = -1:1]
    @test sum(prod_13 .≈ just_13) == 3
end

