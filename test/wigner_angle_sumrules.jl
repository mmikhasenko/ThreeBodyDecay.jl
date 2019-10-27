using ThreeBodyDecay
using Test

let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    tbs = ThreeBodySystem(m0,m1,m2,m3)
    #
    dpp = randomPoint(tbs)
    #
    ζ13_for1 = acos(cosζ13_for1(dpp.σs,tbs.msq))
    ζ21_for1 = acos(cosζ21_for1(dpp.σs,tbs.msq))
    ζ23_for1 = acos(cosζ23_for1(dpp.σs,tbs.msq))
    @test ζ23_for1 ≈ ζ21_for1+ζ13_for1
end
