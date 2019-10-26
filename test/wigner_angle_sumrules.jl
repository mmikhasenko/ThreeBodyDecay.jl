using ThreeBodyDecay
using Test

let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    tbs = ThreeBodySystem(m0,m1,m2,m3)
    #
    τ1 = ((tbs.mthsq[1]+tbs.sthsq[1])/2.0, 0.37, 0.2, 0.7, 0.31)
    τ3  = change_basis_3from1(τ1, tbs)
    τ2  = change_basis_2from3(τ3, tbs)
    dp = DalitzPlotPoint(τ1[1],τ2[1],τ3[1])
    #
    ζ13_for1 = acos(cosζ13_for1(dp,tbs))
    ζ21_for1 = acos(cosζ21_for1(dp,tbs))
    ζ23_for1 = acos(cosζ23_for1(dp,tbs))
    @test ζ23_for1 ≈ ζ21_for1+ζ13_for1
end
