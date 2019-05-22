using ThreeBodyDecay
using Test


let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    tbs = ThreeBodySystem(m0,m1,m2,m3)
    #
    τ1 = ((tbs.mthsq[1]+tbs.sthsq[1])/2.0, 0.37, 0.2, 0.7, 0.31)
    τ3  = change_basis_3from1(τ1, tbs)
    τ2  = change_basis_2from3(τ3, tbs)
    #
    tilde3to1_for1 = acos(cosθtilde3to1_for1(tbs.s,tbs.msq,[τ1[1],τ2[1],τ3[1]]))
    tilde2to3_for1 = acos(cosθtilde2to3_for1(tbs.s,tbs.msq,[τ1[1],τ2[1],τ3[1]]))
    tilde1to2_for1 = acos(cosθtilde3to1_for1(tbs.s,[tbs.msq[1],tbs.msq[3],tbs.msq[2]],[τ1[1],τ3[1],τ2[1]]))
    @test tilde2to3_for1 ≈ tilde1to2_for1+tilde3to1_for1
end
