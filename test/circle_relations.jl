

let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    tbs = ThreeBodySystem(m0,m1,m2,m3)
    #
    τ1 = ((tbs.mthsq[1]+tbs.sthsq[1])/2.0, 0.3, 0.3, 0.3, 0.3)
    τ3  = change_basis_3from1(τ1, tbs)
    # @show [τ1...,        tbs.msq...,    tbs.s]
    τ2  = change_basis_2from3(τ3, tbs)
    τ1p = change_basis_1from2(τ2, tbs)
    @test sum(τ1p .≈ τ1) == 5

end
