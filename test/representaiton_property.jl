
@testset "Representaiton property" begin
    m1 = 0.938; m2 = 0.49367; m3 = 0.13957; m0 = 2.46867
    tbs = ThreeBodySystem(m0,m1,m2,m3)
    #
    τ1 = ((tbs.mthsq[1]+tbs.sthsq[1])/2.0, 0.3, 0.3, 0.3, 0.3)
    τ3  = change_basis_3from1(τ1, tbs)
    τ2  = change_basis_2from3(τ3, tbs)
    #
    dpp = DalitzPlotPoint(τ1[1],gσ2(τ3[1],τ1[1],tbs),τ3[1])
    #
    s = 4
    #
    @test sum(
        sum(wignerD(s,λ,ν,τ1[3],τ1[2],τ1[5])*
            wignerd(s,ν,μ,cosθhat31(dpp,tbs)) for ν=-s:s) ≈ wignerD(s,λ,μ,τ3[3],τ3[2],τ3[5])
        for λ = -s:s, μ = -s:s) == (2s+1)^2
    #
    @test sum(
        sum(wignerD(s,λ,ν,τ1[3],τ1[2],τ1[5])*
                  wignerd(s,ν,μ,cosθhat12(dpp,tbs))* # angle is negative
                  (mod(ν-μ,2)==1 ? -1 : 1) for ν=-s:s) ≈ wignerD(s,λ,μ,τ2[3],τ2[2],τ2[5])
        for λ = -s:s, μ = -s:s) == (2s+1)^2
    #

    # Half-integer spin
    two_s = 3
    #
    @test sum(
        sum(wignerD_doublearg(two_s,two_λ,two_ν,τ1[3],τ1[2],τ1[5])*
            wignerd_doublearg(two_s,two_ν,two_μ,cosθhat31(dpp,tbs)) for two_ν=-two_s:2:two_s) ≈ wignerD_doublearg(two_s,two_λ,two_μ,τ3[3],τ3[2],τ3[5])
        for two_λ = -two_s:2:two_s, two_μ = -two_s:2:two_s) == (two_s+1)^2

    @test sum(
        sum(wignerD_doublearg(two_s,two_λ,two_ν,τ1[3],τ1[2],τ1[5])*
                  wignerd_doublearg(two_s,two_ν,two_μ,cosθhat12(dpp,tbs))*
                  (mod(two_ν-two_μ,4)==2 ? -1 : 1) for two_ν=-two_s:2:s) ≈ wignerD_doublearg(two_s,two_λ,two_μ,τ2[3],τ2[2],τ2[5])
        for two_λ = -two_s:2:two_s, two_μ = -two_s:2:two_s) == (two_s+1)^2

end
