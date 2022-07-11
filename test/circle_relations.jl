using Test
using ThreeBodyDecay

@testset "Circular permutations" begin
    m1 = 0.938
    m2 = 0.49367
    m3 = 0.13957
    m0 = 2.46867
    tbs = ThreeBodySystem(m1, m2, m3, m0=m0)
    #
    τ1 = (sum(lims1(tbs.ms)) / 2.0, 0.3, 0.3, 0.3, 0.3)
    τ3 = change_basis_3from1(τ1, tbs.ms)
    # @show [τ1...,        tbs.ms^2...,    tbs.s]
    τ2 = change_basis_2from3(τ3, tbs.ms)
    τ1p = change_basis_1from2(τ2, tbs.ms)
    @test sum(τ1p .≈ τ1) == 5

    #
    @test τ1[1] ≈ σ1of2(τ2[4], τ2[1], tbs.ms^2)
    @test τ3[1] ≈ σ3of1(τ1[4], τ1[1], tbs.ms^2)
    @test τ2[1] ≈ σ2of3(τ3[4], τ3[1], tbs.ms^2)
    #
    # implementation cosθ12
    σs = (τ1[1], τ2[1], τ3[1])
    @test τ3[4] ≈ cosθ12(σs, tbs.ms^2)
    @test τ1[4] ≈ cosθ23(σs, tbs.ms^2)
    @test τ2[4] ≈ cosθ31(σs, tbs.ms^2)
    #
end
