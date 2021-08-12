
using ThreeBodyDecay
using Test

@testset "correct typearg" begin
    @test typeof(wr(1,1,2)) <: TrivialWR
    @test typeof(wr(2,2))   <: TrivialWR
    # 
    @test typeof(wr(3,1)) <: HatWR{true}
    @test typeof(wr(1,2)) <: HatWR{true}
    # 
    @test typeof(wr(1,2,2)) <: ZetaRepWR{:D,false}
    @test typeof(wr(1,2,1)) <: ZetaRepWR{:S,false}
    @test typeof(wr(1,3,2)) <: ZetaAllWR{false}
    @test typeof(wr(2,1,2)) <: ZetaRepWR{:D,true}
    @test typeof(wr(2,1,1)) <: ZetaRepWR{:S,true}
    @test typeof(wr(1,2,3)) <: ZetaAllWR{true}
end

@testset "Properties of xxx" begin
    @test ispositive(wr(1,1,2)) == true
    @test ispositive(wr(3,3)) == true
    # 
    @test ispositive(wr(1,2,2)) == false
    @test ispositive(wr(2,3,2)) == false
    @test ispositive(wr(2,3,3)) == false
    @test ispositive(wr(3,2,3)) == true
    @test ispositive(wr(3,2,2)) == true
    # 
    @test ispositive(wr(1,2,3)) == true
    @test ispositive(wr(3,2,1)) == false
end

@testset "Consistency with old implementation" begin
    ms = ThreeBodyMasses(m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867)
    #
    σs = randomPoint(ms)

    @test cosζ(wr(3,1,0),σs,ms^2) ≈ cosθhat31(σs,ms^2)
    @test cosζ(wr(1,2,0),σs,ms^2) ≈ cosθhat12(σs,ms^2)
    @test cosζ(wr(2,3,0),σs,ms^2) ≈ cosθhat23(σs,ms^2)
    # 
    @test cosζ(wr(2,3,1),σs,ms^2) ≈ cosζ23_for1(σs,ms^2)
    @test cosζ(wr(1,2,3),σs,ms^2) ≈ cosζ12_for3(σs,ms^2)
    @test cosζ(wr(3,1,2),σs,ms^2) ≈ cosζ31_for2(σs,ms^2)
    # # 
    @test cosζ(wr(2,1,1),σs,ms^2) ≈ cosζ21_for1(σs,ms^2)
    @test cosζ(wr(2,1,2),σs,ms^2) ≈ cosζ21_for2(σs,ms^2)
    @test cosζ(wr(1,3,1),σs,ms^2) ≈ cosζ13_for1(σs,ms^2)
    @test cosζ(wr(1,3,3),σs,ms^2) ≈ cosζ13_for3(σs,ms^2)
    @test cosζ(wr(3,2,3),σs,ms^2) ≈ cosζ32_for3(σs,ms^2)
    @test cosζ(wr(3,2,2),σs,ms^2) ≈ cosζ32_for2(σs,ms^2)
end