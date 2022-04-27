using Test
using ThreeBodyDecay

ms = ThreeBodyMasses(1.1,2.2,3.3;m0=10)

@testset "cosζ 2arg functions" begin
    σs = randomPoint(ms)
    @test cosζ21_for1(σs, ms^2) == cosζ(wr(2,1,1), σs, ms^2)
    @test cosζ21_for2(σs, ms^2) == cosζ(wr(2,1,2), σs, ms^2)
    @test cosζ13_for1(σs, ms^2) == cosζ(wr(1,3,1), σs, ms^2)
    @test cosζ13_for3(σs, ms^2) == cosζ(wr(1,3,3), σs, ms^2)
    @test cosζ32_for3(σs, ms^2) == cosζ(wr(3,2,3), σs, ms^2)
    @test cosζ32_for2(σs, ms^2) == cosζ(wr(3,2,2), σs, ms^2)
end

@testset "cosζ 3arg functions" begin
    σs = randomPoint(ms)
    @test cosζ12_for3(σs, ms^2) == cosζ(wr(1,2,3), σs, ms^2)
    @test cosζ23_for1(σs, ms^2) == cosζ(wr(2,3,1), σs, ms^2)
    @test cosζ31_for2(σs, ms^2) == cosζ(wr(3,1,2), σs, ms^2)
end

@testset "cosζ⁰ functions" begin
    σs = randomPoint(ms)
    @test cosθhat12(σs, ms^2) ≈ cosζ(wr(1,2,0), σs, ms^2)
    @test cosθhat23(σs, ms^2) ≈ cosζ(wr(2,3,0), σs, ms^2)
    @test cosθhat31(σs, ms^2) ≈ cosζ(wr(3,1,0), σs, ms^2)
end
