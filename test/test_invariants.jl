using ThreeBodyDecay
using Test


ms = ThreeBodyMasses(1.1, 3.3, 5.5; m0=20.0)
σs = Invariants(ms; σ3=4.5^2, σ2=10.1^2)

@testset "Invariants structure" begin
    # 
    @test σs.σ1 == σs[1]
    @test σs.σ2 == σs[2]
    @test σs.σ3 == σs[3]
    # 
    @test_throws BoundsError σs[5]
end

@testset "operations and interate" begin
    @test sum(σs) == sum(ms^2)
    @test length(σs) == 3
end

@testset "Creation from two given" begin
    @test σs == Invariants(ms; σ1=σs.σ1, σ2=σs.σ2)
    @test σs == Invariants(ms; σ2=σs.σ2, σ3=σs.σ3)
    @test σs == Invariants(ms; σ3=σs.σ3, σ1=σs.σ1)
    # 
    @test_throws ErrorException Invariants(ms; σ1=1.0, σ2=1.0, σ3=1.0)
    # 
    @test Kibble(σs, ms^2) < 0
end


@testset "consistency of cosθij and σkofj, σkofj functions" begin
    z23 = cosθ23(σs, ms^2)
    @test σs.σ3 ≈ σ3of1(z23, σs.σ1, ms^2)
    @test σs.σ2 ≈ σ2of1(z23, σs.σ1, ms^2)
    # 
    z12 = cosθ12(σs, ms^2)
    @test σs.σ1 ≈ σ1of3(z12, σs.σ3, ms^2)
    @test σs.σ2 ≈ σ2of3(z12, σs.σ3, ms^2)
    # 
    z31 = cosθ31(σs, ms^2)
    @test σs.σ3 ≈ σ3of2(z31, σs.σ2, ms^2)
    @test σs.σ1 ≈ σ1of2(z31, σs.σ2, ms^2)
end


@testset "Random points" begin
    σs = randomPoint(ms)
    @test Kibble(σs, ms^2) < 0
end

tbs = ThreeBodySystem(ms)

@testset "Three-body system" begin
    @test tbs == ThreeBodySystem(ms=ms,
        two_js=ThreeBodySpins(0, 0, 0; two_h0=0))
    @test tbs == ThreeBodySystem(ms,
        ThreeBodySpins(0, 0, 0; two_h0=0))
    @test tbs == ThreeBodySystem(ms.m1, ms.m2, ms.m3; m0=ms.m0,
        two_js=ThreeBodySpins(0, 0, 0; two_h0=0))
end

tbsS = ThreeBodySystem(ms, ThreeBodySpins(1, 3, 5; two_h0=1))
rp = randomPoint(tbsS)

@testset "Random points" begin
    @test Kibble(rp.σs, tbsS.ms^2) < 0
    @test Tuple(rp.two_λs) ∈ collect(itr(tbsS.two_js))
end

