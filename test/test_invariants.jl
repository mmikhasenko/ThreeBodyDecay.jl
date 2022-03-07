using ThreeBodyDecay
using Test

rands = rand(3)
σs = Invariants(rands...)

@testset "Invariants structure" begin
    # 
    @test σs.σ1 == σs[1] == rands[1]
    @test σs.σ2 == σs[2] == rands[2]
    @test σs.σ3 == σs[3] == rands[3]
    # 
    @test_throws ErrorException σs[5]
end

@testset "operations and interate" begin
    @test sum(σs) == sum(rands)
    @test length(σs) == 3
end

ms = ThreeBodyMasses(1.1, 3.3,5.5; m0=20.0)
σs = Invariants(ms; σ3=4.5^2, σ2=10.1^2)

@testset "Creation from two given" begin
    @test sum(σs .≈ Invariants(ms; σ1=σs.σ1, σ2=σs.σ2))==3
    @test sum(σs .≈ Invariants(ms; σ2=σs.σ2, σ3=σs.σ3))==3
    @test sum(σs .≈ Invariants(ms; σ3=σs.σ3, σ1=σs.σ1))==3
    # 
    @test nt(σs).σ1 == σs.σ1
    @test_throws ErrorException Invariants(ms; σ1=1.0, σ2=1.0, σ3=1.0)
    # 
    @test Kibble(σs,ms^2) < 0
end

@testset "Random points" begin
    σs = randomPoint(ms)
    @test Kibble(σs,ms^2) < 0
end

tbs = ThreeBodySystem(ms)

@testset "Three-body system" begin
    @test tbs == ThreeBodySystem(ms=ms,
        two_js = ThreeBodySpins(0,0,0;two_h0=0))
    @test tbs == ThreeBodySystem(ms,
        ThreeBodySpins(0,0,0;two_h0=0))
    @test tbs == ThreeBodySystem(ms.m1, ms.m2, ms.m3; m0=ms.m0,
        two_js = ThreeBodySpins(0,0,0;two_h0=0))    
end

tbsS = ThreeBodySystem(ms,ThreeBodySpins(1,3,5;two_h0=1))
rp = randomPoint(tbsS)

@testset "Random points" begin
    @test Kibble(rp.σs, tbsS.ms^2) < 0
    @test Tuple(rp.two_λs) ∈ collect(itr(tbsS.two_js))
end

