using ThreeBodyDecay
using Test

@testset "LS amplitude integer spin" begin
    tbs = ThreeBodySystem(6, 2, 1, 1.5)
    dpp = randomPoint(tbs)
    #
    dc = decay_chain(1, (s,σ)->1/(4.1^2-σ-0.1im); two_s=2, parity='-', tbs=tbs)
    @show amplitude(dpp, dc)
    @test sum(reim(amplitude(dpp, dc)) .≈ 0.0) == 0
    # testing something else?
end

@testset "LS amplitude integer spin" begin
    tbs = ThreeBodySystem(6, 2, 1, 1.5;
        two_jps = ([1, 0, 0, 1], ['+','-','-','+']))  # 1/2+ 0- 0- 1/2+
    dpp = randomPoint(tbs)
    #
    dc = decay_chain(3, (s,σ)->1/(4.1^2-σ-0.1im); two_s=3, parity='-', tbs=tbs)
    @show amplitude(dpp, dc)
    @test sum(reim(amplitude(dpp, dc)) .≈ 0.0) == 0
    # testing something else?
end
