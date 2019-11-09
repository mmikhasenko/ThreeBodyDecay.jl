using Test
using ThreeBodyDecay

@testset "calculation of the binary decay" begin
    (p1,p2) = four_vectors_in_binary_decay(2rand()-1, π*(2rand()-1);
        m1sq = 5^2, m2sq = 1^2, m0sq = 8^2)
    @test sum(p1+p2 .≈ [0,0,0,8]) == 4
end

@testset "calculation of the binary decay in flight" begin
    m0sq = 4^2; m1sq = 1.2; m2sq = 1.8
    p0 = [1,2,3,sqrt(m0sq+14)]
    p1, p2 = four_vectors_in_binary_decay(p0,0.3,0.5; m1sq=m1sq, m2sq=m2sq)
    @test sum(p1 + p2 .≈ p0) == 4
    @test invmasssq(p1) ≈ m1sq && invmasssq(p2) ≈ m2sq
end
