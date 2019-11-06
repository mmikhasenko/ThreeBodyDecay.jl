using Test
using ThreeBodyDecay

@testset "calculation of the binary decay" begin
    (p1,p2) = four_vectors_in_binary_decay(2rand()-1, π*(2rand()-1);
        m1sq = mD0sq, m2sq = mπ0sq, Msq = mDstar0sq)
    @test invmasssq(p1+p2) ≈ mDstar0sq
end
