using Test
using ThreeBodyDecay

@testset "Line-shape functions" begin
    lsh = BreitWigner(0.77,0.15)
    @test real(amp(0.77^2, lsh)) ≈ 0.0
end
