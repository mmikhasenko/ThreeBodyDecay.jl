using ThreeBodyDecay
using Test

@testset "Wigner angle permutations" begin

    @test phase(1,2,rand(-2:2),rand(-2:2)) == 1.0
    @test phase(2,3,rand(-2:2),rand(-2:2)) == 1.0
    @test phase(3,1,rand(-2:2),rand(-2:2)) == 1.0
    #
    @test phase(2,1,1,-1) == -1.0
    @test phase(3,2,1,-1) == -1.0
    @test phase(1,3,1,-1) == -1.0
    #
    @test phase(1,1,1,-1) == 0.0
    @test phase(2,2,1,1) == 1.0
    @test phase(3,3,3,-1) == 0.0

end
