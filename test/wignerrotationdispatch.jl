
using ThreeBodyDecay
using Test

@testset "correct typearg" begin
    @test typeof(wr(1,1,2)) <: TriavialWignerRotation
    @test typeof(wr(2,2))   <: TriavialWignerRotation
    # 
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
    # 
    @test iseven(wr(2,1,2)) == false
    @test iseven(wr(2,1,1)) == true
    @test iseven(wr(1,2,2)) == false
    @test iseven(wr(1,2,1)) == true
end
