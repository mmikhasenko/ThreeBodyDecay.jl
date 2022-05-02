using Test
using ThreeBodyDecay

@testset "jp_str macro" begin
    @test jp"3+" == str2jp("3+")
    @test jp"3^+" == jp"3+"
    @test jp"1/2+" == jp(1//2, '+')
    @test jp"3/2-" == jp(3//2, '-')
    @test jp"3+" == jp(3, '+')
end


@testset "jp ⊗ jp" begin
    Swave = jp(1//2,'+') ⊗ jp(1,'-')
    Pwave = [sw ⊗ jp(1,'-') for sw in Swave]
    Dwave = [sw ⊗ jp(2,'+') for sw in Swave]
    # 
    @test length(Swave) == 2 
    @test length(vcat(Pwave...)) == 5
    @test length(Set(vcat(Pwave...))) == 3
end

let
    @test length(possible_ls(jp"3/2-", jp"3-"; jp=jp"1/2+")) == 4
    # 
    lsLSv = possible_lsLS(1, jp"1/2+",
        [jp"1+", jp"1/2+", jp"0-", jp"1/2-"])
    @test size(lsLSv) == (1,2)
    #
    lsLSv = possible_lsLS(1, 1, '+',
        ThreeBodySpins(2,1,0; two_h0=1),
        ThreeBodyParities('+', '+', '-'; P0='-'))
    @test size(lsLSv) == (1,2)
end

@testset "x2 over2" begin
    @test x2(1//2) == 1
    @test x2(3) == 6
    @test x2((2,3//2)) == (4,3)
    # 
    @test over2(3) == 3//2
    @test over2(2) == 1
    @test over2((2,3)) == (1,3//2)
end
