using Test
using ThreeBodyDecay



@testset "jp_str macro" begin
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
    @test length(vcat(Pwave)) == 4
    @test length(Set(vcat(Pwave))) == 3
end


@testset "possible ls: binory" begin

@test length(possible_ls(jp"3/2-", jp"3-"; jp=jp"1/2+")) == 4


lsLSv = possible_lsLS(1, jp"1//2+",
    [jp"1+", jp"1/2+", jp"0-", jp"1-"])
@test size(lsLSv) == (1,2)
#
lsLSv = possible_coupling_schemes(1, 1, '+',
    ThreeBodySpins(2,1,0; two_h0=1),
    ThreeBodyParities('+', '+', '-'; P0='-'))
@test size(lsLSv) == (1,2)



# @testset "Coupling schemes" begin
#     # final state particles
#     two_jps = [(1,'+'), (0,'-'), (0,'-')];
#     # decaying particle
#     two_JP_pc = (1,'+'); two_JP_pv = (1,'-')
#     [two_jps...,two_JP_pc]
#     @test length(
#         [coupling_scheme23((2,'-'),[two_jps...,two_JP_pc])...,    # parity conserving
#          coupling_scheme23((2,'-'),[two_jps...,two_JP_pv])...]    # parity violating
#          ) == 4
#     #
#     @test length(
#         [coupling_scheme23((0,'+'),[two_jps...,two_JP_pc])...,    # parity conserving
#          coupling_scheme23((0,'+'),[two_jps...,two_JP_pv])...]    # parity violating
#          ) == 2
# end


# @testset "possible helicities" begin
#     mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
#     tbs = ThreeBodySystem(mJψ,mp,mK,mLb; two_jps = ([2,1,0,1],['-','+','-','+']))
#     dpp = randomPoint(tbs);

#     @test length(possible_helicities(dpp,tbs)) == prod([2,1,0,1].+1)
# end


# TODO:
# test HelicityRecoupling_doublearg
