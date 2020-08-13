using Test
using ThreeBodyDecay

@testset "Coupling schemes" begin
    # final state particles
    two_jps = [(1,'+'), (0,'-'), (0,'-')];
    # decaying particle
    two_JP_pc = (1,'+'); two_JP_pv = (1,'-')
    [two_jps...,two_JP_pc]
    @test length(
        [coupling_scheme23((2,'-'),[two_jps...,two_JP_pc])...,    # parity conserving
         coupling_scheme23((2,'-'),[two_jps...,two_JP_pv])...]    # parity violating
         ) == 4
    #
    @test length(
        [coupling_scheme23((0,'+'),[two_jps...,two_JP_pc])...,    # parity conserving
         coupling_scheme23((0,'+'),[two_jps...,two_JP_pv])...]    # parity violating
         ) == 2
end


@testset "possible helicities" begin
    mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
    tbs = ThreeBodySystem(mLb, mJψ, mp, mK; two_js = [2,1,0,1])

    @test length(possible_helicities(tbs.two_js)) == prod([2,1,0,1].+1)
end

# TODO:
# test HelicityRecoupling_doublearg
