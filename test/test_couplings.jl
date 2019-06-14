
@testset "Coupling schemes" begin
    # final state particles
    two_jps = [(1,"+"), (0,"-"), (0,"-")];
    # decaying particle
    two_JP_pc = (1,"+"); two_JP_pv = (1,"-")
    @test length(
        [coupling_scheme23(two_JP_pc,(2,"-"),two_jps)...,    # parity conserving
         coupling_scheme23(two_JP_pv,(2,"-"),two_jps)...]    # parity violating
         ) == 4
    #
    @test length(
        [coupling_scheme23(two_JP_pc,(0,"+"),two_jps)...,    # parity conserving
         coupling_scheme23(two_JP_pv,(0,"+"),two_jps)...]    # parity violating
         ) == 2
end

# TODO:
# test HelicityRecoupling_doublearg
