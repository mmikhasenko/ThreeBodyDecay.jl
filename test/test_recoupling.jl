using ThreeBodyDecay
using Test
using Parameters


@testset "BasicParityRecoupling" begin
    reaction = jp"1/2+"=>(jp"1/2+",jp"0-")

    H_pc = ParityRecoupling(1,0, jp"1/2+"=>(jp"1/2+",jp"0-"))
    @test H_pc == ParityRecoupling(1,0,'-')

    H_ls = RecoulingsLS(possible_ls(reaction)[1] |> x2, reaction)

    factor = -1 / sqrt(2)
    @test amplitude(H_pc, 1, 0) * factor ≈ amplitude(H_ls, 1, 0)
    @test amplitude(H_pc,-1, 0) * factor ≈ amplitude(H_ls,-1, 0)
end



@testset "DecayChainCouplings" begin
    mΛb = 5.61960
    mπ = 0.14
    # intermediate resonances
    mΣb = 5.81065; ΓΣb = 0.005
    mΣb_x = 5.83032; ΓΣb_x = 0.0094

    mΛb2S = 6.3; # just a peak of the plot

    tbs = ThreeBodySystem(mπ,mΛb,mπ; m0=mΛb2S,
        two_js = ThreeBodySpins(0, 1, 0; two_h0=1))  # (0- 1/2+ 0- | 1/2+)

    dc_pc = DecayChain(
        k=1,
        Xlineshape=σ->BW(σ,mΣb,ΓΣb),
        tbs=tbs, two_s = 1,
        Hij=ParityRecoupling(1,0,'-'),
        HRk=ParityRecoupling(1,0,'-'))
    # 
    σs = randomPoint(tbs.ms)
    # 
    @test amplitude(σs, [0, 1,0, 1], dc_pc) != 0
    @test amplitude(σs, [0, 1,0,-1], dc_pc) != 0
    @test amplitude(σs, [0,-1,0, 1], dc_pc) != 0
    @test amplitude(σs, [0,-1,0,-1], dc_pc) != 0
    # 
    dc_pv = DecayChain(
        k=1,
        Xlineshape=σ->BW(σ,mΣb,ΓΣb),
        tbs=tbs, two_s = 1,
        Hij=ParityRecoupling(1,0,'-'),
        HRk=NoRecoupling(1,0))
    #
    @test amplitude(σs, [0, 1,0, 1], dc_pv) != 0
    @test amplitude(σs, [0,-1,0, 1], dc_pv) != 0
    @test amplitude(σs, [0,-1,0,-1], dc_pv) == 0im
    @test amplitude(σs, [0,-1,0,-1], dc_pv) == 0im
    # 
    dc_pv = DecayChain(
        k=1,
        Xlineshape=σ->BW(σ,mΣb,ΓΣb),
        tbs=tbs, two_s = 1,
        Hij=ParityRecoupling(1,0,'-'),
        HRk=NoRecoupling(1,0))
    #
    @test amplitude(σs, [0, 1,0, 1], dc_pv) != 0
    @test amplitude(σs, [0,-1,0, 1], dc_pv) != 0
    @test amplitude(σs, [0,-1,0,-1], dc_pv) == 0im
    @test amplitude(σs, [0,-1,0,-1], dc_pv) == 0im
end
