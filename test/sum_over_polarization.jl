using Test
using ThreeBodyDecay

@testset "Sum over polarization" begin
    mΛb = 5.61960
    mπ = 0.14
    # intermediate resonances
    mΣb = 5.81065; ΓΣb = 0.005
    mΣb_x = 5.83032; ΓΣb_x = 0.0094

    mΛb2S = 6.3; # just a peak of the plot

    tbs = ThreeBodySystem(mπ,mΛb,mπ,mΛb2S;
        two_jps = ([0, 1, 0, 1], ['-','+','-','+']))  # 1/2+ 0- 0- 1/2+
    dpp = randomPoint(tbs)

    # lineshape
    dc_Σb1 = decay_chain(1, (s,σ) -> BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')
    dc_Σb3 = decay_chain(3, (s,σ) -> BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')

    dc_Σb_x1 = decay_chain(1, (s,σ) -> BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')
    dc_Σb_x3 = decay_chain(3, (s,σ) -> BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')

    full_a(σs,two_λs) = sum(amplitude(σs,two_λs, ch) for ch in [dc_Σb1, dc_Σb_x1, dc_Σb3, dc_Σb_x3])
    full_a(dpp) = full_a(dpp.σs, dpp.two_λs)
    #
    total_I(σs) = sum(abs2(full_a(σs,two_λs)) for two_λs in Iterators.product(
        -tbs.two_js[1]:2:tbs.two_js[1],
        -tbs.two_js[2]:2:tbs.two_js[2],
        -tbs.two_js[3]:2:tbs.two_js[3],
        -tbs.two_js[4]:2:tbs.two_js[4]))

    test_I = summed_over_polarization((σs,two_λs)->abs2(full_a(σs,two_λs)),tbs.two_js)
    @test test_I(dpp.σs) == total_I(dpp.σs)
end
