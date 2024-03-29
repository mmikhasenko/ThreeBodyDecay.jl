using ThreeBodyDecay
using Test
using Parameters

@testset "Infererred type" begin
    tbs = ThreeBodySystem(2.0, 1.0, 1.5; m0=6.0,
        two_js=ThreeBodySpins(1, 0, 0; two_h0=1))  # 1/2+ 0- 0- 1/2+

    dc = DecayChainLS(3, σ -> 1 / (4.1^2 - σ - 0.1im); two_s=3, parity='-', Ps=['+', '-', '-', '+'], tbs=tbs)

    @unpack σs, two_λs = randomPoint(tbs)

    @inferred Complex{Float64} amplitude(dc, σs, two_λs)

    # using InteractiveUtils
    # @code_warntype amplitude(dc, σs, two_λs)
end
