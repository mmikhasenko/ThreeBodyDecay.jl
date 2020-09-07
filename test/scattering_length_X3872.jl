using ThreeBodyDecay
using Test
# using Plots

@testset "Scattering Length Approximation" begin
    mJψ = 3.0969;  mππ = 0.3;
    mD0 = 1.86483; mDstar0 = 2.00685

    cs = (c11 = -0.03, c12 = 0.07, c22=0.0)
    #
    ScA = ScattLenApproximation(
        Ms=(cs.c11,cs.c12,cs.c22),
        ms1=(mD0, mDstar0), ms2=(mJψ,mππ))
    #
    ev = range(-15,10, length=300) # MeV
    sv = ((mD0+mDstar0) .+ ev .* 1e-3).^2 # GeV
    #
    calv = map(x->amp(x, ScA), sv.+1e-7im) # amplitude
    v,i = findmax(imag.(calv))
    #
    # plot(sqrt.(sv), [real.(calv) imag.(calv)])
    # vline!([1.86483+2.00685])
    #
    # peaking structure
    @test abs(-1.705-ev[i]) < 0.001
end
