using ThreeBodyDecay
using Test

@testset "Border and Flat Dalitz" begin
    ms = ThreeBodyMasses(6.0, 1.1, 1.3, 1.5)
    fDP = flatDalitzPlotSample(ms)
    b = border31(ms)
    # make better test!
    #   It would be great to check if all points in `dDP` are inside of `b`
    @test true
end
