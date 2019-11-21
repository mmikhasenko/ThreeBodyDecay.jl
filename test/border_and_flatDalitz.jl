using ThreeBodyDecay
using Test

@testset "Border and Flat Dalitz" begin
    tbs = ThreeBodySystem(1.1,1.3,1.5,6.0)
    fDP = flatDalitzPlotSample31(tbs)
    b = border31(tbs)
    # make better test!
    #   It would be great to check if all points in `dDP` are inside of `b`
    @test true
end
