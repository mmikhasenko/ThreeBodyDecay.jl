using ThreeBodyDecay
using Test

@testset "Border and Flat Dalitz" begin
    ms = ThreeBodyMasses(1.1, 1.3, 1.5; m0=6.0)
    fDP = flatDalitzPlotSample(ms)
    b = border(ms)
    borderKibblemax = maximum(abs.(Kibble.(b, Ref(ms^2))))
    phspKibblemax = maximum(abs.(Kibble.(fDP, Ref(ms^2))))
    # 
    @test borderKibblemax / phspKibblemax < 1e-7
end
