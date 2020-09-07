using ThreeBodyDecay
using Test

@testset "Xic2pKpi" begin
    mp = 0.938; mK = 0.49367; mpi  = 0.13957; mXic = 2.46867
    tbs = ThreeBodySystem(mp,mK,mpi; m0=mXic)
    sample = flatDalitzPlotSample(tbs.ms)
    @test true
end
