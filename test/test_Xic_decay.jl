using ThreeBodyDecay
using Test

@testset "Xic2pKpi" begin

    mp = 0.938; mK = 0.49367; mpi = 0.13957; mXic = 2.46867
    tbs = ThreeBodySystem(mXic,mp,mK,mpi)

    amplitude(two_λ,two_Λ,σ3,σ1,Cls) = amp_b2bzz(two_λ,two_Λ,σ3,σ1,Cls,tbs)

    intensity(σ3,σ1,Cksl) = rate_b2bzz(σ3,σ1,Cksl,tbs)
    intensityCC(σ3,σ1,Csl1,Csl2) = rateCC_b2bzz(σ3,σ1,Csl1,Csl2)

    sample = flatDalitzPlotSample31(tbs)

    Cls = [
        (1,
        BreitWigner(0.77,0.15),
        [twochain(jls(1,2,1),jls(2,2,0))=>1.0im])]
    #
    rate = [intensity(σ3,σ1,Cls) for (σ3,σ1) in zip(sample[1],sample[2])]
    @test sum(rate .> 0) == length(sample[1])
    # return  [sample... rate]
end
