using ThreeBodyDecay
using Test

@testset "Xic2pKpi" begin

    mp = 0.938; mK = 0.49367; mpi  = 0.13957; mXic = 2.46867
    tbs = ThreeBodySystem(mp,mK,mpi,mXic)

    CS = [
        (1,
        BreitWigner(0.77,0.15),
        [twochain(jls(1,2,1),jls(2,2,0))])]

    amplitude(two_λ,two_Λ,dpp,Cs) = amp_b2bzz(two_λ,two_Λ,dpp,CS,tbs,Cs)
    intensity(σ3,σ1,Cs) = rate_b2bzz(DalitzPlotPoint31(σ3,σ1,tbs),CS,tbs,Cs)

    sample = flatDalitzPlotSample31(tbs)

    rate = [intensity(σ3,σ1,[reshape([1,2],1,2)]) for (σ3,σ1) in zip(sample[1],sample[2])]
    @test sum(rate .> 0) == length(sample[1])
    # return  [sample... rate]
end
