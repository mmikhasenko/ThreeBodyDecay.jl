using ThreeBodyDecay

tbs = let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    ThreeBodySystem(m1,m2,m3,m0)
end

using Plots

let
    σ1v = LinRange(tbs.mthsq[1], tbs.sthsq[1],300)
    σ3m = [σ3of1(-1.0,σ,tbs.msq) for σ in σ1v]
    σ3p = [σ3of1( 1.0,σ,tbs.msq) for σ in σ1v]
    plot(σ1v, [σ3m σ3p], lab="")
end

let
    σ3v, σ1v = flatDalitzPlotSample31(tbs; Nev=1000000)
    histogram2d(σ1v, σ3v, lab="", bins=100)
end
