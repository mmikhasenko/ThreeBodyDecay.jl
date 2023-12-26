using ThreeBodyDecay
using Plots

theme(:wong)

tbs = let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    ThreeBodySystem(m1, m2, m3, m0=m0)
end

let
    data = flatDalitzPlotSample(tbs.ms; Nev=1_000_000)
    X, Y = getproperty.(data, :σ1), getproperty.(data, :σ3)
    histogram2d(X, Y, lab="", bins=200)
    plot!(border13(tbs.ms), lab="", lw=5, lc=2)
end
