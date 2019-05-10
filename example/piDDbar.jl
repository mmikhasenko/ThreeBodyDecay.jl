using ThreeBodyDecay

const mDstar = 2.00685 # 2 GeV
const ΓDstar = 30e-6 # 30 keV
tbs = let mD = 1.86483, mπ0 = 0.13498
    ThreeBodySystem(mD+mDstar+0.01,mD,mπ0,mD)
end

using Plots

# amplitude

BWD(σ) = BW(σ,mDstar,100*ΓDstar)
A(σ3,σ1) = BWD(σ1) + BWD(σ3)
A(σ3,σ1,ν) = BWD(σ1) * wignerd(1, ν, 0, cosθ23(gσ2(σ3,σ1,tbs),σ3,tbs)) +
             BWD(σ3) * sum(wignerd(1, ν, λ, cos_plus_θhat3(σ1,gσ2(σ3,σ1,tbs),tbs)) *
                           wignerd(1, λ, 0, cosθ12(σ1,gσ2(σ3,σ1,tbs),tbs)) for λ=-1:1)

# intensity

I(σ3,σ1) = sum(abs2(A(σ3,σ1,ν)) for ν=-1:1)

let
    σ1v = LinRange(tbs.mthsq[1], tbs.sthsq[1],200)
    σ3v = LinRange(tbs.mthsq[3], tbs.sthsq[3],200)
    cal = [Kibble31(σ3,σ1,tbs) < 0.0 ? I(σ3,σ1) : NaN for σ3 in σ3v, σ1 in σ1v]
    heatmap(σ1v, σ3v, cal)
end

# sample density
density2d = getbinned2dDensity(
                (σ3,σ1) -> Kibble31(σ3,σ1,tbs) < 0.0 ? I(σ3,σ1) : 0.0,
                    (tbs.mthsq[1], tbs.sthsq[1]),
                    (tbs.mthsq[3], tbs.sthsq[3]), 100, 100);

let r = hcat([rand(density2d) for _ in 1:100000]...)
    histogram2d(r[2,:], r[1,:], bins=100)
end
