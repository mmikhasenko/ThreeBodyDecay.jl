using ThreeBodyDecay

const mΛb = 5.61960
const mπ = 0.14
# intermediate resonances
const mΣb = 5.81065; const ΓΣb = 0.005
const mΣb_x = 5.83032; const ΓΣb_x = 0.0094

const mΛb2S = 6.3; # just a peak of the plot

tbs = ThreeBodySystem(mπ,mΛb,mπ; m0=mΛb2S,
    two_jps = ([0, 1, 0, 1], ['-','+','-','+']))  # 1/2+ 0- 0- 1/2+
dpp = randomPoint(tbs)

# lineshape
dc_Σb1 = DecayChainLS(1, (s,σ) -> BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')
dc_Σb3 = DecayChainLS(3, (s,σ) -> BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')

dc_Σb_x1 = DecayChainLS(1, (s,σ) -> BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')
dc_Σb_x3 = DecayChainLS(3, (s,σ) -> BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')

dc_f0 = DecayChainLS(2, (s,σ) -> 1.0; tbs=tbs, two_s = 0, parity='+')

full_a(σs,two_λs) = sum(amplitude(σs,two_λs, ch) for ch in [dc_Σb1, dc_Σb_x1, dc_Σb3, dc_Σb_x3])
full_a(dpp) = full_a(dpp.σs, dpp.two_λs)
#
total_I(σs) = sum(abs2(full_a(σs,two_λs)) for two_λs in Iterators.product(
    -tbs.two_js[1]:2:tbs.two_js[1],
    -tbs.two_js[2]:2:tbs.two_js[2],
    -tbs.two_js[3]:2:tbs.two_js[3],
    -tbs.two_js[4]:2:tbs.two_js[4]))

σ3v,σ1v = flatDalitzPlotSample31(tbs)
σ2v = [gσ2(σ3,σ1,tbs.msq) for (σ3,σ1) in zip(σ3v,σ1v)]

cala = [full_a(DalitzPlotPoint31(σ3,σ1,tbs)) for (σ3,σ1) in zip(σ3v,σ1v)]
cali = abs2.(cala)

using Plots
histogram2d(σ3v, σ2v, bins=50)
histogram2d(σ3v, σ2v, weights=cali, bins=50)

let
    σ3v = range(tbs.mthsq[3],tbs.sthsq[3], length=150)
    σ2v = range(tbs.mthsq[2],tbs.sthsq[2], length=150)
    cali = [let dpp = DalitzPlotPoint23(σ2,σ3,tbs)
        Kibble(dpp.σs, tbs.msq) < 0 ? total_I(dpp.σs) : NaN
        end for (σ3,σ2) in Iterators.product(σ3v,σ2v)]
    heatmap(σ2v,σ3v,cali)
end
