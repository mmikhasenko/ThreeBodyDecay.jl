using ThreeBodyDecay
using Plots

const B2psippbar = let mB = 5.27932
    ThreeBodySystem(0.938,3.09,0.938,mB)
end

σ3MC, σ1MC = flatDalitzPlotSample31(B2psippbar)

let
    histogram2d(σ3MC, σ1MC, bins=50)
    hline!([4.312].^2, lab="")
    vline!([4.312].^2, lab="")
end
