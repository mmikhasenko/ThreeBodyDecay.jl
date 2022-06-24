using ThreeBodyDecay
using Plots
using StaticArrays
using Parameters
using BenchmarkTools

@with_kw struct BW
    m::Float64
    Γ::Float64
end

(bw::BW)(σ::Float64) = 1/(bw.m^2-σ-1im*bw.m*bw.Γ)

const mΛb = 5.62
const mJψ = 3.09
const mp=0.938
const mK = 0.49367
# 
const ms = ThreeBodyMasses(mJψ,mp,mK,m0=mΛb)
const tbs = ThreeBodySystem(;ms, two_js=(2,1,0,1))
#
const Ps_pc = SVector('-','+','-','+');
const Ps_pv = SVector('-','+','-','-');

isobars = [
    (key = "Λ1520",  JP="3/2-", lineshape=:BW, m=1.5195,  Γ=0.016),
    (key = "Λ1600",  JP="1/2+", lineshape=:BW, m=1.630,   Γ=0.250),
    (key = "Λ1670",  JP="1/2-", lineshape=:BW, m=1.685,   Γ=0.050),
    (key = "Pc4312", JP="1/2-", lineshape=:BW, m=4.312,  Γ=0.012),
    (key = "Pc4440", JP="1/2-", lineshape=:BW, m=4.440,  Γ=0.008),
    (key = "Pc4457", JP="3/2-", lineshape=:BW, m=4.457,  Γ=0.002)    
]

chains = []
for ξ in isobars
    # 
    @unpack key, JP, lineshape, m, Γ = ξ
    # 
    qn = str2jp(JP)
    k = key[1] == "Λ" ? 1 : 3
    ξf = eval(quote
        $(lineshape)(m=$(m), Γ=$Γ)
    end)
    # 
    dc_pc = DecayChainsLS(k, ξf;
        tbs,
        two_s=ThreeBodyDecay.two_j(qn),
        parity=qn.p,
        Ps=Ps_pc)
    push!(chains, dc_pc...)
    #
    dc_pv = DecayChainsLS(k, ξf;
        tbs,
        two_s=ThreeBodyDecay.two_j(qn),
        parity=qn.p,
        Ps=Ps_pv)
    push!(chains, dc_pv...)
    #
end

cs = rand(length(chains))

I(σsλ, model) = sum(abs2, model.cs .* amplitude.(Ref(σsλ), model.chains))

const σsλ0 = randomPoint(tbs)
# amplitude.(Ref(σsλ0), chains)

const model = (; cs, chains)

I(σsλ0, model)

function timeonN(Nev)
    ph = flatDalitzPlotSample(ms; Nev)
    phλ = DalitzPlotPoint.(ph, Ref(σsλ0.two_λs))
    # 
    return @belapsed I.($(phλ), $(Ref(model)))
end

timeonN(1)
timeonN(10)
timeonN(100)
timeonN(1000)

