using ThreeBodyDecay
using Plots
using StaticArrays
using Parameters
using BenchmarkTools
using PrettyTables

@with_kw struct BW
    m::Float64
    Γ::Float64
end

(bw::BW)(σ::Float64) = 1 / (bw.m^2 - σ - 1im * bw.m * bw.Γ)

const mΛb = 5.62
const mJψ = 3.09
const mp = 0.938
const mK = 0.49367
# 
const ms = ThreeBodyMasses(mJψ, mp, mK, m0=mΛb)
const tbs = ThreeBodySystem(; ms, two_js=(2, 1, 0, 1))
#
const Ps_pc = SVector('-', '+', '-', '+');
const Ps_pv = SVector('-', '+', '-', '-');

isobars = [
    (key="Λ1520", JP="3/2-", lineshape=:BW, m=1.5195, Γ=0.016),
    (key="Λ1600", JP="1/2+", lineshape=:BW, m=1.630, Γ=0.250),
    (key="Λ1670", JP="1/2-", lineshape=:BW, m=1.685, Γ=0.050),
    (key="Pc4312", JP="1/2-", lineshape=:BW, m=4.312, Γ=0.012),
    (key="Pc4440", JP="1/2-", lineshape=:BW, m=4.440, Γ=0.008),
    (key="Pc4457", JP="3/2-", lineshape=:BW, m=4.457, Γ=0.002)
]

chains = []
for ξ in isobars
    # 
    @unpack key, JP, lineshape, m, Γ = ξ
    # 
    qn = str2jp(JP)
    k = key[1] == 'Λ' ? 1 : 3
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
chains = vcat(chains...)
cs = rand(length(chains))

function I(model, σsλ)
    v = 0.0im
    for (c, ch) in zip(model.cs, model.chains)
        v += c * amplitude(ch, σsλ)
    end
    return abs2(v)
end

const σsλ0 = randomPoint(tbs)

const model = (; cs, chains)




function summarystr(dc::DecayChain{BW})
    @unpack m, Γ = dc.Xlineshape
    sm, sΓ = string.(round.(Int, (m, Γ) .* 1e3))
    (kξ, kB) = dc.k == 1 ? ("Λ", "J/ψ") : ("Pc", "K")
    two_ls = dc.HRk.two_ls
    return "[ $(kξ)($(sm)) $(kB) ]_{$(two_ls[2])/2} $(div(two_ls[1],2))-wave"
end

belapsedbychain = [@belapsed amplitude($(ch), $(σsλ0)) for ch in model.chains]
extrema(belapsedbychain .* 1e6)

pretty_table(hcat(summarystr.(model.chains), belapsedbychain .* 1e6),
    header=(["wave name", "@time amplitude(wave)"],
        ["", "[μs]"]), show_omitted_cell_summary=false)
# 
# ┌──────────────────────────────┬───────────────────────┐
# │                    wave name │ @time amplitude(wave) │
# │                              │                  [μs] │
# ├──────────────────────────────┼───────────────────────┤
# │ [ Λ(1520) J/ψ ]_{1/2} 0-wave │                  13.9 │
# │ [ Λ(1520) J/ψ ]_{3/2} 2-wave │                  14.3 │
# │ [ Λ(1520) J/ψ ]_{5/2} 2-wave │                  14.4 │
# │ [ Λ(1520) J/ψ ]_{1/2} 1-wave │                  14.0 │
# │ [ Λ(1520) J/ψ ]_{3/2} 1-wave │                  14.2 │
# │ [ Λ(1520) J/ψ ]_{5/2} 3-wave │                  14.4 │
# │ [ Λ(1630) J/ψ ]_{1/2} 1-wave │                  7.05 │
# │              ⋮               │            ⋮           │
# └──────────────────────────────┴───────────────────────┘
#                                          22 rows omitted



function timeonN(Nev)
    ph = flatDalitzPlotSample(ms; Nev)
    phλ = DalitzPlotPoint.(ph, Ref(σsλ0.two_λs))
    # 
    return @belapsed I.($(Ref(model)), $(phλ))
end


Nev = [1, 10, 100, 1000]
belapsed = timeonN.(Nev)

t = pretty_table(hcat(Nev, belapsed),
    header=(["N", "Time"], ["[ev]", "[s]"]))
# 
# I.(phλ, Ref(model))
# ┌────────┬───────────┐
# │      N │      Time │
# │   [ev] │       [s] │
# ├────────┼───────────┤
# │    1.0 │ 0.0003994 │
# │   10.0 │ 0.0040236 │
# │  100.0 │ 0.0403715 │
# │ 1000.0 │  0.409975 │
# └────────┴───────────┘

# compilation
@profview I(σsλ0, model)
# pure runtime
@profview I(σsλ0, model)

using InteractiveUtils

I(σsλ0, model)
@code_warntype I(σsλ0, model)