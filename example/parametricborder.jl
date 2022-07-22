using ThreeBodyDecay
using Plots
using Polynomials

theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false)
# 

function dalitzcenter(ms)
    thresholds = ((ms.m2 + ms.m3), (ms.m3 + ms.m1), (ms.m1 + ms.m2), 0)
    T0 = sum(abs2, ms) .- sum(abs2, thresholds)
    return T0 / 3 .+ thresholds .^ 2
end

function σs2rθ(σs, ms)

    Δσs′ = Tuple(σs) .- dalitzcenter(ms)[1:3]
    # 
    rcosθ = -Δσs′[1]
    rcosθmπ3 = Δσs′[3]
    rsinθ = 2 / √3 * (rcosθmπ3 - rcosθ / 2)
    # 
    θ = atan(rsinθ, rcosθ)
    r = rcosθ / cos(θ)
    # 
    return (; r, θ)
end

function rlimsdalitz(θ, ms)
    ϕ(σs) = Kibble(σs, ms^2)
    T0 = 3 * dalitzcenter(ms)[4]
    σs = thresholds .^ 2 .+ ThreeBodyDecay.polardalitz2invariants(θ, T0)
    rmax = minimum(filter(x -> x > 0, roots(ϕ(σs))))
    return rmax
end


struct PolarDalitz
    ms::ThreeBodyDecay.MassTuple
end

@recipe function f(pd::PolarDalitz)
    ms = pd.ms
    σsv = border(ms)
    rθv = σs2rθ.(σsv, Ref(ms))
    xv = getproperty.(rθv, :θ)
    yv = getproperty.(rθv, :r)
    projection := :polar
    limits --> (0, :auto)
    # 
    (xv, yv)
end











ms = ThreeBodyMasses(1, 2, 4; m0=13)

σsv = border(ms)
rθv = σs2rθ.(σsv, Ref(ms))
# 
getindex(rθv)
plot(NamedTuple{(:σ1, :σ2)}.(σsv))

let
    θv = range(-π, π, length=200)
    rv = range(0, 1, length=101)
    zv = [rlimsdalitz(θ, ms) for r in rv, θ in θv]
    Plots.heatmap(θv, rv, zv, proj=:polar, title="phase space intensity")
end




let
    ms0 = PolarDalitz(ThreeBodyMasses(1, 2, 4; m0=13))
    plot()
    for α in 1.0:-0.05:0.2
        plot!(PolarDalitz(ThreeBodyMasses(1α, 2α, 4α; m0=13α)), lab="$α",
            legtitle="alpha")
    end
    plot!()
end


let
    ms0 = ThreeBodyMasses(1, 2, 4; m0=13)
    # 
    plot()
    for α in 1.0:-0.02:0.1
        ms = ThreeBodyMasses(1α, 2α, 4α; m0=13α)
        data = NamedTuple{(:σ1, :σ2)}.(border(ms))
        plot!(
            getindex.(data, 1) .+ (dalitzcenter(ms0)[1] - dalitzcenter(ms)[1]),
            getindex.(data, 2) .+ (dalitzcenter(ms0)[2] - dalitzcenter(ms)[2]))
    end
    plot!()
end

let
    ms = ThreeBodyMasses(1, 2, 4; m0=13)
    rθv = σs2rθ.(σsv, Ref(ms))
    plot(
        getproperty.(rθv, :θ),
        getproperty.(rθv, :r), proj=:polar, lims=(0, :auto))
end





function ϕ(α, σs, ms)
    rθ = σs2rθ(σs, ms)
    # 
    !(0 ≤ α ≤ 1) && return NaN
    msα = ThreeBodyMasses(α * ms.m1, α * ms.m2, α * ms.m3; m0=α * ms.m0)
    # 
    σsα = dalitzcenter(msα)[1:3] .+
          map(ThreeBodyDecay.polardalitz2invariants(rθ.θ, 0)) do P
        P(rθ.r)
    end
    return Kibble(σsα, msα^2) / α^2
end

@assert ϕ(1, σsv[1], ms) == Kibble(σsv[1], ms^2)

using NLsolve

function tounitcircle(σs, ms)
    solv = [
        nlsolve(x -> ϕ(atan(x[1]) / π + 0.5, σs, ms),
            [tan(π * (α0 - 0.5))]).zero[1]
        for α0 in 0.1:0.1:0.9]
    sol = maximum(filter(solv) do s
        !isnan(s) && abs(ϕ(atan(s) / π + 0.5, σs, ms)) < 0.0001
    end)
    return atan(sol) / π + 0.5
end
# 
let
    plot()
    for σs in dv[rand(1:100, 3)]
        plot!(α -> ϕ(α, σs, ms), 0, 1)
    end
    hline!([0.0])
    plot!()
end


tounitcircle(dv[2], ms)

dv = flatDalitzPlotSample(ms; Nev=10000)
@time αv = tounitcircle.(dv, Ref(ms))
# histogram(αv, bins=100)
rθdv = σs2rθ.(dv, Ref(ms))


scatter(
    getproperty.(rθdv, :θ),
    αv .^ 2, proj=:polar, lims=(0, :auto), c=2, ms=0.3)

# check
# let
#     thresholds = ((ms.m2 + ms.m3), (ms.m3 + ms.m1), (ms.m1 + ms.m2))
#     T0 = sum(abs2, ms) .- sum(abs2, thresholds)

#     plot(NamedTuple{(:σ1, :σ2)}.(
#         ThreeBodyDecay.MandestamTuple.(
#             map(ThreeBodyDecay.polardalitz2invariants(rθ.θ, T0)) do P
#                 P(rθ.r)
#             end .+ thresholds .^ 2 for rθ in rθv)
#     ))
#     plot!(border12(ms))
# end


