using Plots
using ThreeBodyDecay
using LaTeXStrings
using BenchmarkTools
using StaticArrays

pyplot()

#
ms = (Jψ = 3.09, p=0.938, Bs = 5.366) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.p,  ms.Bs;   # masses m1,m2,m3,m0
            two_jps=([   1,  1//2, 1//2,    0] .|> x2,  # twice spin
                     [ '-',   '+',  '-',  '+'])) # parities
#
# chain 2,3, i.e. (2+3): Λs with the lowest ls, LS
Pc4430_2  = decay_chain(2, (s,σ)->BW(σ, 4.330, 0.02); two_s = 1/2|>x2, parity='+', tbs=tbs)
Pc4430_3  = decay_chain(3, (s,σ)->BW(σ, 4.330, 0.02); two_s = 1/2|>x2, parity='-', tbs=tbs)
Pcs = (Pc4430_2, Pc4430_3)
#
# chain-1, i.e. (1+2): Pentaquarks with the lowest ls, LS
ggS = decay_chain(1, (s,σ)->1.0; two_s = 1|>x2, tbs=tbs, parity='-')
ggL = (ggS,)
#
#
A(σs,two_λs,cs) = sum(c*amplitude(σs,two_λs,dc) for (c, dc) in zip(cs, (Pcs..., ggL...)))
I(σs,cs) = sum(abs2(A(σs,two_λs,cs)) for two_λs in itr(tbs.two_js))
#
const cs0 = [1.0, -1.0, 20.0+0im]
I(σs) = I(σs,cs0) # set the couplings
#
dpp = randomPoint(tbs)
I(dpp.σs)
#
let
    plot(size=(500,450),
        xlab=L"m_{p J/\psi}^2\,(\mathrm{GeV})", ylab=L"m_{p J/\psi}^2\,(\mathrm{GeV})")
    #
    σ3v = range(tbs.mthsq[3], tbs.sthsq[3], length=100)
    σ2v = range(tbs.mthsq[2], tbs.sthsq[2], length=100)
    cal = [Kibble23(σ2,σ3,tbs.msq) < 0.0 ? I([gσ1(σ2,σ3,tbs.msq),σ2,σ3]) : NaN
        for (σ2, σ3) in Iterators.product(σ2v,σ3v)]
    heatmap!(σ2v, σ3v, cal, colorbar=false, c=:viridis, title="unpolarized intensity")
end

function Iz(σ1,z1,cs)
    ((σ1 == tbs.mthsq[1]) || (σ1 == tbs.mthsq[1]) || (z1==1) || (z1==-1)) && return 0.0
    σ3 = σ3of1(z1,σ1,tbs.msq)
    σs = [σ1,gσ2(σ3,σ1,tbs.msq),σ3]
    # @show σ1,z1, Kibble(σs,tbs.msq)
    Kibble(σs,tbs.msq) > 0.0 && return 0.0
    return I(σs,cs)*sqrt(λ(σ1,tbs.msq[2],tbs.msq[3])*λ(tbs.msq[4],σ1,tbs.msq[1]))/σ1
end

let
    plot(size=(500,450),
        xlab=L"m_{p \bar{p}}\,(\mathrm{GeV})", ylab=L"\cos\theta")
    σ1v = range(tbs.mthsq[1], tbs.sthsq[1], length=100)
    z1v = range(-1, 1, length=60)
    cal = [Iz(σ1,z1, cs0) for (z1, σ1) in Iterators.product(z1v,σ1v)]
    heatmap!(σ1v, z1v, cal, colorbar=false, c=:viridis, title="unpolarized intensity")
end

#
@btime Iz(tbs.mthsq[1]+rand()*(tbs.sthsq[1]-tbs.mthsq[1]),2*rand()-1,cs0)
@btime I(dpp.σs, cs0)

@code_warntype amplitude(dpp.σs, dpp.two_λs, ggS)
amplitude(dpp.σs, dpp.two_λs, ggS)

@code_warntype CG_doublearg(1,1,4,2,3,3)
@code_warntype CG_doublearg(1,1,4,2,3,3)

I(dpp.σs, cs0)
Juno.@profile
Juno.@profiler  I(dpp.σs, cs0)
@code_warntype amplitude(dpp.σs,dpp.two_λs,ggS)

@profiler Zksτ(ggS.k,ggS.two_s,0,[0,1,1,0],[0,1,1,0],dpp.σs,tbs)

#        _|                                _|    _|
#    _|_|_|    _|_|    _|_|_|      _|_|_|      _|_|_|_|  _|    _|
#  _|    _|  _|_|_|_|  _|    _|  _|_|      _|    _|      _|    _|
#  _|    _|  _|        _|    _|      _|_|  _|    _|      _|    _|
#    _|_|_|    _|_|_|  _|    _|  _|_|_|    _|      _|_|    _|_|_|
#                                                              _|
#                                                          _|_|

using Random; Random.seed!(1) # hide
using Statistics, Parameters # imported for our implementation
using DynamicHMC
using LogDensityProblems
using TransformVariables
using DynamicHMC.Diagnostics
using Parameters

struct DalitzPosterior{T} # contains the summary statistics
    cs::T
end
#
# calculate summary statistics from a data vector

# define a callable that unpacks parameters, and evaluates the log likelihood
function (problem::DalitzPosterior)(vars)
    @unpack σ1,z1 = vars
    @unpack cs = problem
    v = Iz(σ1,z1,cs)
    v < 0 && error("cannot be σ1,z1,cs = $((σ1,z1,cs))")
    loglikelihood = log(Iz(σ1,z1,cs))
    return loglikelihood
end

problem = DalitzPosterior(cs0)

# problem((xs = [1,3,5],))
#
trans = as((σ1=as(Real,tbs.mthsq[1],tbs.sthsq[1]), z1 = as(Real,-1,1)))
dimension(trans)
ℓ = TransformedLogDensity(trans, problem)
#
LogDensityProblems.dimension(ℓ)
LogDensityProblems.logdensity(ℓ, zeros(2))
#
∇P = ADgradient(:ForwardDiff, ℓ)

@time results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, 5000; reporter = NoProgressReport())
posterior = transform.(trans, results.chain)

let
    # plot(layout=grid(1,2), size=(1000,350))
    histogram2d(getindex.(posterior,1), getindex.(posterior,2), c=:viridis, bins=30)
end
histogram(getindex.(posterior,2), bins=30)

summarize_tree_statistics(results.tree_statistics)
