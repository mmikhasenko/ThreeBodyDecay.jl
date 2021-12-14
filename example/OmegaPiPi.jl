using ThreeBodyDecay
using ThreeBodyDecay.PartialWaveFunctions
using Parameters

using Test
using Plots



struct minusone end
import Base:^
^(x::minusone, n::Int) = isodd(n) ? -1 : 1
macro x_str(s::String)
    minusone()
end


wignerd_sign(j,λ1,λ2, cosθ, ispositive) =
(ispositive ? 1 : x"-1"^*(λ1-λ2)) *
    wignerd(j,λ1,λ2, cosθ)

#
H(j1,λ1,j2,λ2,j,l,s) = 
    clebschgordan(j1,λ1,j2,-λ2,s,λ1-λ2) *
        clebschgordan(l,0,s,λ1-λ2,j,λ1-λ2)

# 
@with_kw struct Chain{T}
    k::Int
    lineshape::T
    jR::Int
    # 
    L::Int
    S::Int
    l::Int
    s::Int
end
# 

function pk(k,σs,msq)
    i,j,_ = ijk(k)
    sqrt(ThreeBodyDecay.λ(σs[k],msq[i],msq[j])/(4σs[k]))
end
qk(k,σs,msq) = sqrt(ThreeBodyDecay.λ(msq[4],σs[k],msq[k])/(4msq[4]))


function O(system, chain::Chain, λs, τ; refζs::Vector{Int}=[1,1,1,1])
    # system parameters
    @unpack js = system
    ms² = system.ms^2
    # kinematics
    @unpack σ1, σ2 = τ
    # chain parameters
    @unpack lineshape, k = chain
    @unpack jR = chain
    @unpack l,s,L,S = chain
    # 
    # evaluation
    i,j,_ = ijk(k)
    σs = Invariants(ms; σ1, σ2)
    #
    cosθ = cosθij(k, σs, ms²)
    p = pk(k,σs,ms²)
    q = qk(k,σs,ms²)
    #
    w0 = wr(k,refζs[4],0); cosζ0 = cosζ(w0, σs, ms²)
    wi = wr(k,refζs[i],i); cosζi = cosζ(wi, σs, ms²)
    wj = wr(k,refζs[j],j); cosζj = cosζ(wj, σs, ms²)
    wk = wr(k,refζs[k],k); cosζk = cosζ(wk, σs, ms²)
    #
    angular = sum(
        wignerd_sign(js[4],λs[4],ν-λs′[k], cosζ0, ispositive(w0)) *
        #
        q^L * H(jR,ν,js[k],λs′[k],js[4],L,S) * x"-1"^(js[k]-λs′[k]) *
        wignerd(jR, ν, λs′[i]-λs′[j], cosθ) *
        p^L * H(js[i],λs′[i],js[j],λs′[j],jR,l,s) * x"-1"^(js[j]-λs′[j]) *
        #
        wignerd_sign(js[i],λs′[i],λs[i], cosζi, ispositive(wi)) *
        wignerd_sign(js[j],λs′[j],λs[j], cosζj, ispositive(wj)) *
        wignerd_sign(js[k],λs′[k],λs[k], cosζk, ispositive(wk)) *
        1.0
        # 
                for ν in -jR:jR,
                    λs′ in Iterators.product(-js[1]:js[1],
                                             -js[2]:js[2],
                                             -js[3]:js[3]))
    # 
    bw = lineshape(σs[k])
    return angular*bw
end

# test
const mω  = 0.78266
const mπ⁻ = 0.13957039
const mπ⁰ = 0.1349768

struct ThreeBodySystem{T,K}
    masses::T
    spins::K
end

ms = ThreeBodyMasses(mπ⁻, mπ⁻, mπ⁰; m0=mω)
js = (j1=0,j2=0,j3=0,j0=1)
system = (;ms, js)

# 
ρBW(σ;p=(m=0.4,Γ=0.05)) = 1 / (p.m^2-σ-1im*p.m*p.Γ)

#
ρs = [Chain(
    k, ρBW, 1, # 
    1, 1, 1, 0) for k in 1:3]
# 
@test O(system, ρ1, (0,0,0,1), randomPoint(system.ms)) != 0.0im
@test O(system, ρ1, (0,0,0,0), randomPoint(system.ms)) == 0im
# 

itr()

# 
A3(σs) = sum(O.(Ref(system), ρs, Ref((0,0,0,1)), Ref(σs)))
I3(σs) = abs2(A3(σs))

plot(system.ms, I3)