using ThreeBodyDecay
using Plots
using Zygote


tbs = let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    ThreeBodySystem(m1,m2,m3,m0)
end

#
# let
#     σ1v = LinRange(tbs.mthsq[1], tbs.sthsq[1],300)
#     σ3m = [σ3of1(σ,-1.0,tbs) for σ in σ1v]
#     σ3p = [σ3of1(σ, 1.0,tbs) for σ in σ1v]
#     plot(σ1v, [σ3m σ3p], lab="")
# end
#
# let
#     σ3v, σ1v = flatDalitzPlotSample31(tbs; Nev=10000)
#     histogram2d(σ1v, σ3v, lab="", bins=100)
# end

const Ks892  = BreitWigner(0.89176, 0.05)
const Λ1520  = BreitWigner(1.5195,  0.0156)
const Δ1232  = BreitWigner(1.232,   0.112)

model = [(1,Ks892),
         (3,Λ1520),
         (2,Δ1232)]

function scalar_amplitude(σs, model, cs)
    mod = [1.1,2.2im, 3.3]
    #
    A = 0.0im
    for (me,c) in zip(model,cs)
        (ch,ξ) = me
    # for (m,c) in zip(mod,cs)
        # A += c * amp(σs[ch],ξ)
        A += c * mod[ch]
    end
    return A
end

const cs0 = [1.2im, 2.1, 2.3+1im]

let σ1 = 1.0, σ3=4.0
    scalar_amplitude([σ1,gσ2(σ3,σ1,tbs),σ3], model, cs0)
end

scalar_amplitude_squared(σs, model, cs) = abs2(scalar_amplitude(σs, model, cs))

# @time let
#     σ3v, σ1v = flatDalitzPlotSample31(tbs; Nev=100000)
#     weights = [scalar_amplitude_squared([σ1,gσ2(σ3,σ1,tbs),σ3], model, cs0) for (σ3, σ1) in zip(σ3v, σ1v)]
#     histogram2d(σ1v, σ3v, lab="", bins=100, weights=weights)
# end

# getbinned2dDensity(g, xlim, ylim,  Nrows, Ncols)

diff_through_pars(cs) = let σ1 = 1.0, σ3=4.0,
    σs = [σ1,gσ2(σ3,σ1,tbs),σ3]
    return gradient(x->scalar_amplitude_squared(σs, model, x), cs)
end

diff_through_pars([1.1im,1.1,1.1])

let
end


# function amp_b2bzz(two_λ,two_Λ,σs,CS,tbs,Cs)
#     σ2 = gσ2(σ3,σ1,tbs)
#     # Wigner rotations
#     D1 = [two_λ==two_λp ? 1.0 : 0.0 for two_λp=-1:2:1]
#     D2 = [((two_λp-two_λ) % 4 == 2 ? -1 : 1) *
#         wignerd_doublearg(1,two_λp,two_λ,cosζ12_for1(σ1,σ2,tbs))
#             for two_λp=-1:2:1]
#     D3 = [wignerd_doublearg(1,two_λp,two_λ,cosζ31_for1(σ3,σ1,tbs))
#             for two_λp=-1:2:1]
#     #
#     σs = (σ1,σ2,σ3)
#     Ds = (D1,D2,D3)
#     #
#     val = zero(Cs[1][1,1]+0.0im)
#     for (cs,W) in zip(CS, Cs) # coupling scheme and couplings
#         (k,ξ,chains) = cs
#         length(chains) == 0 && continue
#         decay_amp = zero(Cs[1][1,1]+0.0im)
#         for two_τ in -two_j(chains[1]):2:two_j(chains[1]), two_λp = -1:2:1
#             Wc = [v1+v2*1im for (v1,v2) in zip(W[:,1],W[:,1])]
#             decay_amp += RkΛsτ(k,two_λp,two_Λ,two_τ,chains, Wc) *
#                     ZkΛsτ(k,two_λp,two_Λ,two_j(chains[1]),two_τ,σ3,σ1,tbs) *
#                     Ds[k][div(two_λp+3,2)]
#         end
#         val += decay_amp * amp(σs[k],ξ)
#     end
#     return val;
# end
