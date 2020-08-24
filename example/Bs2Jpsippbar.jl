using Plots
using ThreeBodyDecay
using LaTeXStrings
pyplot()

#
ms = (Jψ = 3.09, p=0.938, Bs = 5.366) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.p, m0=ms.Bs;   # masses m1,m2,m3,m0
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
# dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
# dpp123 = [dpp.σs[1],dpp.σs[2],dpp.σs[3]]
# dpp132 = [dpp.σs[1],dpp.σs[3],dpp.σs[2]]
#
# I(dpp.σs)
# I(dpp.σs)
#
# I(dpp123), I(dpp132)
#
# sum(abs2(amplitude(dpp123, two_λs, ggS)) for two_λs in itr(tbs.two_js))
# sum(abs2(amplitude(dpp123, two_λs, ggS)) for two_λs in itr(tbs.two_js))
#
let
    plot(size=(500,450),
        xlab=L"m_{p J/\psi}^2\,(\mathrm{GeV})", ylab=L"m_{p J/\psi}^2\,(\mathrm{GeV})")
    #
    σ3v = range(tbs.mthsq[3], tbs.sthsq[3], length=100)
    σ2v = range(tbs.mthsq[2], tbs.sthsq[2], length=100)
    cal = [inphrange((gσ1(σ2,σ3,tbs.msq),σ2,σ3), tbs) ? I([gσ1(σ2,σ3,tbs.msq),σ2,σ3]) : NaN
        for (σ2, σ3) in Iterators.product(σ2v,σ3v)]
    heatmap!(σ2v, σ3v, cal - cal', colorbar=false, c=:viridis, title="unpolarized intensity")
end

let
end



# function diff_I(ϕ; ϕP=0)
#     cs = [2.0, 2.0*cis(ϕP),20*cis(ϕ)]
#     I(dpp123,cs) - I(dpp132,cs)
# end
#
# let
#     plot(ϕ->diff_I(ϕ,ϕP=0), -π, π)
#     # plot(ϕP->diff_I(0.0,ϕP=ϕP), -π, π)
#     vline!([-π/2,π/2,π/4])
# end

# Zksτ(1,2,0, [0,1//2,1//2,0].|>x2, dpp123, tbs)
# Zksτ(1,2,0, [0,1//2,1//2,0].|>x2, dpp132, tbs)
#
# cosθhatk1(1,dpp123,tbs.msq), cosθhatk1(1,dpp132,tbs.msq)
# cosθij(1,dpp123,tbs.msq), cosθij(1,dpp132,tbs.msq)
#
# wignerd_doublearg(2,0,0,0.3)
# wignerd_doublearg(2,0,0,-0.3)
#
#
# cosζk2_for2(1,dpp123,tbs.msq), cosζk2_for2(1,dpp132,tbs.msq)
# cosζk3_for3(1,dpp123,tbs.msq), cosζk3_for3(1,dpp132,tbs.msq)
# phase(2,1,1,-1), phase(3,1,1,-1)
# #
# let k=1,two_s=2,two_τ=0,two_λs=[0,1//2,1//2,0].|>x2,σs=dpp123
#     (i,j) = (2,3)
#     #
#     two_Λ = 0
#     #
#     @show cosθij(1,σs,tbs.msq)
#     #
#     val = 0.0;
#     for two_λs_prime in Iterators.product(0:2:0, -1:2:1, -1:2:1)
#         #
#         val +=
#             wignerd_doublearg(two_s,two_τ,two_λs_prime[2]-two_λs_prime[3],cosθij(1,σs,tbs.msq)) *
#             phase(2,1,two_λs_prime[2],two_λs[2]) *
#             wignerd_doublearg(tbs.two_js[2], two_λs_prime[2], two_λs[2], cosζk2_for2(1,σs,tbs.msq)) *
#             wignerd_doublearg(tbs.two_js[3], two_λs_prime[3], two_λs[3], cosζk3_for3(1,σs,tbs.msq))
#     end
#     val
# end
# phase(2,1,1,-1)
# phase(1,3,1,-1)
#
#
# let k=1,two_s=2,two_τ=0,two_λs=[0,1//2,1//2,0].|>x2,σs=dpp132
#     (i,j) = (2,3)
#     #
#     two_Λ = 0
#     #
#     @show cosθij(1,σs,tbs.msq)
#     #
#     val = 0.0;
#     for two_λs_prime in Iterators.product(0:2:0, -1:2:1, -1:2:1)
#         #
#         val +=
#             wignerd_doublearg(two_s,two_τ,two_λs_prime[2]-two_λs_prime[3],cosθij(1,σs,tbs.msq)) *
#             phase(2,1,two_λs_prime[2],two_λs[2]) *
#             wignerd_doublearg(tbs.two_js[2], two_λs_prime[2], two_λs[2], cosζk2_for2(1,σs,tbs.msq)) *
#             wignerd_doublearg(tbs.two_js[3], two_λs_prime[3], two_λs[3], cosζk3_for3(1,σs,tbs.msq))
#     end
#     val
# end


# phase(k,1,two_Λ,two_τ-two_λs_prime[k]) *
# sqrt(two_J(tbs)+1) * wignerd_doublearg(two_J(tbs),0,two_τ-two_λs_prime[k],cosθhatk1(k,σs,tbs.msq)) *
# sqrt(two_s+1) * wignerd_doublearg(two_s,two_τ,two_λs_prime[i]-two_λs_prime[j],cosθij(k,σs,tbs.msq)) *
#     phase(1,k,two_λs_prime[1],two_λs[1]) * wignerd_doublearg(tbs.two_js[1], two_λs_prime[1], two_λs[1], cosζk1_for1(k,σs,tbs.msq)) *
#     phase(2,k,two_λs_prime[2],two_λs[2]) * wignerd_doublearg(tbs.two_js[2], two_λs_prime[2], two_λs[2], cosζk2_for2(k,σs,tbs.msq)) *
#     phase(3,k,two_λs_prime[3],two_λs[3]) * wignerd_doublearg(tbs.two_js[3], two_λs_prime[3], two_λs[3], cosζk3_for3(k,σs,tbs.msq))
#
