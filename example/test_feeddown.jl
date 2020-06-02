using ThreeBodyDecay
using Parameters
using Plots
using PartialWaveFunctions
theme(:wong)

urand(x1,x2) = x1+(x2-x1)*rand()
const mΞcc = 3.6212
const mΞc  = 2.46794
const mΞc′  = 2.5784
const mπ   = 0.14
const mπ0  = 0.135
const mγ   = 0.0
#
function four_body_decay(vars; ms)
    @unpack s,σ1,cosθ1,cosθ23,ϕ23 = vars
    @unpack m0,m1,m2,m3,m4 = ms
    p0 = [0,0,0,m0]
    p123, p4 = four_vectors_in_binary_decay(p0,1,0; m1sq=s, m2sq=m4^2)
    p23, p1 = four_vectors_in_binary_decay(p123,cosθ1,0; m1sq=σ1, m2sq=m1^2)
    p2,p3 = four_vectors_in_binary_decay(p23,cosθ23,cosθ23; m1sq=m2^2, m2sq=m3^2)
    return p1,p2,p3,p4
end
#
phsp(s,σ;ms) = sqrt(λ(s,σ,ms.m1^2)*λ(ms.m0^2,s,ms.m4^2)*λ(σ,ms.m2^2,ms.m3^3))/(σ*s)
#
function rand_Ξcpipi_phsp(; m_miss=mπ0)
    #
    ms = (m0 = mΞcc, m2 = mΞc, m3 = m_miss, m1 = mπ, m4 = mπ)
    #
    s = urand((ms.m1+ms.m2+ms.m3)^2,(ms.m0-ms.m4)^2)
    σ1 = urand((ms.m2+ms.m3)^2,(√s-ms.m1)^2)
    cosθ1,cosθ23,ϕ23 = urand(-1,1),urand(-1,1),urand(-π,π)
    vars = (s=s,σ1=σ1,cosθ1=cosθ1,cosθ23=cosθ23,ϕ23=ϕ23)
    #
    p1,p2,p3,p4 = four_body_decay(vars; ms=ms)
    (m=sqrt(invmasssq(p1+p2+p4)), w=phsp(s,σ1;ms=ms))
end

function rand_Ξcpipi_Ξc′(;mΞc′=mΞc′)
    ms = (m0 = mΞcc, m2 = mΞc, m3 = mγ, m1 = mπ, m4 = mπ)
    #
    s = urand((ms.m1+ms.m2+ms.m3)^2,(ms.m0-ms.m4)^2)
    σ1 = mΞc′^2
    √s-ms.m1 < √σ1 && return rand_Ξcpipi_Ξc′(;mΞc′=mΞc′)
    cosθ1,cosθ23,ϕ23 = urand(-1,1),urand(-1,1),urand(-π,π)
    vars = (s=s,σ1=σ1,cosθ1=cosθ1,cosθ23=cosθ23,ϕ23=ϕ23)
    #
    p1,p2,p3,p4 = four_body_decay(vars; ms=ms)
    (m=sqrt(invmasssq(p1+p2+p4)), w=phsp(s,σ1;ms=ms))
end
#
function rand_Ξcpipi_ω(;mω=0.78)
    ms = (m0 = mΞcc, m2 = mπ, m3 = mπ, m1 = mπ0, m4 = mΞc)
    #
    s = mω^2
    σ1 = urand((ms.m2+ms.m3)^2,(√s-ms.m1)^2)
    cosθ1,cosθ23,ϕ23 = urand(-1,1),urand(-1,1),urand(-π,π)
    vars = (s=s,σ1=σ1,cosθ1=cosθ1,cosθ23=cosθ23,ϕ23=ϕ23)
    #
    p1,p2,p3,p4 = four_body_decay(vars; ms=ms)
    (m=sqrt(invmasssq(p2+p3+p4)), w=phsp(s,σ1;ms=ms))
end

let Nev = 1e5, bins = range(2.7,3.7,length=200)
    plot(leg=:left, title="feed-down Ξcc")
    calv = [rand_Ξcpipi_Ξc′() for _ in 1:Nev]; w=getindex.(calv,:w)
    stephist!(getindex.(calv,:m), weights=w/sum(w)*Nev , bins=bins, α=0.7, lab="Ξc′ (Ξcγ)")
    calv = [rand_Ξcpipi_ω() for _ in 1:Nev]; w=getindex.(calv,:w)
    stephist!(getindex.(calv,:m), weights=w/sum(w)*Nev, bins=bins, α=0.7, lab="ω (π⁻π⁺π⁰)")
    calv = [rand_Ξcpipi_phsp(m_miss=mπ0) for _ in 1:Nev]; w=getindex.(calv,:w)
    stephist!(getindex.(calv,:m), weights=w/sum(w)*Nev, bins=bins, α=0.7, lab="∼ph.sp.(miss π0)")
    calv = [rand_Ξcpipi_phsp(m_miss=mγ) for _ in 1:Nev]; w=getindex.(calv,:w)
    stephist!(getindex.(calv,:m), weights=w/sum(w)*Nev, bins=bins, α=0.7, lab="∼ph.sp.(miss γ)")
    vline!([mΞcc], c=:black, lab="mΞcc", xlab="m(Ξcππ) (GeV)")
end
