using ThreeBodyDecay
using Parameters

tbs = ThreeBodySystem(1.0,1.5,2.0; m0=6.0, two_js=ThreeBodySpins(0, 1, 0; two_h0=1))
ch = decay_chain(1,(s,σ)->1.0; two_s=1, tbs=tbs, parity='-', Ps=['-','+','-','+'])

let 
    @unpack σs,two_λs = randomPoint(tbs)
    two_λs′ = two_λs
    two_τ = 1
    # 
    @unpack k,two_s = ch
    # 
    Zksτ(k,two_s,two_τ,two_λs,two_λs′,σs,tbs), 
        ThreeBodyDecay.Zksτ_old(k,two_s,two_τ,two_λs,two_λs′,σs,tbs)
end

