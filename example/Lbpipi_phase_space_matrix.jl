using ThreeBodyDecay

const mΛb = 5.61960
const mπ = 0.14
# intermediate resonances
const mΣb = 5.81065; const ΓΣb = 0.005
const mΣb_x = 5.83032; const ΓΣb_x = 0.0094

const mΛb2S = 6.3; # just a peak of the plot

tbs = ThreeBodySystem(mπ,mΛb,mπ; m0=mΛb2S,
    two_jps = ([0, 1, 0, 1], ['-','+','-','+']))  # 1/2+ 0- 0- 1/2+
dpp = randomPoint(tbs)

# lineshape
dc_Σb1 = DecayChainLS(1, σ->BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')
dc_Σb3 = DecayChainLS(3, σ->BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')
full_AΣb(σs,two_λs)   = sum(amplitude(σs, two_λs, ch) for ch in [dc_Σb1, dc_Σb3])

dc_Σb_x1 = DecayChainLS(1, σ->BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')
dc_Σb_x3 = DecayChainLS(3, σ->BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')
full_AΣb_x(σs,two_λs) = sum(amplitude(σs, two_λs, ch) for ch in [dc_Σb_x1, dc_Σb_x3])

dc_f0 = DecayChainLS(2, σ->1.0; tbs=tbs, two_s = 0, parity='+')
full_Af0(σs,two_λs)   = amplitude(σs, two_λs, dc_f0)

####

Φij = let
    As = [full_AΣb, full_AΣb_x, full_Af0]
    [let
        Aij(σs,two_λs)=Ai(σs,two_λs) * conj(Aj(σs,two_λs))
        I = summed_over_polarization(Aij, tbs.two_js)
        three_body_phase_space_integral(I, tbs)
        end for (Ai,Aj) in Iterators.product(As,As)]
end

using LinearAlgebra
Φij_n = let
    d = sqrt.(diag(Φij))
    Φij ./ d ./ d'
end

###
one_AΣb(σs,two_λs)   = amplitude(σs, two_λs, dc_Σb1)
one_AΣb_x(σs,two_λs) = amplitude(σs, two_λs, dc_Σb_x1)

Φ0ij = let
    As = [one_AΣb, one_AΣb_x, full_Af0]
    [let
        Aij(σs,two_λs)=Ai(σs,two_λs) * conj(Aj(σs,two_λs))
        I = summed_over_polarization(Aij, tbs.two_js)
        three_body_phase_space_integral(I, tbs)
        end for (Ai,Aj) in Iterators.product(As,As)]
end

Φ0ij_n = Φ0ij .* [√2; √2; 1] .* [√2 √2 1]
Φ0ij_n
