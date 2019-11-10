using ThreeBodyDecay

const mΛb = 5.61960
const mπ = 0.14
# intermediate resonances
const mΣb = 5.81065; const ΓΣb = 0.005
const mΣb_x = 5.83032; const ΓΣb_x = 0.0094

const mΛb2S = 6.3; # just a peak of the plot

tbs = ThreeBodySystem(mΛb2S,mπ,mΛb,mπ,
    two_jps = ([0, 1, 0, 1], ['-','+','-','+']))  # 1/2+ 0- 0- 1/2+
dpp = randomPoint(tbs)

# lineshape
dc_Σb1 = decay_chain(1, (s,σ) -> BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')
dc_Σb3 = decay_chain(3, (s,σ) -> BW(σ,mΣb,ΓΣb); tbs=tbs, two_s = 1, parity='+')
full_AΣb(σs,two_λs)   = sum(amplitude(σs, two_λs, ch) for ch in [dc_Σb1, dc_Σb3])

dc_Σb_x1 = decay_chain(1, (s,σ) -> BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')
dc_Σb_x3 = decay_chain(3, (s,σ) -> BW(σ,mΣb_x,ΓΣb_x); tbs=tbs, two_s = 3, parity='+')
full_AΣb_x(σs,two_λs) = sum(amplitude(σs, two_λs, ch) for ch in [dc_Σb_x1, dc_Σb_x3])

dc_f0 = decay_chain(2, (s,σ) -> 1.0; tbs=tbs, two_s = 0, parity='+')
full_Af0(σs,two_λs)   = amplitude(σs, two_λs, dc_f0)

####

using Cuba

function three_body_phase_space_integral(function_σs, tbs)
    msq = tbs.msq
    m1sq,m2sq,m3sq,s = msq
    #
    σ1min,σ1max = tbs.mthsq[1],tbs.sthsq[1]
    function integrand(x,f)
        σ1 = σ1min + x[1]*(σ1max-σ1min)
        z = 2*x[2]-1
        σ3 = σ3of1(z,σ1,msq); σ2 = gσ2(σ3,σ1,msq)
        σs = [σ1, σ2, σ3]
        #
        r = function_σs(σs) *
            sqrt(λ(s,σ1, m1sq)*λ(σ1,m2sq,m3sq))/σ1
        f[1],f[2] = reim(r)
    end
    return complex(cuhre(integrand,2,2)[1]...) / s # /(2π*(8π)^2)
end

Φij = let
    As = [full_AΣb, full_AΣb_x, full_Af0]
    [let
        Aij(σs,two_λs)=Ai(σs,two_λs) * conj(Aj(σs,two_λs))
        I = summed_over_polarization(Aij, tbs)
        three_body_phase_space_integral(I, tbs)
        end for (Ai,Aj) in Iterators.product(As,As)]
end

using LinearAlgebra
Φij_n = let
    d = sqrt.(diag(Φij))
    Φij ./ d ./ d'
end
