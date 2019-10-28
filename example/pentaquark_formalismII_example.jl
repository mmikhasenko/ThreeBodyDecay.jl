# Pentaquark decay chain: Lb -> Jψ p K
#   in the s-channel scattering kinematics
#   Lb Jψ -> p K
#
# it is check for an example described in arXiv:1805.02113
# https://arxiv.org/pdf/1805.02113.pdf
# at the page 29

using ThreeBodyDecay

mLb = 5.62; mJψ = 3.09; mp=0.938; mK = 0.49367
const tbs = ThreeBodySystem(mLb,mJψ,mp,mK; two_jps = ([2,1,0,1],['-','+','-','+']))


# helicity couplings
# Lb Jψ
H1(two_λ0,two_λ1;two_LS=(4,3),two_s=3) = ClGd(1,two_λ0,2,-two_λ1,two_LS[2],two_λ0-two_λ1) *
                                         ClGd(two_LS[1],0,two_LS[2],two_λ0-two_λ1,two_s,two_λ0-two_λ1)
# pK
H2(two_λ2;two_l=4,two_s=3) = ClGd(1,two_λ2,0,0,1,two_λ2) *
                             ClGd(two_l,0,1,two_λ2,two_s,two_λ2)

# amplitude
function A(z,two_λs;two_s=3,two_LS=(4,3),two_l=4)
    return H1(two_λs[4],two_λs[1]; two_LS=two_LS, two_s=two_s) *
        H2(two_λs[2]; two_l=two_l, two_s=two_s) *
        # (two_λs[1]!=0 ? 1.0 : 2.0) * # weight of A(λ_J/ψ = 0) could be different with E/m factor
        wignerd_doublearg(two_s,two_λs[4]-two_λs[1],two_λs[2],z)
end

# amplitude squared
Asq_summed(z;two_s=3,two_LS=(4,3),two_l=4) =
    sum(abs2(A(z,two_λs;two_s=3,two_LS=two_LS,two_l=two_l))
        for two_λs in Iterators.product(-tbs.two_js[1]:2:tbs.two_js[1],
                                        -tbs.two_js[2]:2:tbs.two_js[2],
                                        -tbs.two_js[3]:2:tbs.two_js[3],
                                        -tbs.two_js[4]:2:tbs.two_js[4]))
I(σs) = Asq_summed(cosθ23(σs,tbs.msq))

## plotting
using Plots
# S=3/2, D-wave
plot(z->Asq_summed(z),-1,1, lab="", ylim=(0,1))
# P-wave
plot(z->Asq_summed(z;two_s=3,two_LS=(2,3),two_l=4),-1,1, lab="", ylim=(0,1))
