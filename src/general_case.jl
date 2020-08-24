two_J(tbs) = tbs.two_js.two_h0

function Zksτ(k,two_s,two_τ,two_λs,two_λs′,σs,tbs)
    (i,j) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
    #
    two_Λ = two_λs[4]
    #
    val =
        phase(k,1,two_Λ,two_τ-two_λs′[k]) *
            phase(tbs.two_js[k]-two_λs′[k]) * # particle-2 convention
        sqrt(two_J(tbs)+1) * wignerd_doublearg(two_J(tbs),two_Λ,two_τ-two_λs′[k],cosθhatk1(k,σs,tbs.ms^2)) *
            phase(tbs.two_js[j]-two_λs′[j]) * # particle-2 convention
        sqrt(two_s+1) * wignerd_doublearg(two_s,two_τ,two_λs′[i]-two_λs′[j],cosθij(k,σs,tbs.ms^2)) *
            phase(1,k,two_λs′[1],two_λs[1]) * wignerd_doublearg(tbs.two_js[1], two_λs′[1], two_λs[1], cosζk1_for1(k,σs,tbs.ms^2)) *
            phase(2,k,two_λs′[2],two_λs[2]) * wignerd_doublearg(tbs.two_js[2], two_λs′[2], two_λs[2], cosζk2_for2(k,σs,tbs.ms^2)) *
            phase(3,k,two_λs′[3],two_λs[3]) * wignerd_doublearg(tbs.two_js[3], two_λs′[3], two_λs[3], cosζk3_for3(k,σs,tbs.ms^2))
    return val
end

# 
function QTB_mismatch_factor(dc)
    k = dc.k; i,j = ij_from_k(k);
    tbs = dc.tbs
    two_js = tbs.two_js
    two_s = dc.two_s
    #
    avHHsq =
        sum((tbs.two_js[4]!=(two_τ-two_λs[k]) ? 0.0 : 1.0) *
            jls_coupling(two_js[i], two_λs[i], two_js[j], two_λs[j], two_s, dc.two_ls[1], dc.two_ls[2])^2 *
            jls_coupling(two_s, two_τ, two_js[k], two_λs[k], two_js[4], dc.two_LS[1], dc.two_LS[2])^2
            for two_λs in itr(tbs.two_js),
                two_τ in -dc.two_s:2:dc.two_s)
    return (two_J(dc.tbs)+1) * avHHsq
end
