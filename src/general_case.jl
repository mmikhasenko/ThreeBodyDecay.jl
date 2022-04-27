function Zksτ(k, two_s,two_τ, two_λs, two_λs′, σs, tbs; refζs=(1,2,3,1))
    i,j = ij_from_k(k)
    #
    ms² = tbs.ms^2
    # 
    w0 = wr(k,refζs[4],0); cosζ0 = cosζ(w0, σs, ms²)
    wi = wr(k,refζs[i],i); cosζi = cosζ(wi, σs, ms²)
    wj = wr(k,refζs[j],j); cosζj = cosζ(wj, σs, ms²)
    wk = wr(k,refζs[k],k); cosζk = cosζ(wk, σs, ms²)
    # 
    cosθ = cosθij(k,σs,ms²)
    two_λ0 = two_λs[4]
    two_j0 = tbs.two_js[4]
    #
    val = # the constant factors are important for the normalization
        sqrt(two_j0+1) * wignerd_doublearg_sign(two_j0, two_λ0, two_τ-two_λs′[k], cosζ0, ispositive(w0)) *
            # 
            phase(tbs.two_js[k]-two_λs′[k]) * # particle-2 convention
        sqrt(two_s+1) * wignerd_doublearg(two_s, two_τ, two_λs′[i]-two_λs′[j], cosθ) *
            phase(tbs.two_js[j]-two_λs′[j]) * # particle-2 convention
            # 
            wignerd_doublearg_sign(tbs.two_js[i], two_λs′[i], two_λs[i], cosζi, ispositive(wi)) *
            wignerd_doublearg_sign(tbs.two_js[j], two_λs′[j], two_λs[j], cosζj, ispositive(wj)) *
            wignerd_doublearg_sign(tbs.two_js[k], two_λs′[k], two_λs[k], cosζk, ispositive(wk))
    return val
end

wignerd_doublearg_sign(j,λ1,λ2, cosθ, ispositive) =
    (ispositive ? 1 : x"-1"^*(λ1-λ2)) *
        wignerd_doublearg(j,λ1,λ2, cosθ)
