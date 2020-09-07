function Zksτ(k,two_s,two_τ,two_λs,two_λs′,σs,tbs)
    (i,j) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
    #
    two_Λ = two_λs[4]
    two_J = tbs.two_js.two_h0
    #
    val =
        phase(k,1,two_Λ,two_τ-two_λs′[k]) *
            phase(tbs.two_js[k]-two_λs′[k]) * # particle-2 convention
        sqrt(two_J+1) * wignerd_doublearg(two_J,two_Λ,two_τ-two_λs′[k],cosθhatk1(k,σs,tbs.ms^2)) *
            phase(tbs.two_js[j]-two_λs′[j]) * # particle-2 convention
        sqrt(two_s+1) * wignerd_doublearg(two_s,two_τ,two_λs′[i]-two_λs′[j],cosθij(k,σs,tbs.ms^2)) *
            phase(1,k,two_λs′[1],two_λs[1]) * wignerd_doublearg(tbs.two_js[1], two_λs′[1], two_λs[1], cosζk1_for1(k,σs,tbs.ms^2)) *
            phase(2,k,two_λs′[2],two_λs[2]) * wignerd_doublearg(tbs.two_js[2], two_λs′[2], two_λs[2], cosζk2_for2(k,σs,tbs.ms^2)) *
            phase(3,k,two_λs′[3],two_λs[3]) * wignerd_doublearg(tbs.two_js[3], two_λs′[3], two_λs[3], cosζk3_for3(k,σs,tbs.ms^2))
    return val
end
