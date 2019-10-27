
function Zsτ(k,two_s,two_τ,dpp,tbs)
    (i,j) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
    #
    val = 0.0;
    for two_λs_prime in Iterators.product(-tbs.two_js[1]:2:tbs.two_js[1],
                                          -tbs.two_js[2]:2:tbs.two_js[2],
                                          -tbs.two_js[3]:2:tbs.two_js[3])
        #
        val +=
            phase(k,1,two_Λ(dpp),two_τ-two_λs_prime[k]) *
            sqrt(two_J(tbs)+1) * wignerd_doublearg(two_J(tbs),two_Λ(dpp),two_τ-two_λs_prime[k],cosθhatk1(k,dpp.σs,tbs.msq)) *
            sqrt(two_s+1) * wignerd_doublearg(two_s,two_τ,two_λs_prime[i]-two_λs_prime[j],cosθij(k,dpp.σs,tbs.msq)) *
                phase(1,k,two_λs_prime[1],dpp.two_λs[1]) * wignerd_doublearg(tbs.two_js[1], two_λs_prime[1], dpp.two_λs[1], cosζk1_for1(k,dpp.σs,tbs.msq)) *
                phase(2,k,two_λs_prime[2],dpp.two_λs[2]) * wignerd_doublearg(tbs.two_js[2], two_λs_prime[2], dpp.two_λs[2], cosζk2_for2(k,dpp.σs,tbs.msq)) *
                phase(3,k,two_λs_prime[3],dpp.two_λs[3]) * wignerd_doublearg(tbs.two_js[3], two_λs_prime[3], dpp.two_λs[3], cosζk3_for3(k,dpp.σs,tbs.msq))
    end
    return val
end
