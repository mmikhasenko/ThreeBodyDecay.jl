

function Zsτ(k,two_s,two_τ,two_Λ,two_λs,dpp,tbs)
    (i,j) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
    #
    val = 0.0;
    for two_λ1_prime=-tbs.two_js[1]:2:tbs.two_js[1],
        two_λ2_prime=-tbs.two_js[2]:2:tbs.two_js[2],
        two_λ3_prime=-tbs.two_js[3]:2:tbs.two_js[3]
        two_λs = SVector(two_λ1_prime,two_λ2_prime,two_λ3_prime)
        #
        val +=
            phase(k,1,two_λ1_prime,two_λs[1]) *
            sqrt(tbs.two_J+1) * wignerd_doublearg(tbs.two_J,two_Λ,two_τ-two_λs[k],cosθhatk1(k,dpp.σ123,tbs.msq)) *
            sqrt(two_s+1) * wignerd_doublearg(two_s,two_τ,two_λs[i]-two_λs[j],cosθij(k,dpp.σ123,tbs.msq)) *
                phase(1,k,two_λ1_prime,two_λs[1]) * wignerd_doublearg(tbs.two_js[1],two_λ1_prime,two_λs[1],cosζk1_for1(k,dpp.σ123,tbs.msq)) *
                phase(2,k,two_λ1_prime,two_λs[1]) * wignerd_doublearg(tbs.two_js[2],two_λ2_prime,two_λs[2],cosζk2_for2(k,dpp.σ123,tbs.msq)) *
                phase(3,k,two_λ1_prime,two_λs[1]) * wignerd_doublearg(tbs.two_js[3],two_λ3_prime,two_λs[3],cosζk3_for3(k,dpp.σ123,tbs.msq))
    end
    return val
end
