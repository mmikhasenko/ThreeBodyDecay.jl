
function Zksτ(k,two_s,two_τ,two_λs,σs,tbs)
    (i,j) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
    #
    two_Λ = two_λs[4]
    #
    val = 0.0;
    for two_λs_prime in Iterators.product(-tbs.two_js[1]:2:tbs.two_js[1],
                                          -tbs.two_js[2]:2:tbs.two_js[2],
                                          -tbs.two_js[3]:2:tbs.two_js[3])
        #
        val +=
            phase(k,1,two_Λ,two_τ-two_λs_prime[k]) *
            sqrt(two_J(tbs)+1) * wignerd_doublearg(two_J(tbs),two_Λ,two_τ-two_λs_prime[k],cosθhatk1(k,σs,tbs.msq)) *
            sqrt(two_s+1) * wignerd_doublearg(two_s,two_τ,two_λs_prime[i]-two_λs_prime[j],cosθij(k,σs,tbs.msq)) *
                phase(1,k,two_λs_prime[1],two_λs[1]) * wignerd_doublearg(tbs.two_js[1], two_λs_prime[1], two_λs[1], cosζk1_for1(k,σs,tbs.msq)) *
                phase(2,k,two_λs_prime[2],two_λs[2]) * wignerd_doublearg(tbs.two_js[2], two_λs_prime[2], two_λs[2], cosζk2_for2(k,σs,tbs.msq)) *
                phase(3,k,two_λs_prime[3],two_λs[3]) * wignerd_doublearg(tbs.two_js[3], two_λs_prime[3], two_λs[3], cosζk3_for3(k,σs,tbs.msq))
    end
    return val
end
Zksτ(k,two_s,two_τ,dpp,tbs) = Zksτ(k,two_s,two_τ,two_λs,σs,tbs)

function Rksτ(k,two_s,two_τ,coupls,CScheme,dpp,tbs)
    sum(c*clebsch_for_chaink(k, two_s,two_τ, CS, dpp.two_λs, tbs.two_js) for (c,CS) in zip(coupls,CScheme))
end

function Ask(k,two_s,coupls,CScheme,dpp,tbs)
    return sum(Zksτ(k,two_s,two_τ,dpp,tbs) * Rksτ(k,two_s,two_τ,coupls,CScheme,dpp,tbs) for two_τ=-two_s:2:two_τ)
end

function Asummed(setup,dpp,tbs)
    return sum(RZsk(k,two_s,coupls,CScheme,dpp,tbs) for (k,two_s,coupls,CScheme) in setup)
end

function Asq(setup,dpp,tbs)
    return sum(abs2(Asummed(setup,dpp_λs,tbs)) for dpp_λs in possible_helicities(dpp,tbs))
end


#        _|                                                        _|                  _|
#    _|_|_|    _|_|      _|_|_|    _|_|_|  _|    _|        _|_|_|  _|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|_|_|_|  _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#  _|    _|  _|        _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|    _|_|_|        _|_|_|  _|    _|    _|_|_|  _|  _|    _|
#                                                _|
#                                            _|_|

struct decay_chain
    k::Int
    #
    two_s::Int # isobar spin
    #
    two_ls::Tuple # isobar decay ξ->ij
    two_LS::Tuple # 0->ξ k decay
    #
    tbs::ThreeBodySystem # the structure with masses and spins
    #
    Xlineshape::Function # lineshape
end
#     two_λs=error("two_λs is required"),

# function decay_chain(k;
#     two_s::Vector{Int}=error("two_s, i.e. spin of the isobar is required"),
#     two_js::Vector{Int}=error("two_js is required"),
#     Ps::Vector{Char}=['N','N','N','N'], # default values
#     two_ls=(-1,-1), # default values
#     two_LS=(-1,-1) # default values
#     )
#     # do something to figure out minimal LS ls
#     return decay_chain(k, two_s, two_ls, two_LS)
# end

function decay_chain(k, Xlineshape;
    two_s=error("give two_s, i.e. spin of the isobar is required"),
    parity::Char='+',
    tbs=error("give three-body-system structure"),
    two_ls=(-1,-1), # default values
    two_LS=(-1,-1) # default values
    )
    k = k; i,j = ij_from_k(k);
    #
    possible_ls = possibleLS((tbs.two_js[i],tbs.Ps[i]), (tbs.two_js[j],tbs.Ps[j]), (two_s,parity))
    possible_LS = possibleLS((two_s,parity), (tbs.two_js[k],tbs.Ps[k]), (tbs.two_js[4],tbs.Ps[4]))
    #
    length(possible_ls)==0 && error("no ls are possible! check parity")
    length(possible_LS)==0 && error("no LS are possible! check parity")
    #
    used_two_ls = two_ls
    if two_ls[1] < 0
        println("Possible (l,s) for the isobar decay: ",possible_ls,"\nI use the first one: ", possible_ls[1],"\n")
        used_two_ls = possible_ls[1]
    end
    #
    used_two_LS = two_LS
    if two_LS[1] < 0
        println("Possible (L,S) for a decay to isobar-spectator: ",possible_ls,"\nI use the first one: ", possible_ls[1],"\n")
        used_two_LS = possible_LS[1]
    end
    return decay_chain(k, two_s, used_two_ls, used_two_LS, tbs, Xlineshape)
end


function amplitude(σs, two_λs, dc)
    k = dc.k; i,j = ij_from_k(k);
    tbs = dc.tbs
    #
    two_s = dc.two_s
    two_js = tbs.two_js
    #
    f = sum(jls_coupling(two_js[i], two_λs[i], two_js[j], two_λs[j], two_s, dc.two_ls[1], dc.two_ls[2]) *
            Zksτ(k,two_s,two_τ,two_λs,σs,tbs) *
            jls_coupling(two_s, two_τ, two_js[k], two_λs[k], two_js[4], dc.two_LS[1], dc.two_LS[2])
        for two_τ = -two_s:2:two_s)
    return f * dc.Xlineshape(s,σs[k])
    # return Zksτ(k,two_s,two_τ,two_λs,σs,tbs)
end
#
amplitude(dpp, dc) = amplitude(dpp.σs, dpp.two_λs, dc)
