
#        _|                                                        _|                  _|
#    _|_|_|    _|_|      _|_|_|    _|_|_|  _|    _|        _|_|_|  _|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|_|_|_|  _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#  _|    _|  _|        _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|    _|_|_|        _|_|_|  _|    _|    _|_|_|  _|  _|    _|
#                                                _|
#                                            _|_|

@with_kw struct CoulingsLS
    two_ls::Tuple{Int,Int} # isobar decay ξ->ij
    two_LS::Tuple{Int,Int} # 0->ξ k decay
end

function amplitude(cs::CoulingsLS, two_τ, two_s,
    two_λi, two_λj, two_λk,
    two_ji, two_jj, two_jk, two_j0)
    return jls_coupling(two_ji, two_λi, two_jj, two_λj, two_s, cs.two_ls[1], cs.two_ls[2]) *
           jls_coupling(two_s, two_τ, two_jk, two_λk, two_j0, cs.two_LS[1], cs.two_LS[2])
end

@with_kw struct DecayChain{X,V,T}
    k::Int
    #
    two_s::Int # isobar spin
    #
    Xlineshape::X # lineshape
    #
    couplingproduct::V
    #
    tbs::T # the structure with masses and spins
end

function printable_l(two_l)
    l = div(two_l,2)
    waves = ['S', 'P', 'D', 'F', 'G', 'H'];
    return l < 6 ? waves[l+1] : string(l)[1]
end
printable_s(two_s) = two_s//2
printable_ls(two_ls) = (printable_s(two_ls[2]), printable_l(two_ls[1]))

"""
    DecayChainLS(k, Xlineshape;
        two_s, # give two_s, i.e. the spin of the isobar
        parity, # 
        Ps, # need parities: [1,2,3,0]
        tbs) # give three-body-system structure

    Returns the decay chain with the smallest LS, ls
"""
function DecayChainLS(k, Xlineshape;
    tbs=error("give three-body-system structure"),
    two_s=0,
    parity::Char='+',
    Ps=SVector('+','+','+','+'))
    # 
    lsLS = vcat(possible_lsLS(k, two_s, parity, tbs.two_js, Ps)...)
    length(lsLS) == 0 && error("there are no possible LS couplings")
    # 
    lsLS_sorted = sort(lsLS, by=x->x.LS[1])
    @unpack ls, LS = lsLS_sorted[1]
    # 
    return DecayChain(; k, Xlineshape, tbs, two_s,
        couplingproduct=CoulingsLS(two_ls=Int.(2 .* ls), two_LS=Int.(2 .* LS)))
end

"""
    DecayChainsLS(k, Xlineshape;
        two_s, # give two_s, i.e. the spin of the isobar
        parity, # 
        Ps, # need parities: [1,2,3,0]
        tbs) # give three-body-system structure

    Returns an array of the decay chains with all possible couplings
"""
function DecayChainsLS(k, Xlineshape;
    two_s = error("give two_s, i.e. the spin of the isobar, two_s=..."),
    parity::Char = error("give the parity, parity=..."),
    Ps = error("need parities: Ps=[P1,2,3,0]"),
    tbs = error("give three-body-system structure, tbs=..."))
    # 
    LSlsv = possible_lsLS(k, two_s, parity, tbs.two_js, Ps)
    return [DecayChain(;
        k, Xlineshape, tbs, two_s,
        couplingproduct = CoulingsLS(
            two_ls=Int.(2 .* x.ls), two_LS=Int.(2 .* x.LS)))
        for x in LSlsv]
end

function amplitude(σs, two_λs, dc)
    k = dc.k; i,j = ij_from_k(k);
    tbs = dc.tbs
    s = tbs.ms.m0^2
    #
    two_s = dc.two_s
    two_js = tbs.two_js
    #
    itr_two_λs′ = itr(SVector{3}(tbs.two_js[1],tbs.two_js[2],tbs.two_js[3]))
    f = 0.0
    for two_τ = -two_s:2:two_s, two_λs′ in itr_two_λs′
        f += Zksτ(k,two_s,two_τ,two_λs,two_λs′,σs,tbs) *
            amplitude(dc.couplingproduct,
                two_τ, two_s,
                two_λs′[i], two_λs′[j], two_λs′[k],
                two_js[i], two_js[j], two_js[k], two_js[4]
            )
    end
    lineshape = dc.Xlineshape(s,σs[k])
    return f * lineshape
end
#
amplitude(dpp, dc) = amplitude(dpp.σs, dpp.two_λs, dc)
#
summed_over_polarization(fn, two_js) = σs->sum(fn(σs,two_λs) for two_λs in itr(two_js))
#
itr(two_js) = Iterators.ProductIterator(Tuple([-two_j:2:two_j for two_j in two_js]))
