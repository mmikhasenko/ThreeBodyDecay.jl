
#        _|                                                        _|                  _|
#    _|_|_|    _|_|      _|_|_|    _|_|_|  _|    _|        _|_|_|  _|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|_|_|_|  _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#  _|    _|  _|        _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|    _|_|_|        _|_|_|  _|    _|    _|_|_|  _|  _|    _|
#                                                _|
#                                            _|_|



@with_kw struct decay_chain
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

function printable_l(two_l)
    l = div(two_l,2)
    waves = ['S', 'P', 'D', 'F', 'G', 'H'];
    return l < 6 ? waves[l+1] : string(l)[1]
end
printable_s(two_s) = two_s//2
printable_ls(two_ls) = (printable_s(two_ls[2]), printable_l(two_ls[1]))

"""
    decay_chain(k, Xlineshape;
        two_s, # give two_s, i.e. the spin of the isobar
        parity, # 
        Ps, # need parities: [1,2,3,0]
        tbs) # give three-body-system structure

    Returns the decay chain with the smallest LS, ls
"""
function decay_chain(k, Xlineshape;
    two_s=error("give two_s, i.e. the spin of the isobar"),
    tbs=error("give three-body-system structure"),
    parity::Char='+',
    Ps=SVector('+','+','+','+'))
    # 
    (two_LS,two_ls) = first(
            possibleLSls(k; two_s=two_s, Ps=Ps, two_js=tbs.two_js, parity=parity))
    return decay_chain(; k=k, Xlineshape=Xlineshape,
        two_s=two_s, two_ls=Tuple(two_ls), two_LS=Tuple(two_LS),
        tbs=tbs)
end

"""
    decay_chains(k, Xlineshape;
        two_s, # give two_s, i.e. the spin of the isobar
        parity, # 
        Ps, # need parities: [1,2,3,0]
        tbs) # give three-body-system structure

    Returns an array of the decay chains with all possible couplings
"""
function decay_chains(k, Xlineshape;
    two_s=error("give two_s, i.e. the spin of the isobar, two_s=..."),
    parity::Char=error("give the parity, parity=..."),
    Ps=error("need parities: Ps=[P1,2,3,0]"),
    tbs=error("give three-body-system structure, tbs=..."))
    # 
    LSlsv = possibleLSls(k; two_s=two_s, Ps=Ps, two_js=tbs.two_js, parity=parity)
    return [decay_chain(;k=k, Xlineshape=Xlineshape,
        two_s=two_s, two_ls=Tuple(two_ls), two_LS=Tuple(two_LS),
        tbs=tbs) for (two_LS,two_ls) in LSlsv]
end

function possibleLSls(k;
    two_js=error("need spin values: [two_j1,2,3,0]"),
    two_s=error("give two_s, i.e. the spin of the isobar"),
    parity::Char=error("give the parity, parity=..."),
    Ps=error("need parities: [P1,2,3,0]"))
    # 
    k = k; i,j = ij_from_k(k);
    #
    possible_ls = possibleLS((two_js[i],Ps[i]), (two_js[j],Ps[j]), (two_s,parity))
    possible_LS = possibleLS((two_s,parity), (two_js[k],Ps[k]), (two_js[4],Ps[4]))
    #
    return collect(Iterators.product(
        NamedTuple{(:two_L,:two_S)}.(possible_LS),
        NamedTuple{(:two_l,:two_s)}.(possible_ls)))    
end

[(two_ls=Tuple(two_ls), two_LS=Tuple(two_LS)) for (two_LS,two_ls) in collect(Iterators.product(-4:0, 0:4))]

function amplitude(σs, two_λs, dc)
    k = dc.k; i,j = ij_from_k(k);
    tbs = dc.tbs
    s = tbs.ms.m0^2
    #
    two_s = dc.two_s
    two_js = tbs.two_js
    #
    itr_two_λs′ = itr(SVector{3}(tbs.two_js[1],tbs.two_js[2],tbs.two_js[3]))
    f = sum(jls_coupling(two_js[i], two_λs′[i], two_js[j], two_λs′[j], two_s, dc.two_ls[1], dc.two_ls[2]) *
            Zksτ(k,two_s,two_τ,two_λs,two_λs′,σs,tbs) *
            jls_coupling(two_s, two_τ, two_js[k], two_λs′[k], two_js[4], dc.two_LS[1], dc.two_LS[2])
        for two_τ = -two_s:2:two_s, two_λs′ in itr_two_λs′)
    lineshape = dc.Xlineshape(s,σs[k])
    return f * lineshape
end
#
amplitude(dpp, dc) = amplitude(dpp.σs, dpp.two_λs, dc)
#
summed_over_polarization(fn, two_js) = σs->sum(fn(σs,two_λs) for two_λs in itr(two_js))
#
itr(two_js) = Iterators.ProductIterator(Tuple([-two_j:2:two_j for two_j in two_js]))
