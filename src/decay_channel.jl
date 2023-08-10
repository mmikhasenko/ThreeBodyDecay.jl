
#        _|                                                        _|                  _|
#    _|_|_|    _|_|      _|_|_|    _|_|_|  _|    _|        _|_|_|  _|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|_|_|_|  _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#  _|    _|  _|        _|        _|    _|  _|    _|      _|        _|    _|  _|    _|  _|  _|    _|
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|    _|_|_|        _|_|_|  _|    _|    _|_|_|  _|  _|    _|
#                                                _|
#                                            _|_|

abstract type Recoupling end
@with_kw struct NoRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
end

amplitude(cs::NoRecoupling, two_λa, two_λb) =
    (cs.two_λa == two_λa) *
    (cs.two_λb == two_λb)

@with_kw struct ParityRecoupling <: Recoupling
    two_λa::Int
    two_λb::Int
    ηηηphaseisplus::Bool
end
ParityRecoupling(two_λa::Int, two_λb::Int, ηηηphasesign::Char) = ParityRecoupling(two_λa, two_λb, ηηηphasesign == '+')
function ParityRecoupling(two_λa::Int, two_λb::Int, (jp, (jp1, jp2))::Pair{jp,Tuple{jp,jp}})
    ηηη = jp1.p ⊗ jp2.p ⊗ jp.p
    ηηηphase = (2 * (ηηη == '+') - 1) * x"-1"^(jp.j - jp1.j - jp2.j)
    return ParityRecoupling(two_λa, two_λb, ηηηphase == 1)
end

function amplitude(cs::ParityRecoupling, two_λa, two_λb)
    (cs.two_λa == two_λa) * (cs.two_λb == two_λb) && return 1
    (cs.two_λa == -two_λa) * (cs.two_λb == -two_λb) && return 2 * cs.ηηηphaseisplus - 1
    return 0
end

@with_kw struct RecouplingLS <: Recoupling
    two_j::Int
    two_ls::Tuple{Int,Int}
    two_ja::Int
    two_jb::Int
end
RecouplingLS(two_ls, (jp, (jpa, jpb))::Pair{jp,Tuple{jp,jp}}) =
    RecouplingLS(jp.j |> x2, two_ls, jpa.j |> x2, jpb.j |> x2)

amplitude(cs::RecouplingLS, two_λa, two_λb) =
    jls_coupling(cs.two_ja, two_λa, cs.two_jb, two_λb, cs.two_j, cs.two_ls[1], cs.two_ls[2])

@with_kw struct DecayChain{X,V1<:Recoupling,V2<:Recoupling,T}
    k::Int
    #
    two_s::Int # isobar spin
    #
    Xlineshape::X # lineshape
    #
    HRk::V1
    Hij::V2
    #
    tbs::T # the structure with masses and spins
end

function printable_l(two_l)
    l = div(two_l, 2)
    waves = ['S', 'P', 'D', 'F', 'G', 'H']
    return l < 6 ? waves[l+1] : string(l)[1]
end
printable_s(two_s) = two_s // 2
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
    Ps=SVector('+', '+', '+', '+'))
    # 
    lsLS = vcat(possible_lsLS(k, two_s, parity, tbs.two_js, Ps)...)
    length(lsLS) == 0 && error("there are no possible LS couplings")
    # 
    lsLS_sorted = sort(lsLS, by=x -> x.LS[1])
    @unpack ls, LS = lsLS_sorted[1]
    #
    i, j = ij_from_k(k)
    return DecayChain(; k, Xlineshape, tbs, two_s,
        Hij=RecouplingLS(two_s, Int.(2 .* ls), tbs.two_js[i], tbs.two_js[j]),
        HRk=RecouplingLS(tbs.two_js[4], Int.(2 .* LS), two_s, tbs.two_js[k]))
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
    two_s=error("give two_s, i.e. the spin of the isobar, two_s=..."),
    parity::Char=error("give the parity, parity=..."),
    Ps=error("need parities: Ps=[P1,2,3,0]"),
    tbs=error("give three-body-system structure, tbs=..."))
    # 
    i, j = ij_from_k(k)
    LSlsv = possible_lsLS(k, two_s, parity, tbs.two_js, Ps)
    return [DecayChain(;
        k, Xlineshape, tbs, two_s,
        Hij=RecouplingLS(two_s, Int.(2 .* x.ls), tbs.two_js[i], tbs.two_js[j]),
        HRk=RecouplingLS(tbs.two_js[4], Int.(2 .* x.LS), two_s, tbs.two_js[k]))
            for x in LSlsv]
end


wignerd_doublearg_sign(two_j, two_λ1, two_λ2, cosθ, ispositive) =
    (ispositive ? 1 : x"-1"^div(two_λ1 - two_λ2, 2)) *
    wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ)

function amplitude(dc::DecayChain, σs, two_λs; refζs=(1, 2, 3, 1))

    @unpack k, tbs, two_s, HRk, Hij = dc
    # 
    i, j = ij_from_k(k)
    #
    ms² = tbs.ms^2
    two_js = tbs.two_js
    # 
    w0 = wr(k, refζs[4], 0)
    wi = wr(k, refζs[i], i)
    wj = wr(k, refζs[j], j)
    wk = wr(k, refζs[k], k)
    # 
    cosζ0 = cosζ(w0, σs, ms²)
    cosζi = cosζ(wi, σs, ms²)
    cosζj = cosζ(wj, σs, ms²)
    cosζk = cosζ(wk, σs, ms²)
    #
    cosθ = cosθij(k, σs, ms²)
    #
    T = typeof(two_λs[1])
    x1 = one(T)
    # 
    itr_two_λs′ = itr(SVector{4,T}(two_js[1], two_js[2], two_js[3], two_js[4]))
    lineshape = dc.Xlineshape(σs[k])
    f = zero(lineshape)
    for two_λs′ in itr_two_λs′
        # two_λs′[4] = two_τ - two_λs′[k]
        # two_τ = two_λs′[4] + two_λs′[k]
        f +=
            wignerd_doublearg_sign(two_js[4], two_λs[4], two_λs′[4], cosζ0, ispositive(w0)) *
            # 
            amplitude(HRk, two_λs′[4] + two_λs′[k], two_λs′[k]) * phase(two_js[k] - two_λs′[k]) * # particle-2 convention
            sqrt(two_s * x1 + 1) * wignerd_doublearg(two_s, two_λs′[4] + two_λs′[k], two_λs′[i] - two_λs′[j], cosθ) *
            amplitude(Hij, two_λs′[i], two_λs′[j]) * phase(two_js[j] - two_λs′[j]) * # particle-2 convention
            # 
            wignerd_doublearg_sign(two_js[i], two_λs′[i], two_λs[i], cosζi, ispositive(wi)) *
            wignerd_doublearg_sign(two_js[j], two_λs′[j], two_λs[j], cosζj, ispositive(wj)) *
            wignerd_doublearg_sign(two_js[k], two_λs′[k], two_λs[k], cosζk, ispositive(wk))
    end
    return f * lineshape
end
#
amplitude(dc::DecayChain, dpp) = amplitude(dc, dpp.σs, dpp.two_λs)
#
summed_over_polarization(fn, two_js) = σs -> sum(fn(σs, two_λs) for two_λs in itr(two_js))
#
itr(two_js) = Iterators.ProductIterator(Tuple([-two_j:2:two_j for two_j in two_js]))
