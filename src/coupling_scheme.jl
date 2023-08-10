
struct jp{T<:Number}
    j::T
    p::Char
end

jp(v::Tuple{T,Char} where {T<:Number}) = jp(v[1], v[2])

import Base: length
length(jp1::jp) = 0
two_j(jp::jp) = Int(2jp.j)
# 
# dealing with spin 1/2
x2(v) = @. Int(2v)
# 
_over2(two) = mod(two, 2) == 0 ? div(two, 2) : two // 2
over2(two) = _over2.(two)

⊗(p1::Char, p2::Char) = p1 == p2 ? '+' : '-'
⊗(jp1::jp, jp2::jp) = [jp(j, ⊗(jp1.p, jp2.p)) for j in abs(jp1.j - jp2.j):abs(jp1.j + jp2.j)]

function str2jp(pin::String)
    p = filter(x -> x != '^', pin)
    !(contains(p, '/')) && return jp(Meta.parse(p[1:end-1]), p[end])
    p[end-2:end-1] != "/2" && error("the string should be `x/2±`, while it is $(p)")
    two_j = Meta.parse(p[1:end-3])
    !(typeof(two_j) <: Int) && error("the string should be `x/2±`, while it is $(p)")
    return jp(two_j // 2, p[end])
end

macro jp_str(p)
    return str2jp(p)
end

possible_ls((jp, (jp1, jp2))::Pair{A,Tuple{B,C}} where {A<:jp,B<:jp,C<:jp}) = possible_ls(jp1, jp2; jp)

function possible_ls(jp1::jp, jp2::jp; jp::jp)
    ls = Vector{Tuple{Int,Number}}(undef, 0)
    for s in abs(jp1.j - jp2.j):abs(jp1.j + jp2.j)
        for l in Int(abs(jp.j - s)):Int(abs(jp.j + s))
            if jp1.p ⊗ jp2.p ⊗ jp.p == (isodd(l) ? '-' : '+')
                push!(ls, (l, s))
            end
        end
    end
    return sort(ls, by=x -> x[1])
end

function possible_lsLS(k::Int, jpR::jp, jps::Vector{T} where {T<:jp})
    i, j = ij_from_k(k)
    lsv = possible_ls(jps[i], jps[j]; jp=jpR)
    LSv = possible_ls(jpR, jps[k]; jp=jps[4])
    return [(ls=ls, LS=LS) for (ls, LS) in Iterators.product(lsv, LSv)]
end
function possible_lsLS(k::Int, two_s, parity::Char, two_js, Ps)
    jpR = jp(two_s // 2, parity)
    jps = jp.(zip(Tuple(two_js) .// 2, Ps))
    return possible_lsLS(k, jpR, jps)
end

function jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s)
    T1 = one(two_λ1) # type consistency
    return sqrt((two_l * T1 + 1) / (two_j * T1 + 1)) *
           CG_doublearg(two_j1, two_λ1, two_j2, -two_λ2, two_s, two_λ1 - two_λ2) *
           CG_doublearg(two_l, zero(two_λ1 - two_λ2), two_s, two_λ1 - two_λ2, two_j, two_λ1 - two_λ2)
end
