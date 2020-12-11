
vtype = Union{Rational{Int}, Int}
struct jp
    j::vtype
    p::Char
end

jp(v::Tuple{vtype, Char}) = jp(v[1],v[2])

import Base: length
length(jp1::jp) = 0
two_j(jp::jp) = Int(2jp.j)
# 
# dealing with spin 1/2
x2(v) = @. Int(2v)
# 
_over2(two) = mod(two,2)==0 ? div(two,2) : two // 2
over2(two) = _over2.(two)

⊗(p1::Char,p2::Char) = p1==p2 ? '+' : '-'
⊗(jp1::jp, jp2::jp) = [jp(j, ⊗(jp1.p, jp2.p)) for j in abs(jp1.j-jp2.j):abs(jp1.j+jp2.j)]

macro jp_str(p)
    !(contains(p, '/')) && return jp(Meta.parse(p[1:end-1]), p[end])
    p[end-2:end-1] != "/2" && error("the string should be `x/2±`, while it is $(p)")
    two_j = Meta.parse(p[1:end-3])
    !(typeof(two_j) <: Int) && error("the string should be `x/2±`, while it is $(p)")
    return jp(two_j//2, p[end])
end

function possible_ls(jp1::jp, jp2::jp; jp::jp)
    ls = Vector{Tuple{Int,vtype}}(undef,0)
    for s in abs(jp1.j-jp2.j):abs(jp1.j+jp2.j)
        for l in Int(abs(jp.j-s)):Int(abs(jp.j+s))
            if jp1.p ⊗ jp2.p ⊗ jp.p  == (isodd(l) ? '-' : '+')
                push!(ls, (l,s))
            end
        end
    end
    return ls
end

function possible_lsLS(k::Int, jpR::jp, jps::T where T<:Vector{jp})
    i,j = ij_from_k(k)
    lsv = possible_ls(jps[i], jps[j]; jp = jpR)
    LSv = possible_ls(jpR, jps[k]; jp = jps[4])
    return [(ls=ls, LS=LS) for (ls, LS) in Iterators.product(lsv, LSv)]
end
function possible_lsLS(k::Int, two_s::Int, parity::Char, two_js , Ps)
    i,j = ij_from_k(k)
    jpR = jp(two_s//2, parity)
    jps = jp.(zip(two_js .// 2, Ps))
    return possible_lsLS(k, jpR, jps)
end



# function possibleLS(two_jp1,two_jp2,two_jp) # expects tuples
#     two_ls = Vector{Tuple{Int,Int}}(undef,0)
#     for two_s in abs(two_jp1[1]-two_jp2[1]):2:abs(two_jp1[1]+two_jp2[1])
#         for two_l in abs(two_jp[1]-two_s):2:abs(two_jp[1]+two_s)
#             if (two_jp1[2] == '+' ? 1 : -1) *
#                (two_jp2[2] == '+' ? 1 : -1) *
#                (two_jp[2]  == '+' ? 1 : -1) == (two_l % 4 == 2  ? -1 : 1)
#                 push!(two_ls, (two_l,two_s))
#             end
#         end
#     end
#     return two_ls
# end

jls_coupling(two_j1, two_λ1, two_j2, two_λ2, two_j, two_l, two_s) =
    CG_doublearg(two_j1, two_λ1, two_j2, -two_λ2, two_s, two_λ1-two_λ2) *
    CG_doublearg(two_l,  0, two_s, two_λ1-two_λ2, two_j, two_λ1-two_λ2)

# function clebsch_for_chaink(k, two_s_int, two_τ, chain, two_λs, two_js)
#     (i,j) = ij_from_k(k)
#     #
#     v = 1.0;
#     two_λi_λj = two_λs[i]-two_λs[j]
#     v *= CG_doublearg(two_js[i],two_λs[i],two_js[j],-two_λs[j],two_s(chain),two_λi_λj) *
#          CG_doublearg(two_l(chain),0,two_s(chain),two_λi_λj,two_s_int,two_λi_λj)
#     #
#     two_τ_λk = two_τ - two_λs[k]
#     v *= CG_doublearg(two_s_int,two_τ,two_js[k],-two_λs[k],two_S(chain),two_τ_λk) *
#          CG_doublearg(two_L(chain),0,two_S(chain),two_τ_λk,two_J(chain),two_τ_λk)
#     return v
# end

# #(j1λ1j2λ2,JLS)
# function HelicityRecoupling_doublearg(HLSpairs)
#     v = 1.0;
#     for p in HLSpairs
#         (two_j1,two_λ1,two_j2,two_λ2) = p[1]
#         (two_J,two_L,two_S) = p[2]
#         two_λ1_λ2 = two_λ1-two_λ2
#         v *= CG_doublearg(two_j1,two_λ1,two_j2,-two_λ2,two_S,two_λ1_λ2) *
#              CG_doublearg(two_L,0,two_S,two_λ1_λ2,two_J,two_λ1_λ2)
#     end
#     return  v;
# end
# function HelicityRecoupling(HLSpairs)
#     v = 1.0;
#     for p in HLSpairs
#         (j1,λ1,j2,λ2) = p[1]
#         (J,L,S) = p[2]
#         v *= CG_doublearg(2*j1,2*λ1,2*j2,-2*λ2,2*S,2*(λ1-λ2)) *
#              CG_doublearg(2*L,0,2*S,2*(λ1-λ2),2*J,2*(λ1-λ2))
#     end
#     return  v;
# end
