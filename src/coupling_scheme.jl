
struct jls
    two_j::Int
    two_l::Int
    two_s::Int
end

struct twochain
    first::jls
    second::jls
end

two_J(ch::twochain) = ch.first.two_j
two_L(ch::twochain) = ch.first.two_l
two_S(ch::twochain) = ch.first.two_s
two_j(ch::twochain) = ch.second.two_j
two_l(ch::twochain) = ch.second.two_l
two_s(ch::twochain) = ch.second.two_s

function posibleLS(two_jp,two_jp1,two_jp2)
    two_ls = Vector{Tuple{Int,Int}}(undef,0)
    for two_s in abs(two_jp1[1]-two_jp2[1]):2:abs(two_jp1[1]+two_jp2[1])
        for two_l in abs(two_jp[1]-two_s):2:abs(two_jp[1]+two_s)
            if (two_jp1[2] == "+" ? 1 : -1) *
               (two_jp2[2] == "+" ? 1 : -1) *
               (two_jp[2]  == "+" ? 1 : -1) == (two_l % 4 == 2  ? -1 : 1)
                push!(two_ls, (two_l,two_s))
            end
        end
    end
    return two_ls
end

function coupling_scheme23(two_JP,two_jp,two_jps)
    pls = posibleLS(two_jp,two_jps[2],two_jps[3])
    length(pls) != 1 && error("length(pls) = $(length(pls)) != 1 for the second decay: $(two_jp...)=>$(two_jps[2]...),$(two_jps[3]...). Check quantum numbers!");
    return [twochain(
                jls(two_JP[1],two_LS...),  # JLS
                jls(two_jp[1],pls[1]...))  # jls
                    for two_LS in posibleLS(two_JP,two_jp,two_jps[1])]
end
coupling_scheme12(two_JP,two_jp,two_jps) = coupling_scheme23(two_JP,two_jp,[two_jps[3],two_jps[1],two_jps[2]])
coupling_scheme31(two_JP,two_jp,two_jps) = coupling_scheme23(two_JP,two_jp,[two_jps[2],two_jps[3],two_jps[1]])
