using ThreeBodyDecay
using Parameters
using StaticArrays
using PartialWaveFunctions
using Plots
using Plots.PlotMeasures: mm
using LaTeXStrings

struct AngularBusiness
    k::Int
    two_js::Union{ThreeBodySpins, SVector{4, Int}}
    two_ls::Union{NamedTuple,Tuple}
    two_LS::Union{NamedTuple,Tuple}
    two_j::Int
end

function angular_ampitude(z, two_λs, ab::AngularBusiness)
    k = ab.k
    #
    @unpack k, two_ls, two_LS, two_js, two_j = ab
    i,j = ij_from_k(k)
    #
    two_l, two_s = ab.two_ls
    two_L, two_S = ab.two_LS
    #
    hR = jls_coupling(two_js[i], two_λs[i], two_js[j], two_λs[j], two_j, two_l, two_s)
    h0 = jls_coupling(two_j, two_λs[4]+two_λs[k], two_js[k], two_λs[k], two_js[4], two_L, two_S)
    #
    return h0*wignerd_doublearg(two_j,two_λs[4]+two_λs[k],two_λs[i]-two_λs[j], z)*hR
end

angular_intensity(z,ab::AngularBusiness) = 
    sum(abs2, angular_ampitude(z,two_λs,ab) for two_λs in itr(ab.two_js))
#
#
function explore(k, jp::jp, jps::Vector)
    i,j = ij_from_k(k)
    lsv = possible_ls(jps[i], jps[j]; jp = jp)
    LSv = possible_ls(jp,  jps[k]; jp = jps[4])
    ABs = [AngularBusiness(k, SVector{4}(two_j.(jps)), x2.(ls), x2.(LS), two_j(jp))
        for ls in lsv, LS in LSv]
    return ABs
end
# 
#
function vector_of_vector_of_plots(vov)
    l = max(length.(vov)...)
    plot(layout=grid(3,l), size=(300*l,300*3))    
    for (i,ab_all) in enumerate(vov)
        for k in 1:l
            sp = (i-1)*l + k
            if k > length(ab_all)
                plot!(sp=sp, xaxis=false, yaxis=false, grid=false)
                continue
            end
            ab = ab_all[k]
            plot!(sp=sp, ab)
        end        
    end
    plot!()
end


@recipe function f(vov::Vector{Vector{T}} where T)
    l = max(length.(vov)...)
    # 
    size --> (300*l,300*length(vov))
    layout := grid(length(vov), l)
    # 
    for (i,objects) in enumerate(vov)
        for k in 1:l
            sp = (i-1)*l + k
            if k > length(objects)
                @series begin
                    label := ""
                    subplot := sp
                    xaxis := false
                    yaxis := false
                    grid := false
                    ()
                end
                continue
            end
            ab = objects[k]
            @series begin
                subplot := sp
                ab
            end
        end        
    end
end

@recipe function f(ab::AngularBusiness)
    xv = range(-1,1, length=100)
    calv = angular_intensity.(xv, Ref(ab))

    framestyle --> :origin
    title --> "l=$(over2(ab.two_ls[1])), L=$(over2(ab.two_LS[1]))"

    xv, calv
end


# Ωb⁻ → Ξc⁺π⁻K⁻
let
    jps_pc = [jp"1/2+", jp"0-", jp"0-", jp"1/2+"]
    jps_pv = [jp"1/2+", jp"0-", jp"0-", jp"1/2-"]

    jpRv = [jp"1/2+", jp"3/2+", jp"5/2+"]

    ABs_all = []
    for jp in jpRv
        ABm = explore.(3, Ref(jp), [jps_pc, jps_pv])
        push!(ABs_all, hcat(ABm...))
    end
    plot([vcat(v...) for v in ABs_all], lab="", xlab=L"\cos\,\theta")
end

# B⁻ → Λc⁻Λc⁺K⁻
let
    jps_pc = [jp"1/2-", jp"1/2+", jp"0-", jp"0+"]
    jps_pv = [jp"1/2-", jp"1/2+", jp"0-", jp"0-"]

    jpRv = [jp"1/2+", jp"3/2+", jp"5/2+"]

    ABs_all = []
    for jp in jpRv
        ABm = explore.(1, Ref(jp), [jps_pc, jps_pv])
        push!(ABs_all, hcat(ABm...))
    end
    plot([vcat(v...) for v in ABs_all], lab="", xlab=L"\cos\,\theta")
end

# Λb⁰ → J/ψ p K⁻
let
    jps_pc = [jp"1-", jp"1/2+", jp"0-", jp"1/2+"]
    jps_pv = [jp"1-", jp"1/2+", jp"0-", jp"1/2-"]

    jpRv = [jp"1/2+", jp"3/2+", jp"5/2+"]

    ABs_all = []
    for jp in jpRv
        ABm = explore.(3, Ref(jp), [jps_pc, jps_pv])
        push!(ABs_all, hcat(ABm...))
    end
    plot([vcat(v...) for v in ABs_all], lab="", xlab=L"\cos\,\theta")
end


# Ωb⁻ → Ωc⁰π⁻π⁰
let
    jps_pc = [jp"1/2+", jp"0-", jp"0-", jp"1/2+"]
    jps_pv = [jp"1/2+", jp"0-", jp"0-", jp"1/2-"]

    jpRv = [jp"1-"]

    ABs_all = []
    for jp in jpRv
        ABm = explore.(1, Ref(jp), [jps_pc, jps_pv])
        push!(ABs_all, hcat(ABm...))
    end
    plot([vcat(v...) for v in ABs_all], lab="", xlab=L"\cos\,\theta", bottom_margin=5mm)
end

