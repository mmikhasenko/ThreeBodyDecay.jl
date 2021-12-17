
abstract type AbstractWignerRotation end
struct TriavialWignerRotation <: AbstractWignerRotation
    k::Int
end
struct WignerRotation{N} <: AbstractWignerRotation
    k::Int
    ispositive::Bool
    iseven::Bool
end

const Arg0WignerRotation = WignerRotation{0}
const Arg2WignerRotation = WignerRotation{2}
const Arg3WignerRotation = WignerRotation{3}

WignerRotation{N}(k::Int, ispositive::Bool) where N =
    WignerRotation{N}(k, ispositive, true)

ispositive(wr::TriavialWignerRotation) = true
ispositive(wr::WignerRotation) = wr.ispositive

import Base:iseven
iseven(wr::TriavialWignerRotation) = true
iseven(wr::WignerRotation{0}) = true
iseven(wr::WignerRotation{3}) = true
iseven(wr::WignerRotation{2}) = wr.iseven


ijk(k::Int) = (k+1,k+2,k) |> x->mod.(x,Ref(Base.OneTo(3)))
ijk(wr::AbstractWignerRotation) = ijk(wr.k)

issequential(i,j) = (j-i) ∈ (1,-2)


"""
    wr(system_a, reference_b, particle_c=0)

Create a WignerRotation object of the right type based on provided indices.
The daughter particles are numbered 1,2,3, the mother particle is 0.
 - `system_a` tells which isobar is considered, and
 - `reference_b` tell which system is used as a reference.
For `system_a` and `reference_b` the spectator notations are used, i.e.
1 for the system (2,3), 2 for the system (3,1), and 3 for the system (1,2).
"""
function wr(system_a, reference_b, particle_c=0)
    system_a == reference_b && return TriavialWignerRotation(particle_c)
    S = issequential(system_a, reference_b)
    A,B = S ? (system_a, reference_b) : (reference_b, system_a)
    # 
    particle_c == 0 && return Arg0WignerRotation(A,S)
    #
    particle_c ∉ (system_a, reference_b) && return Arg3WignerRotation(particle_c,S)
    #
    T = (particle_c == A)
    return Arg2WignerRotation(particle_c,!(S),T)
end



cosζ(wr::TriavialWignerRotation,σs,msq) = 1.0

function cosζ(wr::Arg0WignerRotation,σs,msq)
    i,j,k = ijk(wr)
    # 
    s = msq[4]
    EE4s = (s+msq[i]-σs[i])*(s+msq[k]-σs[k])
    pp4s = sqrt(Kallen(s,msq[i],σs[i])*Kallen(s,msq[k],σs[k]))
    rest = σs[j]-msq[i]-msq[k]
    return (EE4s-2s*rest)/pp4s
end

function cosζ(wr::Arg2WignerRotation,σs,msq)
    i,j,k = ijk(wr)
    # 
    if !(iseven(wr))
        i,j = j,i 
    end
    # 
    msq[k] ≈ 0 && return 1.0
    # 
    s = msq[4]
    EE4mksq = (s+msq[k]-σs[k])*(σs[i]-msq[k]-msq[j])
    pp4mksq = sqrt(Kallen(s,msq[k],σs[k])*Kallen(msq[k],msq[j],σs[i]))
    rest = σs[j]-s-msq[j]
    return (2msq[k]*rest+EE4mksq)/pp4mksq
end

# 
function cosζ(wr::Arg3WignerRotation,σs,msq)
    i,j,k = ijk(wr)
    # 
    msq[k] ≈ 0 && return 1.0
    # 
    s = msq[4]
    EE4m1sq = (σs[i]-msq[j]-msq[k])*(σs[j]-msq[k]-msq[i])
    pp4m1sq = sqrt(Kallen(σs[i],msq[j],msq[k])*Kallen(σs[j],msq[k],msq[i]))
    rest = msq[i]+msq[j]-σs[k]
    return (2msq[k]*rest+EE4m1sq)/pp4m1sq
end
