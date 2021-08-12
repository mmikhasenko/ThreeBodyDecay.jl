
abstract type IndexWignerRotation end
struct TrivialWR <: IndexWignerRotation
    k::Int
end
struct HatWR{S} <: IndexWignerRotation
    k::Int
end
struct ZetaRepWR{T,S} <: IndexWignerRotation
    k::Int
end
struct ZetaAllWR{S} <: IndexWignerRotation
    k::Int
end


ispositive(wr::TrivialWR) = true
ispositive(wr::HatWR{S}) where S = S
ispositive(wr::ZetaRepWR{T,S}) where {T,S} = S
ispositive(wr::ZetaAllWR{S}) where S = S

sameparticlereference(wr::ZetaRepWR{T,S}) where {T,S} = T==:S

ijk(k::Int) = (k+1,k+2,k) |> x->mod.(x,Ref(Base.OneTo(3)))
ijk(wr::IndexWignerRotation) = ijk(wr.k)

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
    system_a == reference_b && return TrivialWR(particle_c)
    S = issequential(system_a, reference_b)
    A,B = S ? (system_a, reference_b) : (reference_b, system_a)
    # 
    particle_c == 0 && return HatWR{S}(A)
    #
    particle_c ∉ (system_a, reference_b) && return ZetaAllWR{S}(particle_c)
    #
    T = (particle_c == A) ? :S : :D
    return ZetaRepWR{T,!(S)}(particle_c)
end



cosζ(wr::TrivialWR,σs,msq) = 1.0

function cosζ(wr::HatWR,σs,msq)
    i,j,k = ijk(wr)
    # 
    s = msq[4]
    EE4s = (s+msq[i]-σs[i])*(s+msq[k]-σs[k])
    pp4s = sqrt(λ(s,msq[i],σs[i])*λ(s,msq[k],σs[k]))
    rest = σs[j]-msq[i]-msq[k]
    return (EE4s-2s*rest)/pp4s
end

function cosζ(wr::ZetaRepWR,σs,msq)
    i,j,k = ijk(wr)
    # 
    if !(sameparticlereference(wr))
        i,j = j,i 
    end
    # 
    msq[k] ≈ 0 && return 1.0
    # 
    s = msq[4]
    EE4mksq = (s+msq[k]-σs[k])*(σs[i]-msq[k]-msq[j])
    pp4mksq = sqrt(λ(s,msq[k],σs[k])*λ(msq[k],msq[j],σs[i]))
    rest = σs[j]-s-msq[j]
    return (2msq[k]*rest+EE4mksq)/pp4mksq
end

# 
function cosζ(wr::ZetaAllWR,σs,msq)
    i,j,k = ijk(wr)
    # 
    msq[k] ≈ 0 && return 1.0
    # 
    s = msq[4]
    EE4m1sq = (σs[i]-msq[j]-msq[k])*(σs[j]-msq[k]-msq[i])
    pp4m1sq = sqrt(λ(σs[i],msq[j],msq[k])*λ(σs[j],msq[k],msq[i]))
    rest = msq[i]+msq[j]-σs[k]
    return (2msq[k]*rest+EE4m1sq)/pp4m1sq
end
