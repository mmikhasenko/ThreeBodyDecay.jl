
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

WignerRotation{N}(k::Int, ispositive::Bool) where {N} =
    WignerRotation{N}(k, ispositive, true)

ispositive(wr::TriavialWignerRotation) = true
ispositive(wr::WignerRotation) = wr.ispositive

import Base: iseven
iseven(wr::TriavialWignerRotation) = true
iseven(wr::WignerRotation{0}) = true
iseven(wr::WignerRotation{3}) = true
iseven(wr::WignerRotation{2}) = wr.iseven


ijk(k::Int) = (k + 1, k + 2, k) |> x -> mod.(x, Ref(Base.OneTo(3)))
ijk(wr::AbstractWignerRotation) = ijk(wr.k)
# 
ij_from_k(k::Int) = ijk(k)

issequential(i, j) = (j - i) ∈ (1, -2)


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
    A, B = S ? (system_a, reference_b) : (reference_b, system_a)
    # 
    particle_c == 0 && return Arg0WignerRotation(A, S)
    #
    particle_c ∉ (system_a, reference_b) && return Arg3WignerRotation(particle_c, S)
    #
    T = (particle_c == A)
    return Arg2WignerRotation(particle_c, !(S), T)
end



cosζ(wr::TriavialWignerRotation, σs, msq) = one(σs[1])

function cosζ(wr::Arg0WignerRotation, σs, msq)
    i, j, k = ijk(wr)
    # 
    s = msq[4]
    EE4s = (s + msq[i] - σs[i]) * (s + msq[k] - σs[k])
    pp4s = sqrt(Kallen(s, msq[i], σs[i]) * Kallen(s, msq[k], σs[k]))
    rest = σs[j] - msq[i] - msq[k]
    return (EE4s - 2s * rest) / pp4s
end

function cosζ(wr::Arg2WignerRotation, σs, msq)
    i, j, k = ijk(wr)
    # 
    if !(iseven(wr))
        i, j = j, i
    end
    # 
    msq[k] ≈ 0 && return one(σs[1])
    # 
    s = msq[4]
    EE4mksq = (s + msq[k] - σs[k]) * (σs[i] - msq[k] - msq[j])
    pp4mksq = sqrt(Kallen(s, msq[k], σs[k]) * Kallen(msq[k], msq[j], σs[i]))
    rest = σs[j] - s - msq[j]
    return (2msq[k] * rest + EE4mksq) / pp4mksq
end

# 
function cosζ(wr::Arg3WignerRotation, σs, msq)
    i, j, k = ijk(wr)
    # 
    msq[k] ≈ 0 && return one(σs[1])
    # 
    s = msq[4]
    EE4m1sq = (σs[i] - msq[j] - msq[k]) * (σs[j] - msq[k] - msq[i])
    pp4m1sq = sqrt(Kallen(σs[i], msq[j], msq[k]) * Kallen(σs[j], msq[k], msq[i]))
    rest = msq[i] + msq[j] - σs[k]
    return (2msq[k] * rest + EE4m1sq) / pp4m1sq
end

# explicit
cosζ21_for1(σs, ms²) = cosζ(wr(2, 1, 1), σs, ms²)
cosζ21_for2(σs, ms²) = cosζ(wr(2, 1, 2), σs, ms²)
cosζ13_for1(σs, ms²) = cosζ(wr(1, 3, 1), σs, ms²)
cosζ13_for3(σs, ms²) = cosζ(wr(1, 3, 3), σs, ms²)
cosζ32_for3(σs, ms²) = cosζ(wr(3, 2, 3), σs, ms²)
cosζ32_for2(σs, ms²) = cosζ(wr(3, 2, 2), σs, ms²)

cosζ12_for3(σs, ms²) = cosζ(wr(1, 2, 3), σs, ms²)
cosζ23_for1(σs, ms²) = cosζ(wr(2, 3, 1), σs, ms²)
cosζ31_for2(σs, ms²) = cosζ(wr(3, 1, 2), σs, ms²)

cosζ12_for0(σs, ms²) = cosζ(wr(1, 2, 0), σs, ms²)
cosζ23_for0(σs, ms²) = cosζ(wr(2, 3, 0), σs, ms²)
cosζ31_for0(σs, ms²) = cosζ(wr(3, 1, 0), σs, ms²)

# 
cosζk1_for1(k, σs, ms²) = cosζ(wr(k, 1, 1), σs, ms²)
cosζk2_for2(k, σs, ms²) = cosζ(wr(k, 2, 2), σs, ms²)
cosζk3_for3(k, σs, ms²) = cosζ(wr(k, 3, 3), σs, ms²)

cosθhatk1(k, σs, ms²) = cosζ(wr(k, 1, 0), σs, ms²)
cosθhatk2(k, σs, ms²) = cosζ(wr(k, 2, 0), σs, ms²)
cosθhatk3(k, σs, ms²) = cosζ(wr(k, 3, 0), σs, ms²)


"""
    Phase for wigner d-functions for clockwise rotations
"""
phase(two_λ1_minus_λ2) = (abs(two_λ1_minus_λ2) % 4 == 2 ? -one(two_λ1_minus_λ2) : one(two_λ1_minus_λ2))
phase(two_λ1, two_λ2) = phase(two_λ1 - two_λ2)

