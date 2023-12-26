
Kallen(x, y, z) = x^2 + y^2 + z^2 - 2x * y - 2y * z - 2z * x
sqrtKallenFact(a, b, c) = sqrt(a - (b + c)) * sqrt(a - (b - c)) * sqrt(a + (b + c)) * sqrt(a + (b - c))
Kibble(σs, msq) = Kallen(Kallen(msq[4], msq[1], σs[1]),
    Kallen(msq[4], msq[2], σs[2]),
    Kallen(msq[4], msq[3], σs[3]))
#
"""
    σjofk(k,z,σi,msq)

Computes invariant σj = (p0-pj)² from
the scattering angle z=cosθij in the rest from of (i,j),
given the mass of the system m(i,j)² = σk

Explicit forms: `σ3of1`, `σ1of2`, `σ2of3`.

See also `σkofj(k,z,σj,msq)`
"""
function σjofk(k, z, σk, msq)
    i, j, _ = ijk(k)
    #
    s = msq[4]
    # σj = (p0-pj)² in the rest frame of (i,j)
    EE4σ = (σk + msq[j] - msq[i]) * (σk + s - msq[k])
    p²q²4σ = Kallen(σk, msq[i], msq[j]) * Kallen(s, σk, msq[k])
    p²q²4σ = (p²q²4σ < 0) ? 0.0 : p²q²4σ # for numerical errors
    σi = s + msq[j] - (EE4σ - sqrt(p²q²4σ) * z) / (2σk)
    return σi
end
σ3of1(z, σ1, msq) = σjofk(1, z, σ1, msq)
σ1of2(z, σ2, msq) = σjofk(2, z, σ2, msq)
σ2of3(z, σ3, msq) = σjofk(3, z, σ3, msq)


"""
    σiofk(k,z,σj,msq)

Computes invariant σi = (p0 - pi)² from
the scattering angle z=cosθij in the rest from of (i,j),
given the mass of the system m(i,j)² = σk

Explicit forms: `σ3of2`, `σ1of3`, `σ2of1`.
 
See also `σjofj(k,z,σk,msq)`
"""
function σiofk(k, z, σk, msq)
    σj = σjofk(k, z, σk, msq)
    sum(msq) - σj - σk
end
σ3of2(z, σ2, msq) = σiofk(2, z, σ2, msq)
σ1of3(z, σ3, msq) = σiofk(3, z, σ3, msq)
σ2of1(z, σ1, msq) = σiofk(1, z, σ1, msq)


# Scattering angle

"""
    cosθij(k,σs,msq)

Isobar decay angle for the chain-k, i.e. 
an angle of between vectors pi and -pk in the (ij) rest frame.

Explicit forms: `cosθ23`, `cosθ31`, `cosθ12`.
"""
function cosθij(k, σs, msq)
    (i, j) = ij_from_k(k)
    #
    s = msq[4]
    EE4σ = (σs[k] + msq[i] - msq[j]) * (s - σs[k] - msq[k])
    pp4σ = sqrt(Kallen(σs[k], msq[i], msq[j]) * Kallen(s, σs[k], msq[k]))
    rest = σs[j] - msq[k] - msq[i]
    return (2σs[k] * rest - EE4σ) / pp4σ
end
cosθ23(σs, msq) = cosθij(1, σs, msq)
cosθ31(σs, msq) = cosθij(2, σs, msq)
cosθ12(σs, msq) = cosθij(3, σs, msq)
