"""
    x,y = squaredalitz(k,σs,ms)

    Calculate normalized values in square coordinates, -1 ≤ x,y ≤ 1
"""
function squaredalitz(k, σs, ms::MassTuple)
    cθ = cosθij(k, σs, ms^2)
    cθ = cθ ≥ 1.0 ? 1.0 : (cθ ≤ -1.0 ? -1.0 : cθ)
    y = acos(cθ) / π
    mlimsk = sqrt.(lims(k, ms))
    xn = 2 * (sqrt(σs[k]) - mlimsk[1]) / (mlimsk[2] - mlimsk[1]) - 1
    x = acos(xn) / π
    return (x, y)
end

squaredalitz1(σs, ms) = squaredalitz(1, σs, ms)
squaredalitz2(σs, ms) = squaredalitz(2, σs, ms)
squaredalitz3(σs, ms) = squaredalitz(3, σs, ms)

"""
    invsquaredalitz(k,x,y,ms)

    Inverse transformation from square coordinates
"""
function invsquaredalitz(k, x, y, ms::MassTuple)
    mlimsk = sqrt.(lims(k, ms))
    # 
    _mk = mlimsk[1] + (mlimsk[2] - mlimsk[1]) * (cos(x * π) + 1) / 2
    _σk = _mk^2
    _cosθij = cos(y * π)
    # 
    k == 1 && return Invariants(ms, σ1=_σk, σ3=σ3of1(_cosθij, _σk, ms^2))
    k == 2 && return Invariants(ms, σ2=_σk, σ1=σ1of2(_cosθij, _σk, ms^2))
    k != 3 && error("k should be 1,2, or 3, while it is $k")
    return Invariants(ms, σ3=_σk, σ2=σ2of3(_cosθij, _σk, ms^2))
end

invsquaredalitz1(x, y, ms) = invsquaredalitz(1, x, y, ms)
invsquaredalitz2(x, y, ms) = invsquaredalitz(2, x, y, ms)
invsquaredalitz3(x, y, ms) = invsquaredalitz(3, x, y, ms)

"""
    jacobian_squaredalitz(k,σs,ms)

calculates jacobian of transformation to square coordinates
"""
function jacobean_squaredalitz(k, σs, ms::MassTuple)
    (i, j) = ij_from_k(k)
    #
    cθ = cosθij(k, σs, ms^2)
    cθ = cθ ≥ 1.0 ? 1.0 : (cθ ≤ -1.0 ? -1.0 : cθ)
    dydcθ = 1 / sqrt(1 - cθ^2) / π
    #
    m = sqrt(σs[k])
    sqrtlimsk = sqrt.(lims(k, ms))
    xn = 2 * (m - sqrtlimsk[1]) / (sqrtlimsk[2] - sqrtlimsk[1]) - 1
    dxdm = 2 / (sqrtlimsk[2] - sqrtlimsk[1]) / sqrt(1 - xn^2) / π
    ρρm = sqrt(Kallen(σs[k], ms[i]^2, ms[j]^2) * Kallen(ms[4]^2, σs[k], ms[k]^2)) / (σs[k] * ms[4]^2) * m
    return ρρm / dydcθ / dxdm #
end

jacobean_squaredalitz1(σs, ms) = jacobean_squaredalitz(1, σs, ms)
jacobean_squaredalitz2(σs, ms) = jacobean_squaredalitz(2, σs, ms)
jacobean_squaredalitz3(σs, ms) = jacobean_squaredalitz(3, σs, ms)
