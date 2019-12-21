### four vectors
################################################################################

function four_vectors_in_binary_decay(cosθ, ϕ;
        m1sq::Real=error("give mass 1 squared"),
        m2sq::Real=error("give mass 2 squared"),
        m0sq::Real=error("give decay mass squared"))
    E1 = (m0sq+m1sq-m2sq)/(2sqrt(m0sq));
    E2 = (m0sq-m1sq+m2sq)/(2sqrt(m0sq));
    p = sqrt(λ(m0sq,m1sq,m2sq))/(2sqrt(m0sq));
    #
    sinθ = sqrt(1-cosθ^2)
    p1 = MVector( p*sinθ*cos(ϕ),  p*sinθ*sin(ϕ),  p*cosθ, E1)
    p2 = MVector(-p*sinθ*cos(ϕ), -p*sinθ*sin(ϕ), -p*cosθ, E2)
    return (p1, p2)
end

function four_vectors_in_binary_decay(p0,cosθ,ϕ;
        m1sq::Real=error("give mass 1 squared"),
        m2sq::Real=error("give mass 2 squared"))
    psq = sum(abs2,p0[1:3])
    m0sq = p0[4]^2-psq
    #
    γ = p0[4]/sqrt(m0sq)
    cosθ0 = p0[3]/sqrt(psq)
    ϕ0 = atan(p0[2],p0[1])
    #
    p1, p2 = four_vectors_in_binary_decay(cosθ,ϕ; m1sq=m1sq, m2sq=m2sq, m0sq=m0sq)
    boostz!(p1,γ); roty_cos!(p1,cosθ0); rotz!(p1, ϕ0)
    boostz!(p2,γ); roty_cos!(p2,cosθ0); rotz!(p2, ϕ0)
    #
    return (p1, p2)
end

invmasssq(p) = p[4]^2-sum(abs2, p[1:3])

function rotz!(p,θ)
    c, s = cos(θ), sin(θ)
    p[1], p[2] = [c -s; s c]*[p[1], p[2]]
    return
end
function roty!(p,θ)
    c, s = cos(θ), sin(θ)
    p[1], p[3] = [c s; -s c]*[p[1], p[3]]
    return
end
function roty_cos!(p,cosθ)
    c, s = cosθ, sqrt(1-cosθ^2)
    p[1], p[3] = [c s; -s c]*[p[1], p[3]]
    return
end
function roty_cos_inv!(p,cosθ)
    c, s = cosθ, sqrt(1-cosθ^2)
    p[1], p[3] = [c -s; s c]*[p[1], p[3]]
    return
end
function boostz!(p,γ)
    γ, βγ = γ*sign(γ), sqrt(γ^2-1)*sign(γ)
    p[3], p[4] = [γ βγ; βγ γ]*[p[3], p[4]]
    return
end
