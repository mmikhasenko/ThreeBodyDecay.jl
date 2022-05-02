
Kallen(x,y,z) = x^2+y^2+z^2 - 2x*y - 2y*z - 2z*x
KallenFact(a²,b²,c²) = ((a,b,c)=sqrt.((a²,b²,c²)); (a-(b+c))*(a-(b-c))*(a+(b+c))*(a+(b-c)))
Kibble(σs,msq) = Kallen(Kallen(msq[4],msq[1],σs[1]),
                        Kallen(msq[4],msq[2],σs[2]),
                        Kallen(msq[4],msq[3],σs[3]))
#
"""
    σkofi(k,z,σi,msq)

Computes invariant σk = (pi + pj)² from
the scattering angle z=cosθjk in the rest from of (j,k),
given the mass of the system m(j,k)² = σi

Explicit forms: `σ3of1`, `σ1of2`, `σ2of3`.

See also `σkofj(k,z,σj,msq)`
"""
function σkofi(k,z,σi,msq)
    i,j,_ = ijk(k)
    #
    s = msq[4]
    EE4σ = (σi+msq[j]-msq[k])*(s-σi-msq[i])
    p²q²4σ = Kallen(σi,msq[j],msq[k])*Kallen(s,σi,msq[i])
    p²q²4σ = (p²q²4σ < 0) ? 0.0 : p²q²4σ # for numerical errors
    σk = msq[i]+msq[j] + (EE4σ + sqrt(p²q²4σ)*z) / (2σi)
    return σk
end
σ3of1(z,σ1,msq) = σkofi(3,z,σ1,msq)
σ1of2(z,σ2,msq) = σkofi(1,z,σ2,msq)
σ2of3(z,σ3,msq) = σkofi(2,z,σ3,msq)

"""
    σkofj(k,z,σj,msq)

Computes invariant σk = (pi + pj)² from
the scattering angle z=cosθki in the rest from of (k,i),
given the mass of the system m(k,i)² = σj

Explicit forms: `σ3of2`, `σ1of3`, `σ2of1`.
 
See also `σkofi(k,z,σi,msq)`
"""
function σkofj(k,z,σj,msq)
    σi = σkofi(k,z,σj,msq)
    sum(msq)-σi-σj
end
σ3of2(z,σ2,msq) = σkofj(3,z,σ2,msq)
σ1of3(z,σ3,msq) = σkofj(1,z,σ3,msq)
σ2of1(z,σ1,msq) = σkofj(2,z,σ1,msq)

# Scattering angle

"""
    cosθij(k,σs,msq)

Isobar decay angle for the chain-k, i.e. 
an angle of between vectors pi and -pk in the (ij) rest frame.

Explicit forms: `cosθ23`, `cosθ31`, `cosθ12`.
"""
function cosθij(k,σs,msq)
    (i,j) = ij_from_k(k)
    #
    s = msq[4]
    EE4σ = (σs[k]+msq[i]-msq[j])*(s-σs[k]-msq[k])
    pp4σ = sqrt(Kallen(σs[k],msq[i],msq[j])*Kallen(s,σs[k],msq[k]))
    rest = σs[j]-msq[k]-msq[i]
    return (2σs[k]*rest-EE4σ)/pp4σ
end
cosθ23(σs,msq) = cosθij(1,σs,msq)
cosθ31(σs,msq) = cosθij(2,σs,msq)
cosθ12(σs,msq) = cosθij(3,σs,msq)
