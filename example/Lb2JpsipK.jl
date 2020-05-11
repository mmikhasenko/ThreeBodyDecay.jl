using ThreeBodyDecay
using Plots


H1(two_τ,two_j,two_L) = CG_doublearg(two_L,0,two_j,two_τ,1,two_τ)
H2(two_μ,two_λ,two_j,two_s,two_l) = CG_doublearg(2,two_μ,1,         -two_λ,two_s,two_μ-two_λ)*
                              CG_doublearg(two_l,0,two_s,two_μ-two_λ,two_j,two_μ-two_λ)
#
I(cosθ,two_j,two_s,two_l,two_L) = sum(abs2,
    H1(two_τ,two_j,two_L)*
    wignerd_doublearg(two_j,two_τ,two_μ-two_λ,cosθ)*
    H2(two_μ,two_λ,two_j,two_s,two_l) for two_μ=-2:2:2, two_λ=-1:2:1, two_τ=(-two_j):2:two_j)

let jp = (3,"-")
    ( :LS=>possibleLS((1,"+"),jp,(0,"-")),
      :ls=>possibleLS(jp,(2,"-"),(1,"+")) )
end

let
    two_j,two_s,two_l,two_L = 3,3,2,2
    plot(cosθ->I(cosθ,two_j,two_s,two_l,two_L),-1,1)
end
############################################################

const Lb2JpK = let mLb = 5.62, mJψ = 3.09, mp=0.938, mK = 0.49367
    ThreeBodySystem(mJψ,mp,mK,mLb)
end

function plot_dalitz_with_projections(f=(σ3,σ1)->1.0)
    # layout = @layout [a{0.8h}; grid(1,2)]
    layout = @layout [a{0.65w,0.7h} b; c _]
    plot(layout=layout, size=(800,600), link=:both)
    #
    σ1v = LinRange(Lb2JpK.mthsq[1],Lb2JpK.sthsq[1], 202)
    σ3v = LinRange(Lb2JpK.mthsq[3],Lb2JpK.sthsq[3], 200)
    cal = [Kibble31(σ3,σ1,Lb2JpK) < 0 ? f(σ3,σ1) : NaN for σ3 in σ3v, σ1 in σ1v]
    heatmap!(σ1v, σ3v, cal, c=:viridis, colorbar=false, sp=1,
        ylab="m[J/psi p] (GeV^2)", xlab="m[p K] (GeV^2)")
    #
    calz = map(z->isnan(z) ? 0.0 : z, cal)
    plot!(sum(calz, dims=2)[:,1], σ3v, lab="", xaxis=false, l=(2,:black), sp=2)
    plot!(σ1v, sum(calz, dims=1)[1,:], lab="", yaxis=false, l=(2,:black), sp=3)
end


############################################################

function A(σ3,σ1, CS, Cs)
    σ = [σ1, 0.0, σ3]
    return sum(c*amp(σ[ch[1]], ch[2]) for (ch,c) in zip(CS,Cs))
end

I(σ3,σ1, CS, Cs) = abs2(A(σ3,σ1, CS, Cs))

Λ1405  = BreitWigner(1.405,   0.090)
Λ1520  = BreitWigner(1.5195,  0.0156)
Λ1690  = BreitWigner(1.685,   0.050)
Λ1810  = BreitWigner(1.80,    0.090)

model = [(1,Λ1405),
 (1,Λ1520),
 (1,BreitWigner(1.6,0.2)),
 (1,Λ1690),
 (1,Λ1810),
 (3,BreitWigner(4.45,0.06))]


σ3v, σ1v = flatDalitzPlotSample31(tbs; Nev=1000)


I(σ3v[1],σ1v[1],
    model,
    [0.9, 1.0, 0.7, 0.2, 0.3, 0.2im])

# let
#     plot_dalitz_with_projections((σ3,σ1)->I(σ3,σ1,
#         model,
#         [0.9, 1.0, 0.7, 0.2, 0.3, 0.2im]))
# end
