
struct binned1dDensity
    grid::Array{Float64,1}
    cumarr::Array{Float64}
    density
end

function getbinned1dDensity(g, lim, Nbins)
    grid = [σi for σi=LinRange(lim[1], lim[2], Nbins)]
    #
    weights = [g(v...) for v in (grid[2:end] .+ grid[1:end-1]) ./ 2]
    weights ./= sum(weights)
    cumarr = [0;cumsum(vcat(weights...), dims=1)]
    return binned1dDensity(grid,cumarr,g)
end

import Base:rand
function rand(bD::binned1dDensity)
    # get
    binind = findfirst(bD.cumarr .> rand())-1
    σl, σr = bD.grid[binind], bD.grid[binind+1]
    σ = σl+rand()*(σr-σl)
    bD.density(σ) == 0.0 && return rand(bD)
    return σ
end

#################

struct binned2dDensity
    grid::Matrix{Array{Float64,1}}
    cumarr::Array{Float64}
    density
end

arg1_lims(binned2dDensity) = (binned2dDensity.grid[1,1][1],binned2dDensity.grid[end,end][1])
arg2_lims(binned2dDensity) = (binned2dDensity.grid[1,1][2],binned2dDensity.grid[end,end][2])

function getbinned2dDensity(g, arg1_lims, arg2_lims,  N1, N2)
    grid = hcat([[[x1,x2]
        for x2=range(arg2_lims[1], arg2_lims[2], length=N1)]
        for x1=range(arg1_lims[1], arg1_lims[2], length=N2)]...)
    #
    weights = [g(v...) for v in (grid[2:end,2:end] .+ grid[1:end-1,1:end-1]) ./ 2]
    # normalize by the cell size
    for i=1:size(grid,1)-1, j=1:size(grid,2)-1
        weights[i,j] *= (grid[i+1,j][2] - grid[i,j][2])*(grid[i,j+1][1] - grid[i,j][1])
    end
    weights ./= sum(weights)
    cumarr = [0;cumsum(vcat(weights...), dims=1)]
    return binned2dDensity(grid,cumarr, g)
end

import Base:rand
function rand(bD::binned2dDensity)
    # get
    binind = findfirst(bD.cumarr .> rand())-2
    # pars back
    Nrows, Ncols = size(bD.grid)
    indCol = div(binind, Nrows-1) + 1
    indRow = mod(binind, Nrows-1) + 1
    indRow+1 > Nrows && error("indRow+1 > Nrows: binind = $binind")
    indCol+1 > Ncols && error("indCol+1 > Ncols: binind = $binind")
    s1, σ1 = bD.grid[indRow,indCol]
    s2, σ2 = bD.grid[indRow+1,indCol+1]
#         s1 != s2 && error("Something is wrong!")
    s = s1+rand()*(s2-s1)
    σ = σ1+rand()*(σ2-σ1)
    bD.density(s,σ) == 0.0 && return rand(bD)
    return [s, σ]
end

function gridded_density_function(x, bD::binned2dDensity)
    (s,σ) = x
    (bD.grid[1,1][1] > s || s > bD.grid[1,end][1]) && return 0.0
    (bD.grid[1,1][2] > σ || σ > bD.grid[end,1][2]) && return 0.0
    #
    mp_s = map(x->x[1], bD.grid[1,:])
    i0 = findfirst(mp_s .> s)
    s0 = (mp_s[i0]+mp_s[i0-1])/2
    #
    mp_σ = map(x->x[2], bD.grid[:,1])
    j0 = findfirst(mp_σ .> σ)
    σ0 = (mp_σ[j0]+mp_σ[j0-1])/2
    #
    return bD.density(s0,σ0)!=0.0 ?  bD.density(s0,σ0) : 0.0;
end
