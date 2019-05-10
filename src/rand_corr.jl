
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

xlim(binned2dDensity) = (binned2dDensity.grid[1,1][1],binned2dDensity.grid[end,end][1])
ylim(binned2dDensity) = (binned2dDensity.grid[1,1][2],binned2dDensity.grid[end,end][2])

function getbinned2dDensity(g, xlim, ylim,  Nrows, Ncols)
    grid = hcat([[[si,σi] for σi=LinRange(xlim[1], xlim[2], Nrows)] for si=LinRange(ylim[1], ylim[2], Ncols)]...)
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
