using Test
using ThreeBodyDecay
using Statistics

using Random
Random.seed!(1234)

ms = ThreeBodyMasses(1.1, 1.2, 1.3; m0=4.2)

let
    σs = randomPoint(ms)
    # 
    x, y = squaredalitz1(σs, ms)
    @test -1 ≤ x ≤ 1
    @test -1 ≤ y ≤ 1
    x, y = squaredalitz2(σs, ms)
    @test -1 ≤ x ≤ 1
    @test -1 ≤ y ≤ 1
    x, y = squaredalitz3(σs, ms)
    @test -1 ≤ x ≤ 1
    @test -1 ≤ y ≤ 1
end

let
    x, y = rand(2)
    # 
    σs = invsquaredalitz1(x, y, ms)
    @test prod((x, y) .≈ squaredalitz1(σs, ms))
    # 
    σs = invsquaredalitz2(x, y, ms)
    @test prod((x, y) .≈ squaredalitz2(σs, ms))
    # 
    σs = invsquaredalitz3(x, y, ms)
    @test prod((x, y) .≈ squaredalitz3(σs, ms))
end


let Nd = 100_000, Nb = 10
    ms = ThreeBodyMasses(1.1, 1.2, 1.3; m0=4.2)

    xv, yv = rand(Nd), rand(Nd)

    σsv = [invsquaredalitz1(x, y, ms) for (x, y) in zip(xv, yv)]
    ws = [jacobean_squaredalitz1(σs, ms) for σs in σsv]


    # this should be consistent with flat
    # histogram2d(getproperty.(σsv,:σ1), getproperty.(σsv,:σ3), weights=ws, bins=30)

    # test in the grid: count cells that are entirely inside of Dalitz
    σ1v = range(lims1(ms)..., length=Nb)
    σ3v = range(lims1(ms)..., length=Nb)
    samples = [sum(x -> (r1[1] < x[1].σ1 < r1[2] && r3[1] < x[1].σ3 < r3[2]) * x[2], zip(σsv, ws))
               for (r1, r3) in Iterators.product(zip(σ1v[1:end-1], σ1v[2:end]), zip(σ3v[1:end-1], σ3v[2:end]))
               if (Kibble(Invariants(ms; σ1=r1[1], σ3=r3[1]), ms^2) < 0) *
               (Kibble(Invariants(ms; σ1=r1[2], σ3=r3[1]), ms^2) < 0) *
               (Kibble(Invariants(ms; σ1=r1[2], σ3=r3[2]), ms^2) < 0) *
               (Kibble(Invariants(ms; σ1=r1[1], σ3=r3[2]), ms^2) < 0)]
    nev = vcat(samples...)
    @test mean(nev) > 10 * sqrt(cov(nev))
end
