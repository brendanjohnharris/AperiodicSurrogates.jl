using AperiodicSurrogates
using Test
using TimeseriesSurrogates, CairoMakie
using DelimitedFiles

@testset "AperiodicSurrogates.jl" begin
    x = readdlm(joinpath(@__DIR__, "test.csv")) |> vec
    fs = 1250.0
    s = surrogate(x, AP(fs))
    fig = surroplot(x, s; resolution=(720, 1080), cx=:cornflowerblue, cs=:crimson)
    xlims!(fig.content[1], (0, 2000))
    save(joinpath(@__DIR__, "test.png"), fig)
end
@testset "Fourier Surrogates" begin
    x = readdlm(joinpath(@__DIR__, "test.csv")) |> vec
    s = surrogate(x, FT())
    fig = surroplot(x, s; resolution=(720, 1080), cx=:cornflowerblue, cs=:crimson)
    xlims!(fig.content[1], (0, 2000))
    save(joinpath(@__DIR__, "test_FT.png"), fig)
end
