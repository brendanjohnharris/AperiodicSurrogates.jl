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
    xlims!(fig.content[3], (0, 0.2))
    ylims!(fig.content[3], (1e-9, nothing))
    save(joinpath(@__DIR__, "test.png"), fig)
end
