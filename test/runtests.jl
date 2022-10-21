using AperiodicSurrogates
using Test
using TimeseriesSurrogates, CairoMakie
using DelimitedFiles

@testset "AperiodicSurrogates.jl" begin
    x = readdlm("test/test.csv") |> vec
    fs = 1250.0
    s = surrogate(x, AP(fs))
    fig = surroplot(x, s; resolution=(720, 1080), cx=:cornflowerblue, cs=:crimson)
    xlims!(fig.content[1], (0, 2000))
    xlims!(fig.content[3], (0, 0.2))
    ylims!(fig.content[3], (1e-9, nothing))
    save("test/test.png", fig)
end
