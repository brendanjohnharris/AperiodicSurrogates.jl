using AperiodicSurrogates
using Documenter

DocMeta.setdocmeta!(AperiodicSurrogates, :DocTestSetup, :(using AperiodicSurrogates); recursive=true)

makedocs(;
    modules=[AperiodicSurrogates],
    authors="brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
    repo="https://github.com/brendanjohnharris/AperiodicSurrogates.jl/blob/{commit}{path}#{line}",
    sitename="AperiodicSurrogates.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
