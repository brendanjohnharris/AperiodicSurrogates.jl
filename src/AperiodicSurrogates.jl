module AperiodicSurrogates
using Reexport
using TimeseriesSurrogates
using StatsBase
using TimeseriesSurrogates.AbstractFFTs
using TimeseriesSurrogates.Distributions
using IntervalSets
using PyFOOOF
using PythonCall

@reexport using TimeseriesSurrogates
export Aperiodic, AP

function fooof(f, res::AbstractVector; freqrange=[1.0, 300.0])
    freqrange = pylist([(first(freqrange)), (last(freqrange))])
    f = collect(f)
    spectrum = collect(res)
    fm = PyFOOOF.FOOOF(peak_width_limits=pylist([0.5, 20.0]), max_n_peaks=4, aperiodic_mode="knee") # Maybe relax these in the future?
    fm.add_data(Py(f).to_numpy(), Py(spectrum).to_numpy(), freqrange)
    fm.fit()
    # if fm.aperiodic_params_[2] < 0.0
    #     fm = PyFOOOF.FOOOF(peak_width_limits=pylist([0.5, 20.0]), max_n_peaks=4, aperiodic_mode="fixed")
    #     fm.add_data(f, spectrum, freqrange)
    #     fm.fit()
    # end
    return fm
end

function aperiodicfit(f, res::AbstractVector; kwargs...)
    fm = fooof(f, res; kwargs...)
    # * The aperiodic model, as described in doi.org/10.1038/s41593-020-00744-x
    b, k, χ = pyconvert.((Float64,), fm.aperiodic_params_)
    k = max(k, 0.01)
    # b, k, χ = length(ps) ==3 ? ps : (ps[1], 0.0, ps[2])
    L = f -> 10.0 .^ (b - log10(k + (f)^χ))
end

"""
    Aperiodic(freqrange=[10, 300])
Produce a surrogate with the random phases and the same power spectrum as the aperiodic FOOOF fit of the input time series.
`fs` is the sampling frequency of the input time series.
`freqrange` determines the range of frequency values to use for the FOOOF fit.
"""
struct Aperiodic <: Surrogate
    fs::Number
    freqrange
end
Aperiodic(fs) = Aperiodic(fs, (1, 300))
Aperiodic(fs, freqrange::AbstractInterval, args...) = Aperiodic(fs, freqrange |> extrema)
const AP = Aperiodic

function TimeseriesSurrogates.surrogenerator(x::AbstractVector, ap::Aperiodic, rng=Random.default_rng())
    # * First, measure the power spectrum of the input signal
    m = mean(x)
    forward = plan_rfft(x)
    f = rfftfreq(length(x), ap.fs)
    𝓕 = forward * (x .- m)
    inverse = plan_irfft(𝓕, length(x))
    shuffled𝓕 = zero(𝓕)
    s = similar(x)
    n = length(𝓕)
    _r = abs.(𝓕)
    coeffs = zero(_r)

    # * Then, fit the aperiodic model to the power spectrum and calculate its Fourier coefficients
    fm = aperiodicfit(f, _r .^ 2, freqrange=ap.freqrange)
    r = f .|> fm .|> sqrt
    r[f.<minimum(ap.freqrange)] = _r[f.<minimum(ap.freqrange)] # Ensures no funky behaviour outside of the fit bounds

    # * Scale the coefficients so the surrogate has the same energy as the original time series
    r = r .* sum(_r) ./ sum(r)

    init = (; inverse, m, coeffs, n, r, shuffled𝓕)
    return TimeseriesSurrogates.SurrogateGenerator(ap, x, s, init, rng)
end

function (sg::TimeseriesSurrogates.SurrogateGenerator{<:Aperiodic})()
    inverse, m, coeffs, n, r, shuffled𝓕 =
        getfield.(Ref(sg.init),
            (:inverse, :m, :coeffs, :n, :r, :shuffled𝓕))
    s, rng = sg.s, sg.rng

    coeffs .= rand(rng, Uniform(0, 2π), n)
    shuffled𝓕 .= r .* exp.(coeffs .* 1im)
    s .= inverse * shuffled𝓕 .+ m
    return s
end


end
