module AperiodicSurrogates
using Reexport
using TimeseriesSurrogates
using TimeseriesSurrogates.StatsBase
using TimeseriesSurrogates.AbstractFFTs
using IntervalSets
using PyFOOOF
using PyFOOOF.PyCall

@reexport using TimeseriesSurrogates
export Aperiodic, AP

function fooof(f, res::AbstractVector; freqrange=[10.0, 300.0])
    freqrange = py"[$(first(freqrange)), $(last(freqrange))]"o
    f = collect(f)
    spectrum = collect(res)
    fm = PyFOOOF.FOOOF(peak_width_limits=py"[0.5, 20.0]"o, max_n_peaks=4, aperiodic_mode="knee") # Maybe relax these in the future?
    fm.add_data(f, spectrum, freqrange)
    fm.fit()
    return fm
end

function aperiodicfit(f, res::AbstractVector; kwargs...)
    fm = fooof(f, res; kwargs...)
    # * The aperiodic model, as described in doi.org/10.1038/s41593-020-00744-x
    b, k, χ = fm.aperiodic_params_
    L = f -> 10.0.^(b - log10(k + (f)^χ))
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
Aperiodic(fs) = Aperiodic(fs, (10, 300))
Aperiodic(fs, freqrange::AbstractInterval, args...) = Aperiodic(fs, freqrange|>extrema)
const AP = Aperiodic

function TimeseriesSurrogates.surrogenerator(x::AbstractVector, ap::Aperiodic, rng=Random.default_rng())
    # * First, measure the power spectrum of the input signal
    m = mean(x)
    forward = plan_rfft(x)
    f = rfftfreq(length(x), ap.fs)
    𝓕 = forward*(x .- m)
    inverse = plan_irfft(𝓕, length(x))
    shuffled𝓕 = zero(𝓕)
    s = similar(x)
    n = length(𝓕)
    r = abs.(𝓕)
    coeffs = zero(r)

    # * Then, fit the aperiodic model to the power spectrum and calculate its Fourier coefficients
    fm = aperiodicfit(f, r.^2, freqrange=ap.freqrange)
    r = f .|> fm .|> sqrt

    init = (inverse = inverse, m = m, coeffs = coeffs, n = n, r = r,
            ϕ = ϕ, shuffled𝓕 = shuffled𝓕)
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