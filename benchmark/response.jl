# Compare the computation of response between Seis.jl,
# SAC and Obspy

using DSP: unwrap
using Seis
using Plots
using PyCall
using StationXML

struct PolesZerosResponse{T}
    zeros::Vector{Complex{T}}
    poles::Vector{Complex{T}}
    scale::T
end

function PolesZerosResponse(nzeros::Integer, poles, scale)
    zeros = fill(zero(eltype(poles)), nzeros)
    PolesZerosResponse(zeros, poles, scale)
end

Base.broadcastable(paz::PolesZerosResponse) = Ref(paz)

Inv = pyimport("obspy.signal.invsim")
Signal = pyimport("scipy.signal")

"""
    frequency_response(zeros, poles, scale, f) -> resp

Compute the complex response `resp` at a frequency `f` in Hz for
an analog system described by complex `poles` and `zeros`,
plus a real `scale` factor.

The value `resp` or ``G`` is given by Eq. (10.11) of Scherbaum (2001),

``G(\omega) = a_0 \prod_{j=1}^M (i\omega - z_j) / \prod_{k=1}^N (i\omega - p_k)``

where ``\omega`` is the angular frequency ``\omega = 2 \pi f``,
`z_j` is the ``j``th zero, ``p_k`` is the ``k``th pole, and there
are ``M`` zeros and ``N`` poles.  ``a_0`` is the `scale` factor.
"""
function frequency_response(zeros, poles, scale, f)
    ω = 2π*f
    iω = im*ω
    response = scale*prod(z -> iω - z, zeros; init=one(iω)) /
                     prod(p -> iω - p, poles; init=one(iω))
end
frequency_response(paz::PolesZerosResponse, f) =
    frequency_response(paz.zeros, paz.poles, paz.scale, f)

function frequency_response2(zeros, poles, scale, f)
    Z = zeros
    P = poles

    cf = 2π*im*f
    γ = zero(f)
    n = scale*one(complex(f))
    d = zero(complex(f))
    ϵ = eps(f)^2

    for z in Z
        n *= cf - z
    end
    for p in P
        d = conj(cf - p)
        n *= d
        n /= max(ϵ, abs2(d))
    end
    g = abs2(n)
    if g > γ
        γ = g
    end
    n
end
frequency_response2(paz::PolesZerosResponse, f) =
    frequency_response2(paz.zeros, paz.poles, paz.scale, f)

_amp_phase(resp) = abs.(resp), angle.(resp)

function _amp_phase(paz::PolesZerosResponse, freqs)
    resp = frequency_response.(paz, freqs)
    _amp_phase(resp)
end

function bode_plot!(p, paz::PolesZerosResponse, freqs; label="", line=1)
    amp, phase = _amp_phase(paz, freqs)
    plot!(p[1], freqs, amp; label=label, line)
    plot!(p[2], freqs, phase; label=label, line)
    p
end

function bode_plot!(p, freqs, resp; label="", line=1)
    amp, phase = _amp_phase(resp)
    plot!(p[1], freqs, amp; label=label, line)
    plot!(p[2], freqs, phase; label=label, line)
    p
end   

function bode_plot(paz::PolesZerosResponse, freqs; label="", line=1, size=(600,700))
    amp, phase = _amp_phase(paz, freqs)
    xlim = (iszero(first(freqs)) ? freqs[begin+1] : first(freqs), freqs[end])
    p = plot(
        plot(ylabel="Amplitude", xticks=nothing, xlim=xlim,
            yaxis=:log10, legend=:bottomright),
        plot(xlabel="Frequency / Hz", ylabel="Phase", xlim=xlim, legend=false),
        xaxis=:log10, layout=(2,1), framestyle=:box, fontfamily="Helvetica",
        grid=false, size=size
    )
    bode_plot!(p, paz, freqs; label, line)
end


paz_igu16 = PolesZerosResponse(
    2,
    [-22.211059 + 22.217768im, -22.211059 - 22.217768im],
    76.7
)

paz_flat = PolesZerosResponse(ComplexF64[2/(1 - 1im)], ComplexF64[1 + 1im], 1.0)

JSA_file = joinpath(dirname(pathof(StationXML)), "..", "test", "data", "JSA.xml")
paz_JSA = let sxml = StationXML.read(JSA_file),
        cha = sxml.network[1].station[1].channel[1]
    poles_zeros = cha.response.stage[1].poles_zeros
    poles = [complex(pole.real.value, pole.imaginary.value) for pole in poles_zeros.pole]
    zeros = [complex(zero.real.value, zero.imaginary.value) for zero in poles_zeros.zero]
    gain = prod(x -> x.stage_gain.value, cha.response.stage)
    normalization = cha.response.stage[1].poles_zeros.normalization_factor
    gain = gain*normalization
    PolesZerosResponse(zeros, poles, gain)
end


freqs = 10 .^ (-2:0.1:3)
resp = frequency_response.(paz_igu16, freqs)
resp2 = frequency_response2.(paz_igu16, freqs)

# Scipy response
freqs_py, resp_py = let p = paz_igu16, f = freqs
    fs, resp = Signal.freqs_zpk(p.zeros, p.poles, p.scale, 2π.*f)
    fs./2π, resp
end

# Obspy response
resp_obspy, freqs_obspy = let p = paz_igu16, f = freqs
    delta = 1/2maximum(f)
    nfft = 2*length(f)
    Inv.paz_to_freq_resp(p.poles, p.zeros, p.scale,
        delta, nfft, freq=true)
end

p = begin
    p = bode_plot(paz_igu16, freqs, label="Seis", line=5, size=(500,600))
    bode_plot!(p, freqs, resp2, label="Seis 2", line=(4,:dash))
    bode_plot!(p, freqs_py, resp_py, label="Scipy", line=3)
    bode_plot!(p, freqs_obspy, resp_obspy, label="Obspy", line=(2,:dash))
    ylims!(p[1], 1e-4, 100)
    p
end
