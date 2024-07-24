using DSP
using Makie, GLMakie
GLMakie.activate!()

function make_signal(;
    delta = 0.1,
    npts = 1_000_000,
    fs = (0.000_1, 0.001, 0.01, 0.1),
    ϕs = (deg2rad(-90), 0, deg2rad(30), deg2rad(60)),
    amps = (0.1, 2, 1, 4)
)
    sample_times = (0:(npts - 1)).*delta
    long_signal = zeros(Float32, npts)
    for i in eachindex(long_signal, sample_times)
        t = sample_times[i]
        for (f, ϕ, amp) in zip(fs, ϕs, amps)
            ω = 2π*f
            long_signal[i] += amp*sin(ω*t + ϕ)
        end
    end
    sample_times, long_signal
end

begin
    resample_rate = 0.0121

    long_times, long_signal = make_signal(; npts=100_000)
    npts = length(long_signal)

    fig, ax, _ = Makie.lines(long_times, long_signal; label="Input",
        axis=(; limits=(9980, 10010, nothing, nothing)))

    # Filter in one go
    resampled_long = DSP.resample(long_signal, resample_rate)
    resampled_times = (0:(length(resampled_long) - 1)).*(step(long_times)/resample_rate)

    Makie.lines!(resampled_times, resampled_long, label="All in one")

    # Chunk it up
    nchunks = 13
    chunk_len = floor(Int, npts/nchunks)
    short_signals = Iterators.partition(long_signal, chunk_len)
    short_times = Iterators.partition(long_times, chunk_len)

    # Filter statefully
    filter = DSP.resample_filter(resample_rate)
    firfilter = DSP.FIRFilter(filter, resample_rate)
    τ = DSP.timedelay(firfilter)
    DSP.setphase!(firfilter, τ)
    resampled_chunked = reduce(vcat, map(enumerate(short_signals)) do (ichunk, chunk)
        outlen = outputlength(firfilter, length(chunk))
        if outlen < 0
            error("chunks too short to downsample at rate $resample_rate (buffer length is $outlen points)")
        end

        if ichunk == length(short_signals)

            outLen = ceil(Int, length(chunk)*resample_rate)
            reqInlen = inputlength(firfilter, outLen)
            reqZerosLen = reqInlen - length(chunk)
            chunk_padded = [chunk; zeros(eltype(chunk), reqZerosLen)]

            @show τ
            chunk_padded = [chunk; zeros(eltype(chunk), ceil(Int, τ))]
            DSP.filt(firfilter, chunk_padded)
        else
            DSP.filt(firfilter, chunk)
        end
    end)
    resampled_times_short = (0:(length(resampled_chunked) - 1)).*(step(long_times)/resample_rate)

    length(resampled_long) == length(resampled_chunked) || @warn("Trace lengths not the same")

    Makie.lines!(resampled_times_short, resampled_chunked, label="Chunked", linestyle=:dash)
    Makie.vlines!([first.(short_times); last(long_times)])
    Makie.axislegend(ax)
    display(fig)
end
