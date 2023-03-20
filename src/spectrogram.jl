"""
    spectrogram(trace::AbstractTrace; length=(endtime(t) - starttime(t))/20, overlap=0.5, window=nothing, pad=nothing) -> spec

Calculate the spectrogram for the data in `trace`.  Each segment is `length` s
long, and overlaps with the next by a fraction of `overlap` (i.e., `overlap`
must be 0 s or greater, and less than 1).  Note that these values are rounded
to an integer number of samples in both cases.

Windows are by default padded to the next length which allows for fast FFT
computation.  Setting `pad` to a number larger than 1 pads the trace with
zeroes to the nearest integer length which is `pad` times the trace length.
This allows for finer sampling in frequency space.

By default, no windowing of each segment (of `length` s) is performed;
pass a windowing function or vector of amplitudes to `window` to window each
segment before the periodogram is computed.  See `DSP.periodogram`
for more details of the kind of windowing which is possible.  (The spectrogram
calculation is performed by `DSP.spectrogram` and returns a
`DSP.Periodograms.Spectrogram` object.)

# Examples
```
julia> t = sample_data();

julia> spec = spectrogram(t);

julia> spec.time
52.90999984741211:0.25:62.40999984741211
```
"""
function DSP.spectrogram(t::AbstractTrace;
        length=max((endtime(t) - starttime(t))/20, 2t.delta), overlap=0.5, window=nothing,
        pad=nothing)
    length > 0 || throw(ArgumentError("`length` must be greater than 0 s"))
    length >= t.delta || throw(ArgumentError(
        "`length` must at least the sampling interval"))
    0 <= overlap < 1 || throw(ArgumentError(
        "`overlap` must be greater than or equal to 0, and less than 1"))
    pad === nothing || 1 <= pad || throw(ArgumentError("`pad` must be 1 or greater"))
    delta = t.delta
    n = round(Int, length/delta)
    noverlap = floor(Int, overlap*n)
    nfft = (pad === nothing || pad == 1) ? DSP.nextfastfft(n) : round(Int, pad*n)
    # Compute spectrogram
    spec = DSP.spectrogram(trace(t), n, noverlap;
        fs=1/delta, window=window, nfft=nfft)
    # Update times in spectrogram to reflect sample times of trace
    typeof(spec)(spec.power, spec.freq, spec.time .+ starttime(t))
end
