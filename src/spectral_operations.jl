# Operations on Fourier-domain traces

"""
    differentiate!(t::AbstractFourierTrace) -> t

Differentiate the frequency-domain trace `t` in the Fourier domain by
multiplying the cofficients ``Y_k`` by ``2\\pi f_k i k``, where ``k <= N/2``
and ``N`` is the number of samples in the time-domain trace `ifft(t)`
and ``f_k`` is the frequency corresponding to that coefficient.
The operation is performed in-place.

To return an updated copy, use the out-of-place form
[`differentiate(::AbstractFourierTrace)`](@ref).
"""
function differentiate!(t::AbstractFourierTrace)
    # Follow Algorithm 1 of Steven Johnson at
    # https://math.mit.edu/~stevenj/fft-deriv.pdf
    # but for a one-sided FFT
    data = trace(t)
    data .*= (2*π*im).*frequencies(t)
    # Remove Nyquist component
    if iseven(nsamples(t))
        data[end] = zero(eltype(data))
    end
    t
end

"""
    differentiate(t::AbstractFourierTrace) -> t′

Differentiate the frequency-domain trace `t` in the Fourier domain by
multiplying the cofficients ``Y_k`` by ``2\\pi f_k i k``, where ``k <= N/2``
and ``N`` is the number of samples in the time-domain trace `ifft(t)`
and ``f_k`` is the frequency corresponding to that coefficient.
The operation is performed on a copy of `t`, which is returned.

See also: [`differentiate!(::AbstractFourierTrace)`](@ref).
"""
differentiate(t::AbstractFourierTrace) = differentiate!(deepcopy(t))

"""
    integrate!(t::AbstractFourierTrace) -> t

Integrate the frequency-domain trace `t` in the Fourier domain by
dividing the cofficients ``Y_k`` by ``2\\pi f_k i k``, where ``k <= N/2``
and ``N`` is the number of samples in the time-domain trace `ifft(t)`
and ``f_k`` is the frequency corresponding to that coefficient.
The operation is performed in-place.

Note that in this implementation, the zero-frequency (DC) term is set to
zero, such that differentiating the result of this function may produce
a different trace to that which was originally integrated.

To return an updated copy, use the out-of-place form
[`integrate(::AbstractFourierTrace)`](@ref).
"""
function integrate!(t::AbstractFourierTrace)
    data = trace(t)
    data ./= (2*π*im).*frequencies(t)
    data[begin] = zero(eltype(data))
    t
end

"""
    integrate(t::AbstractFourierTrace) -> t′

Integrate the frequency-domain trace `t` in the Fourier domain by
dividing the cofficients ``Y_k`` by ``2\\pi f_k i k``, where ``k <= N/2``
and ``N`` is the number of samples in the time-domain trace `ifft(t)`
and ``f_k`` is the frequency corresponding to that coefficient.
The operation is performed on a copy of `t`, which is returned.

Note that in this implementation, the zero-frequency (DC) term is set to
zero, such that differentiating the result of this function may produce
a different trace to that which was originally integrated.

See also: [`integrate!(::AbstractFourierTrace)`](@ref).
"""
integrate(t::AbstractFourierTrace) = integrate!(deepcopy(t))
