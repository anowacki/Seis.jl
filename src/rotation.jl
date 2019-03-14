# Trace rotations

"""
    rotate_through!(t1::Trace, t2::Trace, phi)

For two traces `t1` an `t2` which are horizontal and orthgonal, rotate
them *clockwise* by `phi`° about the vertical axis.

This is a reference frame transformation (passive rotation) and hence particle motion
will appear to rotate anti-clockwise.
"""
function rotate_through!(t1::AbstractTrace, t2::AbstractTrace, phi)
    T = promote_type(eltype(trace(t1)), eltype(trace(t2)))
    if !traces_are_orthogonal(t1, t2)
        throw(ArgumentError("traces must be orthogonal"))
    elseif t1.sta.inc != t2.sta.inc != 90
        throw(ArgumentError("traces must be horizontal"))
    elseif nsamples(t1) != nsamples(t2)
        throw(ArgumentError("traces must be same length"))
    elseif t1.delta != t2.delta
        throw(ArgumentError("traces must have same delta"))
    end
    sinp, cosp = sincos(deg2rad(T(phi)))
    @inbounds for i = 1:nsamples(t1)
        t1.t[i], t2.t[i] = cosp*t1.t[i] - sinp*t2.t[i], sinp*t1.t[i] + cosp*t2.t[i]
    end
    for t in (t1, t2)
        t.sta.azi = mod(t.sta.azi + phi, 360)
        t.sta.cha = string(round(t.sta.azi, digits=2, base=10))
    end
    t1, t2
end
rotate_through(t1::AbstractTrace, t2::AbstractTrace, phi) =
    rotate_through!(deepcopy(t1), deepcopy(t2), phi)
@doc (@doc rotate_through!) rotate_through

"""
    rotate_to_gcp!(t1, t2; reverse=false) -> t1, t2

Rotate the pair of traces `t1` and `t2` in place so that t1
points along the radial direction (the backazimuth plus 180°),
and `t2` is 90° clockwise from that.

If `reverse` is `true`, then `t2` is rotated to be 90° anticlockwise
from `t1`, so that the polarity is reversed.

The component names of the radial and transverse traces are updated
to be 'R', and either 'T' or '-T' respectively for normal and reverse
polarity.

Traces must be orthogonal and horizontal.
"""
function rotate_to_gcp!(t1, t2; reverse::Bool=false)
    all(!ismissing, (t1.sta.azi, t1.sta.inc, t2.sta.azi, t2.sta.inc)) ||
        throw(ArgumentError("Traces must have sta.inc and sta.azi defined"))
    t2_is_clockwise_of_t1 = angle_difference(t1.sta.azi, t2.sta.azi) > 0
    t2_is_clockwise_of_t1 || ((t1, t2) = (t2, t1))
    β = backazimuth(t1)
    β ≈ backazimuth(t2) ||
        throw(ArgumentError("Backazimuth is not the same for both traces"))
    ϕ = mod(β + 180 - t1.sta.azi, 360)
    rotate_through!(t2, t1, ϕ)
    reverse && flip_component!(t2)
    t1.sta.cha = "R"
    @assert t1.sta.azi ≈ mod(β + 180, 360)
    t2.sta.cha = reverse ? "-T" : "T"
    t1, t2
end

"""
    rotate_to_gcp!(t::AbstractArray{T); kwargs...)

Rotate pairs of traces in the array `t` so they are in the order
`[R1, T1, R2, T2, ...]`.
"""
function rotate_to_gcp!(t::AbstractArray{T}; kwargs...) where T<:AbstractTrace
    length(t) % 2 == 0 ||
        throw(ArgumentError("Array of traces must be a multiple of two"))
    for i in 1:length(t)÷2
        i1 = 2i - 1
        i2 = 2i
        rotate_to_gcp!(t[i1], t[i2]; kwargs...)
        if t[i1].sta.cha == "T"
            # FIXME: This shouldn't be necessary, but is probably due to swapping
            # of the components in rotate_to_gcp!(t1, t2) above.
            t[i1], t[i2] = t[i2], t[i1]
        end
    end
    t
end

"""
    rotate_to_gcp(t1, t2; reverse=false) -> R, T
    rotate_to_gcp(t::AbstractArray{<:AbstractTrace}; reverse=false) -> t′

Copying version of `rotate_to_gcp!` which returns the radial `R` and
transverse `T` traces in the first form, or pairs of radial and
transverse traces in the modified array `t′`
"""
rotate_to_gcp(args...; kwargs...) = rotate_to_gcp!(deepcopy.(args)...; kwargs...)
