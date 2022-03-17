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
polarity, unless the component code is a valid SEED identifier which
seems rotatable and matches for the traces; then the correct component name is used.
(E.g., `"BHE"` and `"BHN"` become `"BHR"` and `"BHT"`.)

Traces must be orthogonal and horizontal.
"""
function rotate_to_gcp!(t1, t2; reverse::Bool=false)
    all(!ismissing, (t1.sta.azi, t1.sta.inc, t2.sta.azi, t2.sta.inc)) ||
        throw(ArgumentError("Traces must have sta.inc and sta.azi defined"))
    t2_is_clockwise_of_t1 = angle_difference(t1.sta.azi, t2.sta.azi) > 0
    t2_is_clockwise_of_t1 || ((t1, t2) = (t2, t1))
    # Keep band and instrument codes only if the channel name is a valid SEED identifier and seems rotatable,
    # and we aren't reversing the transverse
    code = first(t1.sta.cha, 2)
    keep_code = _is_rotatable_seed_channel_name(t1.sta.cha) &&
        _is_rotatable_seed_channel_name(t2.sta.cha) && !reverse && code == first(t2.sta.cha, 2)
    β = backazimuth(t1)
    β ≈ backazimuth(t2) ||
        throw(ArgumentError("Backazimuth is not the same for both traces"))
    ϕ = mod(β + 180 - t1.sta.azi, 360)
    rotate_through!(t2, t1, ϕ)
    @assert t1.sta.azi ≈ mod(β + 180, 360)
    reverse && flip!(t2)
    # Channel code
    rcha = "R"
    tcha = reverse ? "-T" : "T"
    if keep_code
        t1.sta.cha = code * rcha
        t2.sta.cha = code * tcha
    else
        t1.sta.cha = rcha
        t2.sta.cha = tcha
    end
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
        if last(t[i1].sta.cha) == 'T'
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

"""
    rotate_to_azimuth_incidence!(t1, t2, t3, azimuth, incidence[; tol]) -> x, y, z

Rotate a mutually-orthogonal set of traces `t1`, `t2` and `t3` such that the trace
`x` points along the direction defined by `azimuth`, `y` points perpendicular to `x`
and has components only in the direction along `azimuth` and the vertical;
`z` is perpendicular to both and lies in the horizontal plane.
`x`, `y` and `z` form a right-handed set and are commonly known as L, Q and T
respectively.  (See [`rotate_to_lqt!`](@ref).)

Note that in this form, the underlying data and metadata in the input traces is
altered, and the returned traces point to the same data as the input traces, but
possibly in a different order.  Use [`rotate_to_azimuth_incidence`](@ref) to return
copies of the traces and leave the original traces unmodified.

See also: [`rotate_to_azimuth_incidence`](@ref), [`rotate_to_gcp!`](@ref),
[`rotate_to_lqt!`](@ref), [`rotate_through`](@ref)
"""
function rotate_to_azimuth_incidence!(t1::AbstractTrace, t2::AbstractTrace, t3::AbstractTrace,
        azimuth=backazimuth(t1)+180, incidence=incidence(t1);
        tol=_angle_tol(t1, t2, t3))
    # Check angles
    0 <= incidence <= 180 ||
        throw(ArgumentError("incidence must be in the range 0 to 180"))
    # Check traces
    if any(x -> x.sta.azi === missing || x.sta.inc === missing, (t1, t2, t3))
        throw(ArgumentError("traces must contain station orientation information"))
    elseif !are_orthogonal(t1, t2, t3)
        throw(ArgumentError("traces must be orthogonal"))
    elseif !(axes(trace(t1)) == axes(trace(t2)) == axes(trace(t3)))
        throw(ArgumentError("traces must all be the same length"))
    elseif !(eltype(t1) == eltype(t2) == eltype(t3))
        throw(ArgumentError("traces must all have the same element type"))
    elseif !(times(t1) == times(t2) == times(t3))
        throw(ArgumentError("traces must all have the same start time and sampling interval"))
    end
    # Sort traces so that we start with a right-handed set
    x, y, z, perm = sort_traces_right_handed(t1, t2, t3)

    # Component directions
    uˣ, uʸ, uᶻ = _direction_vector.((x, y, z))

    # New x direction is like L
    azˣ, incˣ = mod(azimuth, 360), incidence
    uˣ′ = _direction_vector(azˣ, incˣ)

    # New y is 'Q' component, accounting for flip of direction
    azʸ, incʸ = if incidence >= 90
        azimuth, incidence - 90
    else
        mod(azimuth + 180, 360), 90 - incidence
    end
    uʸ′ = _direction_vector(azʸ, incʸ) # Like Q

    # New z is 'T' component
    azᶻ = mod(azimuth + 90, 360)
    incᶻ = 90
    uᶻ′ = uˣ′ × uʸ′

    rotation_matrix = _construct_rotation_matrix(uˣ, uʸ, uᶻ, uˣ′, uʸ′, uᶻ′)

    # Do the rotation
    tx = trace(x)
    ty = trace(y)
    tz = trace(z)
    @inbounds for i in eachindex(tx)
        tx[i], ty[i], tz[i] = rotation_matrix*StaticArrays.@SVector[tx[i], ty[i], tz[i]]
    end

    # Set headers
    x.sta.azi, x.sta.inc = azˣ, incˣ
    y.sta.azi, y.sta.inc = azʸ, incʸ
    z.sta.azi, z.sta.inc = azᶻ, incᶻ

    for (i, t) in enumerate((x, y, z))
        code = if isapprox(t.sta.inc, 0, atol=tol)
            "Z"
        elseif isapprox(t.sta.inc, 90, atol=tol) && isapprox(t.sta.azi, 0, atol=tol)
            "N"
        elseif isapprox(t.sta.inc, 90, atol=tol) && isapprox(t.sta.azi, 90, atol=tol)
            "E"
        else
            string(i)
        end

        if _is_rotatable_seed_channel_name(t.sta.cha)
            cha = t.sta.cha
            t.sta.cha = cha[1:2] * code
        else
            t.sta.cha = code
        end
    end

    x, y, z
end

"""
    rotate_to_lqt!(t1, t2, t3[, azimuth, [incidence]]; [tol]) -> L, Q, T

Rotate the three traces `t1`, `t2` and `t3` in place, and return them in the
order `L`, `Q` and `T`, forming a right-handed set.

The `L` trace records positive motion along the event-receiver direction,
defined by the local `azimuth` at the receiver (i.e., backazimuth + 180°)
and `incidence` angle (measured downwards away from the positive upwards direction).

The `T` direction is transverse to `L`, such that when looking along the
direction towards the station, `T` is on the right and lies in the horizontal plane.

The `Q` direction is perpendicular to both, lying in the saggital plane, with
some component in the vertical direction.

If `azimuth` and `inclination` are not passed in explicitly, then they are
determined using the event-station geometry.  This is only possible for
`Trace`s in Cartesian geometry, since `inclination` is not well-determined in
a geographic system.  For the Cartesian case, we assume straight-line paths
between event and receiver.

`tol` specifies the angle in ° by which the traces must be orthogonal; see
[`traces_are_orthogonal`](@ref)

!!! note
    If `incidence` is 0° or 90° (the direction is vertical), then the choice
    of `Q` and `T` is arbitrary and simply chosen such that `Q` points along
    the local `azimuth`.

# Diagram
## Plan view
```
 ⋆ (event)       North
  `                ↑
    `              |
      `            |
        `
          ∇ (station)
         ⊙ `
      -  Q   `
    ↙          ↘
   T            L
```

## Side view
```
                 Q
                 ↑
 Up               |
 ↑                 |
 |                  |      __ → L
 |      (station)  ∇ ⊙__---
             .__--   T
        .__--
    __--
  ⋆ (event)
```

# Examples
```
julia> e, n, z = sample_data(:regional)[1:3];

julia> azi = azimuth(e) + 180
216.85858118935954

julia> inc = 30; # Determined from other calculation

julia> l, q, t = rotate_to_lqt!(e, n, z, azi, inc)
(Seis.Trace(.ELK..1: delta=0.025, b=-5.000015, nsamples=12000), Seis.Trace(.ELK..2: delta=0.025, b=-5.000015, nsamples=12000), Seis.Trace(.ELK..T: delta=0.025, b=-5.000015, nsamples=12000))

julia> l.sta.azi, l.sta.inc, l.sta.cha
(216.85858f0, 30.0f0, "1")

julia> t.sta.cha
"T"

julia> l == e, q == n, t == z # Original traces are modified and returned in a possibly different order
(true, true, true)
```
"""
function rotate_to_lqt!(t1, t2, t3, azimuth=backazimuth(t1)+180, incidence=incidence(t1);
        tol=_angle_tol(t1, t2, t3))
    l, q, t = rotate_to_azimuth_incidence!(t1, t2, t3, azimuth, incidence; tol=tol)
    if _is_rotatable_seed_channel_name(l.sta.cha)
        cha = l.sta.cha
        t.sta.cha = cha[1:2] * "T"
    else
        t.sta.cha = "T"
    end
    l, q, t
end

"""
    rotate_to_lqt(t1, t2, t3[, azimuth[, incidence]; tol) -> l, q, t

Copying version of [`rotate_to_lqt!`](@ref).
"""
rotate_to_lqt(t1, t2, t3, args...; kwargs...) =
    rotate_to_lqt!(deepcopy.((t1, t2, t3))..., args...; kwargs...)

"""
    rotate_to_enz!(t1, t2, t3[; tol]) -> e, n, z

Rotate three orthogonal traces `t1`, `t2` and `t3` in place so that they point
east (`e`), north (`n`) and vertically (`z`).

`e`, `n` and `z` are bindings to the same data as `t1`, `t2` and `t3`,
but they may not be returned in the same order as passed in if the original
set of traces are not given as a right-handed set.
"""
function rotate_to_enz!(t1, t2, t3; tol=_angle_tol(t1, t2, t3))
    n, z, e = rotate_to_azimuth_incidence!(t1, t2, t3, 0, 90; tol=tol)
    for (t, cha) in zip((e, n, z), ("E", "N", "Z"))
        if _is_rotatable_seed_channel_name(t.sta.cha)
            t.sta.cha = t.sta.cha[1:2] * cha
        else
            t.sta.cha = cha
        end
    end
    e, n, z
end

"""
    rotate_to_enz(t1, t2, t3[; tol]) -> e, n, z

Copying version of [`rotate_to_enz!`](@ref).
"""
rotate_to_enz(t1, t2, t3; kwargs...) = rotate_to_enz!(deepcopy.((t1, t2, t3))...; kwargs...)

"""
    sort_traces_right_handed(t1, t2, t3) -> x, y, z, perm

Sort the traces `t1`, `t2` and `t3` such that they form a right-handed set;
requiring the traces to be mutually orthogonal.  The direction of `x`, `y` and `z`
are arbitrary.  `perm` is a length-three tuple containing the indices (from 1 to 3)
of the new order of traces, such that `x` is `(t1, t2, t3)[perm[1]]` and so on.

# Example
```
julia> e, n, z = sample_data(:regional)[1:3];

julia> e.sta.cha, n.sta.cha, z.sta.cha # Confirm the order
("e", "n", "z")

julia> u, v, w, perm = Seis.sort_traces_right_handed(n, e, z); # A left-handed arrangement

julia> [u, v, w].sta.cha
3-element Vector{String}:
 "n"
 "z"
 "e"

julia> [n, e, z][[perm...]] # Indexing by perm gives us the new order
3-element Vector{Trace{Float32, Vector{Float32}, Seis.Geographic{Float32}}}:
 Seis.Trace(.ELK..n: delta=0.025, b=-5.000015, nsamples=12000)
 Seis.Trace(.ELK..z: delta=0.025, b=-5.000015, nsamples=12000)
 Seis.Trace(.ELK..e: delta=0.025, b=-5.000015, nsamples=12000)
```
"""
function sort_traces_right_handed(t1::AbstractTrace, t2::AbstractTrace, t3::AbstractTrace;
        tol=_angle_tol(t1, t2, t3))
    are_orthogonal(t1, t2, t3; tol=tol) ||
        throw(ArgumentError("traces must be orthogonal"))
    # If t1, t2, t3 form a right-handed set and 𝐮₁, 𝐮₂ and 𝐮₃ are the unit vectors
    # pointing along the channels, then 𝐮₁×𝐮₂ = 𝐮₃.  Otherwise, 𝐮₁×𝐮₂ = –𝐮₃.
    u1, u2, u3 = _direction_vector.((t1, t2, t3))
    u1_cross_u2 = u1 × u2
    # Only compare to the nearest degree
    if _directions_are_parallel(u1_cross_u2, u3, 1)
        return t1, t2, t3, (1, 2, 3)
    else
        return t1, t3, t2, (1, 3, 2)
    end
end

"""
    _is_rotatable_seed_channel_name(s) -> ::Bool

Return `true` if `s` is a string conforming to SEED channel naming conventions which
describes a channel for which it makes sense to rotate a matching pair of these channels.

Returns `false` if the above does not hold or `s` is `missing`.

#### References
- https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
"""
_is_rotatable_seed_channel_name(s) = occursin(r"^[FGDCESHBMLVURPTQAO][HLGMNAFJPSXYZ][ZNEABCTR123UVW]$", s)
_is_rotatable_seed_channel_name(::Missing) = false

"""
    _construct_rotation_matrix(uˣ, uʸ, uᶻ, uˣ′, uʸ′, uᶻ′) -> R

Construct the rotation matrix `R` which maps the basis formed by `(uˣ, uʸ, uᶻ)`
to the new basis formed by `(uˣ′, uʸ′, uᶻ′)`, where each `uⁱ` is a single unit
vector.
"""
function _construct_rotation_matrix(uˣ, uʸ, uᶻ, uˣ′, uʸ′, uᶻ′)
    # If `R¹` is formed by placing `uˣ`, `uʸ` and `uᶻ` in the columns of a matrix, and
    # likewise `R²` for the second basis set, and `R` is the rotation matrix,
    # then `R * R¹ = R²` and `R = inv(R¹) * R²`.  Thankfully `inv(R¹)` is just `(R¹)ᵀ`
    # because the vectors are orthogonal, thus we avoid doing `R¹\R²`.
    R¹ = StaticArrays.@SMatrix[uˣ[1] uʸ[1] uᶻ[1]
                               uˣ[2] uʸ[2] uᶻ[2]
                               uˣ[3] uʸ[3] uᶻ[3]]
    R² = StaticArrays.@SMatrix[uˣ′[1] uʸ′[1] uᶻ′[1]
                               uˣ′[2] uʸ′[2] uᶻ′[2]
                               uˣ′[3] uʸ′[3] uᶻ′[3]]
    transpose(R²)*R¹
end

"""
    _rotate_by_vector(v, k, θ) -> v′

Rotate the vector `v` about the vector `k` by `θ` radians,
following the right-hand rule.  `k` need not be a unit vector.

This implements Rodrigues's formula.

# References
- https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
"""
function _rotate_by_vector(v, k, θ)
    k̂ = normalize(k)
    cosθ = cos(θ)
    u′ = v*cosθ + (k̂×v)*sin(θ) + k̂*(k̂⋅v)*(1 - cosθ)
end
