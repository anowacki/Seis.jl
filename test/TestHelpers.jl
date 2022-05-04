"""
Module `TestHelpers` contains helper functions for the package's tests.

If you are running individual tests, include this file and do `using .TestHelpers`.
"""
module TestHelpers

using LinearAlgebra: ×, ⋅
using Seis
import Rotations
using StaticArrays: SVector, @SVector

export
    azimuth_inclination,
    compare_copy_modify_func,
    random_basis_traces,
    random_unit_vector,
    unit_vector,
    vector_and_deviated_vector

"Make sure that copying and in-place versions respectively do not and do modify the trace"
function compare_copy_modify_func(f, f!, args...; kwargs...)
    t = Trace(100rand(), rand(), rand(100))
    t′ = deepcopy(t)
    f(t, args...; kwargs...) == f!(t′, args...; kwargs...) && t != t′
end

"Randomly oriented unit vector"
random_unit_vector(T=Float64) = normalize(rand(SVector{3,T}) .- one(T)/2)

"Return the azimuth and inclination in the Seis frame from a unit vector"
azimuth_inclination(v) = atand(v[1], v[2]), acosd(v[3])

"Get a unit vector from an azimuth and inclination"
unit_vector(azi, inc) = @SVector[sind(azi)*sind(inc), cosd(azi)*sind(inc), cosd(inc)]

"""
Create a randomly-oriented set of traces making up a right-handed set,
with 
"""
function random_basis_traces(T=Float64, P=Seis.Geographic{T})
    x, y, z = [Trace{T, Vector{T}, P}(0, 1, 0) for _ in 1:3]
    x.sta.cha, y.sta.cha, z.sta.cha = "X", "Y", "Z"
    # Do vector creation and calculation in Float64
    # A random unit vector gives the first trace's direction
    ux = random_unit_vector(Float64)
    # Then find an orthogonal vector
    v = random_unit_vector(Float64)
    uy = normalize(ux × v)
    # Third direction gives right-handed set x, y, z
    uz = ux × uy
    # Only round from Float64 to lower precision at this point during
    # auto-conversion to T
    x.sta.azi, x.sta.inc = azimuth_inclination(ux)
    y.sta.azi, y.sta.inc = azimuth_inclination(uy)
    z.sta.azi, z.sta.inc = azimuth_inclination(uz)
    x, y, z
end
random_basis_traces(::Type{Trace{T}}) where T = random_basis_traces(T, Seis.Geographic{T})
random_basis_traces(::Type{CartTrace{T}}) where T = random_basis_traces(T, Seis.Cartesian{T})

"""
Fill the traces `t1`, `t2` and `t3` (an orthogonal set) with `data_l`, `data_q` and
`data_t` (which must all be the same length), and which are defined to be
respectively along the L, Q and T directions defined by `azimuth` and `inclination`.
"""
function project_onto_traces!(t1, t2, t3, azimuth, inclination, data_l, data_q, data_t)
    length(data_l) == length(data_q) == length(data_t) ||
        throw(ArgumentError("all data vectors must be the same length"))
    npts = length(data_l)

    # Get unit vector for each of L, Q and T
    l⃗ = unit_vector(azimuth, inclination)
    q⃗ = if inclination >= 90
        unit_vector(azimuth, inclination - 90)
    else
        unit_vector(azimuth - 180, 90 - inclination)
    end
    t⃗ = unit_vector(azimuth + 90, 90)

    for t in (t1, t2, t3)
        resize!(trace(t), npts)
        trace(t) .= 0
        # Get unit vector for this trace
        v = unit_vector(t.sta.azi, t.sta.inc)
        # Project new data onto existing traces
        for (data_vec, data) in zip((l⃗, q⃗, t⃗), (data_l, data_q, data_t))
            for i in 1:npts
                trace(t)[i] += data[i]*v⋅data_vec
            end
        end
    end
    t1, t2, t3
end

"Return a vector and one randomly deviated away from it by `θ`°."
function vector_and_deviated_vector(T, θ)
    v = random_unit_vector()
    w = random_unit_vector()
    # Vector normal to v (and w)
    u = normalize(v × w)
    # Rotate about this new vector u
    rot = Rotations.AngleAxis(deg2rad(θ), u...)
    v′ = rot*v
    v, v′
end
vector_and_deviated_vector(θ) = vector_and_deviated_vector(Float64, θ)

end # module
