"""
Module `TestHelpers` contains helper functions for the package's tests.

If you are running individual tests, include this file and do `using .TestHelpers`.
"""
module TestHelpers

using LinearAlgebra: ×
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

"Create a randomly-oriented set of traces making up a right-handed set"
function random_basis_traces(T=Float64, P=Seis.Geographic{T})
    x, y, z = [Trace{T, Vector{T}, P}(0, 1, 0) for _ in 1:3]
    x.sta.cha, y.sta.cha, z.sta.cha = "X", "Y", "Z"
    # Do vector creation and calculation in Float64
    # A random unit vector give the first trace's direction
    ux = random_unit_vector(Float64)
    # Then find an orthogonal vector
    v = random_unit_vector(Float64)
    uy = normalize(ux × v)
    # Third direction gives right-handed set x, y, z
    uz = ux × uy
    # Only round to lower precision at this point
    x.sta.azi, x.sta.inc = azimuth_inclination(ux)
    y.sta.azi, y.sta.inc = azimuth_inclination(uy)
    z.sta.azi, z.sta.inc = azimuth_inclination(uz)
    x, y, z
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
