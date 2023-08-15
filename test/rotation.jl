# Test rotation to great circle path
using Test
using Seis
using LinearAlgebra: ×
using StaticArrays: SVector, @SVector
import Geodesics

import .TestHelpers

trace_permutations(x, y, z) = (x, y, z), (x, z, y), (y, x, z), (y, z, x), (z, x, y), (z, y, x)

angles_are_same(a, b, tol=Seis._angle_tol(typeof(float(a)), typeof(float(b)))) =
    abs(Seis.angle_difference(a, b)) < tol

@testset "Trace rotation" begin
    @testset "rotate_through" begin
        @testset "Argument checking" begin
            @testset "Not orthogonal" begin
                @testset "Default tolerance" begin
                    @test_throws ArgumentError rotate_through(sample_data(), sample_data(), 10)
                end
                @testset "Custom tolerance" begin
                    e, n, z = sample_data(:local)[1:3]
                    e.sta.azi += 1e-5
                    @test_throws ArgumentError rotate_through(e, n, 10, tol=1e-6)
                end
            end
            @testset "Different samples" begin
                e, n = sample_data(:local)[1:2]
                pop!(trace(e))
                @test_throws ArgumentError rotate_through(e, n, 1)
            end
            @testset "Different delta" begin
                e, n = sample_data(:local)[1:2]
                e.delta *= 2
                @test_throws ArgumentError rotate_through(e, n, 1)
            end
        end

        atol = 1e-6
        @testset "Polarisation $pol" for pol in rand(0:359, 3)
            @testset "Rotation $rot" for rot in rand(-180:180, 2)
                @testset "Horizontals" begin
                    n = Trace(0, 1, 1)
                    e = deepcopy(n)
                    n.sta.azi, e.sta.azi = 0, 90
                    n.sta.inc, e.sta.inc = 90, 90
                    # Initial polarisation
                    @testset "N to E" begin
                        # Try rotating from n -> e
                        trace(n)[1], trace(e)[1] = cosd(pol), sind(pol)
                        rotate_through!(deepcopy.((n, e))..., rot) == rotate_through(n, e, rot)
                        rotate_through!(deepcopy.((n, e))..., rot) != (n, e)
                        n′, e′ = rotate_through(n, e, rot)
                        # Polarisation should now be as if the polarisation was
                        # pol - rot
                        @test all(isapprox.(
                            (trace(n′)[1], trace(e′)[1]), (cosd(pol - rot), sind(pol - rot)), atol=atol))
                        @test angles_are_same(n′.sta.azi, mod(rot, 360))
                        @test angles_are_same(e′.sta.azi, mod(rot + 90, 360))
                    end

                    @testset "E to N" begin
                        e′, n′ = rotate_through(e, n, rot)
                        # Polarisation now as if pol + rot
                        @test all(isapprox.(
                            (trace(n′)[1], trace(e′)[1]), (cosd(pol + rot), sind(pol + rot)), atol=atol))
                        @test angles_are_same(n′.sta.azi, mod(-rot, 360), atol)
                        @test angles_are_same(e′.sta.azi, mod(90 - rot, 360), atol)
                    end

                    @testset "Reversible" begin
                        @testset "Negative rotation" begin
                            n″, e″ = rotate_through(
                                rotate_through(n, e, rot)..., -rot)
                            @test trace(n) ≈ trace(n″)
                            @test trace(e) ≈ trace(e″)
                            @test angles_are_same(n.sta.azi, n″.sta.azi, atol)
                            @test n.sta.inc ≈ n″.sta.inc atol=atol
                            @test angles_are_same(e.sta.azi, e″.sta.azi, atol)
                            @test e.sta.inc ≈ e″.sta.inc atol=atol
                        end

                        @testset "Flip order" begin
                            n″, e″ = rotate_through(
                                rotate_through(n, e, rot)[[2,1]]..., rot)[[2,1]]
                            @test trace(n) ≈ trace(n″)
                            @test trace(e) ≈ trace(e″)
                            @test angles_are_same(n.sta.azi, n″.sta.azi, atol)
                            @test n.sta.inc ≈ n″.sta.inc atol=atol
                            @test angles_are_same(e.sta.azi, e″.sta.azi, atol)
                            @test e.sta.inc ≈ e″.sta.inc atol=atol
                        end
                    end

                    @testset "Arbitrary" begin
                        t1 = Trace(0, 1, 1)
                        t2 = deepcopy(t1)
                        t1.sta.azi, t2.sta.azi = (45, 45)
                        t1.sta.inc, t2.sta.inc = (45, 135)
                    end
                end
            end
        end

        @testset "In-place v copying" begin
            e, n = sample_data(:local)[1:2]
            e′, n′ = deepcopy.((e, n))
            @test rotate_through(e, n, 10) == rotate_through!(e′, n′, 10)
            @test e′ != e
        end
    end

    @testset "rotate_to_gcp" begin
        # The azimuth and distance of the event
        az = rand(0:359) # degrees, approximate for now
        d = 0.1 # degrees

        # Pair of traces representing a spike arriving on the radial component
        b = 0
        npts = 10
        delta = 1
        s = Trace.([b,b], [delta,delta], [npts,npts])
        s.sta.lon = 0.0
        s.sta.lat = 0.0
        s.evt.lat = -d*cosd(az)
        s.evt.lon = -d*sind(az)

        # Make azimuth more accurate as using great circle calculation
        az = mod(backazimuth(s[1]) + 180, 360)
        for ss in s
            ss.t .= zeros(npts)
        end
        # Randomly swap the order of the components, which are in random orientations
        ne = rand(Bool)
        az1 = rand(0:359)
        s.sta.azi = ne ? [az1, az1+90] : [az1+90, az1]
        s.sta.cha = ne ? ["1", "2"] : ["2", "1"]
        s.sta.inc = 90.0
        imax = npts÷2
        s[1].t[imax], s[2].t[imax] = ne ? (cosd(az-az1), sind(az-az1)) :
                                          (sind(az-az1), cosd(az-az1))

        # Create all possible combinations of rotations
        (r, t) = rotate_to_gcp(s[1], s[2])
        (r′, t′) = rotate_to_gcp(s[2], s[1])
        rt = rotate_to_gcp(s)
        rt′ = rotate_to_gcp(s[[2,1]])
        (R, T) = rotate_to_gcp!(deepcopy(s[1]), deepcopy(s[2]))
        (R′, T′) = rotate_to_gcp!(deepcopy(s[2]), deepcopy(s[1]))
        RT = rotate_to_gcp!(deepcopy(s))
        RT′ = rotate_to_gcp!(deepcopy(s[[2,1]]))

        list_of_radials = [r, r′, rt[1], rt′[1], R, R′, RT[1], RT′[1]]
        list_of_transverses = [t, t′, rt[2], rt′[2], T, T′, RT[2], RT′[2]]

        # Test rotation has been done correctly
        atol = sqrt(eps(Float64))
        for (r1,t1) in zip(list_of_radials, list_of_transverses)
            @test r1.t[imax] ≈ 1
            @test all(isapprox.([r1.t[1:imax-1]; r1.t[imax+1:end]], 0, atol=atol))
            @test all(isapprox.(t1.t, 0, atol=atol))
        end

        # Test order of arguments doesn't stop returns being in order r, t
        for r1 in list_of_radials, r2 in list_of_radials
            @test r1 == r2
        end
        for t1 in list_of_transverses, t2 in list_of_transverses
            @test t1 == t2
        end

        # Renaming of channels
        s[1].sta.cha = "LXN"
        s[2].sta.cha = "LXE"
        r, t = rotate_to_gcp(s[1], s[2])
        @test r.sta.cha == "LXR"
        @test t.sta.cha == "LXT"
        s[1].sta.cha, s[2].sta.cha = "WeirdA", "WeirdB"
        r, t = rotate_to_gcp(s[1], s[2])
        @test r.sta.cha == "R"
        @test t.sta.cha == "T"
        r, t = rotate_to_gcp(s[1], s[2], reverse=true)
        @test r.sta.cha == "R"
        @test t.sta.cha == "-T"
    end

    @testset "sort_traces_right_handed" begin
        x, y, z = random_basis_traces()

        @testset "$([traces...].sta.cha)" for traces in trace_permutations(x, y, z)
            @test are_orthogonal(traces..., tol=1)
            x′, y′, z′ = sort_traces_right_handed(traces...)
            ux′ = TestHelpers.unit_vector(x′.sta.azi, x′.sta.inc)
            uy′ = TestHelpers.unit_vector(y′.sta.azi, y′.sta.inc)
            uz′ = TestHelpers.unit_vector(z′.sta.azi, z′.sta.inc)
            @test ux′ × uy′ ≈ uz′
        end
    end

    @testset "rotate_to_azimuth_incidence" begin
        @testset "Arguments" begin
            e, n, z = sample_data(:regional)[1:3]

            @testset "Incidence" begin
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e, n, z, 0, -1)
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e, n, z, 0, 181)
            end

            @testset "Orthogonal" begin
                e′ = deepcopy(e)
                e′.sta.azi += 5
                @testset "Default tolerance" begin
                    @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
                end
                @testset "Custom tolerance" begin
                    @test length(rotate_to_azimuth_incidence(e′, n, z, 0, 0, tol=90)) == 3
                end
            end

            @testset "Trace data" begin
                e′ = deepcopy(e)
                pop!(trace(e′))
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
            end

            @testset "Header info" begin
                @testset "Missing $field" for field in (:azi, :inc)
                    e′ = deepcopy(e)
                    setproperty!(e′.sta, field, missing)
                    @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
                end
            end

            @testset "Eltype" begin
                e′ = convert(Trace{Float32, Vector{Float64}, Seis.Geographic{Float32}}, e)
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
            end

            @testset "Sampling interval" begin
                e′ = deepcopy(e)
                e′.delta *= 2
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
                e′.delta = e.delta
                e′.b += 1
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
            end
        end

        @testset "Simple rotations" begin
            @testset "$T" for T in (Float32, Float64)
                e = Trace{T, Vector{Float64}}(0, 1, [0])
                e.sta.azi, e.sta.inc = 90, 90
                e.sta.cha = "E"
                n = Trace{T, Vector{Float64}}(0, 1, [0])
                n.sta.azi, n.sta.inc = 0, 90
                n.sta.cha = "N"
                z = Trace{T, Vector{Float64}}(0, 1, [1])
                z.sta.azi, z.sta.inc = 0, 0
                z.sta.cha = "Z"

                atol = √eps(T)

                @testset "Vertical" begin
                    azimuth = 0
                    incidence = 0
                    @testset "$([traces...].sta.cha)" for traces in trace_permutations(e, n, z)
                        l, q, t = rotate_to_azimuth_incidence!(deepcopy.(traces)...,
                            azimuth, incidence)
                        @test trace(l)[1] ≈ 1 atol=atol
                        @test trace(q)[1] ≈ 0 atol=atol
                        @test trace(t)[1] ≈ 0 atol=atol
                        @test l.sta.azi == azimuth
                        @test l.sta.inc == incidence
                        @test q.sta.azi ≈ 180
                        @test q.sta.inc ≈ 90
                        @test t.sta.azi ≈ 90
                        @test t.sta.azi == 90
                    end
                end

                @testset "East" begin
                    azimuth = 90
                    incidence = 90
                    @testset "$([traces...].sta.cha)" for traces in trace_permutations(e, n, z)
                        l, q, t = rotate_to_azimuth_incidence!(deepcopy.(traces)...,
                            azimuth, incidence)
                        @test trace(l)[1] ≈ 0
                        @test trace(q)[1] ≈ 1
                        @test trace(t)[1] ≈ 0
                    end
                end

                @testset "North" begin
                    azimuth = 0
                    incidence = 90
                    @testset "$([traces...].sta.cha)" for traces in trace_permutations(e, n, z)
                        l, q, t = rotate_to_azimuth_incidence!(deepcopy.((e, n, z))...,
                            azimuth, incidence)
                        @test trace(l)[1] ≈ 0
                        @test trace(q)[1] ≈ 1
                        @test trace(t)[1] ≈ 0
                    end
                end

                @testset "Complex polarisation" begin
                    azimuth = 360rand(T)
                    incidence = 180rand(T)
                    trace(e)[1] = sind(azimuth)*sind(incidence)
                    trace(n)[1] = cosd(azimuth)*sind(incidence)
                    trace(z)[1] = cosd(incidence)

                    @testset "In line" begin
                        @testset "$([traces...].sta.cha)" for traces in trace_permutations(e, n, z)
                            l, q, t = rotate_to_azimuth_incidence!(deepcopy.(traces)...,
                                azimuth, incidence)
                            @test trace(l)[1] ≈ 1 atol=atol
                            @test trace(q)[1] ≈ 0 atol=atol
                            @test trace(t)[1] ≈ 0 atol=atol
                            @test l.sta.azi == azimuth
                            @test l.sta.inc == incidence
                            @test q.sta.azi ≈ (incidence < 90 ? mod(azimuth + 180, 360) : azimuth)
                            @test q.sta.inc ≈ (incidence >= 90 ? incidence - 90 : 90 - incidence)
                            @test t.sta.azi ≈ mod(azimuth + 90, 360)
                            @test t.sta.inc == 90
                        end
                    end

                    @testset "Perpendicular on Q" begin
                        trace(e)[1] = sind(azimuth)*sind(incidence - 90)
                        trace(n)[1] = cosd(azimuth)*sind(incidence - 90)
                        trace(z)[1] = cosd(incidence - 90)

                        @testset "$([traces...].sta.cha)" for traces in trace_permutations(e, n, z)
                            l, q, t = rotate_to_azimuth_incidence!(deepcopy.(traces)...,
                                azimuth, incidence)
                            @test trace(l)[1] ≈ 0 atol=atol
                            @test trace(q)[1] ≈ 1 atol=atol
                            @test trace(t)[1] ≈ 0 atol=atol
                        end
                    end

                    @testset "Perpendicular on -T" begin
                        trace(e)[1] = sind(azimuth - 90)*sind(90)
                        trace(n)[1] = cosd(azimuth - 90)*sind(90)
                        trace(z)[1] = cosd(90)

                        @testset "$([traces...].sta.cha)" for traces in trace_permutations(e, n, z)
                            l, q, t = rotate_to_azimuth_incidence!(deepcopy.(traces)...,
                                azimuth, incidence)
                            @test trace(l)[1] ≈ 0 atol=atol
                            @test trace(q)[1] ≈ 0 atol=atol
                            @test trace(t)[1] ≈ -1 atol=atol
                        end
                    end
                end
            end

            @testset "In-place v copying" begin
                e, n, z = sample_data(:local)[1:3]
                e′, n′, z′ = deepcopy.((e, n, z))
                az, inc = 45, 10
                @test rotate_to_azimuth_incidence(e, n, z, az, inc) ==
                    rotate_to_azimuth_incidence!(e′, n′, z′, az, inc)
                @test e′ != e
                @test n′ != n
                @test z′ != z
            end
        end
    end

    @testset "rotate_to_lqt" begin
        @testset "Random orientations" begin
            @testset "$P" for P in (Seis.Geographic, Seis.Cartesian)
                @testset "$T" for T in (Float32, Float64)
                    x, y, z = TestHelpers.random_basis_traces(T, P{T})
                    azi, inc = 45, 90
                    TestHelpers.project_onto_traces!(x, y, z, azi, inc, [1], [2], [3])
                    l, q, t = rotate_to_lqt(x, y, z, azi, inc)
                    @test trace(l)[1] ≈ 1
                    @test trace(q)[1] ≈ 2
                    @test trace(t)[1] ≈ 3
                    @test angles_are_same(l.sta.azi, azi)
                    @test angles_are_same(l.sta.inc, 90)
                    # Q is vertical so azimuth not important, but it will be azi
                    @test angles_are_same(q.sta.azi, azi)
                    @test angles_are_same(q.sta.inc, 0)
                    @test angles_are_same(t.sta.azi, azi + 90)
                    @test angles_are_same(t.sta.inc, 90)
                end
            end
        end

        @testset "Using trace headers" begin
            @testset "$P" for P in (Seis.Geographic, Seis.Cartesian)
                @testset "$T" for T in (Float32, Float64)
                    x, y, z = TestHelpers.random_basis_traces(T, P{T})
                    azi = 360rand(T)
                    inc = 180rand(T)
                    TestHelpers.project_onto_traces!(x, y, z, azi, inc, [1], [2], [3])
                    # Helpful array of traces
                    t = [x, y, z]
                    # Define headers appropriately for geometry
                    if P <: Seis.Geographic
                        t.sta.lon = t.sta.lat = t.sta.elev = 0
                        # Approximate location of event
                        Δ = 0.01
                        lon, lat = Geodesics.angular_step(x.sta.lon, x.sta.lat, azi + 180, Δ)
                        t.evt.lon = lon
                        t.evt.lat = lat
                        t.evt.dep = 20
                        # Refined azimuth at station
                        azi = backazimuth(x) + 180
                        l, q, t = rotate_to_lqt(x, y, z, inc)
                    else
                        t.sta.x = t.sta.y = t.sta.z = 0
                        # Compute position which gives correct azimuth and incidence
                        t.evt.x = 0 - sind(azi)*sind(inc)
                        t.evt.y = 0 - cosd(azi)*sind(inc)
                        t.evt.z = 0 - cosd(inc)
                        l, q, t = rotate_to_lqt(x, y, z)
                    end
                    # Reference traces, possibly using updated azimuth above
                    lqt_test = rotate_to_lqt(x, y, z, azi, inc)

                    @test angles_are_same(l.sta.azi, lqt_test[1].sta.azi)
                    @test angles_are_same(l.sta.inc, lqt_test[1].sta.inc)
                    @test angles_are_same(q.sta.azi, lqt_test[2].sta.azi)
                    @test angles_are_same(q.sta.inc, lqt_test[2].sta.inc)
                    @test angles_are_same(t.sta.azi, lqt_test[3].sta.azi)
                    @test angles_are_same(t.sta.inc, lqt_test[3].sta.inc)
                    @test trace(l)[1] ≈ 1
                    @test trace(q)[1] ≈ 2
                    @test trace(t)[1] ≈ 3
                end
            end
        end

        @testset "In-place v copying" begin
            @testset "$T" for T in (Float32, Float64)
                x, y, z = TestHelpers.random_basis_traces(T)
                for t in (x, y, z)
                    push!(trace(t), rand())
                end
                x′, y′, z′ = deepcopy.((x, y, z))
                azi = 360rand()
                inc = 180rand()
                @test rotate_to_lqt!(x′, y′, z′, azi, inc) ==
                    rotate_to_lqt(x, y, z, azi, inc)
                @test x != x′
                @test y != y′
                @test z != z′
            end
        end
    end

    @testset "rotate_to_enz" begin
        @testset "ENZ" begin
            @testset "$T" for T in (Trace, CartTrace)
                @testset "$(T{F})" for F in (Float16, Float32, Float64)
                    e, n, z = t = [T{F}(0, 1, 0) for _ in 1:3]
                    t.sta.azi = 90, 0, 0
                    t.sta.inc = 90, 90, 0
                    t.sta.cha = "E", "N", "Z"
                    push!.(trace.(t), 1:3)
                    @test rotate_to_enz(e, n, z) == (e, n, z)
                    @test rotate_to_enz(n, e, z) == (e, n, z)
                    @test rotate_to_enz(z, n, e) == (e, n, z)
                    @test rotate_to_enz(z, e, n) == (e, n, z)
                    e′, n′, z′ = t′ = deepcopy(t)
                    @test rotate_to_enz(e′, n′, z′) == (e, n, z)
                    @test rotate_to_enz(n′, e′, z′) == (e, n, z)
                    @test rotate_to_enz(z′, n′, e′) == (e, n, z)
                    @test rotate_to_enz(z′, e′, n′) == (e, n, z)
                end
            end
        end

        @testset "Random" begin
            @testset "$TR" for TR in (Trace, CartTrace)
                @testset "$T" for T in (Float32, Float64)
                    t1, t2, t3 = TestHelpers.random_basis_traces(TR{T})
                    t1.sta.cha = "HH1"
                    t2.sta.cha = "HH2"
                    t3.sta.cha = "HH3"
                    azi = 360rand()
                    inc = 90rand()
                    # With these angles, 1->Z, 2->-N, 3->E and rotate_to_enz
                    # returns in the order E, N, Z, so we assign numbers which
                    # will end up as E=1, N=2, Z=3.
                    TestHelpers.project_onto_traces!(t1, t2, t3, 0, 0, [3], [-2], [1])
                    e, n, z = rotate_to_enz(t1, t2, t3)

                    @testset "Station angles" begin
                        @test angles_are_same(e.sta.azi, 90)
                        @test angles_are_same(e.sta.inc, 90)
                        @test angles_are_same(n.sta.azi, 0)
                        @test angles_are_same(n.sta.inc, 90)
                        @test angles_are_same(z.sta.azi, 0)
                        @test angles_are_same(z.sta.inc, 0)
                    end

                    @testset "Data" begin
                        @test trace(e)[1] ≈ 1
                        @test trace(n)[1] ≈ 2
                        @test trace(z)[1] ≈ 3
                    end

                    @testset "Channel names" begin
                        @test e.sta.cha == "HHE"
                        @test n.sta.cha == "HHN"
                        @test z.sta.cha == "HHZ"                        
                    end

                    @testset "In-place versus copying" begin
                        t1′, t2′, t3′ = deepcopy.((t1, t2, t3))
                        e′, n′, z′ = rotate_to_enz!(t1′, t2′, t3′)
                        @test (e, n, z) == (e′, n′, z′)
                        @test t1 != t1′
                        @test t2 != t2′
                        @test t3 != t3′
                    end
                end
            end
        end

    end

    @testset "_is_rotatable_seed_channel_name" begin
        @test Seis._is_rotatable_seed_channel_name("BH1") == true
        @test Seis._is_rotatable_seed_channel_name("LXU") == true
        @test Seis._is_rotatable_seed_channel_name("BHL") == false
        @test Seis._is_rotatable_seed_channel_name(missing) == false
    end

    @testset "_rotate_by_vector" begin
        @test Seis._rotate_by_vector([0, 0, 1], [1, 0, 0], π/2) ≈ [0, -1, 0]
        @test Seis._rotate_by_vector([0, 0, 1], [0, 1, 0], π/2) ≈ [1, 0, 0]
        # No need to normalise the rotation vector
        @test Seis._rotate_by_vector([0, 0, 1], [0, -10, 0], π/2) ≈ [-1, 0, 0]
        v′ = Seis._rotate_by_vector(@SVector[0, 0, 1], @SVector[0, -10, 0], π/2)
        @test v′ isa SVector
        @test v′ ≈ @SVector[-1, 0, 0]
    end
end
