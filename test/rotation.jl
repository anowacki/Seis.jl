# Test rotation to great circle path
using Test
using Seis
using LinearAlgebra: ×
using StaticArrays: SVector, @SVector

import .TestHelpers

trace_permutations(x, y, z) = (x, y, z), (x, z, y), (y, x, z), (y, z, x), (z, x, y), (z, y, x)

@testset "Trace rotation" begin
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
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
            end

            @testset "Trace data" begin
                e′ = deepcopy(e)
                pop!(trace(e′))
                @test_throws ArgumentError rotate_to_azimuth_incidence!(e′, n, z, 0, 0)
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
        end

        @testset "Random trace orientation" begin
            T = Float64
            x, y, z = TestHelpers.random_basis_traces(T)
            # Add a single data point which we will modify
            push!(trace(x), 0)
            push!(trace(y), 0)
            push!(trace(z), 0)

            @testset "On L" begin
                v = TestHelpers.random_unit_vector(T)
                azimuth, incidence = TestHelpers.azimuth_inclination(v)
                @info "TODO: Add tests for trace rotation"
            end


        end
    end

    @testset "rotate_to_lqt" begin
        @testset "$T" for T in (Float32, Float64)
            x, y, z = TestHelpers.random_basis_traces(T)
            @info "TODO: Add tests for `rotate_to_lqt` and implement for `CartTrace`"
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
                    x, y, z = TestHelpers.random_basis_traces(TR{T})
                    @info "TODO: Add tests for `rotate_to_enz`"
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
