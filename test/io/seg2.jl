using Dates: DateTime
using Test
using Seis

# Path to test file
seg2_file_path(file) = joinpath(dirname(pathof(Seis)), "..", "test", "test_data", "io", file)

# List of SEG-2 files for testing
seg2_files = ("geode.seg2", "seistronix.seg2", "vipa_15.seg2", "smartseis_20bit.seg2")

@testset "SEG2" begin
    @testset "Read from disk" begin
        for file in seg2_files
            @test read_seg2(seg2_file_path(file)) isa AbstractArray{<:Seis.AbstractTrace}
        end
    end

    @testset "Headers" begin
        @testset "Geode" begin
            t = read_seg2(seg2_file_path("geode.seg2"))
            @test length(t) == 48

            @test all(x -> x == 0.00025, t.delta)
            @test all(iszero, t.b)

            @test t.sta.sta == string.(1:48)
            @test all(x -> x === missing, t.sta.net)
            @test all(x -> x === missing, t.sta.loc)
            @test all(x -> x === missing, t.sta.cha)

            @test all(t.sta.x .== (0:47).*3)
            @test all(x -> x === missing, t.sta.y)
            @test all(x -> x === missing, t.sta.z)

            @test all(t.evt.time .== DateTime("2025-06-11T12:00:07"))
            @test all(t.evt.id .== "2003")
            @test all(t.evt.x .== 141)
            @test all(x -> x === missing, t.evt.y)
            @test all(x -> x === missing, t.evt.z)

            @test all(x -> haskey(x.meta, :seg2), t)
            @test all(x -> x.meta.seg2.digital_low_cut_filter == "0 0", t)
            @test all(x -> x.meta.seg2.descaling_factor == "1.698500E-004", t)
            @test all(
                x -> x.meta.seg2.file_descriptor["INSTRUMENT"] == "GEOMETRICS SEISMODULES CONTROLLER 0000",
                t
            )
            @test all(
                x -> x.meta.seg2.file_descriptor["NOTE"] == Dict(
                      "SHOT_INCREMENT"  => "0.00",
                      "BASE_INTERVAL"   => "3.00",
                      "PHONE_INCREMENT" => "0.00",
                      "DISPLAY_FILTERS" => "0 0",
                      "AGC_WINDOW"      => "0",
                ),
                t
            )

            @test trace(t[1])[1:10] == Float32[
                -0.016963765, -0.0171823, -0.017255036, -0.017831378, -0.017445045, -0.017813183, -0.017661417, -0.017147293, -0.017137721, -0.016328318
            ]
            @test trace(t[end])[end-9:end] == Float32[
                -0.027062118, -0.00016504187, 0.03143936, 0.065828055, 0.095990926, 0.11593198, 0.124560975, 0.120522626, 0.10570981, 0.084255785
            ]
        end

        @testset "Seistronix" begin
            t = read_seg2(seg2_file_path("seistronix.seg2"))
            @test length(t) == 48

            @test all(x -> x == 0.0005, t.delta)
            @test all(iszero, t.b)

            @test t.sta.sta == string.(" ", 1:48)
            @test all(x -> x === missing, t.sta.net)
            @test all(x -> x === missing, t.sta.loc)
            @test all(x -> x === missing, t.sta.cha)

            @test all(iszero, t.sta.x)
            @test all(iszero, t.sta.y)
            @test all(iszero, t.sta.z)

            @test all(x -> x === missing, t.evt.time)
            @test all(t.evt.id .== " 30")
            @test all(t.evt.x .== 1.5)
            @test all(x -> x === missing, t.evt.y)
            @test all(x -> x === missing, t.evt.z)

            @test all(x -> haskey(x.meta, :seg2), t)
            @test all(x -> x.meta.seg2.alias_filter == " 825 78", t)
            @test all(x -> x.meta.seg2.fixed_gain == "24", t)
            @test all(
                x -> x.meta.seg2.file_descriptor["INSTRUMENT"] == "Seistronix RAS-24",
                t
            )

            @test trace(t[1])[1:10] == Float32[
                -0.0005939894, -0.0011653507, -0.00046387746,
                0.0024551563, 0.0017140837, 0.00061661756, 0.0006618739,
                -0.000118797885, -0.0006901591, -0.0009220978
            ]
            @test trace(t[end])[end-9:end] == Float32[
                0.007257985, 0.0075634653, 0.0076370067, 0.007823689,
                0.0081913965, 0.008151798, 0.00808957, 0.008106541,
                0.007942487, 0.0081178555
            ]

        end

        @testset "Vipa" begin
            t = read_seg2(seg2_file_path("vipa_15.seg2"))
            @test length(t) == 3

            @test all(x -> x == 0.001, t.delta)
            @test all(iszero, t.b)

            @test t.sta.sta == string.(1:3)
            @test all(x -> x === missing, t.sta.net)
            @test all(x -> x === missing, t.sta.loc)
            @test all(x -> x === missing, t.sta.cha)

            @test all(x -> x === missing, t.sta.x)
            @test all(x -> x === missing, t.sta.y)
            @test all(x -> x === missing, t.sta.z)

            @test all(t.evt.time .== DateTime("2013-01-07T09:30:41"))
            @test all(x -> x === missing, t.evt.id)
            @test all(x -> x === missing, t.evt.x)
            @test all(x -> x === missing, t.evt.y)
            @test all(x -> x === missing, t.evt.z)

            @test all(x -> haskey(x.meta, :seg2), t)
            @test all(x -> x.meta.seg2.sensor_type_name == "DMT-3D/DIN", t)
            @test all(x -> x.meta.seg2.scale_unit == "mm/s", t)
            @test all(t.meta.seg2.registration_direction .== ["X", "Y", "Z"])
            @test all(
                x -> x.meta.seg2.file_descriptor["INSTRUMENT"] == "DMT_VIPA_01-0000143912a3",
                t
            )

            @test trace(t[1])[1:10] == Float32[
                -0.0002391158, -0.0002825914, -0.0004782316, -0.0003912804, -0.0002391158, -0.0001521646, -0.0005869206, -0.0001956402, -0.0003043292, 0.0003478048
            ]
            @test trace(t[end])[end-9:end] == Float32[
                -8.5926f-5, -0.0002362965, -8.5926f-5, -0.0001933335, 2.14815f-5, 0.0001933335, 8.5926f-5, 0.000128889, -6.44445f-5, -0.0001503705
            ]
        end

        @testset "SmartSeis (20-bit float data)" begin
            # File obtained from https://github.com/obspy/obspy/raw/refs/heads/master/obspy/io/seg2/tests/data/20180307_031245000.0.seg2
            # at commit 4f8ff29fe3455877df1dc18c5ff495d8f59a5673
            ts = read_seg2(seg2_file_path("smartseis_20bit.seg2"))
            @test length(ts) == 1

            t = only(ts)
            @test length(trace(t)) == 2048

            @test t.delta == 0.000125

            @test t.sta.sta == "1"

            @test t.sta.x == 1004.0
            @test t.sta.y === missing
            @test t.sta.z === missing

            @test t.evt.time == DateTime("2018-03-07T03:12:45")
            @test t.evt.x == 1000.0
            @test t.evt.y === missing
            @test t.evt.z === missing

            @test t.meta.seg2.file_descriptor["INSTRUMENT"] == "GEOMETRICS SmartSeis 0000"


            @test trace(t)[1:10] == Float32[
                -0.0029975, -0.00329725, -0.004046625, -0.004796, -0.00569525, -0.005245625, -0.00629475, -0.007044125, -0.007343875, -0.007943375
            ]
        end

        @testset "Precision" begin
            @testset "$T (no warning)" for T in (Float16, Float32, Float64)
                for file in seg2_files
                    @test eltype(trace(first(read_seg2(
                        seg2_file_path(file),
                        CartTrace{Float64,Vector{T}};
                        warn=false
                    )))) == T
                end
            end

            @testset "Float16 (warning)" begin
                @test_logs(
                    (:warn, "chosen trace eltype (Float16) is smaller than file trace eltype (20-bit floating point (SEG-D)); precision will be lost"),
                    read_seg2(
                        seg2_file_path("smartseis_20bit.seg2"),
                        CartTrace{Float64,Vector{Float16}}
                    )
                )
            end
        end
    end
end
