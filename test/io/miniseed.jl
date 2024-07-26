# Miniseed reading and writing

using Test, Seis
import Dates

io_data_path(file) = joinpath(@__DIR__, "..", "test_data", "io", file)

@testset "Miniseed" begin
    @testset "Reading" begin
        let file = "miniseed_GB.CWF.single_sample_gaps.mseed", path = io_data_path(file)
            @testset "Single file" begin
                t = read_mseed(path)
                @test length(t) == 2
                @test channel_code.(t) == "GB.CWF.." .* ["BHZ", "HHZ"]
                @test nsamples.(t) == [3000, 6000]
                @test trace(t[1])[1:10] == Float32[1466.0, 1466.0, 1453.0, 1449.0,
                    1449.0, 1443.0, 1441.0, 1443.0, 1444.0, 1439.0]
                @test trace(t[2])[1:10] == Float32[1469.0, 1469.0, 1463.0, 1465.0,
                    1447.0, 1449.0, 1457.0, 1450.0, 1447.0, 1446.0]
                @test t.delta == [0.02, 0.01]
                @test t[1].meta.mseed_file == t[2].meta.mseed_file == path
            end

            @testset "Globbing" begin
                @testset "Simple globbing" begin
                    @test read_mseed(path) == read_mseed(
                        "miniseed_GB.CWF.si??le_sample_gaps.mseed", io_data_path(""))
                end

                @testset "Globing with trace type" begin
                    @test read_mseed(path, CartTrace{Float64, Vector{Float64}}) ==
                        read_mseed("miniseed_GB.CWF.*_gaps.mseed", io_data_path(""),
                            CartTrace{Float64, Vector{Float64}})
                end

                @testset "No matches" begin
                    @test read_mseed("pattern which should not match anything",
                        io_data_path("")) == []
                end
            end

            @testset "Trace type Trace{$T}" for T in (Float32, Float64)
                @testset "$V" for V in (Vector{Float32}, Vector{Float64})
                    @testset "$P" for P in (Seis.Geographic{T}, Seis.Cartesian{T})
                        typ = Trace{T,V,P}
                        t = read_mseed(path, typ)
                        @test eltype(t) == typ
                    end
                end                
            end
        end

        @testset "Data" begin
            let path = io_data_path("miniseed_GB.CWF.single_sample_gaps.mseed"),
                    data = read(path)
                t = read_mseed(path)
                t.meta.mseed_file = missing
                @test t == read_mseed(data)
            end
        end

        @testset "Gapped data" begin
            path = io_data_path("miniseed_2sample_gap.mseed")

            @testset "Default maximum_gap" begin
                t = read_mseed(path)
                @test length(t) == 2
                @test trace(t[1]) == [1, 2, 3]
                @test trace(t[2]) == [4, 5, 6]
                @test all(x -> channel_code(x) == "AN.STA2..HHZ", t)
                @test startdate(t[2]) - enddate(t[1]) == Dates.Second(2)
            end

            @testset "No gap allowed" begin
                @test read_mseed(path, maximum_gap=0) == read_mseed(path)
            end

            @testset "Single segment from larger gap" begin
                t = read_mseed(path, maximum_gap=2)
                @test length(t) == 1
                @test nsamples(t[1]) == 6
                @test trace(t[1]) == [1, 2, 3, 4, 5, 6]
                @test channel_code(t[1]) == "AN.STA2..HHZ"
            end
        end

        @testset "Header only" begin
            path = io_data_path("miniseed_submillisecond_startdate.mseed")

            @testset "Default" begin
                @test read_mseed(path) == read_mseed(path; header_only=false)
            end

            @testset "Meta" begin                
                t = only(read_mseed(path))
                t2 = only(read_mseed(path; header_only=true))
                @test isempty(trace(t2))
                @test t2.meta.mseed_nsamples == nsamples(t)
                @test starttime(t) == starttime(t2)
                @test startdate(t) == startdate(t2)
                @test endtime(t) == t2.meta.mseed_endtime
                @test enddate(t) == t2.meta.mseed_enddate
            end

            @testset "Data" begin
                data = read(path)
                t2 = only(read_mseed(path; header_only=true))
                t3 = only(read_mseed(data; header_only=true))
                t3.meta.mseed_file = t2.meta.mseed_file
                @test t2 == t3
            end
        end

        @testset "Wrong format" begin
            @test_throws ErrorException read_mseed(io_data_path("../seis.sac"))
        end
    end

    @testset "Writing" begin
        mktempdir() do dir
            file = joinpath(dir, "test.mseed")

            @testset "Single trace" begin
                @testset "Eltype $T" for T in (Float32, Float64)
                    b = 10
                    delta = 0.1
                    t = Trace{T,Vector{T}}(b, delta, [1, 2, 3])
                    otime = Dates.DateTime(2001)
                    origin_time!(t, otime)
                    t.sta.net = "AN"
                    t.sta.sta = "ABCDE"
                    t.sta.loc = ""
                    t.sta.cha = "XYZ"
                    write_mseed(file, t)
                    t′ = only(read_mseed(file, Trace{T,Vector{T},Seis.Geographic{T}}))
                    @test t′.b == 0
                    @test t.sta == t′.sta
                    @test trace(t) == trace(t′)
                    @test starttime(t′) == 0
                    @test startdate(t′) == otime + Dates.Second(b)
                end
            end

            @testset "Multiple traces" begin
                Ts = (Float32, Float64)
                Vs = collect(Vector{T} for T in Ts)
                @testset "Output T=$Tout" for Tout in Ts
                    @testset "Output V=$Vout" for Vout in Vs
                        @testset "Input T=$Tin" for Tin in Ts
                            @testset "Output V=$Vin" for Vin in Vs
                                tout = Trace{Tout,Vout,Seis.Geographic{Tout}}
                                tin = Trace{Tin,Vin,Seis.Geographic{Tin}}
                                sta = Seis.Station{Tout}(
                                    net="A", sta="B", loc="C", cha="DEF")
                                otime = Dates.DateTime(1990)

                                @testset "Traces joined up" begin
                                    t = tout[Trace(i, 0.01, rand(100)) for i in 0:2]
                                    t.sta = sta
                                    origin_time!.(t, otime)
                                    write_mseed(file, t)
                                    t′ = read_mseed(file, tin)
                                    @test length(t′) == 1
                                    @test nsamples(t′[1]) == 300
                                    @test trace(t′[1]) ≈
                                        Tin.(reduce(vcat, trace.(t)))
                                    @test startdate(t′[1]) == otime
                                end

                                @testset "Traces separate" begin
                                    t = tout[Trace(2i, 0.01, rand(100)) for i in 0:2]
                                    t.sta = sta
                                    origin_time!.(t, otime)
                                    write_mseed(file, t)
                                    t′ = read_mseed(file, tin)
                                    @test length(t′) == 3
                                    for (tt, tt′) in zip(t, t′)
                                        @test trace(tt′) ≈ Tin.(trace(tt′))
                                        @test tt′.b == 0
                                        @test startdate(tt′) ==
                                            otime + Dates.Second(Int(tt.b))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Sub-millisecond start date" begin
        @testset "Reading" begin
            file = io_data_path("miniseed_submillisecond_startdate.mseed")
            t = only(read_mseed(file))
            @test t.b ≈ 0.000123
            @test startdate(t) == Dates.DateTime(2000)
        end

        @testset "Writing" begin
            mktemp() do path, io
                t = Trace{Float64,Vector{Float32}}(1.000456789, 1, [0.f0])
                origin_time!(t, Dates.DateTime(1999, 12, 31, 23, 59, 59, 123))
                t.sta.net, t.sta.sta, t.sta.loc, t.sta.cha =
                    "AB", "CD", "EF", "HG2"
                write_mseed(path, t)
                t′ = only(read_mseed(path, Trace{Float64, Vector{Float32}}))
                # Writing in miniSEED version 2 format only retains µs precision,
                # truncating any further digits.
                @test t′.b ≈ 0.000456 atol = 1e-9
                @test startdate(t′) == Dates.DateTime(2000) + Dates.Millisecond(123)
            end
        end
    end
end
