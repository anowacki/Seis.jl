using Test
import Dates
using Seis
import SAC

@testset "IO" begin
    @testset "SAC" begin
        let dir = joinpath(@__DIR__(), "test_data"), file = "seis.sac",
                filepath = joinpath(dir, file),
                t = read_sac(filepath),
                s = SAC.read(filepath)
            # Reading single file
            @test t.t == s.t
            @test t.evt.id == "K8108838"
            @test t.sta.sta == "CDV"
            @test t.meta.SAC_lpspol   == true
            @test t.meta.SAC_nevid    == 0
            @test t.meta.SAC_iftype   == 1
            @test t.meta.file         == joinpath(dir, file)
            @test t.meta.SAC_idep     == 50
            @test t.meta.SAC_iztype   == 9
            @test t.meta.SAC_lcalda   == true
            @test t.meta.SAC_unused18 == false
            @test t.meta.SAC_lovrok   == true
            @test t.meta.SAC_norid    == 0
            @test t.meta.SAC_ievtyp   == 42

            # Reading with globbing pattern
            t = read_sac(file, dir, echo=false)
            @test length(t) == 1
            t′, fmatch = SAC.read_wild(file, dir, echo=false)
            @test length(t′) == 1
            @test t[1].t == t′[1].t
            @test t.meta.file == fmatch

            # Writing
            mktempdir() do tempdir
                t = read_sac(filepath)
                tempfile = joinpath(tempdir, "file.sac")
                write_sac(t, tempfile)
                t′ = read_sac(tempfile)
                delete!(t.meta, :file)
                delete!(t′.meta, :file)
                @test t == t′
            end
        end
    end

    @testset "miniSEED" begin
        let file = joinpath(@__DIR__(), "test_data", "test.mseed")
            tt = read_mseed(file)
            @test tt isa Vector{Trace{Seis.DEFAULT_FLOAT,Vector{Seis.DEFAULT_FLOAT},Seis.DEFAULT_STRING}}
            @test length(tt) == 1
            t = tt[1]
            @test channel_code(t) == "NL.HGN.00.BHZ"
            @test t.b == 0
            @test t.delta == 0.025
            @test t.evt.time == Dates.DateTime(2003, 05, 29, 02, 13, 22, 43)
            @test nsamples(t) == 11947
            @test t.meta.MSEED_encoding == "STEIM2"
            @test t.meta.MSEED_dataquality == "R"
            @test t.meta.MSEED_filesize == 8192
            @test t.meta.MSEED_record_length == 4096
            @test t.meta.MSEED_number_of_records == 2
            @test t.meta.MSEED_byteorder == ">"
            @test t.meta.file == file
            t.meta.file = missing
            @test tt == parse_mseed(read(file))
        end
    end

    @testset "Channel code" begin
        let t = Trace(-10, 0.0001, [0,1])
            t.sta.net = "AB"
            t.sta.sta = "CDEF"
            t.sta.loc = "GH"
            t.sta.cha = "IJK"
            @test channel_code(t) == "AB.CDEF.GH.IJK"
            @test channel_code(t.sta) == "AB.CDEF.GH.IJK"
        end
    end
end
