using Test
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

            # Writing and reading picks
            mktempdir() do tempdir
                tempfile = joinpath(tempdir, "file.sac")
                t = Trace(0, 1, rand(10))
                t.picks.A = (1, "Pick")
                t.picks.f = 2
                t.picks.T0 = 3
                t.picks.T9 = (4, "pickname")
                t.picks.pickone = (5, "xyz")
                t.picks.picktwo = 6
                write_sac(t, tempfile)
                t′ = read_sac(tempfile)
                @test t′.picks.A.time == 1
                @test t′.picks.A.name == "Pick"
                @test t′.picks.F.time == 2
                @test ismissing(t′.picks.F.name)
                @test t′.picks.T0.time == 3
                @test t′.picks.T9 == Seis.Pick{Float32}((4, "pickname"))
                @test ismissing(t′.picks.pickone)
                @test ismissing(t′.picks.picktwo)
            end
        end

        # Removing null bytes from strings on read
        mktempdir() do dir
            file = joinpath(dir, "file.sac")
            t = sample_data()
            write_sac(t, file)
            data = read(file)
            # Modify KSTNM
            data[441:448] .= [UInt8(c) for c in "CDV\0\0\0\0\0"]
            write(file, data)
            t = read_sac(file)
            @test t.sta.sta == "CDV"
        end

        # Endianness
        mktempdir() do dir
            file = joinpath(dir, "file.sac")
            t = sample_data()
            # Write bigendian by default
            write_sac(t, file)
            @test SAC.machine_is_little_endian != SAC.file_is_native_endian(file)
            # Write littleendian
            write_sac(t, file, littleendian=true)
            @test SAC.machine_is_little_endian == SAC.file_is_native_endian(file)
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
