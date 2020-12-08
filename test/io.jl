using Test
using Seis
import Seis.SAC

@testset "IO" begin
    @testset "SAC" begin
        let dir = joinpath(@__DIR__(), "test_data"), file = "seis.sac",
                filepath = joinpath(dir, file),
                t = read_sac(filepath),
                s = SAC.read(filepath)
            # Reading single file
            @test t isa Trace{Float32, Vector{Float32}, Seis.Geographic{Float32}}
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
            # Origin time shifting
            @test t.b ≈ s.b - s.o
            @test t.picks.A.time ≈ s.a - s.o

            # Reading with globbing pattern
            t = read_sac(file, dir, echo=false)
            @test length(t) == 1
            t′, fmatch = SAC.read_wild(file, dir, echo=false)
            @test length(t′) == 1
            @test t[1].t == t′[1].t
            @test t.meta.file == fmatch
            t = read_sac("pattern which should match nothing", dir)
            @test isempty(t)
            @test eltype(t) == Trace{Float32, Vector{Float32}, Seis.Geographic{Float32}}

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

            # Setting file type on write
            mktemp() do file, io
                t = Trace(0, 1, rand(2))
                write_sac(t, file)
                t′ = read_sac(file)
                @test t′.meta.SAC_iftype == SAC.SAC_ITIME
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

            # Reading only headers
            let t = read_sac(filepath, header_only=true),
                    t2 = read_sac(filepath)
                empty!(trace(t2))
                @test isempty(trace(t))
                @test t == t2
            end
            let ts = read_sac("*.sac", dirname(filepath),
                        header_only=true, echo=false),
                    t2s = read_sac("*.sac", dirname(filepath), echo=false)
                empty!(trace(t2s[1]))
                @test length(ts) == 1
                @test isempty(trace(ts[1]))
                @test ts == t2s
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
            @test SAC.MACHINE_IS_LITTLE_ENDIAN != SAC.file_is_native_endian(file)
            # Write littleendian
            write_sac(t, file, littleendian=true)
            @test SAC.MACHINE_IS_LITTLE_ENDIAN == SAC.file_is_native_endian(file)
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
