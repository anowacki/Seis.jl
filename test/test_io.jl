using Compat.Test
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
end