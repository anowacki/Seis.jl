using Test
using Seis
import Seis.SAC
using Statistics: mean

sample_data_path = joinpath(dirname(pathof(Seis)), "..", "data", "seis.sac")

@testset "SAC" begin
    @testset "Construction" begin
        # Construction from array
        let npts = 1000, delta = 1, b = -5
            local t = rand(SAC.SACFloat, npts)
            local s = SAC.SACTrace(t, delta)
            @test s.t == t
            @test s.delta == delta
            @test s.npts == npts

            s = SAC.SACTrace(t, delta, b)
            @test s.b == b
            @test s.e ≈ b + (npts - 1)*delta
            @test s.depmax == maximum(t)
            @test s.depmin == minimum(t)
            @test s.depmen == mean(t)
        end

        # Construction without array
        let npts = rand(1:1000), b = 2, delta = 0.1
            @test SAC.SACTrace(delta, npts).npts == npts
            @test SAC.SACTrace(delta, npts).b ≈ 0.0
            @test SAC.SACTrace(delta, npts, b).b ≈ b
        end

        ## Field modification and access
        # Single traces
        let a = rand(SAC.SACFloat, 100), delta = rand(), s = SAC.SACTrace(a, delta)
            s[:kevnm] = "ABCD"
            @test s.kevnm == "ABCD"
            @test s[:kevnm] == "ABCD"
            s[:user0] = 1
            @test s.user0 ≈ 1
            @test s[:user0] ≈ 1
            s[:t] = 1:100
            @test s.t ≈ 1:100
            @test s.t ≈ s[:t]
            s[:t] = 2:101
            @test s.t ≈ 2:101
            @test s[:depmin] ≈ 2
            @test s[:depmax] ≈ 101
            @test s[:depmen] ≈ mean(2:101)
            local k = 1
            s[:t] = k
            @test s[:depmin] ≈ k
            @test s[:depmax] ≈ k
            @test all(s[:t] .== SAC.SACFloat(k))
            k = 1:100
            s[:t] = k
            @test all(s[:t] .== SAC.SACFloat.(k))
        end
    end

    @testset "Great circles" begin
        let s = SAC.SACTrace(1, 2)
            s[:evlo], s[:evla], s[:stlo], s[:stla] = 0, 0, 12, 15
            @test s[:gcarc] ≈ 19.047431084635505
            @test s[:az] ≈ 37.992268575139384
            @test s[:baz] ≈ 219.57789440455088
        end
    end

    @testset "IO" begin
        mktemp() do tempfile, tempio

            # Reading
            let s = SAC.SACTrace(sample_data())
                # Reals
                @test s.t[1:100] ≈ Float32[-0.09728,-0.09728,-0.09856,-0.09856,-0.09728,-0.096,-0.09472,-0.09344,-0.09344,-0.09344,-0.09344,-0.09344,-0.09472,-0.09472,-0.09344,-0.09344,-0.09216,-0.09216,-0.09216,-0.09216,-0.09088,-0.09088,-0.09216,-0.09344,-0.09472,-0.09472,-0.09472,-0.09472,-0.09472,-0.09472,-0.09344,-0.09344,-0.09216,-0.09088,-0.09088,-0.09216,-0.09216,-0.09216,-0.09344,-0.09472,-0.096,-0.09856,-0.09856,-0.09856,-0.09728,-0.09728,-0.09856,-0.09984,-0.09984,-0.09984,-0.09984,-0.09984,-0.10112,-0.10112,-0.10112,-0.10112,-0.1024,-0.1024,-0.10368,-0.1024,-0.10496,-0.10496,-0.10624,-0.10368,-0.10368,-0.1024,-0.10368,-0.10368,-0.10368,-0.10368,-0.10496,-0.10624,-0.10624,-0.10496,-0.10368,-0.10368,-0.10496,-0.10624,-0.10496,-0.10496,-0.10496,-0.10624,-0.1088,-0.1088,-0.1088,-0.1088,-0.1088,-0.10624,-0.10496,-0.10496,-0.10496,-0.10624,-0.10624,-0.10624,-0.10368,-0.10368,-0.10368,-0.10368,-0.1024,-0.10112]
                # Integers
                @test s.nzyear == 1981
                # Strings
                @test s.kevnm == "K8108838"
                @test s.kstnm == "CDV"
                @test s.khole == SAC.SAC_CNULL
                # Logicals
                @test s.leven == true
                # Origin time (without shifting done in Seis.read_sac)
                @test SAC.read(sample_data_path)[:o] == -43.2f0
            end

            # Writing
            let s1 = SAC.SACTrace(sample_data())
                @test begin
                    SAC.write([s1], [tempfile])
                    s2 = SAC.read(tempfile)
                    s1 == s2
                end
                @test begin
                    s1 == SAC.SACTrace(sample_data())
                    SAC.write(s1, tempfile)
                    s2 = SAC.read(tempfile)
                    s1 == s2
                end
            end

            let s = SAC.SACTrace(sample_data())
                s.kevnm = SAC.SAC_CNULL
                SAC.write(s, tempfile)
                local d = open(tempfile, "r") do f
                    read(f)
                end
                @test String(d[449:(449+15)]) == "-12345  -12345  "
            end

            # Handling spaces in kevnm correctly
            let s = SAC.SACTrace(0.1, 2), kevnm = "1234567 910"
                s.kevnm = kevnm
                SAC.write(s, tempfile)
                local s2 = SAC.read(tempfile)
                @test s2.kevnm == kevnm
            end

            # Cutting
            let
                # Shift O time to be 0
                SAC.write(SAC.SACTrace(sample_data()), tempfile)
                local s1 = SAC.read(tempfile)
                s1 = SAC.SACTrace(cut!(Trace(s1), 54, 55))
                local s2 = SAC.read_cut(tempfile, 54, 55)
                @test s1 == s2
            end
        end
    end
end
