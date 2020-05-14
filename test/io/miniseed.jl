# Miniseed reading and writing

using Test, Seis

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
                @test t.delta == Float32[0.02, 0.01]
                @test t[1].meta.mseed_file == t[2].meta.mseed_file == path
            end

            @testset "Data" begin
                let path = io_data_path("miniseed_GB.CWF.single_sample_gaps.mseed"),
                        data = read(path)
                    t = read_mseed(path)
                    t.meta.mseed_file = missing
                    @test t == read_mseed(data)
                end
            end

            @testset "Globbing" begin
                @test read_mseed(path) == read_mseed(
                    "miniseed_GB.CWF.si??le_sample_gaps.mseed", io_data_path(""))
                @test read_mseed("pattern which should not match anything",
                    io_data_path("")) == []
            end

            @testset "Gapped data" begin
                @test length(read_mseed(path)) == 2
                t = read_mseed(path, maximum_gap=0)
                @test length(t) == 4
                @test count(t.sta.cha .== "HHZ") == 3
                @test sum(nsamples.(t[t.sta.cha .== "HHZ"])) == 6000
                @test read_mseed(path, maximum_offset=0) == read_mseed(path)
            end

            @testset "Wrong format" begin
                @test_throws ErrorException read_mseed(
                    joinpath(@__DIR__, "..", "test_data", "seis.sac"))
            end
        end
    end

    @testset "Writing" begin
        
    end
end
