using Test
using Seis

@testset "Sample data" begin
    @test_throws ArgumentError sample_data(:nonsense_kind_of_data)
    let t1 = sample_data(), t2 = read_sac(joinpath(@__DIR__(), "..", "data", "seis.sac"))
        t1.meta.file = t2.meta.file = ""
        @test t1 == t2
    end
    let t = sample_data(:regional)
        @test length(t) == 12
        @test t.sta.sta == ["ELK", "ELK", "ELK", "KNB", "KNB", "KNB",
                            "LAC", "LAC", "LAC", "MNV", "MNV", "MNV"]
        @test t.sta.cha == repeat(["e", "n", "z"], 4)
    end
    @test length(sample_data(:array)) == 60
end
