# Periodograms and spectrograms
using Test
using Seis

using DSP: DSP

@testset "Periodograms and spectrograms" begin
    @testset "spectrogram" begin
        @testset "ArgumentErrors" begin
            let t = sample_data()
                @test_throws ArgumentError spectrogram(t, length=-1)
                @test_throws ArgumentError spectrogram(t, overlap=-1)
                @test_throws ArgumentError spectrogram(t, overlap=1)
                @test_throws ArgumentError spectrogram(t, pad=0.5)
                @test_throws ArgumentError spectrogram(t, length=0.00001)
            end
        end

        @testset "Types" begin
            spec = spectrogram(sample_data())
            @test hasproperty(spec, :time)
            @test hasproperty(spec, :freq)
            @test hasproperty(spec, :power)
        end

        @testset "Array sizes" begin
            t = Trace(0, 0.01, rand(1000))
            len = 0.5
            overlap = 0.75
            spec = spectrogram(t, length=len, overlap=overlap, window=DSP.hamming)
            @test extrema(spec.freq) == (0, 1/2t.delta)
            @test first(spec.time) ≈ starttime(t) + len/2
        end

        @testset "Defaults" begin
            let t = sample_data()
                compare_specs(s1, s2) =
                    all(getfield(s1, f) == getfield(s2, f)
                        for f in (:time, :freq, :power))
                @test compare_specs(spectrogram(t), spectrogram(t, overlap=0.5))
                @test compare_specs(spectrogram(t),
                    spectrogram(t, length=(endtime(t) - starttime(t))/20))
                @test compare_specs(spectrogram(t), spectrogram(t, pad=1))
            end
        end
    end
end
