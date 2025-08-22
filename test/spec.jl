# Fourier traces

using Seis
using Test

@testset "FourierTrace" begin
    @testset "fft/ifft" begin
        t = sample_data()
        f = fft(t)
        t′ = ifft(f)

        @testset "Field references" begin
            @testset "Trace to FourierTrace" begin
                @test t.b == f.b
                @test t.evt === f.evt
                @test t.sta === f.sta
                @test t.meta === f.meta
            end

            @testset "Trace to ifft(FourierTrace)" begin
                @test t.evt === t′.evt
                @test t.sta === t′.sta
                @test t.meta === t′.meta
            end
        end

        @testset "nsamples $odd_even" for odd_even in (:odd, :even)
            t2, t2′ = if odd_even == :odd
                t2 = deepcopy(t)
                pop!(trace(t2))
                t2, ifft(fft(t2))
            else
                t, t′
            end

            @testset "Round-trip" begin
                for p in propertynames(t2)
                    v, v′ = getfield.((t2, t2′), p)
                    if p == :t
                        @test v ≈ v′
                    else
                        @test isequal(v, v′)
                    end
                end
            end
        end
    end
end
