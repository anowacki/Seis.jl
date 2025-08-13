# Spectral operations on Fourier traces

using Test
using Seis

@testset "Spectral operations" begin
    @testset "Differentiation" begin
        @testset "In-place v out-of-place" begin
            T = fft(sample_data())
            T2 = deepcopy(T)
            @test differentiate(T) == differentiate!(T2)
            @test T == fft(sample_data())
            @test T2 == differentiate(T)
        end
        
        @testset "Sine" begin
            f = 0.5
            ω = 2π*f
            npts = 512
            b = -0.5
            delta = 0.01
            taper_width = 0.5
            t_untapered = Seis.Synth.sines(b, delta, npts, [f])
            # Single sines have lots of edge effects, so taper to avoid
            t = taper(t_untapered, taper_width)
            # Analytic solution
            t′_analytic = deepcopy(t_untapered)
            trace(t′_analytic) .= ω.*cos.(ω.*times(t′_analytic))
            taper!(t′_analytic, taper_width)
            # Time-domain solution
            t′_time = differentiate(t; points=3)
            # Frequency-domain solution
            t′_freq = ifft(differentiate(fft(t)))

            @test maximum(abs, trace(t′_freq)[begin+1:end-1] .- trace(t′_time)) < 0.001
        end
    end

    @testset "Integration" begin
        @testset "In-place v out-of-place" begin
            T = fft(sample_data())
            T2 = deepcopy(T)
            @test integrate(T) == integrate!(T2)
            @test T == fft(sample_data())
            @test T2 == integrate(T)
        end

        @testset "Sine" begin
            f = 0.5
            ω = 2π*f
            npts = 512
            b = -0.5
            delta = 0.01
            t = taper(Seis.Synth.sines(b, delta, npts, [f]), 0.5)
            # Time-domain solution, with mean removed
            ∫t_time = remove_mean(integrate(t))
            # Frequency-domain solution
            ∫t_freq = ifft(integrate(fft(t)))

            @test maximum(abs, trace(∫t_freq)[begin+1:end] .- trace(∫t_time)) < 0.005
        end
    end
end