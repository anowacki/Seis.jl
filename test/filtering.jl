# Filtering
using Test
import DSP, FFTW
using Seis

function compare_traces(test_file, f, args...; rtol=1e-4, kwargs...)
    impulse = read_sac(joinpath(@__DIR__, "test_data", "filtering", "impulse.sac"))
    reference = read_sac(joinpath(@__DIR__, "test_data", "filtering", test_file))
    isapprox(f(impulse, args...; kwargs...).t, reference.t, rtol=rtol)
end

@testset "Filtering" begin
    @testset "Bandpass" begin
        @test compare_traces("impulse_bp_bu_c_0.01-0.1_npoles_4_passes_2.sac",
            bandpass, 0.01, 0.1; poles=4, twopass=true)
    end

    @testset "Bandstop" begin
        @test compare_traces("impulse_br_bu_c_0.01-0.1_npoles_6_passes_1.sac",
            bandstop, 0.01, 0.1; poles=6)
    end

    @testset "Highpass" begin
        @test compare_traces("impulse_hp_bu_c_0.01_npoles_5_passes_1.sac",
            highpass, 0.01; poles=5)
    end

    @testset "Lowpass" begin
        @test compare_traces("impulse_lp_bu_c_0.01_npoles_3_passes_1.sac",
            lowpass, 0.01; poles=3)
    end
end
