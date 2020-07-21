using Test
using Seis
using Seis.Plot

@testset "Plotting" begin
    @testset "_intercept" begin
        @test Seis.Plot._intercept(1, 1, 2, 2, 1.5) == 1.5
        @test Seis.Plot._intercept(2, 2, 1, 1, 1.5) == 1.5
    end

    @testset "_below_above" begin
        @test isequal(Seis.Plot._below_above(0:2, [1, 3, 5], 3),
                    ([0, 1, NaN], [1, 3, NaN], [1, 2], [3, 5]))
        @test isequal(Seis.Plot._below_above([-1, 2, 4], [-1, -1, 5], -1),
            ([], [], [2,4], [-1, 5]))
        @test isequal(Seis.Plot._below_above([2, 4, 6], [0, 1, 0], 0.5),
            ([2, 3, NaN, 5, 6], [0, 0.5, NaN, 0.5, 0], [3, 4, 5, NaN], [0.5, 1, 0.5, NaN]))
    end
end
