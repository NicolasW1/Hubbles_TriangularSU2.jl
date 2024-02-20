using Hubbles_TriangularSU2
using Test, StaticArrays

@testset "Hubbles_TriangularSU2.jl" begin
    # Write your tests here.

    @testset "Parameter" begin
        f_params_1 = Parameter(1., 0., 0., 2.)
        f_params_2 = Parameter(1.1, 0.1, -0.25, 2.14)

        @test saveDict(f_params_1) == Dict("t1" => 1., "t2" => 0., "t3" => 0., "mu" => 2.)
        @test saveDict(f_params_2) == Dict("t1" => 1.1, "t2" => 0.1, "t3" => -0.25, "mu" => 2.14)

        @testset "t₁ Parameter" begin
            t₁_params_1 = t₁Params(2.0)
            t₁_params_2 = t₁Params(1.89)

            @test t₁_params_1 == t₁Params(2)
            @test Parameter(t₁_params_1) == f_params_1
            @test Parameter(t₁_params_2) == Parameter(1., 0., 0., 1.89)

            @test saveDict(t₁_params_1) == Dict("t1" => 1., "t2" => 0., "t3" => 0., "mu" => 2.)
            @test saveDict(t₁_params_2) == Dict("t1" => 1., "t2" => 0., "t3" => 0., "mu" => 1.89)

            if Hubbles_TriangularSU2.format_param_precision == 3
                @test parameterString(t₁_params_1) == "t1_1.000_mu_2.000"
                @test parameterString(t₁_params_2) == "t1_1.000_mu_1.890"
            end
        end

        @testset "t₂ Parameter" begin
            t₂_params_1 = t₂Params(0.)
            t₂_params_2 = t₂Params(0.11)

            @test t₂_params_1 == t₂Params(0)
            @test Parameter(t₂_params_1) == Parameter(1., 0., 0., 2.)
            @test Parameter(t₂_params_2) == Parameter(1., 0.11, 0., 2.22)

            @test saveDict(t₂_params_1) == Dict("t1" => 1., "t2" => 0., "t3" => 0., "mu" => 2.)
            @test saveDict(t₂_params_2) == Dict("t1" => 1., "t2" => 0.11, "t3" => 0., "mu" => 2.22)

            if Hubbles_TriangularSU2.format_param_precision == 3
                @test parameterString(t₂_params_1) == "t2_0.000_mu_2.000"
                @test parameterString(t₂_params_2) == "t2_0.110_mu_2.220"
            end
        end
    end

    @testset "FormFactors" begin
        num_ff = [1, 19, 61, 127, 217]

        for s = 0:4
            ffs = filtered_formfactors(s)
            @test length(ffs) == 1 + (num_ff[s+1] - 1) / 2
            calc_vals = rand(ComplexF64, length(ffs))

            all = Vector{ComplexF64}(undef, num_ff[s+1])
            restore_formfactors!(all, i->calc_vals[i])

            @test all[1] == calc_vals[1]
            @test all[2:2:end] == calc_vals[2:end]
            @test all[3:2:end] == conj.(calc_vals[2:end])
        end
    end

    @testset "Output" begin
        momentum_1 = SVector(0.1 * pi, -1.25)
        momentum_2 = SVector(-0.1, 15.38)
        momentum_3 = SVector(-pi, -pi)
        momentum_4 = SVector(1., 1.)

        t₁_params_1 = t₁Params(2.0)
        ext = ExternalInput(momentum_1, nothing, 8.57, t₁_params_1)

        if Hubbles_TriangularSU2.format_momentum_precision == 6
            @test momentumString(momentum_1) == "qx_0.314159_qy_-1.250000"
            @test momentumString(momentum_2) == "qx_-0.100000_qy_15.380000"
            @test momentumString(momentum_3) == "qx_-3.141593_qy_-3.141593"
            @test momentumString(momentum_4) == "qx_1.000000_qy_1.000000"

            # This generate the output folder
            # @show outputFolder(ext)

            @test outputFileName(ParticleHole, ext) == "ParticleHole_qx_0.314159_qy_-1.250000"
        end
    end
end
