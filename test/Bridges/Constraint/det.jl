using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MathOptInterface.Test
const MOIU = MathOptInterface.Utilities
const MOIB = MathOptInterface.Bridges

include("../utilities.jl")

mock = MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()))
config = MOIT.TestConfig()

@testset "LogDet" begin
    bridged_mock = MOIB.Constraint.LogDet{Float64}(mock)

    @testset "logdet1test" begin
        var_primal = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1]
        mock.optimize! = (mock::MOIU.MockOptimizer) -> MOIU.mock_optimize!(mock, var_primal)

        MOIT.logdett1vtest(bridged_mock, config)
        MOIT.logdett1ftest(bridged_mock, config)

        # set primal/dual start is not yet implemented for LogDet bridge
        ci = first(MOI.get(bridged_mock, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, MOI.LogDetConeTriangle}()))
        test_delete_bridge(bridged_mock, ci, 5, (
            (MOI.VectorAffineFunction{Float64}, MOI.ExponentialCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle, 0),
            (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}, 0)
        ))
    end

    @testset "logdet2test" begin
        var_primal = [log(5)]
        mock.optimize! = (mock::MOIU.MockOptimizer) -> MOIU.mock_optimize!(mock, var_primal)

        MOIT.logdett1vtest(bridged_mock, config)
        MOIT.logdett1ftest(bridged_mock, config)

        # set primal/dual start is not yet implemented for LogDet bridge
        ci = first(MOI.get(bridged_mock, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, MOI.LogDetConeTriangle}()))
        test_delete_bridge(bridged_mock, ci, 5, (
            (MOI.VectorAffineFunction{Float64}, MOI.ExponentialCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle, 0),
            (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}, 0)
        ))

    end
end

@testset "RootDet" begin
    bridged_mock = MOIB.Constraint.RootDet{Float64}(mock)

    @testset "rootdet1test" begin
        var_primal = [1, 1, 0, 1, 1, 0, 1]
        mock.optimize! = (mock::MOIU.MockOptimizer) -> MOIU.mock_optimize!(mock, var_primal)

        MOIT.rootdett1vtest(bridged_mock, config)
        MOIT.rootdett1ftest(bridged_mock, config)

        # set primal/dual start is not yet implemented for RootDet bridge
        ci = first(MOI.get(bridged_mock, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, MOI.RootDetConeTriangle}()))
        test_delete_bridge(bridged_mock, ci, 4, (
            (MOI.VectorAffineFunction{Float64}, MOI.RotatedSecondOrderCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.GeometricMeanCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle, 0)
        ))
    end

    @testset "rootdet2test" begin
        var_primal = [5 ^ inv(3)]
        mock.optimize! = (mock::MOIU.MockOptimizer) -> MOIU.mock_optimize!(mock, var_primal)

        MOIT.rootdett1vtest(bridged_mock, config)
        MOIT.rootdett1ftest(bridged_mock, config)

        # set primal/dual start is not yet implemented for RootDet bridge
        ci = first(MOI.get(bridged_mock, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, MOI.RootDetConeTriangle}()))
        test_delete_bridge(bridged_mock, ci, 4, (
            (MOI.VectorAffineFunction{Float64}, MOI.RotatedSecondOrderCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.GeometricMeanCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle, 0)
        ))
    end
end
