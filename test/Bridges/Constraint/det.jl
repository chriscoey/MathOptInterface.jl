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


    # opt.x = [1.609437887625386, 2.999999981363997, 1.9999999845086225, 0.6666666696605013, 0.9999999922779824, 0.3333333325096238, 2.4999999868779224, 1.0986122771846563, -0.40546510530752866, 0.9162907216367258]
    # opt.s = [-5.8884671951948826e-9, 1.0986122771846563, 0.9999999994264411, 2.9999999834775517, -0.40546510530752866, 0.9999999994264411, 0.6666666717740563, 0.9162907216367258, 0.9999999994264411, 2.4999999889914775, 2.9999999941142654, 1.9999999941142654, 1.9999999970571327, 0.9999999970571327, 0.9999999970571327, 2.9999999941142654, 2.9999999796785226, 1.9999999845086225, 0.9999999922779825, 2.9999999826213903, 0.0, 0.666666667975027, 0.3333333325096239, 0.0, 0.6666666709178946, 0.0, 0.0, 2.4999999851924484, 0.0, 0.0, 2.4999999881353157]
    # opt.y = Float64[]
    # opt.z = [-1.000000000069177,
    # -0.9999999994369781, 0.09861229209769833, 0.3333333362978019,
    # -0.9999999994194586, -1.4054650936159223, 1.4999999969291422,
    # -0.9999999994612729, -0.0837092604550414, 0.4000000016715615,

    # 0.9999999979988198, -0.9999999994818665, 1.5999999980602015, -2.601097313972902e-10, -0.1999999983277935, 0.39999999742064857, -0.33333333111240304, 3.315616913751688e-11, -1.2768924444409872e-11, 0.33333333276122146, 0.9999999981404057, -1.499999996538671, -1.731699247129028e-11, -1.9040558981184468e-10, 1.5000000034702832, 2.6154790380704195e-10, 0.19999999815264463, -0.39999999713014683, 1.2155316527092914e-11, 3.2402157386207174e-11, 0.3999999993375813]

    @testset "logdet2test" begin
        var_primal = [log(5), 3, 2, 2/3, 1, 1/3, 2.5, log(3), log(2/3), log(2.5)]
        exp_duals = [[-1, log(3) - 1, 1/3], [-1, log(2/3) - 1, 3/2], [-1, log(2.5) - 1, 0.4]]
        psd_dual = [1, -1, 1.6, 0, -0.2, 0.4, -1/3, 0, 0, 1/3, 1, -1.5, 0, 0, 1.5, 0, 0.2, -0.4, 0, 0, 0.4]
        mock.optimize! = (mock::MOIU.MockOptimizer) -> MOIU.mock_optimize!(mock, var_primal
            (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}) => [-1.0],
            (MOI.VectorAffineFunction{Float64}, MOI.ExponentialCone) => exp_duals,
            (MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle) => [psd_dual])

        MOIT.logdett1vtest(bridged_mock, config)
        MOIT.logdett1ftest(bridged_mock, config)

        # set primal/dual start is not yet implemented for LogDet bridge
        ci = first(MOI.get(bridged_mock, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, MOI.LogDetConeTriangle}()))
        test_delete_bridge(bridged_mock, ci, 5, (
            (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}, 0)
            (MOI.VectorAffineFunction{Float64}, MOI.ExponentialCone, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle, 0),
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
