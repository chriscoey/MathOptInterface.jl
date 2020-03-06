using Test
import MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

@testset "Default objective sense" begin
    MOIT.default_objective_test(MOIU.MockOptimizer(MOIU.Model{Float64}()))
end

@testset "Default statuses" begin
    model = MOIU.MockOptimizer(MOIU.Model{Float64}())
    MOIT.default_status_test(model)
    MOI.empty!(model)
    MOIT.default_status_test(model)
end

@testset "Name test" begin
    MOIT.nametest(MOIU.MockOptimizer(MOIU.Model{Float64}()))
end

@testset "Optimizer attributes" begin
    optimizer = MOIU.MockOptimizer(MOIU.Model{Float64}())
    @test MOI.supports(optimizer, MOIU.MockModelAttribute())
    MOI.set(optimizer, MOIU.MockModelAttribute(), 10)
    @test MOI.get(optimizer, MOIU.MockModelAttribute()) == 10

    v1 = MOI.add_variable(optimizer)
    @test MOI.supports(optimizer, MOIU.MockVariableAttribute(), typeof(v1))
    MOI.set(optimizer, MOIU.MockVariableAttribute(), v1, 11)
    @test MOI.get(optimizer, MOIU.MockVariableAttribute(), v1) == 11
    MOI.set(optimizer, MOIU.MockVariableAttribute(), [v1], [-11])
    @test MOI.get(optimizer, MOIU.MockVariableAttribute(), [v1]) == [-11]

    @test MOI.supports_constraint(optimizer, MOI.SingleVariable, MOI.GreaterThan{Float64})
    c1 = MOI.add_constraint(optimizer, MOI.SingleVariable(v1), MOI.GreaterThan(1.0))
    @test MOI.supports(optimizer, MOIU.MockConstraintAttribute(), typeof(c1))
    MOI.set(optimizer, MOIU.MockConstraintAttribute(), c1, 12)
    @test MOI.get(optimizer, MOIU.MockConstraintAttribute(), c1) == 12
    MOI.set(optimizer, MOIU.MockConstraintAttribute(), [c1], [-12])
    @test MOI.get(optimizer, MOIU.MockConstraintAttribute(), [c1]) == [-12]
end

@testset "Optimizer solve no result" begin
    optimizer = MOIU.MockOptimizer(MOIU.Model{Float64}())

    v1 = MOI.add_variable(optimizer)

    # Load fake solution
    MOI.set(optimizer, MOI.TerminationStatus(), MOI.INFEASIBLE)

    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.INFEASIBLE
end

@testset "Optimizer solve with result" begin
    optimizer = MOIU.MockOptimizer(MOIU.Model{Float64}(),
                                   eval_objective_value=false,
                                   eval_variable_constraint_dual=false)

    v = MOI.add_variables(optimizer, 2)
    c1 = MOI.add_constraint(optimizer, MOI.SingleVariable(v[1]), MOI.GreaterThan(1.0))
    soc = MOI.add_constraint(optimizer, MOI.VectorOfVariables(v), MOI.SecondOrderCone(2))
    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(v[1]))
    MOI.set(optimizer, MOI.ResultCount(), 1)
    @test_throws(
        ErrorException("No mock primal is set for variable `$(v[1])`."),
        MOI.get(optimizer, MOI.VariablePrimal(), v[1])
    )
    @test_throws(
        MOI.InvalidIndex(MOI.VariableIndex(-1)),
        MOI.get(optimizer, MOI.VariablePrimal(), MOI.VariableIndex(-1))
    )
    MOI.set(optimizer, MOI.ResultCount(), 0)
    err = MOI.ResultIndexBoundsError(MOI.VariablePrimal(1), 0)
    @test_throws err MOI.get(optimizer, MOI.VariablePrimal(), v[1])
    err = MOI.ResultIndexBoundsError(MOI.ConstraintDual(1), 0)
    @test_throws err MOI.get(optimizer, MOI.ConstraintDual(), c1)

    # Load fake solution
    # TODO: Provide a more compact API for this.
    MOI.set(optimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(optimizer, MOI.ObjectiveValue(), 1.0)
    MOI.set(optimizer, MOI.DualObjectiveValue(2), 2.0)
    MOI.set(optimizer, MOI.ResultCount(), 2)
    MOI.set(optimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(optimizer, MOI.DualStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(optimizer, MOI.VariablePrimal(), v, [1.0, 2.0])
    MOI.set(optimizer, MOI.VariablePrimal(), v[1], 3.0)
    MOI.set(optimizer, MOI.ConstraintDual(2), c1, 5.9)
    MOI.set(optimizer, MOI.ConstraintDual(2), soc, [1.0,2.0])

    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.ResultCount()) == 2
    @test MOI.get(optimizer, MOI.ObjectiveValue()) == 1.0
    @test isnan(MOI.get(optimizer, MOI.ObjectiveValue(2)))
    optimizer.eval_objective_value = true
    @test MOI.get(optimizer, MOI.ObjectiveValue()) == 3.0
    @test_throws(
        ErrorException("No mock primal is set for variable `$(v[1])` at result index `2`."),
        MOI.get(optimizer, MOI.ObjectiveValue(2))
    )
    @test_throws(
        ErrorException("No mock dual is set for constraint `$c1` at result index `1`."),
        MOI.get(optimizer, MOI.DualObjectiveValue())
    )
    @test MOI.get(optimizer, MOI.DualObjectiveValue(2)) == 5.9
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.DualStatus(2)) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.VariablePrimal(), v) == [3.0, 2.0]
    @test MOI.get(optimizer, MOI.VariablePrimal(), v[1]) == 3.0
    @test_throws(
        ErrorException("No mock primal is set for variable `$(v[1])` at result index `2`."),
        MOI.get(optimizer, MOI.VariablePrimal(2), v[1])
    )
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), c1) == 3.0
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), soc) == [3.0, 2.0]
    @test_throws(
        ErrorException("No mock primal is set for variable `$(v[1])` at result index `2`."),
        @show MOI.get(optimizer, MOI.ConstraintPrimal(2), c1)
    )
    @test_throws(
        ErrorException("No mock primal is set for variable `$(v[1])` at result index `2`."),
        MOI.get(optimizer, MOI.ConstraintPrimal(2), soc)
    )
    @test MOI.get(optimizer, MOI.ConstraintDual(2), c1) == 5.9
    @test MOI.get(optimizer, MOI.ConstraintDual(2), soc) == [1.0,2.0]
    @test_throws(
        ErrorException("No mock dual is set for constraint `$c1` at result index `1`."),
        MOI.get(optimizer, MOI.ConstraintDual(1), c1)
    )
end

@testset "Delete" begin
    mock = MOIU.MockOptimizer(MOIU.Model{Float64}())
    MOIT.delete_test(mock)
end
