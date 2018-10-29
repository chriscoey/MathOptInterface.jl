```@meta
CurrentModule = MathOptInterface
```

# API Reference

[Some introduction to API. List basic standalone methods.]

## Attributes

List of attribute categories.

```@docs
AbstractOptimizerAttribute
AbstractModelAttribute
AbstractVariableAttribute
AbstractConstraintAttribute
```

Attributes can be set in different ways:

* it is either set when the model is created like [`SolverName`](@ref) and
  [`RawSolver`](@ref),
* or explicitly when the model is copied like [`ObjectiveSense`](@ref),
* or implicitly, e.g., [`NumberOfVariables`](@ref) is implicitly set by
  [`add_variable`](@ref) and [`ConstraintFunction`](@ref) is implicitly set by
  [`add_constraint`](@ref).
* or it is set to contain the result of the optimization during
  [`optimize!`](@ref) like [`VariablePrimal`](@ref).

The following functions allow to distinguish between some of these different
categories:
```@docs
is_set_by_optimize
is_copyable
```

Functions for getting and setting attributes.

```@docs
get
get!
set
supports
```

## Model Interface

```@docs
ModelLike
isempty
empty!
write_to_file
read_from_file
```

Copying

```@docs
copy_to
```

List of model attributes

```@docs
Name
ObjectiveSense
NumberOfVariables
ListOfVariableIndices
ListOfConstraints
NumberOfConstraints
ListOfConstraintIndices
ListOfOptimizerAttributesSet
ListOfModelAttributesSet
ListOfVariableAttributesSet
ListOfConstraintAttributesSet
```


## Optimizers

```@docs
AbstractOptimizer
optimize!
```

List of attributes optimizers attributes

```@docs
SolverName
```

List of attributes useful for optimizers


```@docs
RawSolver
ResultCount
ObjectiveFunction
ObjectiveFunctionType
ObjectiveValue
ObjectiveBound
RelativeGap
SolveTime
SimplexIterations
BarrierIterations
NodeCount
TerminationStatus
PrimalStatus
DualStatus
```

### Termination Status

The `TerminationStatus` attribute indicates why the optimizer stopped executing.
The value of the attribute is of type `TerminationStatusCode`.

```@docs
TerminationStatusCode
```

### Result Status

The `PrimalStatus` and `DualStatus` attributes indicate how to interpret the result returned by the solver.
The value of the attribute is of type `ResultStatusCode`.

```@docs
ResultStatusCode
```

## Variables and Constraints

### Basis Status

The `BasisStatus` attribute of a variable or constraint describes its status with respect to a basis, if one is known.
The value of the attribute is of type `BasisStatusCode`.

```@docs
BasisStatusCode
```

### Index types

```@docs
VariableIndex
ConstraintIndex
is_valid
delete(::ModelLike, ::Index)
```

### Variables

Functions for adding variables. For deleting, see index types section.

```@docs
add_variables
add_variable
```

List of attributes associated with variables. [category AbstractVariableAttribute]
Calls to `get` and `set` should include as an argument a single `VariableIndex` or a vector of `VariableIndex` objects.

```@docs
VariableName
VariablePrimalStart
VariablePrimal
VariableBasisStatus
```

### Constraints

Functions for adding and modifying constraints.

```@docs
is_valid(::ModelLike,::ConstraintIndex)
add_constraint
add_constraints
transform
supports_constraint
```

List of attributes associated with constraints. [category AbstractConstraintAttribute]
Calls to `get` and `set` should include as an argument a single `ConstraintIndex` or a vector of `ConstraintIndex{F,S}` objects.

```@docs
ConstraintName
ConstraintPrimalStart
ConstraintDualStart
ConstraintPrimal
ConstraintDual
ConstraintBasisStatus
ConstraintFunction
ConstraintSet
```

## Functions and function modifications

List of recognized functions.
```@docs
AbstractFunction
SingleVariable
VectorOfVariables
ScalarAffineTerm
ScalarAffineFunction
VectorAffineTerm
VectorAffineFunction
ScalarQuadraticTerm
ScalarQuadraticFunction
VectorQuadraticTerm
VectorQuadraticFunction
```

Functions for getting and setting properties of sets.

```@docs
output_dimension
```

## Sets

List of recognized sets.

```@docs
AbstractSet
Reals
Zeros
Nonnegatives
Nonpositives
GreaterThan
LessThan
EqualTo
Interval
SecondOrderCone
RotatedSecondOrderCone
GeometricMeanCone
ExponentialCone
DualExponentialCone
PowerCone
DualPowerCone
PositiveSemidefiniteConeTriangle
PositiveSemidefiniteConeSquare
LogDetConeTriangle
LogDetConeSquare
RootDetConeTriangle
RootDetConeSquare
Integer
ZeroOne
Semicontinuous
Semiinteger
SOS1
SOS2
```

Functions for getting and setting properties of sets.

```@docs
dimension
```

## Modifications

Functions for modifying objective and constraint functions.

```@docs
modify
AbstractFunctionModification
ScalarConstantChange
VectorConstantChange
ScalarCoefficientChange
MultirowChange
```

## Nonlinear programming (NLP)

### Attributes

```@docs
NLPBlock
NLPBoundsPair
NLPBlockData
NLPBlockDual
NLPBlockDualStart
```
### NLP evaluator methods

```@docs
AbstractNLPEvaluator
initialize
features_available
eval_objective
eval_constraint
eval_objective_gradient
jacobian_structure
hessian_lagrangian_structure
eval_constraint_jacobian
eval_constraint_jacobian_product
eval_constraint_jacobian_transpose_product
eval_hessian_lagrangian
eval_hessian_lagrangian_product
objective_expr
constraint_expr
```

## Errors

When an MOI call fails on a model, precise errors should be thrown when possible
instead of simply calling `error` with a message. The docstrings for the
respective methods describe the errors that the implementation should thrown in
certain situations. This error-reporting system allows code to distinguish
between internal errors (that should be shown to the user) and unsupported
operations which may have automatic workarounds.

When an invalid index is used in an MOI call, an [`InvalidIndex`](@ref) should
be thrown:
```@docs
InvalidIndex
```

The rest of the errors defined in MOI fall in two categories represented by the
following two abstract types:
```@docs
UnsupportedError
NotAllowedError
```

The different [`UnsupportedError`](@ref) and [`NotAllowedError`](@ref) are the
following errors:
```@docs
UnsupportedAttribute
SetAttributeNotAllowed
AddVariableNotAllowed
UnsupportedConstraint
AddConstraintNotAllowed
ModifyConstraintNotAllowed
ModifyObjectiveNotAllowed
DeleteNotAllowed
```

## Bridges

Bridges can be used for automatic reformulation of a certain constraint type into equivalent constraints.
```@docs
Bridges.AbstractBridge
Bridges.AbstractBridgeOptimizer
Bridges.SingleBridgeOptimizer
Bridges.LazyBridgeOptimizer
Bridges.add_bridge
```

Below is the list of bridges implemented in this package.
```@docs
Bridges.SplitIntervalBridge
Bridges.RSOCBridge
Bridges.GeoMeanBridge
Bridges.SquarePSDBridge
Bridges.RootDetBridge
Bridges.LogDetBridge
Bridges.SOCtoPSDBridge
Bridges.RSOCtoPSDBridge
```
For each bridge defined in this package, a corresponding bridge optimizer is available with the same name without the "Bridge" suffix, e.g., `SplitInterval` is an `SingleBridgeOptimizer` for the `SplitIntervalBridge`.

## Allocate-Load API

The Allocate-Load API allows to implement [`copy_to`](@ref) in a way that still
allows transformations to be applied in the copy between the cache and the
model if the transformations are implemented as MOI layers implemented the
Allocate-Load API.

Allocate-Load Interface: 2-pass copy of a MathOptInterface model
Some solver wrappers (e.g. SCS, ECOS, SDOI) do not supporting copying an optimization model using `MOI.add_constraints`, `MOI.add_variables` and `MOI.set`
as they first need to figure out some information about a model before being able to pass the problem data to the solver.

During the first pass (called allocate) : the model collects the relevant information about the problem so that
on the second pass (called load), the constraints can be loaded directly to the solver (in case of SDOI) or written directly into the matrix of constraints (in case of SCS and ECOS).

To support `MOI.copy_to` using this 2-pass mechanism, implement the allocate-load interface defined below and do:
MOI.copy_to(dest::ModelType, src::MOI.ModelLike) = MOIU.allocate_load(dest, src)
In the implementation of the allocate-load interface, it can be assumed that the different functions will the called in the following order:
1) `allocate_variables`
2) `allocate` and `allocate_constraint`
3) `load_variables` and `allocate_constraint`
4) `load` and `load_constraint`
The interface is not meant to be used to create new constraints with `allocate_constraint` followed by `load_constraint` after a solve, it is only meant for being used in this order to implement `MOI.copy_to`.


```@docs
allocate_variables
allocate
allocate_constraint
load_variables
load
load_constraint
allocate_load
```
