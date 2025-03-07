"""
    AbstractOperator

An abstract type representing a general operator. This serves as a base type for defining various specific operators in the context of the Green Function Monte Carlo project.
"""
abstract type AbstractOperator end

"""
    DiagonalOperator

An abstract type representing a diagonal operator in the context of Green Function Monte Carlo simulations.
A diagonal operator is special in the sense that it will not change the configuration of the system when applied to it and will only return a number.
# Interface: 
- (D::DiagonalOperator)(x) : return the value of the operator applied to the configuration `x`
- (D::DiagonalOperator)(x, params): return the value of the operator applied to the configuration `x` with parameters `params`, which can be used to store buffers.
"""
abstract type DiagonalOperator end

"""
    get_moves(O::AbstractOperator)
Retrieves the possible moves for a given `AbstractOperator` object `O`. For general local operators, the length should be proportional to the number of sites in the system.
"""
get_moves(O::AbstractOperator) = throw(MethodError(get_moves, (O,)))

"""
    get_weights(O::AbstractOperator)

Return the weights associated with an `AbstractOperator` object `O`. 
"""
get_weights(O::AbstractOperator) = throw(MethodError(get_weights, (O,)))

(D::DiagonalOperator)(x) = throw(MethodError(D, (x,)))
(D::DiagonalOperator)(x, params) = D(x)