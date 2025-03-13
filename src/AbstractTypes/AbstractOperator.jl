"""
    AbstractOperator

An abstract type representing a general operator. This serves as a base type for defining various specific operators in the context of the Green Function Monte Carlo project.
"""
abstract type AbstractOperator end
abstract type OffdiagonalOperator end

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
function get_moves end

function (D::DiagonalOperator) end

"""
    AbstractSignFreeOperator <: AbstractOperator

An abstract type representing a sign-free operator in the context of Green Function Monte Carlo simulations. 
# Interface:
- get_diagonal(O::AbstractSignFreeOperator): return the diagonal operator associated with the sign-free operator `O`
- get_offdiagonal_elements(O::AbstractSignFreeOperator): return the weights associated with the off-diagonal operator `O`
"""
abstract type AbstractSignFreeOperator <: AbstractOperator end
"""
    get_offdiagonal_elements(O::AbstractSignFreeOperator)

Return the weights associated with an `OffdiagonalOperator` object `O`. 
"""
function get_offdiagonal_elements end

function get_diagonal end