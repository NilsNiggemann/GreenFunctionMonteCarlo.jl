"""
    struct EqualWeightSuperposition <: AbstractGuidingFunction

Represents a guiding function that models an equal weight superposition of states, i.e. ψ(x) = 1 for all x, provided that the configuration x is a valid configuration satisfying all constraints.
This structure is a subtype of `AbstractGuidingFunction` and is used in the context
of variational calculations within the Green Function Monte Carlo framework.

# See Also
- `AbstractGuidingFunction`: The abstract type that this struct extends.
"""
struct EqualWeightSuperposition <: AbstractGuidingFunction end
(psi::EqualWeightSuperposition)(x::AbstractArray) = 1.

guidingfunc_name(F::EqualWeightSuperposition) = "EqualWeightSuperposition"
get_params(ψG::EqualWeightSuperposition) = zeros(0)
allocate_GWF_buffer(logψ::EqualWeightSuperposition,conf) = EmptyGWFBuffer()

function log_psi_diff(x::AbstractArray, move::AbstractMove, logψ::EqualWeightSuperposition, Buffer::AbstractGuidingFunctionBuffer, Hilbert::AbstractHilbertSpace)
    isapplicable(x,move, Hilbert) || return -Inf
    return 0.
end