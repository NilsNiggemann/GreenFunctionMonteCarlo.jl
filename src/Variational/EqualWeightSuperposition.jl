struct EqualWeightSuperposition <: AbstractGuidingFunction end
(psi::EqualWeightSuperposition)(x::Any) = 1.

guidingfunc_name(F::EqualWeightSuperposition) = "EqualWeightSuperposition"
get_params(ψG::EqualWeightSuperposition) = zeros(0)
function log_psi_diff(x, dx::AbstractMove,logψ::EqualWeightSuperposition, Buffer::AbstractGuidingFunctionBuffer,Hilbert::AbstractHilbertSpace)
    isapplicable(dx,x) || return -Inf
    return 0.
end
allocate_GWF_buffer(ψG::EqualWeightSuperposition,S::StencilSpinConfig) = EmptyGuidingFunctionBuffer()