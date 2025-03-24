"""provides a naive wrapper for a guiding function which does not use any buffer. Useful for debugging and testing"""
struct NaiveFunction{T<: AbstractGuidingFunction} <: AbstractGuidingFunction
    logpsi::T
end
(N::NaiveFunction)(x::Any) = N.logpsi(x)

guidingfunc_name(F::NaiveFunction) = "NaiveFunction"
get_params(ψG::NaiveFunction) = get_params(ψG.logpsi)
allocate_GWF_buffer(logψ::NaiveFunction,conf) = NotImplementedBuffer()