"""provides a naive wrapper for a guiding function which does not use any buffer. Useful for debugging and testing"""
struct NaiveFunction{T} <: AbstractGuidingFunction
    logpsi::T
end
(N::NaiveFunction)(x::Any) = N.logpsi(x)

guidingfunc_name(F::NaiveFunction) = "Naive($(F.logpsi))"
get_params(ψG::NaiveFunction) = get_params(ψG.logpsi)
allocate_GWF_buffer(logψ::NaiveFunction,conf) = NotImplementedBuffer()

struct ParametrizedFunction{T,T2} <: AbstractGuidingFunction
    logpsi::T
    params::T2
end
(N::ParametrizedFunction)(x::Any) = N.logpsi(x,N.params)

guidingfunc_name(F::ParametrizedFunction) = "Parametrized($(F.logpsi), $(F.params))"
get_params(ψG::ParametrizedFunction) = ψG.params

allocate_GWF_buffer(logψ::ParametrizedFunction,conf) = NotImplementedBuffer()
Adapt.adapt_structure(to, x::ParametrizedFunction) = ParametrizedFunction(Adapt.adapt(to, x.logpsi), Adapt.adapt(to, x.params))