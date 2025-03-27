
struct CorrelationFunction{T1,PlanType<:FFTW.FFTWPlan} <: AbstractObservable
    x::T1
    Sq::T1
    plan::PlanType
end
"""
Given the dimension `dims`, allocates a functor that computes the FFT for a given spin configuration of size dims.
"""
function CorrelationFunction(dims)
    Sq = zeros(ComplexF32,dims)
    x = zeros(ComplexF32,dims)

    plan = FFTW.plan_fft(x)

    return CorrelationFunction(x,Sq,plan)
end

obs(Corr::CorrelationFunction) = Corr.Sq
Base.copy(Corr::CorrelationFunction) = CorrelationFunction(copy(Corr.x),copy(Corr.Sq),Corr.plan)

function (Corr::CorrelationFunction{T,PlanType})(out::T,Conf::PlanType) where {T,PlanType}
    FFTW.mul!(out, Corr.plan, Conf)
    Nsites = length(Conf)
    Nsites⁻¹ = 1/Nsites
    prefac = Nsites⁻¹
    @inbounds for i in eachindex(out)
        out[i] = abs2(out[i])*prefac
    end
    out
end

function compute_corr!(out,Conf::T,Corr::CorrelationFunction{T}) where {T}
    FFTW.mul!(out, Corr.plan, Conf)
    Nsites = length(Conf)
    Nsites⁻¹ = 1/Nsites
    prefac = Nsites⁻¹
    @inbounds for i in eachindex(out)
        out[i] = abs2(out[i])*prefac
    end
    out
end

function (Corr::CorrelationFunction)(out,Conf)
    copyto!(Corr.x,Conf)
    compute_corr!(out,Corr.x,Corr)
end

@views function getCorrelationFunction(res,m_values)
    pMax = maximum(m_values)
    Gnp = precomputeNormalizedAccWeight(res.TotalWeights,2pMax)

    Conf = res.SaveConfigs[:,begin,begin]

    SqFunc = CorrelationFunction(size(Conf))
    SaveConfs = res.SaveConfigs
    reconfTable = res.reconfigurationTable
    res_m = getObs_diagonal(Gnp,SaveConfs,reconfTable,SqFunc,m_values)
    newRes_m = [similar(real(res),size(res)) for res in res_m]

    return newRes_m
end