function getObs_diagonal(Gnps::AbstractMatrix,AllConfigs::AbstractArray{T,3},reconfigurationTable,ObsFunc::AbstractObservable,m_values) where T
    N = lastindex(AllConfigs,ndims(AllConfigs))

    Nw = lastindex(AllConfigs,ndims(AllConfigs)-1)
    Nw⁻¹ = 1/Nw
    m_max = maximum(m_values)

    PopulationMatrix = similar(reconfigurationTable)
    m_values = 0:m_max-1
    Obs = obs(ObsFunc)
    
    num_m = [zeros(size(Obs)) for _ in m_values]
    denom_m = zeros(length(m_values))

    ObsBuffer = [similar(Obs) for α in 1:Nw, n in 1:(m_max)]
        
    wrap_idx(n) = (n-1) % (m_max) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    _fill_obs_buffer!(ObsBuffer,1:m_max,ObsFunc,AllConfigs,m_max)

    for n in 1:N
        getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max-1)
        _fill_obs_buffer!(ObsBuffer,n,ObsFunc,AllConfigs,m_max)
        for m_index in eachindex(m_values)
            m = m_values[m_index]
            Gnp = Gnps[n,1+2m]
            Gnp == 0 && continue

            denom_m[m_index] += Gnp

            num_m_i = num_m[m_index]
            for α in 1:Nw
                mult = PopulationMatrix[α,m_index]
                mult == 0 && continue
                mult *= Nw⁻¹
                O = obsBuffer(α,n-m)
                Gnpmult = Gnp*mult
                LoopVectorization.@turbo num_m_i .+= O.*Gnpmult
            end
        end
    end
    for i in eachindex(num_m)
        num_m[i] ./= denom_m[i]
    end
    return num_m
    
end
function getObs_diagonal(Observer,ObsFunc::AbstractObservable,m_values)
    Obs = getObs(Observer)
    Gnps = precomputeNormalizedAccWeight(Obs.TotalWeights,2*maximum(m_values))
    SaveConfigs = Obs.SaveConfigs
    reconfigurationTable = Obs.reconfigurationTable
    return getObs_diagonal(Gnps,SaveConfigs,reconfigurationTable,ObsFunc,m_values)
end

function getPopulationMatrix!(PopulationMatrix,reconfigurationTable::AbstractMatrix,n,projectionLength)
    nMax = size(reconfigurationTable,2)
    PopulationMatrixParent = parent(PopulationMatrix)
    # fill!(PopulationMatrixParent,0)
    for i_m in 0:min(n-1,projectionLength)
        PopulationMatrixParent[:,i_m+1] .= 0
        for α in axes(reconfigurationTable,1)
            if i_m == 0
                pop = 1
            else
                pop = PopulationMatrixParent[α,i_m]
            end
            pop == 0 && continue
            α´ = reconfigurationTable[α,n-i_m]
            α´ == 0 && continue
            PopulationMatrixParent[α´,i_m+1] += pop
        end
    end
    return PopulationMatrix
end

function _fill_obs_buffer!(ObsBuffer,nRange,ObsFunc!,AllConfigs::AbstractArray{T,3},pMax) where T
    wrap_idx(n) = (n-1) % (pMax) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    for n in nRange, α in axes(AllConfigs,2)
        conf = @view AllConfigs[:,α,n]
        ObsFunc!(obsBuffer(α,n),conf)
    end
    return
end
