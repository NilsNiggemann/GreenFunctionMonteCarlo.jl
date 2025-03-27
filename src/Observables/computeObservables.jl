function getObs_diagonal(Gnp,AllConfigs::AbstractArray{T,3},reconfigurationTable,ObsFunc::AbstractObservable,m_values) where T
    N = lastindex(AllConfigs,ndims(AllConfigs))

    pMax = maximum(m_values)

    Obs = obs(ObsFunc)
    num_m = [zeros(size(Obs)) for _ in m_values]
    denom = 0.

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    WalkerMultiplicities = zeros(Int,Nw)
    ObsBuffer = [similar(Obs) for α in 1:Nw, n in 1:(pMax)]
    
    wrap_idx(n) = (n-1) % (pMax) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    _fill_obs_buffer!(ObsBuffer,1:pMax,ObsFunc,AllConfigs,pMax)

    for n in pMax+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        _fill_obs_buffer!(ObsBuffer,n-1,ObsFunc,AllConfigs,pMax)
        for (i_m,m) in enumerate(m_values)
            WalkerMultiplicities .= 0
            for α in 1:Nw
                α´ = α
                for i_m in 1:m
                    α´ = reconfigurationTable[α´,n-i_m]
                end
                WalkerMultiplicities[α´] += 1
            end


            for α in 1:Nw
                mult = WalkerMultiplicities[α]
                mult == 0 && continue

                O = obsBuffer(α,n-m)
                @. num_m[i_m] += O*Gn*mult
            end
        end
    end
    for i in eachindex(num_m)
        num_m[i] ./= denom
    end
    return num_m
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
