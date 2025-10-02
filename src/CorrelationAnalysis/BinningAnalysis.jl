"""
    data_bunch(data::AbstractVector{T}) where T

Groups or "bunches" the input data vector into bins for statistical autocorrelation analysis. 
    
# Arguments
- `data::AbstractVector{T}`: The input data vector to be binned.

# Returns
- A binned version of the input data, typically as a vector or collection of bins, depending on the implementation.

# Example


# See Also
- [`bin_data`](@ref)
- [`autocorrelation`](@ref)
"""
function data_bunch(data::AbstractVector{T}) where T
    Base.require_one_based_indexing(data)
    N = length(data) ÷ 2
    Σ = zero(T)
    Σ′ = zero(T)
    newData = zeros(Float64, N)
    for i in 1:N
        Σ += data[2i-1] + data[2i]
        Σ′ += data[2i-1] ^2 + data[2i]^2
        newData[i] = (data[2i-1] + data[2i]) / 2
    end
    mean = Σ / (2N)
    arg = (Σ′ / (2N) - mean^2)
    if arg >0
        sig_error = sqrt(arg / 2N)
    elseif arg >= -1e-30
        sig_error = 0.0
    else
        sig_error = NaN
    end
    return mean, sig_error, newData
end

function bunching_errors_recursive(x)
    errors_bunching = Float64[]
    while length(x) > 1
        x_avg, x_err, x = data_bunch(x)
        push!(errors_bunching, x_err)
    end

    return errors_bunching

end