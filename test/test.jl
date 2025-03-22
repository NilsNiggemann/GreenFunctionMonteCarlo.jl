struct B 
    x::Ref{Union{Nothing,Int}}
    y::Int
end

struct B2{T}
    x::Ref{Union{Nothing,Int}}
    y::T
end

b = B(nothing,1)# works, returns B(Base.RefValue{Union{Nothing, Int64}}(nothing), 1)

b2 = B2(nothing,1)# ERROR: MethodError: no method matching B2(::Nothing, ::Int64)
b2 = B2(convert(Ref{Union{Nothing,Int}},nothing),1) #works, returns B2{Int64}(Base.RefValue{Union{Nothing, Int64}}(nothing), 1)
b2 = B2{Int}(nothing, 1) #works, returns B2{Int64}(Base.RefValue{Union{Nothing, Int64}}(nothing), 1)
