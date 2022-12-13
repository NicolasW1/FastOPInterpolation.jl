abstract type AbstractBasisGeometry{T,N} end

eltype(::Type{<:AbstractBasisGeometry{T, N}}) where {T, N} = T

order(::Type{<:AbstractBasisGeometry{T, N}}) where {T, N} = N
order(x::AbstractBasisGeometry) = order(typeof(x))

dim(x::AbstractBasisGeometry) = dim(typeof(x))
dof(x::AbstractBasisGeometry) = dof(typeof(x))

RecurrenceBuffer(x::AbstractBasisGeometry) = RecurrenceBuffer(typeof(x))