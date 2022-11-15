include("nodes/disk_nodes.jl")

export DiskRecurrenceParameter, DiskGeometry, DiskElement
export mapToReference!, nodes

# Recurrence Parameter
struct DiskRecurrenceParameter{T}
	a::T
end


# Geometry
struct DiskGeometry{T}
	center::Tuple{T,T}
	radius::T
end

@fastmath @inbounds @inline function mapToReference!(x₀::AbstractVector{T}, x::AbstractVector{T}, disk::DiskGeometry{T})::Nothing where {T}
	x₀[1] = (x[1] - disk.center[1]) / disk.radius
	x₀[2] = (x[2] - disk.center[2]) / disk.radius

	return nothing
end
@fastmath @inbounds @inline function mapFromReference!(x::AbstractVector{T}, x₀::AbstractVector{T}, disk::DiskGeometry{T})::Nothing where {T}
	x[1] = (disk.radius * x₀[1]) + disk.center[1]
	x[2] = (disk.radius * x₀[2]) + disk.center[2]

	return nothing
end


# Base Struct
struct DiskElement{T,N} <: AbstractBasisGeometry{T,N}
	geometry::DiskGeometry{T}
	params::DiskRecurrenceParameter{T}
end


function DiskElement(n::Integer,  geometry::DiskGeometry{T}, params::DiskRecurrenceParameter{T}) where {T}
	DiskElement{T, n}(geometry, params)
end

dim(::Type{DiskElement{T, N}}) where {T, N} = 2
dof(::Type{DiskElement{T, N}}) where {T, N} = ((N+1) * (N+2)) ÷ 2


@inbounds function nodes(disk::DiskElement{T,N}) where {T,N}
	n = order(disk)
	ref_nodes = getRefDiskNodes(T, n)
	element_nodes = Matrix{T}(undef, dim(disk), dof(disk))

	for i=1:dof(disk)
		mapFromReference!((@view element_nodes[:,i]), (@view ref_nodes[i,:]), disk.geometry)
	end

	element_nodes
end

# Recurrence
mutable struct DiskRecurrenceBuffer{T, N} <: AbstractRecurrenceBuffer
	x::MVector{2,T}

	i::Int
	j::Int

	Pᵢ::T
	Pᵢ₋₁::T
	Pᵢ₋₂::T
	Pⱼ::T
	Pⱼ₋₁::T
	Pⱼ₋₂::T
end

function RecurrenceBuffer(::Type{DiskElement{T,N}}) where {T,N}
	type_int = typeof(N)
	DiskRecurrenceBuffer{T,N}(zeros(MVector{2,T}), zero(type_int), zero(type_int), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
end

@inline function resetBuffer!(buffer::DiskRecurrenceBuffer{T,N})::Nothing where {T,N}
	buffer.i = N + 1
	buffer.j = -one(N)

	buffer.Pᵢ = zero(T)
	buffer.Pᵢ₋₁ = zero(T)
	buffer.Pᵢ₋₂ = zero(T)
	buffer.Pⱼ = zero(T)
	buffer.Pⱼ₋₁ = zero(T)
	buffer.Pⱼ₋₂ = zero(T)

	return nothing
end

@inline @inbounds @fastmath function getBasisElement!(buffer::DiskRecurrenceBuffer{T, N}, params::DiskRecurrenceParameter{T}) where {T,N}
	#j=k
	#i=n-k
	x₀::T, y₀::T = buffer.x[1], buffer.x[2]
	ytilde::T = 1-x₀ ≈ zero(T) ? zero(T) : (1-x₀^2)^(-one(T)/2) * y₀

	@inline function update_i()
		buffer.Pᵢ = jacobiRecurrenceRelation(buffer.Pᵢ₋₁, buffer.Pᵢ₋₂, buffer.i, params.a + buffer.j + one(T)/2, params.a + buffer.j + one(T)/2, x₀)
		buffer.Pᵢ₋₂ = buffer.Pᵢ₋₁
		buffer.Pᵢ₋₁ = buffer.Pᵢ
	end

	@inline function update_j()
		buffer.Pⱼ = jacobiRecurrenceRelation(buffer.Pⱼ₋₁, buffer.Pⱼ₋₂, buffer.j, params.a, params.a, ytilde)
		buffer.Pⱼ₋₂ = buffer.Pⱼ₋₁
		buffer.Pⱼ₋₁ = buffer.Pⱼ
	end

	if buffer.i + buffer.j == N
		buffer.j += 1
		buffer.i = zero(N)

		buffer.Pᵢ₋₁ = zero(T)
		buffer.Pᵢ₋₂ = zero(T)

		update_i()
		update_j()

	else
		buffer.i += 1

		update_i()
	end


	return buffer.Pᵢ * buffer.Pⱼ * (1-x₀^2)^(buffer.j / 2)
end