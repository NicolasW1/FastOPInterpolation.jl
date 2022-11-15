include("nodes/triangle_nodes.jl")

export TriangleJacobiParameter, TriangleGeometry, TriangleElement
export mapToReference!, nodes

# Recurrence Parameter
struct TriangleJacobiParameter{T}
	a::T
	b::T
	c::T
end


# Geometry
struct TriangleGeometry{T}
	vertex_1::Tuple{T,T}
	vertex_2::Tuple{T,T}
	vertex_3::Tuple{T,T}

	prefactor::T
end

function TriangleGeometry(vertex_1::Tuple{T,T}, vertex_2::Tuple{T,T}, vertex_3::Tuple{T,T}) where {T}
	prefactor = 1 / (vertex_1[2]*vertex_2[1] - vertex_1[1]*vertex_2[2] - vertex_1[2]*vertex_3[1] + vertex_2[2]*vertex_3[1] + vertex_1[1]*vertex_3[2] - vertex_2[1]*vertex_3[2])

	TriangleGeometry(vertex_1, vertex_2, vertex_3, prefactor)
end

@fastmath @inbounds @inline function mapToReference!(x₀::AbstractVector{T}, x::AbstractVector{T}, tri::TriangleGeometry{T}) where {T}
	x₀[1] = tri.prefactor * ((tri.vertex_1[2] - tri.vertex_3[2]) * x[1] + (tri.vertex_3[1] - tri.vertex_1[1]) * x[2] + tri.vertex_1[1]*tri.vertex_3[2] - tri.vertex_1[2]*tri.vertex_3[1])
	x₀[2] = tri.prefactor * ((tri.vertex_2[2] - tri.vertex_1[2]) * x[1] + (tri.vertex_1[1] - tri.vertex_2[1]) * x[2] + tri.vertex_1[2]*tri.vertex_2[1] - tri.vertex_1[1]*tri.vertex_2[2])

	return nothing
end


# Base Struct
struct TriangleElement{T,N} <: AbstractBasisGeometry{T,N}
	geometry::TriangleGeometry{T}
	params::TriangleJacobiParameter{T}
end


function TriangleElement(n::Integer,  geometry::TriangleGeometry{T}, params::TriangleJacobiParameter{T}) where {T}
	TriangleElement{T, n}(geometry, params)
end

dim(::Type{TriangleElement{T, N}}) where {T, N} = 2
dof(::Type{TriangleElement{T, N}}) where {T, N} = ((N+1) * (N+2)) ÷ 2


@inbounds function nodes(triangle::TriangleElement{T,N}) where {T,N}
	n = order(triangle)
	ref_nodes = getRefTriangleNodes_Barycentric(T, n)
	element_nodes = Matrix{T}(undef, dim(triangle), dof(triangle))

	for i=1:dof(triangle)
		for j=1:dim(triangle)
			element_nodes[j,i] = ref_nodes[i,1] * triangle.geometry.vertex_1[j] + ref_nodes[i,2] * triangle.geometry.vertex_2[j] + ref_nodes[i,3] * triangle.geometry.vertex_3[j]
		end
	end

	element_nodes
end

# Recurrence
mutable struct TriangleRecurrenceBuffer{T, N} <: AbstractRecurrenceBuffer
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

function RecurrenceBuffer(::Type{TriangleElement{T,N}}) where {T,N}
	type_int = typeof(N)
	TriangleRecurrenceBuffer{T,N}(zeros(MVector{2,T}), zero(type_int), zero(type_int), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
end

@inline function resetBuffer!(buffer::TriangleRecurrenceBuffer{T,N}) where {T,N}
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

@inline @inbounds @fastmath function getBasisElement!(buffer::TriangleRecurrenceBuffer{T, N}, params::TriangleJacobiParameter{T}) where {T,N}
	#j=k
	#i=n-k
	x₀::T, y₀::T = buffer.x[1], buffer.x[2]
	xtilde::T = one(T) - 2 * x₀
	ytilde::T = 1-x₀ ≈ zero(T) ? zero(T) : 2 * (y₀/(1-x₀)) - one(T)

	@inline function update_i()
		buffer.Pᵢ = jacobiRecurrenceRelation(buffer.Pᵢ₋₁, buffer.Pᵢ₋₂, buffer.i, params.a, 2*buffer.j+one(T) + params.b + params.c, xtilde)
		buffer.Pᵢ₋₂ = buffer.Pᵢ₋₁
		buffer.Pᵢ₋₁ = buffer.Pᵢ
	end

	@inline function update_j()
		buffer.Pⱼ = jacobiRecurrenceRelation(buffer.Pⱼ₋₁, buffer.Pⱼ₋₂, buffer.j, params.b, params.c, ytilde)
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


	return buffer.Pᵢ * buffer.Pⱼ * (1-x₀)^buffer.j
end