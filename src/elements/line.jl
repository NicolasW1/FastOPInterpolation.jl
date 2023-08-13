# Recurrence Parameter
struct LineJacobiParameter{T}
    α::T
    β::T
end

# Geometry
struct LineGeometry{T}
    vertex_1::Tuple{T}
    vertex_2::Tuple{T}
end

function LineGeometry(lower_bound::T, upper_bound::T) where {T}
    LineGeometry((lower_bound,), (upper_bound,))
end


# """
#   mapToReference!(x₀, x, line)
#   Maps the point `x` (as vector) form the geometric object `line` to the point `x₀` (as vector) on the reference element `[-1,1]`
# """
@inline function mapToReference!(x₀::AbstractVector{T}, x::AbstractVector{T}, line::LineGeometry{T})::Nothing where {T}
    @inbounds x₀[1] = (2 * x[1] - (line.vertex_1[1] + line.vertex_2[1])) / (line.vertex_2[1] - line.vertex_1[1])

    return nothing
end
@inline function mapFromReference!(x::AbstractVector{T}, x₀::AbstractVector{T}, line::LineGeometry{T})::Nothing where {T}
    @inbounds x[1] = line.vertex_1[1] + (x₀[1] + 1) / 2 * (line.vertex_2[1] - line.vertex_1[1])

    return nothing
end


"""
    LineElement

Basic struct holding orthogonal polynomials on the line.
"""
struct LineElement{T,N} <: AbstractBasisGeometry{T,N}
    geometry::LineGeometry{T}
    params::LineJacobiParameter{T}
end

"""
    LineElement(n, geometry, params)

Construct interpolation object on a line with order `n`, some `geometry` and parameters `params`.
"""
function LineElement(n::Integer, geometry::LineGeometry{T}, params::LineJacobiParameter{T}) where {T}
    LineElement{T, n}(geometry, params)
end

dim(::Type{LineElement{T, N}}) where {T, N} = 1
dof(::Type{LineElement{T, N}}) where {T, N} = N + 1


function nodes(line::LineElement{T,N}) where {T,N}
    n = order(line)
    node_vector = n>zero(n) ? -cos.(one(T) * π / n * collect(0:n)) : [zero(T)]
    ref_nodes = collect(transpose(node_vector))
    element_nodes = similar(ref_nodes)

    for i=1:dof(line)
        @inbounds mapFromReference!((@view element_nodes[:,i]), (@view ref_nodes[:,i]), line.geometry)
    end

    element_nodes
end

mutable struct LineRecurrenceBuffer{T, I} <: AbstractRecurrenceBuffer
    x::MVector{1,T}

    i::I

    Pᵢ::T
    Pᵢ₋₁::T
    Pᵢ₋₂::T
end

function RecurrenceBuffer(::Type{LineElement{T,N}}) where {T,N}
    LineRecurrenceBuffer(zeros(MVector{1,T}), 0, zero(T), zero(T), zero(T))
end

@inline function resetBuffer!(buffer::LineRecurrenceBuffer{T,I})::Nothing where {T,I}
    buffer.i = -one(I)

    buffer.Pᵢ = zero(T)
    buffer.Pᵢ₋₁ = zero(T)
    buffer.Pᵢ₋₂ = zero(T)

    return nothing
end

@inline function getBasisElement!(buffer::LineRecurrenceBuffer{T, I}, params::LineJacobiParameter{T}) where {T,I}
    buffer.i += 1

    @inbounds buffer.Pᵢ = jacobiRecurrenceRelation(buffer.Pᵢ₋₁, buffer.Pᵢ₋₂, buffer.i, params.α, params.β, buffer.x[1])
    buffer.Pᵢ₋₂ = buffer.Pᵢ₋₁
    buffer.Pᵢ₋₁ = buffer.Pᵢ

    return buffer.Pᵢ
end


@inline function recursive_evaluate(line::LineElement{T,N},
    coeff_i::A,
    x::S
    ;
    state= ( zero(promote_type(S,T,A)) ,
    0,
    zero(promote_type(S,T,A)) ,
    zero( promote_type(S,T,A) )  
    ) 
    )    where {S<:Number,T,N,A<:Number}
    
    
    result=state[1]
    order=state[2] 
    
    (result,order, Pᵢ₋₁ ,Pᵢ₋₂) =state 
    
    Pᵢ = jacobiRecurrenceRelation(Pᵢ₋₁, Pᵢ₋₂, order, line.params.α, line.params.β, x)
    result = muladd( Pᵢ,coeff_i,result)
 

    return ( result , order+1 ,Pᵢ,Pᵢ₋₁ )
end