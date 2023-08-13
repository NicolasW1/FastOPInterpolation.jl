abstract type AbstractBasisGeometry{T,N} end

eltype(::Type{<:AbstractBasisGeometry{T, N}}) where {T, N} = T

order(::Type{<:AbstractBasisGeometry{T, N}}) where {T, N} = N
order(x::AbstractBasisGeometry) = order(typeof(x))

dim(x::AbstractBasisGeometry) = dim(typeof(x))
dof(x::AbstractBasisGeometry) = dof(typeof(x))

RecurrenceBuffer(x::AbstractBasisGeometry) = RecurrenceBuffer(typeof(x))


abstract type Abstract_Element end

abstract type Concrete_Element <: Abstract_Element end 

abstract type Non_Concrete_Element <: Abstract_Element end 

struct Composite_Element{A <:Concrete_Element,B<:Non_Concrete_Element}  <:Non_Concrete_Element
    first::A 
    last::B 
end

Composite_Element(a,rest...)= Composite_Element(a,Composite_Element(rest...))