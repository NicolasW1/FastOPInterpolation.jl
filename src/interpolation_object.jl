mutable struct InterpolationElement{T, D, X, Y}
    modes::Array{T, D}
    invV::NTuple{D, Matrix{T}}

    elements::X

    dims::NTuple{D, Int}
    orders::NTuple{D, Int}
    dofs::NTuple{D, Int}
    arg_ranges::NTuple{D, UnitRange{Int}}

    recurrenceBuffers::Y
    evalBuffer::MVector{D, T}
end

dim(element::InterpolationElement) = sum(element.dims)
dof(element::InterpolationElement) = prod(element.dofs)

function nodes(intElement::InterpolationElement)
    Tuple(nodes(x) for x in intElement.elements)
end

function singleVandermonde(buffer::AbstractRecurrenceBuffer, element::AbstractBasisGeometry{T,N}, nodes::AbstractMatrix{T}) where {T,N}
    V = Matrix{T}(undef, dof(element), dof(element))

    for i in 1:dof(element)
        resetBuffer!(buffer)
        mapToReference!(buffer.x, nodes[:,i], element.geometry)
        for j in 1:dof(element)
            V[i,j] = getBasisElement!(buffer, element.params)
        end
    end

    inv(V)
end

"""
    setVandermonde(intElement::InterpolationElement{T,D,X,Y}, nodes::NTuple{D, AbstractMatrix{T}})

Update the Vandermonde matrix of an interpolation element for set of nodes.
"""
function setVandermonde(intElement::InterpolationElement{T,D,X,Y}, nodes::NTuple{D, AbstractMatrix{T}})::Nothing where {T,D,X,Y}
    intElement.invV = ntuple(i -> singleVandermonde(intElement.recurrenceBuffers[i], intElement.elements[i], nodes[i]) , D)

    return nothing
end

"""
    InterpolationElement(elements::X, constructNodes::Bool = true)

Initalize an interpolation element from a collection of tuple of simple elements.
The optional second argument `constructNodes` determines whether or not predefined nodes should be used.
"""
function InterpolationElement(elements::X, constructNodes::Bool = true) where {X}
    T = eltype(typeof(elements[1]))
    D = length(elements)

    dims = map(a->dim(a), elements)
    orders = map(a->order(a), elements)
    dofs = map(a->dof(a), elements)

    arg_ranges = generateRanges(dims)

    empty_invV = Tuple([zeros(T, i, i) for i in dofs])
    recurrenceBuffers = Tuple([RecurrenceBuffer(x) for x in elements])
    evalBuffer = zeros(fieldtype(InterpolationElement{T, D, X}, :evalBuffer))

    interpObj = InterpolationElement(zeros(T, dofs), empty_invV, elements, dims, orders, dofs, arg_ranges, recurrenceBuffers, evalBuffer)

    if constructNodes
        interpNodes = nodes(interpObj)
        setVandermonde(interpObj, interpNodes)

        return interpObj, interpNodes
    end

    return interpObj
end

