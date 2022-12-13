@inline function initBuffer!(buffer, element, range, x)
    resetBuffer!(buffer)
    @inbounds mapToReference!(buffer.x, (@view x[range]), element.geometry)
end

# this circumvents directly looping over the tuple, which is not type stable as the tuple gets accessed with a runtime index
@inline function innerEvaluationLoop(tupleElement::Tuple, tupleBuffer::Tuple, i::Int, intElement::InterpolationElement{T,D,X,Y}, ind::CartesianIndex) where {T,D,X,Y}
    if length(tupleElement) == 0
        return nothing
    end

    @inbounds element = tupleElement[begin]
    @inbounds buffer = tupleBuffer[begin]

    @inbounds if ind[i]==one(i)
        resetBuffer!(buffer)
        @inbounds intElement.evalBuffer[i] = getBasisElement!(buffer, element.params)
        @inbounds innerEvaluationLoop(tupleElement[2:end], tupleBuffer[2:end], i+1, intElement, ind)
    else
        @inbounds intElement.evalBuffer[i] = getBasisElement!(buffer, element.params)
        return nothing
    end
end

@fastmath function (intElement::InterpolationElement{T,D,X,Y})(x) where {T,D,X,Y}
    result = zero(T)

    broadcast(((a,b,c)->initBuffer!(a,b,c,x)), intElement.recurrenceBuffers, intElement.elements, intElement.arg_ranges)

    for ind in CartesianIndices(intElement.dofs)
        innerEvaluationLoop(intElement.elements, intElement.recurrenceBuffers, 1, intElement, ind)
        @inbounds result += intElement.modes[ind] * prod(intElement.evalBuffer)
    end

    return result
end