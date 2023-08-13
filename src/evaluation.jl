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







@inline function evaluate(elm::T,coeff::A,x) where {T<:Union{LineElement,TriangleElement,DiskElement},A<:AbstractVector}
    
    state=recursive_evaluate(elm, first(coeff),x)
   

    while state[2]+1 <=length(coeff)
        index=state[2] +1
        state =recursive_evaluate(elm, coeff[index],x;state=state)
    end

    return first(state)
    
end


@inline function evaluate(line::T,coeff::A,x::Tuple{S}) where {T<:Union{LineElement,TriangleElement,DiskElement},A<:AbstractVector,S}
    evaluate(line,coeff,first(x))
end


@inline function evaluate(elem::A,coeff,x) where {A}
    
    #we pick the first that is the outer loop 
    x_first =first(x)
    x_rest=Base.tail(x)
    
    nfirst=first(size(coeff))
    n_rest =Base.tail(size(coeff))
   
    #we create a tuple of ranges to properly handly the view  
    restranges=ntuple(i-> 1:n_rest[i],length(n_rest))
    
    #we do the first step 
    c=evaluate(elem.last,view(coeff,1,restranges...),x_rest)
    state=recursive_evaluate(elem.first,c,x_first)
    index=1
    
    #we do the the other step 
    while state[2]+1 <=nfirst
        index=index +1
        c=evaluate(elem.last,view(coeff,index,restranges...),x_rest)
        state =recursive_evaluate(elem.first, c,x_first;state=state)
    end

    return first(state)


end 