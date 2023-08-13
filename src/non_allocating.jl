
abstract type Element end

abstract type Concrete_Element <: Element end 

abstract type Non_Concrete_Element <: Element end 

struct Composite_Element{A,B}  <:Non_Concrete_Element
    first::A
    last::B 
end

Composite_Element(a,rest...)= Composite_Element(a,Composite_Element(rest...))


@fastmath @inline function jacobiRecurrenceRelation(Pₙ₋₁::A, Pₙ₋₂::B, n::Integer, α::C, β::D, x::E) where {A,B,C,D,E} 
    type=promote_type(A,B,C,D,E)
    if n==0
        one(type)
    elseif n==1
        (α - β + x * (2 + α + β))/2*one(type)
    else
        1/(2*n*(n + α + β)*(2*n + α + β - 2)) * ((2*n + α + β - 1) * ((2*n + α + β) * (2*n + α + β - 2) * x + α^2 - β^2) * Pₙ₋₁ - (2*(n + α - 1)*(n + β - 1)*(2*n + α + β)) * Pₙ₋₂)*one(type)
    end
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
    Pᵢ₋₂ = Pᵢ₋₁
    Pᵢ₋₁ = Pᵢ

    return ( result , order+1 ,Pᵢ₋₁,Pᵢ₋₂ )
end 

@inline function recursive_evaluate(triangle::TriangleElement{T,N},
    coeff_i::A,
    x::Tuple{S,S}
    ;
    state= 
    ( zero(promote_type(S,T,A)) ,
    0, N + 1, -1 , 
    zero(promote_type(S,T,A)) ,
    zero( promote_type(S,T,A) ),
    zero(promote_type(S,T,A)) ,
    zero( promote_type(S,T,A) ) 
    ) 
    )    where {S<:Number,T,N,A<:Number}
    
    x₀, y₀ = x
    xtilde = one(S) - 2 * x₀
    ytilde = 1-x₀ ≈ zero(S) ? zero(S) : 2 * (y₀/(1-x₀)) - one(S)
    
    result=state[1]
    order=state[2]
    
    (result,order,i,j, Pᵢ₋₁ ,Pᵢ₋₂,Pⱼ₋₁,Pⱼ₋₂) =state 

    if i + j == N
        j += 1
        i = zero(N)

        Pᵢ₋₁ = zero(T)
        Pᵢ₋₂ = zero(T)

        Pᵢ = jacobiRecurrenceRelation(Pᵢ₋₁, Pᵢ₋₂, i, triangle.params.a, 2*j+one(T) + triangle.params.b + triangle.params.c, xtilde)
        Pᵢ₋₂ = Pᵢ₋₁
        Pᵢ₋₁ = Pᵢ

        Pⱼ = jacobiRecurrenceRelation(Pⱼ₋₁, Pⱼ₋₂, j, triangle.params.b, triangle.params.c, ytilde)
        Pⱼ₋₂ = Pⱼ₋₁
        Pⱼ₋₁ = Pⱼ

        totelm=Pᵢ * Pⱼ * (1-x₀)^j
        result = muladd( totelm,coeff_i,result)

        return (result,order+1,i,j, Pᵢ₋₁ ,Pᵢ₋₂,Pⱼ₋₁,Pⱼ₋₂)

    else
        i += 1

        Pᵢ = jacobiRecurrenceRelation(Pᵢ₋₁, Pᵢ₋₂, i, triangle.params.a, 2*j+one(T) + triangle.params.b + triangle.params.c, xtilde)
        Pᵢ₋₂ = Pᵢ₋₁
        Pᵢ₋₁ = Pᵢ
        
        Pⱼ= Pⱼ₋₁

        totelm=Pᵢ * Pⱼ * (1-x₀)^j
        result = muladd( totelm,coeff_i,result)

        return (result,order+1,i,j, Pᵢ₋₁ ,Pᵢ₋₂,Pⱼ₋₁,Pⱼ₋₂)


    end


end 


@inline function recursive_evaluate(disk::DiskElement{T,N},
    coeff_i::A,
    x::Tuple{S,S}
    ;
    state= 
    ( zero(promote_type(S,T,A)) ,
    0, N + 1, -1 , 
    zero(promote_type(S,T,A)) ,
    zero( promote_type(S,T,A) ),
    zero(promote_type(S,T,A)) ,
    zero( promote_type(S,T,A) ) 
    ) 
    )    where {S<:Number,T,N,A<:Number}
    
    x₀, y₀ = x
    ytilde = 1-abs(x₀) ≈ zero(T) ? zero(T) : (1-x₀^2)^(-one(T)/2) * y₀
    
    result=state[1]
    order=state[2]
    
    (result,order,i,j, Pᵢ₋₁ ,Pᵢ₋₂,Pⱼ₋₁,Pⱼ₋₂) =state 

    if i + j == N
        j += 1
        i = zero(N)

        Pᵢ₋₁ = zero(T)
        Pᵢ₋₂ = zero(T)

        Pᵢ = jacobiRecurrenceRelation(Pᵢ₋₁, Pᵢ₋₂, i, a + j + one(T)/2, disk.params.a + disk.buffer.j + one(T)/2, x₀)
        Pᵢ₋₂ = Pᵢ₋₁
        Pᵢ₋₁ = Pᵢ

        Pⱼ = jacobiRecurrenceRelation(Pⱼ₋₁, Pⱼ₋₂, j, disk.params.a, disk.params.a, ytilde)
        Pⱼ₋₂ = Pⱼ₋₁
        Pⱼ₋₁ = Pⱼ

        totelm=Pᵢ * Pⱼ * (1-x₀^2)^(j / 2)
        result = muladd( totelm,coeff_i,result)

        return (result,order+1,i,j, Pᵢ₋₁ ,Pᵢ₋₂,Pⱼ₋₁,Pⱼ₋₂)

    else
        i += 1

        Pᵢ = jacobiRecurrenceRelation(Pᵢ₋₁, Pᵢ₋₂, i, a + j + one(T)/2, disk.params.a + disk.buffer.j + one(T)/2, x₀)
        Pᵢ₋₂ = Pᵢ₋₁
        Pᵢ₋₁ = Pᵢ
        
        Pⱼ= Pⱼ₋₁

        totelm=Pᵢ * Pⱼ * (1-x₀^2)^(j / 2)
        result = muladd( totelm,coeff_i,result)

        return (result,order+1,i,j, Pᵢ₋₁ ,Pᵢ₋₂,Pⱼ₋₁,Pⱼ₋₂)


    end


end



@inline function evaluate(elm::T,coeff::A,x) where {T<:Union{LineElement,TriangleElement,DiskElement},A<:AbstractVector}
    
    state=recursive_evaluate(elm, first(coeff),x)
    index=1

    while state[2]+1 <=length(coeff)
        index=index +1
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