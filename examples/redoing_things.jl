parentDirectory(path, n=1) = joinpath(splitpath(normpath(path))[1:end-n])

push!(LOAD_PATH, joinpath(parentDirectory(@__DIR__,2)));
using FastOPInterpolation

# the syntax is Element(order, geometry, JacobiParameter)
line1 = LineElement(51, LineGeometry(-1.5,3.5), LineJacobiParameter(1.2,3.5))
triangle1 = TriangleElement(18, TriangleGeometry((0.0,0.0),(1.0,0.0), (0.0,1.0)), TriangleJacobiParameter(0.,0.,0.))
disk1 = DiskElement(14, DiskGeometry((0.0,0.0), 1.0), DiskRecurrenceParameter(0.))

# domain = line ⊗ triangle
intP1, t_nodes1 = InterpolationElement((line1, triangle1));
# use a disk instead of a triangle as domain
intP2, t_nodes2 = InterpolationElement((line1, disk1));

# the funciton we wish to interpolate as a test
f1(x) = sin(x[2] *  x[1] + x[3])

# auto generate nodes, please note that you can also provide your own nodes
f_nodes1 = generateTensorNodes(t_nodes1);
f_nodes2 = generateTensorNodes(t_nodes2);

# generate values on nodes
t_vals1 = [f1(x) for x in eachcol(f_nodes1)];
t_vals2 = [f1(x) for x in eachcol(f_nodes2)];

# update interpolations
updateInterpolation!(intP1, t_vals1)
updateInterpolation!(intP2, t_vals2)

f1([0.23, 0.15, 0.31])
intP1([0.23, 0.15, 0.31])
intP2([0.23, 0.15, 0.31])

# in both cases the error is at the level of machine precision, due to the high order interpolation
f1([0.23, 0.15, 0.31]) - intP1([0.23, 0.15, 0.31])
f1([0.23, 0.15, 0.31]) - intP2([0.23, 0.15, 0.31])



#### Test new things 

line1 = LineElement(51, LineGeometry(-1.,1.), LineJacobiParameter(1.2,3.5))
intP1, t_nodes1 = InterpolationElement((line1,));
intP1

t_nodes1

f_nodes1 = generateTensorNodes(t_nodes1)

t_vals1 = [f1(x) for x in eachcol(f_nodes1)];


line1

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






@inline function evaluate(line::T,coeff::A,x) where {T<:Union{LineElement,TriangleElement},A<:AbstractVector}
    
    state=recursive_evaluate(line, first(coeff),x)
    index=1

    while state[2]+1 <=length(coeff)
        index=index +1
        state =recursive_evaluate(line, coeff[index],x;state=state)
    end

    return first(state)
    
end


@inline function evaluate(line::T,coeff::A,x::Tuple{S}) where {T<:Union{LineElement,TriangleElement,DiskElement},A<:AbstractVector,S}
    evaluate(line,coeff,first(x))
end

coeff=[1. for i in 1:10 ]


evaluate(line1,coeff,(1,))

using BenchmarkTools



@benchmark evaluate($line1,$coeff,1)



abstract type Element end

abstract type Concrete_Element <: Element end 

abstract type Non_Concrete_Element <: Element end 

struct Composite_Element{A,B}  <:Non_Concrete_Element
    first::A
    last::B 
end

Composite_Element(a,rest...)= Composite_Element(a,Composite_Element(rest...))



@inline function evaluate(elem::A,coeff,x) where {A}
    
    #we pick the first that is the outer loop 
    x_first =first(x)
    x_rest=Base.tail(x)
    
    nfirst=first(size(coeff))
    n_rest =Base.tail(size(coeff))
   
    
    restranges=ntuple(i-> 1:n_rest[i],length(n_rest))
    
    c=evaluate(elem.last,view(coeff,1,restranges...),x_rest)
    state=recursive_evaluate(elem.first,c,x_first)
    index=1
    while state[2]+1 <=nfirst
        index=index +1
        c=evaluate(elem.last,view(coeff,index,restranges...),x_rest)
        state =recursive_evaluate(elem.first, c,x_first;state=state)
    end

    return first(state)


end 

composite= Composite_Element(line1,line1)

using Tullio

@tullio coeff2[a,b]:= coeff[a]*coeff[b]


evaluate(composite,coeff2,(1,1,))




@code_warntype evaluate(composite,coeff2,(1,1))


@benchmark evaluate($composite,$coeff2,(1,1))



@tullio coeff3[a,b,c]:= coeff[a]*coeff[b]*coeff[c]

composite= Composite_Element(line1,line1,line1)

evaluate(composite,coeff3,(1,1,1))

@benchmark evaluate($composite,$coeff3,(1,1,1))

###### now the triangle

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





evaluate(triangle1,coeff,(0.5,0.5))
@code_warntype evaluate(triangle1,coeff,(0.5,0.5))

@benchmark evaluate($triangle1,$coeff,(0.5,0.5))

prova=((1,),(2,3))



composite=Composite_Element(line1,triangle1)

evaluate(composite,coeff2,(0.5,(0.5,0.5)))



@code_warntype evaluate(composite,coeff2,(0.5,(0.5,0.5)))

@benchmark evaluate($composite,$coeff2,(0.5,(0.5,0.5)))



composite=Composite_Element(triangle1,line1,triangle1)

evaluate(composite,coeff3,((0.5,0.5),0.5,(0.5,0.5)))



@code_warntype evaluate(composite,coeff3,((0.5,0.5),0.5,(0.5,0.5)))

@benchmark evaluate($composite,$coeff3,((0.5,0.5),0.5,(0.5,0.5)))






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