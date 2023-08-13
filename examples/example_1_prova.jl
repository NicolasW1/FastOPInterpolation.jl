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


@inline function three_term_recurrence(Pᵢ₋₁, Pᵢ₋₂, order, line::LineElement{T,N},x) where {T,N}

    jacobiRecurrenceRelation(Pᵢ₋₁, Pᵢ₋₂, order, line.params.α, line.params.β, x)
end 




@inline function evaluate(line::LineElement{T,N},coeff::A,x::S) where {S<:Number,T,N,A<:AbstractVector}
    type_promoted=promote_type(S,T,eltype(A))
    #Pᵢ=zero(x)
    Pᵢ₋₁=zero(type_promoted)
    Pᵢ₋₂=zero(type_promoted)
    result=zero(type_promoted)
    for (i,c) in enumerate(coeff)
        Pᵢ = three_term_recurrence(Pᵢ₋₁, Pᵢ₋₂, i-1, line, x)
        result = Pᵢ*c
        Pᵢ₋₂ = Pᵢ₋₁
        Pᵢ₋₁ = Pᵢ
    end 
    return result 
    
end


function recursive_evaluate(line::LineElement{T,N},coeff::A,x::S;state= ( zero(promote_type(S,T,eltype(A))),0,zero(promote_type(S,T,eltype(A))),zero(promote_type(S,T,eltype(A)))  ) )    where {S<:Number,T,N,A<:AbstractVector}
    
    c= coeff[order+1]
    (result,order, Pᵢ₋₁ ,Pᵢ₋₂) =state 
    Pᵢ = three_term_recurrence(Pᵢ₋₁, Pᵢ₋₂, order, line, x)
    result = muladd( Pᵢ,c,result)
    Pᵢ₋₂ = Pᵢ₋₁
    Pᵢ₋₁ = Pᵢ

    return ( result , order ,Pᵢ₋₁,Pᵢ₋₂ )
end 




@inline function evaluate(line::LineElement{T,N},coeff::A,x::Tuple{S}) where {S<:Number,T,N,A<:AbstractVector}
    evaluate(line,coeff,first(x)) 
    
end




coeff=[1. for i in 1:10 ]


evaluate(1.,line1,coeff)

using BenchmarkTools

@code_warntype evaluate(1.,line1,coeff)

@benchmark evaluate(1,$line1,$coeff)



struct Composite_Element{A,B}  
    first::A
    last::B 
end

Composite_Element(a,rest...)= Composite_Element(a,Composite_Element(rest...))



@inline function evaluate(elem::A,coeff,x) where {A<:Composite_Element}
    
    #we pick the first that is the outer loop 
    x_first,x_rest... =x 
    nfirst,rest... =size(coeff)
    restranges=ntuple(i-> 1:rest[i],length(rest))

    
    Pᵢ₋₁=zero(eltype(coeff))
    Pᵢ₋₂=zero(eltype(coeff))
    result=zero(eltype(coeff))
    
    for i in Base.OneTo(nfirst)
        
        c=evaluate(x_rest,elem.last,view(coeff,i,restranges...))
        Pᵢ = three_term_recurrence(Pᵢ₋₁, Pᵢ₋₂, i-1, elem.first, x_first)
        result = Pᵢ*c
        Pᵢ₋₂ = Pᵢ₋₁
        Pᵢ₋₁ = Pᵢ
    end 
    return result 


end 

composite= Composite_Element(line1,line1)

using Tullio

@tullio coeff2[a,b]:= coeff[a]*coeff[b]


evaluate(composite,coeff2,(1,1,))




@code_warntype evaluate(composite,coeff2,(1,1,))

@benchmark evaluate($composite,$coeff2,(1,1,))



@tullio coeff3[a,b,c]:= coeff[a]*coeff[b]*coeff[c]

composite= Composite_Element(line1,line1,line1)

evaluate(composite,coeff3,(1,1,))

@benchmark evaluate($composite,$coeff3,(1,1,1))