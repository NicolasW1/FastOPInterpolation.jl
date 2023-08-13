abstract type AbstractRecurrenceBuffer end

"""
    initalizeRecurrence!(buffer, x)

Triggers a reset of the buffer and sets the spatial variable of the buffer to `x`.
"""
@inline function initalizeRecurrence!(buffer::AbstractRecurrenceBuffer, x::AbstractVector)::Nothing
    resetBuffer!(buffer)
    buffer.x .= x

    return nothing
end

"""
    jacobiRecurrenceRelation(Pₙ₋₁, Pₙ₋₂, n, α, β, x)

Classical three terms recurrence relation for Jacobi-Polynomials defined on ``(-1,1)``.
Requires the two pervious terms `Pₙ₋₁`, `Pₙ₋₂`, the order `n`, the parameters `α>-1`, `β>-1` and the location `x` as input.
The value of the previous terms is irrelevant for `n=0` and `n=1`.
"""
@fastmath @inline function jacobiRecurrenceRelation(Pₙ₋₁::A, Pₙ₋₂::B, n::Integer, α::C, β::D, x::E) where {A,B,C,D,E} 
    type=promote_type(A,B,C,D,E)
    if n==0
        one(type)
    elseif n==1
        (α - β + x * (2 + α + β))/2*one(type)
    else
        one(type)/(2*n*(n + α + β)*(2*n + α + β - 2)) * ((2*n + α + β - 1) * ((2*n + α + β) * (2*n + α + β - 2) * x + α^2 - β^2) * Pₙ₋₁ - (2*(n + α - 1)*(n + β - 1)*(2*n + α + β)) * Pₙ₋₂)
    end
end