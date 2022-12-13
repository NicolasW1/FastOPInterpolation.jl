"""
    updateInterpolation!(intElement::InterpolationElement, node_values::AbstractVector)

Update the modes of the interpolation from a function values.
The structure of the node values `fᵢ=f(xᵢ)` must be the same as the one used for the `xᵢ` when initalizing the Vandermonde matrix.
"""
@inline function updateInterpolation!(intElement::InterpolationElement, node_values::AbstractVector)
    mul!(vec(intElement.modes), kronecker(reverse(intElement.invV)...), node_values)

    return nothing
end