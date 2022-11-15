export updateInterpolation!

@inline @inbounds @fastmath function updateInterpolation!(intElement::InterpolationElement, node_values::AbstractVector)
	# works, but the reverse is kinda curious
	mul!(vec(intElement.modes), kronecker(reverse(intElement.invV)...), node_values)

	return nothing
end