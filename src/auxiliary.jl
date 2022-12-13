
"""
    generateRanges(dimTuple)

Returns a Tuple of ranges which span `1:sum(dimTuple)` with the property `length(res[i]) == dimTuple[i]`.
Used to associate dimensions with rows in matrices.
"""
function generateRanges(dimTuple::NTuple{D,Int})::NTuple{D, UnitRange{Int}} where {D}
    Tuple( (k>1 ? sum(dimTuple[1:k-1])+1 : 1):sum(dimTuple[1:k])  for (k,l) in enumerate(dimTuple) )
end