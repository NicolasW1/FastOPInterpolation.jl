"""
    generateTensorNodes!(res, tuple_nodes)

Generates a tensor product of input nodes `tuple_nodes`. The input is a tuple of node matrices.
Columns of the node input store coordinates and rows store different nodes.
The output is written into `res`. If the input nodes `x_i` have dimensions `size(d_i, n_i)` the output `size(res) = (sum(d_i), prod(n_i))`.
The nodes are iterated first to last, i.e.
`([e1_1 e1_2], [e2_1 e2_2]) -> [[e1_1, e2_1] [e1_2, e2_1] [e1_1, e2_2] [e1_2, e2_2]]`
"""
function generateTensorNodes!(res::AbstractMatrix{T}, tuple_nodes::NTuple{D, AbstractMatrix{T}})::Nothing where {D, T}
    @inbounds space_dims::NTuple{D, Int} = map(x->size(x)[1], tuple_nodes)
    @inbounds num_points::NTuple{D, Int} = map(x->size(x)[2], tuple_nodes)

    space_ranges = generateRanges(space_dims)

    cart_inds = CartesianIndices(num_points)
    lin_inds = LinearIndices(num_points)

    for ind in cart_inds
        for i=1:D
            @inbounds res[space_ranges[i], lin_inds[ind]] .= tuple_nodes[i][:, ind[i]]
        end
    end

    return nothing
end

generateTensorNodes!(res::AbstractMatrix{T}, tuple_nodes::AbstractMatrix{T}...) where {T} = generateTensorNodes!(res, tuple_nodes)

"""
    generateTensorNodes(tuple_nodes)

Allocating version of `generateTensorNodes!`
"""
function generateTensorNodes(tuple_nodes::NTuple{D, AbstractMatrix{T}}) where {D, T}
    @inbounds space_dims::NTuple{D, Int} = map(x->size(x)[1], tuple_nodes)
    @inbounds num_points::NTuple{D, Int} = map(x->size(x)[2], tuple_nodes)

    dim_res::Int = sum(space_dims)
    num_pts_res::Int = prod(num_points)

    res = Matrix{T}(undef, dim_res, num_pts_res)
    generateTensorNodes!(res, tuple_nodes)

    return res
end

"""
    generateTensorNodes(tuple_nodes)

Allocating version of `generateTensorNodes!`
"""
generateTensorNodes(tuple_nodes::AbstractMatrix{T}...) where {T} = generateTensorNodes(tuple_nodes)