module FastOPInterpolation

using StaticArrays
using LinearAlgebra
using Kronecker

include("auxiliary.jl")

include("general_nodes.jl")
include("general_recurrence.jl")
include("general_elements.jl")

include("elements/line.jl")
include("elements/triangle.jl")
include("elements/disk.jl")

include("interpolation_object.jl")
include("evaluation.jl")
include("update_modes.jl")

end