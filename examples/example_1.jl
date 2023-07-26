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