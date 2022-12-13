using FastOPInterpolation
using Test

@testset "FastOPInterpolation.jl" begin
    f1(x) = sin(x[2] *  x[1] + x[3])
    f2(x) = 1 / (sum(x .^ 2) + 25)

    line1 = LineElement(51, LineGeometry(-1.5, 3.5), LineJacobiParameter(1.2, 3.5))
    triangle1 = TriangleElement(18, TriangleGeometry((0.0,0.0),(1.0,0.0), (0.0,1.0)), TriangleJacobiParameter(0.,0.,0.))
    disk1 = DiskElement(18, DiskGeometry((0.5, 0.5), 2.0), DiskRecurrenceParameter(0.))

end
