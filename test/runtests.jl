using FastOPInterpolation
using Test

@testset "general" begin
    f1(x) = sin(x[5] *  x[1] + x[3] - x[4] * x[2])
    f2(x) = 1 / (sum(x .^ 2) + 25)
    x₀ = [0.23, 0.4, 0.1, 1., 1.]
    x₁ = [-1.4, 0.05, 0.8, 0., -1.]

    line = LineElement(51, LineGeometry(-1.5, 3.5), LineJacobiParameter(1.2, 3.5))
    triangle = TriangleElement(18, TriangleGeometry((0.0,0.0),(1.0,0.0), (0.0,1.0)), TriangleJacobiParameter(0.,0.,0.))
    disk = DiskElement(30, DiskGeometry((0.5, 0.5), 2.0), DiskRecurrenceParameter(0.))

    intP, t_nodes = InterpolationElement((line, triangle, disk))
    f_nodes = generateTensorNodes(t_nodes)

    t_vals1 = [f1(x) for x in eachcol(f_nodes)]
    updateInterpolation!(intP, t_vals1)

    @test f1(x₀) ≈ intP(x₀)
    @test f1(x₁) ≈ intP(x₁)

    t_vals2 = [f2(x) for x in eachcol(f_nodes)]
    updateInterpolation!(intP, t_vals2)

    @test f2(x₀) ≈ intP(x₀)
    @test f2(x₁) ≈ intP(x₁)
end
