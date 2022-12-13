# FastOPInterpolation

Fast multidimensional interpolation based on orthogonal polynomials. Supports domains build from tensor products, e.g. a prism as product of a line and a triangle.

Currently implemented combinations of Domains & Polynomials are

| Geometry               | Line     | Triangle       | Disk           |
|------------------------|----------|----------------|----------------|
| Dimension              | 1d       | 2d             | 2d             |
| Polynomials            | Jacobi   | Koornwinder IV | Koornwinder II |
| Max order predef nodes | $\infty$ | 18             | 30             |


You can specify the set of nodes yourself or use the predefined ones, their maximum order is also stated in the table above.

A complete working example can be found [the example folder](examples/).

For a given set of orthogonal polynomials (with their recurrence relation) it is easy to extend the code to other geometries. Please do not hesitate to get in contact with us (via mail or github issue) to request further geometries or features. You are of course also always invited to implement it yourself and create a pull request.

Current limitations include:
 - no shared-memory parallelization (due to buffers)
 - autodiff (same reason as shared-memory parallelization)
 - only scalar functions (matching type of coordinates and return type of the function to be interpolated)

All three limitations are partially related to the heavy use internal buffers to maximize performance for repeated evaluations.
We are working on a simulation to circument these limitations.

Similar packages (interpolation with orthogonal polynomials) with a different focus:
 - [FastChebInterp](https://github.com/stevengj/FastChebInterp.jl)
 - [ClassicalOrthogonalPolynomials](https://github.com/JuliaApproximation/ClassicalOrthogonalPolynomials.jl)

The package was developed by Nicolas Wink and Eduardo Grossi. If you should use it in your scientific work please cite this repository (there is a button on the GitHub page).