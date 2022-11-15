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

Unfortunatley, this package is currently not compatible with shared-memory parallelization or autodiff due to the extensive usage of buffers to maximize performance for repeated evaluations. We are aware of the issue and currently working on a solution.

The package was developed by Nicolas Wink and Eduardo Grossi. If you should use it in your scientific work please cite this repository (there is a button on the GitHub page).