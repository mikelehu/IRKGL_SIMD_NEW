# IRKGL_SIMD.jl

- SIMD-vectorized implementation of high order IRK integrators

We present a preliminary version of a SIMD-vectorized implementation of the sixteenth order 8-stage implicit Runge-Kutta integrator  IRKGL16 implemented in the Julia package IRKGaussLegendre.jl. For numerical integrations of typical non-stiff problems performed in double precision, we show that a vectorized implementation of IRKGL16 that exploits the SIMD-based parallelism can outperform high order explicit Runge-Kutta schemes available in the standard package DifferentialEquations.jl.


## Description


The solver IRKGL16 is an implicit Runge-Kutta integrator of collocation type based on the Gauss-Legendre quadrature formula of 8 nodes. It is intended for high precision numerical integration of non-stiff systems of ordinary differential equations. In its sequential implementation, the scheme has interesting properties (symplecticness and time-symmetry) that make it particularly useful for long-term integrations of conservative problems. Such properties are also very useful for Scientific Machine Learning applications, as gradients can be exactly calculated by  integrating backward in time the adjoint equations.

For numerical integration of typical non-stiff problems with very high accuracy beyond the precision offered by double precision (i.e., standard IEEE binary64 floating precision) arithmetic our sequential implementation of IRKGL16 is more efficient than high order explicit Runge-Kutta schemes implemented in the standard package DifferentialEquations.jl. However,  our sequential implementation of IRKGL16 is generally unable to outperform them in double precision arithmetic.

We show that a vectorized implementation of IRKGL16 that exploits the SIMD-based parallelism offered by modern processor can be more efficient than high order explicit Runge-Kutta methods even for double precision computations. We demonstrate that by comparing our vectorized implementation of IRKGL16 with a 9th order explicit Runge-Kutta method (Vern9 from DifferentialEquations.jl) for different benchmark problems.

Our current implementation (https://github.com/mikelehu/IRKGL_SIMD.jl) depends on the Julia package SIMD.jl to efficiently perform computations on vectors with eight Float64 numbers. The right-hand side of the system of ODEs to be integrated has to be implemented as a generic function defined in terms of the arithmetic operations and elementary functions implemented for vectors in the package SIMD.jl. The state variables must be collected in an array of Float64 or Float32 floating point numbers. The SIMD-based vectorization process is performed automatically under the hood.

