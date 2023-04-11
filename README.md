# IRKGL_SIMD.jl

- SIMD-vectorized implementation of high order IRK integrators

We present a preliminary version of a SIMD-vectorized implementation of the sixteenth order 8-stage implicit Runge-Kutta integrator  IRKGL16 implemented in the Julia package IRKGaussLegendre.jl. For numerical integrations of typical non-stiff problems performed in double precision, we show that a vectorized implementation of IRKGL16 that exploits the SIMD-based parallelism can outperform high order explicit Runge-Kutta schemes available in the standard package DifferentialEquations.jl.


## Description


The solver IRKGL16 is an implicit Runge-Kutta integrator of collocation type based on the Gauss-Legendre quadrature formula of 8 nodes. It is intended for high precision numerical integration of non-stiff systems of ordinary differential equations. In its sequential implementation, the scheme has interesting properties (symplecticness and time-symmetry) that make it particularly useful for long-term integrations of conservative problems. Such properties are also very useful for Scientific Machine Learning applications, as gradients can be exactly calculated by  integrating backward in time the adjoint equations.

For numerical integration of typical non-stiff problems with very high accuracy beyond the precision offered by double precision (i.e., standard IEEE binary64 floating precision) arithmetic our sequential implementation of IRKGL16 is more efficient than high order explicit Runge-Kutta schemes implemented in the standard package DifferentialEquations.jl. However,  our sequential implementation of IRKGL16 is generally unable to outperform them in double precision arithmetic.

We show that a vectorized implementation of IRKGL16 that exploits the SIMD-based parallelism offered by modern processor can be more efficient than high order explicit Runge-Kutta methods even for double precision computations. We demonstrate that by comparing our vectorized implementation of IRKGL16 with a 9th order explicit Runge-Kutta method (Vern9 from DifferentialEquations.jl) for different benchmark problems.

Our current implementation (https://github.com/mikelehu/IRKGL_SIMD.jl) depends on the Julia package SIMD.jl to efficiently perform computations on vectors with eight Float64 numbers. The right-hand side of the system of ODEs to be integrated has to be implemented as a generic function defined in terms of the arithmetic operations and elementary functions implemented for vectors in the package SIMD.jl. The state variables must be collected in an array of Float64 or Float32 floating point numbers. The SIMD-based vectorization process is performed automatically under the hood.


## Installation

```julia
julia>using Pkg
julia>Pkg.add(path="https://github.com/mikelehu/IRKGL_SIMD")
julia>using IRKGL_SIMD
```

## Solver options

### Available common arguments

- dt: stepsize
- save_everystep: default is true
- adaptive: false  (soon  we will include the adaptive step size strategy)
- maxiters: maximum number of iterations before stopping


### No common arguments

- s: number of stages of IRKGL, two options s=8 or s=4. Default 8.
- partitioned: false=General First Order ODE, true=Second Order ODE. Default false
- initial_interp: initialization method for stages.
        - =0  simplest initialization
        - =-1 interpolating from the stage values of previous step
- dim: the user must assign ```julia length(size(u0))``` value
- m: output saved at every m steps. Default 1.
- myoutputs: default false

## Example: Rigid-body problem

```julia
using LinearAlgebra,Plots
using DiffEqDevTools,BenchmarkTools
using OrdinaryDiffEq, IRKGL_SIMD, Parameters
using JLD2, FileIO
```

### Step 1: Defining the problem

```julia
function RigidBody!(du,u,p,t)
  I1=p[1]
  I2=p[2]
  I3=p[3]      
  du[1]  = I1*u[2]*u[3]
  du[2]  = I2*u[1]*u[3]
  du[3]  = I3*u[1]*u[2] +  0.25*sin(t)^2
end

p = [-2.0,1.25,-0.5]
u0 = [1.0;0.0;0.9]

dim=length(size(u0))
prob = ODEProblem(RigidBody!,u0,tspan,p);
```

### Step 2: Solving the problem

```julia
t0=0.
tF=100.0
dt=0.02
tspan=(t0,tF)

s=8
sol = solve(prob,IRKGL_simd(s=s, dim=dim,initial_interp=0),dt=dt)
plot(sol, xlabel="t", title="Rigid Body problem")
```

![Rigid-body](/Benchmarks/Rigid-body-Example1.png)


### Step 3: Benchmark

#### Test-solution

```julia
#sol =solve(prob_B,Vern9(),save_everystep=false, abstol=1e-30,reltol=1e-30)
#@save "./Data/Rigid_Body_sol100.jld2" sol
#test_sol = TestSolution(sol)

@load "./Data/Rigid_Body_sol100.jld2" sol
test_sol = TestSolution(sol)
final_state=sol.u[end]
```

#### Working-precision diagrams

```julia
tols=abstols=reltols=1.0 ./ 10.0 .^ (6:13)

dts8= 3.5*[0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05]
dtsVern=dts8/3;


nruns=10
s=8

setups = [ Dict(:alg=>Vern9(),:adaptive=>false,:dts=>dtsVern)
          Dict(:alg=>IRKGL_simd(s=s, dim=dim,initial_interp=0),:adaptive=>false,:dts=>dts8)
]
solnames = ["Vern9","IRKGL16-simd"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;appxsol=test_sol,save_everystep=false,numruns=nruns,names=solnames);

plot(wp)
```
![Rigid-body benchmark](/Benchmarks/Rigid-body-Example2.png)


## More Examples

[Benchmarks](https://github.com/mikelehu/IRKGL_SIMD.jl/tree/master/Benchmarks)
