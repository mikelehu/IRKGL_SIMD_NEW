__precompile__()

module IRKGL_SIMD2

using Reexport
@reexport using DiffEqBase

using LinearAlgebra
using Parameters
using OrdinaryDiffEq
using SIMD

export  IRKGL_simd_2, IRKNGL_simd_2
export  VecArray, floatType

include("../src/IRKGL_Coefficients.jl")
include("IRKGL_SIMD_Solver.jl")
include("IRKGL_SIMD_Step_Functions.jl")
include("IRKNGL_SIMD_Solver.jl")
include("IRKNGL_SIMD_Step_Functions.jl")


end # module
