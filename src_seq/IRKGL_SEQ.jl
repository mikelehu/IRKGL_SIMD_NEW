__precompile__()

module IRKGL_SEQ

using Reexport
@reexport using DiffEqBase

using LinearAlgebra
using Parameters
using OrdinaryDiffEq
using RecursiveArrayTools
using BenchmarkTools
using SIMD

export  IRKGL_Seq, IRKAlgorithm
export  IRKNGL_Seq
export  IRKGLCoefficients_adap


include("IRKGL_Coefficients_adap.jl")
include("IRKGL_Seq_Solver.jl")
include("IRKGL_Step_Functions.jl")
include("IRKNGL_Seq_Solver.jl")
include("IRKNGL_Step_Functions.jl")

end # module
