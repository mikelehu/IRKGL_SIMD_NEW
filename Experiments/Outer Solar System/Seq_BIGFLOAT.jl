using LinearAlgebra
using Plots
using IRKGaussLegendre
using OrdinaryDiffEq
#using BenchmarkTools


PATH_SRC="../src_seq1/"

include(string(PATH_SRC,"IRKGL_SEQ1.jl"))
using .IRKGL_SEQ1   

PATH_SRC="../src_seq2/"

include(string(PATH_SRC,"IRKGL_SEQ2.jl"))
using .IRKGL_SEQ2 


PATH_ODES="../ODEProblems/"

include(string(PATH_ODES,"Initial5Body.jl"))
include(string(PATH_ODES,"Nbody.jl"))
include(string(PATH_ODES,"Nbody2nd.jl"))


u0, Gm, bodylist = Initial5Body(BigFloat)
q0=u0[:,:,1]
v0=u0[:,:,2]
dim=length(size(u0))

N = length(Gm)

show(bodylist)
E0=NbodyEnergy(u0,Gm)


t0 = BigFloat(0.)
dt = BigFloat(500.)  # 500.
tF=  BigFloat(1e6)
tF = BigFloat(10*dt)

#n = 100
#n = 5

#m = convert(Int64,ceil(abs(tF-t0)/(n*dt)))
#n = convert(Int64,ceil(abs(tF-t0)/(m*dt))) # Number of macro-steps (Output is saved for n+1 time values)
#dt = (tF-t0)/(n*m)

m=1
println("dt = $dt, m=$m")

prob = ODEProblem(NbodyODE!, u0, (t0,tF), Gm)

# tol=1e-6

# Formulazio berria
alg=IRKGL_Seq2(s=8,partitioned=false, initial_interp=0, m=m,myoutputs=true)
sol1,iters1=solve(prob,alg,dt=dt, adaptive=false); 