
include("IRKGL_SIMD_Coefficients.jl")
include("VecArray_def.jl")


struct IRKNGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}
    odef::fType # function defining the ODE system
    p::pType # parameters and so
    b::Vec{s_,floatType}
    c::Vec{s_,floatType}
    mu::VecArray{s_,floatType,2}
    nu::VecArray{s_,floatType,2}
    nu1::VecArray{s_,floatType,2}
    nu2::VecArray{s_,floatType,2}
    theta::VecArray{s_,floatType,2}
    omega::VecArray{s_,floatType,2}
    alpha::VecArray{s_,floatType,2}
    g::Vec{s_,floatType}
    d::Vec{s_,floatType}
    s_beta::Vec{s_,floatType}
    U::VecArray{s_,floatType,dim_}
    U_::VecArray{s_,floatType,dim_}
    L::VecArray{s_,floatType,dim_}
    dU::VecArray{s_,floatType,dim_}
    F::VecArray{s_,floatType,dim_}
    Dmin::Array{floatType,dim}
    maxiters::Int64
    step_number::Array{Int64,0}
    tau::Array{floatType,0}
    first3iters::Array{Int64,1}
    initial_interp::Array{Int64,0}
    length_u::Int64
    length_q::Int64
    nrmdigits::Array{Int64, 0}
    E_weights::Array{floatType,dim}
end


function Rdigits(x::Real,r::Real)
    mx=r*x
    mxx=mx+x
    return mxx-mx
end

function Rdigits(x::Vec,r::Real)
    mx=r*x
    mxx=mx+x
    return mxx-mx
end


abstract type IRKAlgorithm{s,initial_interp, dim,floatType,m,myoutputs,nrmbits} <: OrdinaryDiffEqAlgorithm end
struct IRKNGL_simd{s, initial_interp, dim,floatType,m,myoutputs,nrmbits} <: IRKAlgorithm{s, initial_interp, dim,floatType,m,myoutputs,nrmbits} end
IRKNGL_simd(;s=8, initial_interp=1, dim=1,floatType=Float64,m=1,myoutputs=false,nrmbits=0)=IRKNGL_simd{s, initial_interp, dim,floatType,m,myoutputs,nrmbits}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
    alg::IRKNGL_simd{s,initial_interp, dim,floatType, m,myoutputs,nrmbits}, args...;
    dt,
    saveat=[],
    gamma=[],
    save_everystep=true,
    adaptive=false,
    maxiters=100,
    reltol=1e-6,
    abstol=1e-6,
    kwargs...) where {floatType<: Union{Float32,Float64},uType,tType,isinplace,dim,s,m,initial_interp,myoutputs,nrmbits}

    trace=false
   
    destats = DiffEqBase.DEStats(0)
    @unpack f,u0,tspan,p,kwargs=prob

    tType2=eltype(tspan)

    utype = Vector{floatType}
    ttype = floatType

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_interp
    first3iters=[0,0,0]

    dts = Array{tType2}(undef, 1)
    uiType=eltype(u0)

    dtprev=zero(tType2)
    signdt=sign(tspan[2]-tspan[1]) 
               
    if signdt==1 
        t0=prob.tspan[1]
        tf=prob.tspan[2]   # forward
     else
        t0=prob.tspan[2]
        tf=prob.tspan[1]   # backward
     end

    dts=[dt,dtprev,signdt]

    Treltol=convert(uiType,reltol)
    Tabstol=convert(uiType,abstol)

    (b_, c_, mu_, nu_, nu1_, nu2_,theta_, omega_, g_, d_) = IRKGLCoefficients(s,dt)
    length_u = length(u0)
    length_q=div(length_u,2)
    dims = size(u0)
    indices=1:length_u
    indices2 = (length_q+1):length_u

    c = vload(Vec{s,floatType}, c_, 1)
    b = vload(Vec{s,floatType}, b_, 1)
    nu=VecArray{s,floatType,2}(nu_)
    nu1=VecArray{s,floatType,2}(nu1_)
    nu2=VecArray{s,floatType,2}(nu2_)
    mu=VecArray{s,floatType,2}(mu_)
    theta=VecArray{s,floatType,2}(theta_)
    omega=VecArray{s,floatType,2}(omega_)
    alpha=deepcopy(theta)
    g = vload(Vec{s,floatType}, g_, 1)
    d = vload(Vec{s,floatType}, d_, 1)
    s_beta= deepcopy(g)

    ej=zero(u0)

    u0type=typeof(u0)
    uu = u0type[]
    tt = ttype[]

    zz=zeros(Float64, s, dims...)
    U=VecArray{s,Float64,length(dims)+1}(zz)
    U_=deepcopy(U)
    L=deepcopy(U)
    dU=deepcopy(U)
    F=deepcopy(U)

    Dmin=Array{uiType}(undef,length_q)
    for i in 1:length_q Dmin[i]=0 end

    E_weights=Array{uiType}(undef,length_u)
    if isempty(gamma)
        E_weights.=one(uiType)
     else
#        s_gamma=sum(gamma)
        for k in indices
#             E_weights[k]=gamma[k]/sum(gamma) 
             E_weights[k]=gamma[k]
        end
     end


    tau=Array{tType2,0}(undef)
    tau[]=0

    nrmdig = Array{Int64, 0}(undef)
    if (nrmbits > 0)
        nrmdig[] = 2^nrmbits
    else
        nrmdig[] = 0
    end

    irkngl_cache = IRKNGL_SIMD_Cache(f,p,b,c,mu,nu,nu1,nu2,theta,omega,alpha,g,d,s_beta,
                                     U, U_, L, dU,
                                     F,Dmin,maxiters,step_number,tau, first3iters,
                                     init_interp,length_u,length_q,nrmdig, E_weights)

    if adaptive

        iters = Float64[]

        push!(uu,copy(u0))
        push!(tt,t0)
        push!(iters,0.)

        tj = [t0, zero(t0)]
        uj = copy(u0)

        j=0

        if isempty(saveat)             
            tstops=[tf]        
        else
            tstops=saveat
            save_everystep=false
#           m=1
        end

        for tj1 in tstops

             trace ? println("Urratsa tj=", tj[1], ",tj1=", tj1) : nothing

             cont=true

             while cont

                tit=0
                k=0

                for i in 1:m
                    j+=1
                    k+=1

                    irkngl_cache.step_number[] += 1

                    if (irkngl_cache.step_number[] ==1)

                        (status,j_iter) = IRKNGLstep_SIMD_fixed!(tj,tj1, uj,ej,prob,dts,irkngl_cache)
               
                        E2=zero(uiType)
                        for k in indices2
                            Fk=getindex_(F,k) 
                            E=sum(g*Fk)   # E = dot(g, Fk)
#                            E=E_weights[k]*sum(g*Fk)   # E = dot(g, Fk)
                            E2=muladd(E,E,E2)
                        end
      
                        irkngl_cache.tau[]=(dt*sqrt(E2))^(1/(s-1))   
                        trace ? println("tau=", tau[]) : nothing

                    else
                        (status,j_iter) = IRKNGLstep_SIMD_adap!(tj,tj1, uj,ej,prob,dts,irkngl_cache)
                    end

                    if (status=="Failure")
                        println("Fail")
                        sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                        if (myoutputs==true)
                             return(sol,iters)
                        else
                             return(sol)
                        end
                    end

                    tit+= j_iter

                    if (dts[1]==0)
                        cont=false
                        break
                    end
                end
myoutputs
                if save_everystep !=false || (cont==false)
                    push!(uu,uj+ej)
                    push!(tt,tj[1]+tj[2])
                    push!(iters, tit/k)
                end

            end #end while

            dts[1]=dt
            dts[2]=0
            irkngl_cache.step_number[] = 0


        end # for tstops

        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)

        if (myoutputs==true)
            return(sol,iters)
        else
            return(sol)
        end


    else # adaptive=false

        iters = Float64[]

        push!(uu,copy(u0))
        push!(tt,t0)
        push!(iters,0.)
                                   
        tj = [t0, zero(t0)]
        uj = copy(u0)
                                    
        j=0

        if isempty(saveat)             
            tstops=[tf]        
        else
            tstops=saveat
            save_everystep=false
#           m=1
        end

        for tj1 in tstops

           cont=true
                               
           while cont

                tit=0
                k=0
  
                for i in 1:m
                  j+=1
                  k+=1
  
                  irkngl_cache.step_number[] += 1
                 (status,j_iter) = IRKNGLstep_SIMD_fixed!(tj,tj1,uj,ej,prob,dts,irkngl_cache)
  
                  if (status=="Failure")
                    println("Fail")
                    sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= ReturnCode.Failure)                   
                    return(sol)
                  end
  
                  tit+= j_iter
  
                  if (dts[1]==0)
                       cont=false
                       break
                  end
                end
  
                if save_everystep !=false || (cont==false)
                   push!(uu,uj+ej)
                   push!(tt,tj[1]+tj[2])
                   push!(iters, tit/k)
                end

            end # while

            dts[1]=dt
            dts[2]=0

        end # for tstops 

        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)
  
        if (myoutputs==true)
            return(sol,iters)
        else
            return(sol)
        end
 
    end # adaptive-else

end

