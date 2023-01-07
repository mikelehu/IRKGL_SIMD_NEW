
include("IRKGL_Coefficients_adap.jl")


struct IRKNGL_Cache{uT,tT,fT,pT}
    odef::fT # function defining the ODE system
    p::pT # parameters and so
    b::Vector{tT}
    c::Vector{tT}
    mu::Array{tT,2}
    nu::Array{tT,2}
    nu1::Array{tT,2}
    nu2::Array{tT,2}
    theta::Array{tT,2}
    omega::Array{tT,2}
    alpha::Array{tT,2}
    g::Vector{tT}
    d::Vector{tT}
    U::Vector{uT}
    U_::Vector{uT}
    L::Vector{uT}
    dU::Vector{uT}
    F::Vector{uT}
    Dmin::Array{tT}
    maxiters::Int64
    step_number::Array{Int64,0}
    tau::Array{tT,0}
    first3iters::Array{Int64,1}
    initial_interp::Array{Int64,0}
    length_u::Int64
    length_q::Int64
    nrmdigits::Array{Int64, 0}
end


function Rdigits(x::Real,r::Real)
    mx=r*x
    mxx=mx+x
    return mxx-mx
end


abstract type IRKAlgorithm{s,initial_interp, m,myoutputs,nrmbits,backward} <: OrdinaryDiffEqAlgorithm end
struct IRKNGL_Seq{s,initial_interp, m,myoutputs,nrmbits,backward} <: IRKAlgorithm{s, initial_interp, m,myoutputs,nrmbits,backward} end
IRKNGL_Seq(;s=8, initial_interp=1, m=1,myoutputs=false,nrmbits=0,backward=false)=IRKNGL_Seq{s, initial_interp, m,myoutputs,nrmbits,backward}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
    alg::IRKNGL_Seq{s, initial_interp, m,myoutputs,nrmbits,backward}, args...;
    dt,
    saveat=[],
    save_everystep=true,
    adaptive=false,
    maxiters=100,
    reltol=1e-6,
    abstol=1e-6,
    kwargs...) where {uType,tType,isinplace,s,m,initial_interp,myoutputs,nrmbits,backward}

    trace=false

 #   println("Integrazio berria dt=",dt)

    destats = DiffEqBase.DEStats(0)
    @unpack f,u0,tspan,p,kwargs=prob
    t0=tspan[1]
    tf=tspan[2]
    tType2=eltype(tspan)


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

    (b, c, mu, nu, nu1, nu2, theta,omega, g, d) = IRKGLCoefficients_adap(s,dt)
    alpha=similar(theta)
    length_u = length(u0)
    length_q=div(length_u,2)
    dims = size(u0)
    indices=1:length_u
    indices1 = 1: length_q
    indices2 = (length_q+1):length_u

     ej=zero(u0)

     u0type=typeof(u0)
     uu = u0type[]
     tt = tType2[]

     U=Array{uType}(undef, s)
     for i in 1:s
         U[i]=zero(u0)
     end

     U_=deepcopy(U)
     L=deepcopy(U)
     dU=deepcopy(U)
     F=deepcopy(U)

     Dmin=Array{uiType}(undef,length_q)
     for i in 1:length_q Dmin[i]=0 end

     tau=Array{tType2,0}(undef)
     tau[]=0
 
     nrmdig = Array{Int64, 0}(undef)
     if (nrmbits > 0)
          nrmdig[] = 2^nrmbits
     else
         nrmdig[] = 0
     end
 
 
     irkngl_cache = IRKNGL_Cache(f,p,b,c,mu,nu,nu1,nu2,theta,omega,alpha,g,d,
                                 U, U_, L, dU, 
                                 F,Dmin,maxiters,step_number,tau,first3iters,
                                 init_interp,length_u,length_q,nrmdig)



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

                    j_iter1=0
                    (status,j_iter) = IRKNGLstep_fixed!(tj,tj1,uj,ej,prob,dts,irkngl_cache)

                    E2=zero(uiType)
                    for k in indices2   
                        E= g[1]*F[1][k]
                        for js in 2:s
                            E= muladd(g[js], F[js][k],  E)
                        end
                        E2=muladd(E,E,E2)
                    end

                    irkngl_cache.tau[]=(dt*sqrt(E2))^(1/(s-1))

                    trace ? println("tau=", tau[]) : nothing

                else
                    (status,j_iter) = IRKNGLstep_adap!(tj,tj1,uj,ej,prob,dts,irkngl_cache)
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

            if save_everystep !=false || (cont==false)
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
                push!(iters, tit/k)
            end

        end # end while

        dts[1]=dt
        dts[2]=0
        irkngl_cache.step_number[] = 0
    
    end # for stops

    sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode=ReturnCode.Success)

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
                (status,j_iter) = IRKNGLstep_fixed!(tj,tj1,uj,ej,prob,dts,irkngl_cache)
  
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

        end # end while

        dts[1]=dt
        dts[2]=0

    end # for tstops 

    sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)

    if backward==false

        if (myoutputs==true)
              return(sol,iters)
        else
              return(sol)
        end

      else

        dtprev=zero(tType2)
        t0=tspan[2]
        tf=tspan[1]
        signdt=sign(tf-t0) 
        dts=[dt,dtprev,signdt]
 
        uu = u0type[]
        tt = tType2[]
        
        iters_back = Float64[]
 
        push!(uu,uj+ej)
        push!(tt,tj[1]+tj[2])
        push!(iters_back,0.)
    
        cont=true
        j=0
                                        
        while cont
    
            tit=0
            k=0
      
            for i in 1:m
               j+=1
               k+=1
      
               irkngl_cache.step_number[] += 1
               (status,j_iter) = IRKNGLstep_fixed!(tj,tf,uj,ej,prob,dts,irkngl_cache)
      
               if (status=="Failure")
                  println("Fail")
                  sol_back=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= ReturnCode.Failure)
                  return(sol_back)
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
                push!(iters_back, tit/k)
            end
        end # end while
    
        sol_back=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)
    
        if (myoutputs==true)
            return(sol,iters,sol_back, iters_back)
        else
            return(sol,solback)
        end
    
    end # else backward
        
  end #else adaptive
  
  
  end


