

struct IRKNGL_Cache{uT,tT,fT,pT}
    odef::fT # function defining the ODE system
    p::pT # parameters and so
    b::Vector{tT}
    c::Vector{tT}
    mu::Array{tT,2}
    nu::Array{tT,2}
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
    Dmin::Array{tT,1}
    maxiters::Int64
    step_number::Array{Int64,0}
    tau::Array{tT,0}
    initial_interp::Array{Int64,0}
    length_u::Int64
    length_q::Int64
    nrmdigits::Array{Int64, 0}
    E_weights::Array{tT,1}
end


function Rdigits(x::Real,r::Real)
    mx=r*x
    mxx=mx+x
    return mxx-mx
end


abstract type IRKNAlgorithm{s,initial_interp, m,myoutputs,nrmbits} <: OrdinaryDiffEqAlgorithm end
struct IRKNGL_Seq{s,initial_interp, m,myoutputs,nrmbits} <: IRKNAlgorithm{s, initial_interp, m,myoutputs,nrmbits} end
IRKNGL_Seq(;s=8, initial_interp=1, m=1,myoutputs=false,nrmbits=0)=IRKNGL_Seq{s, initial_interp, m,myoutputs,nrmbits}()

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
    alg::IRKNGL_Seq{s, initial_interp, m,myoutputs,nrmbits}, args...;
    dt,
    saveat=[],
    gamma=[],
    save_everystep=true,
    adaptive=false,
    maxiters=100,
    reltol=1e-6,
    abstol=1e-6,
    kwargs...) where {uType,tType,isinplace,s,m,initial_interp,myoutputs,nrmbits}

    trace=false

    destats = DiffEqBase.DEStats(0)
    @unpack f,u0,tspan,p,kwargs=prob
    f= SciMLBase.unwrapped_f(prob.f) 

    tType2=eltype(tspan)

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_interp

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

    (b, c, mu, nu, theta,omega, g, d) = IRKGLCoefficients_adap(s,dt)
    alpha=similar(theta)
    length_u = length(u0)
    length_q=div(length_u,2)
    dims = size(u0)
    indices=1:length_u
    indices1 = 1: length_q
    indices2 = (length_q+1):length_u


     uu = uType[]
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
     Dmin.=zero(uiType)

     E_weights=Array{uiType}(undef,length_u)
     if isempty(gamma)
         E_weights.=one(uiType)
     else
         for k in indices
              E_weights[k]=gamma[k]
         end
     end

     tau=Array{uiType,0}(undef)
     tau[]=zero(uiType)
 
     nrmdig = Array{Int64, 0}(undef)
     if (nrmbits > 0)
          nrmdig[] = 2^nrmbits
     else
         nrmdig[] = 0
     end
  
     irkngl_cache = IRKNGL_Cache(f,p,b,c,mu,nu,theta,omega,alpha,g,d,
                                 U, U_, L, dU, 
                                 F,Dmin,maxiters,step_number,tau,
                                 init_interp,length_u,length_q,nrmdig, E_weights)


    iters = Float64[]
                                
    push!(uu,copy(u0))
    push!(tt,t0)
    push!(iters,0.)
                                                             
    tj = [t0, zero(t0)]
    uj = copy(u0)  
    ej=zero(u0)  

    if isempty(saveat)     

        tstops=[Inf] 
#        save_everystep=true

    else

        if saveat isa Number 

            if t0<tf
                tstops=vcat(Inf, abs.(reverse(t0:saveat:tf)))
            else
                tstops=vcat(Inf, abs.(tf:-saveat:t0))
            end
            
        else
            tstops=vcat(Inf,abs.(reverse(saveat)))
        end

#       m=1
        save_everystep=false
        tj_=similar(tj)
        uj_=similar(uj)
        ej_=similar(ej)
        L_=deepcopy(L)
        dts_=similar(dts)
        if tstops[end]==t0 pop!(tstops) end 

    end

    tout=pop!(tstops)
    save_step=false


    if adaptive

        cont=true

        while cont

            tit=0
            j=0

            trace ? println("Urratsa j=",step_number[],",tj=", tj[1]) : nothing

            for i in 1:m

                j+=1
                step_number[] += 1

                if step_number[] ==1

                    (status,j_iter) = IRKNGLstep_fixed!(tj,tf,uj,ej,prob,dts,irkngl_cache)

                    E2=zero(uiType)
                    for k in indices2   
                        Ek= g[1]*F[1][k]
                        for js in 2:s
                            Ek= muladd(g[js], F[js][k],  Ek)
                        end
                        Ek=dt*Ek
                        wEk= E_weights[k]*Ek
                        E2=muladd(Ek,wEk,E2)
                    end

                    tau[]=E2^(1/(2s-2)) 

                    trace ? println("tau=", tau[]) : nothing

                else
                    (status,j_iter) = IRKNGLstep_adap!(tj,tf,uj,ej,prob,dts,irkngl_cache)
                end

                if (status=="Failure")
                    println("Fail")
                    sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                    if (myoutputs==true)
                        return(sol,iters, step_number[])
                    else
                        return(sol)
                    end
                end

                tit+= j_iter

                if tout<tj[1]+tj[2]

                    tj_.=tj
                    
                    for k in indices
                        uj_[k]=uj[k]
                        ej_[k]=ej[k]
                    end

                    for is in 1:s
                        for k in indices
                            L_[is][k]=L[is][k]
                        end
                    end

                    init_interp_=init_interp

                    dts_[1]=abs(tj[1]+tj[2]-tout)
                    dts_[2]=zero(dts[2])
                    dts_[3]=-dts[3]

                    init_interp=0
                    (status_,j_iter_) = IRKNGLstep_fixed!(tj_,tf, uj_,ej_,prob,dts_,irkngl_cache)
                    
                    if (status=="Failure")
                        println("Fail- Computing tout")
                        sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                        if (myoutputs==true)
                          return(sol,iters, step_number[])
                        else
                          return(sol)
                        end
                    end

                    for is in 1:s
                        for k in indices
                            L[is][k]=L_[is][k]
                        end
                    end

                    init_interp=init_interp_

                    push!(uu,uj_+ej_)
                    push!(tt,tj_[1]+tj_[2])

                    tout=pop!(tstops)

                elseif tout==tj[1]+tj[2]
                    save_step=true
                    tout=pop!(tstops)

                end         

                if (dts[1]==0)
                    cont=false
                    break
                end

            end # for 

            if save_step==true || save_everystep !=false || (cont==false)
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
                push!(iters, tit/j)
                save_step=false
            end

        end # end while

        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode=ReturnCode.Success)

        if (myoutputs==true)
           return(sol,iters, step_number[])
        else
           return(sol)
        end
  
    else # adaptive=false
                                     
        cont=true
        
        while cont

            tit=0
            j=0
  
            for i in 1:m

                j+=1
                step_number[] += 1

                (status,j_iter) = IRKNGLstep_fixed!(tj,tf,uj,ej,prob,dts,irkngl_cache)
  
                if (status=="Failure")
                    println("Fail")
                    sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= ReturnCode.Failure)
                    if (myoutputs==true)
                        return(sol,iters, step_number[])
                    else
                        return(sol)
                    end                    
                end
  
                tit+= j_iter

                if tout<tj[1]+tj[2]

                    tj_.=tj

                    for k in indices
                        uj_[k]=uj[k]
                        ej_[k]=ej[k]
                    end

                    for is in 1:s
                        for k in indices
                            L_[is][k]=L[is][k]
                        end
                    end

                    init_interp_=init_interp

                    dts_[1]=abs(tj[1]+tj[2]-tout)
                    dts_[2]=zero(dts[2])
                    dts_[3]=-dts[3]

                    init_interp=0
                    (status_,j_iter_) = IRKNGLstep_fixed!(tj_,tf, uj_,ej_,prob,dts_,irkngl_cache)
                    
                    if (status=="Failure")
                        println("Fail- Computing tout")
                        sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                        if (myoutputs==true)
                          return(sol,iters, step_number[])
                        else
                          return(sol)
                        end
                    end

                    for is in 1:s
                        for k in indices
                            L[is][k]=L_[is][k]
                        end
                    end

                    init_interp=init_interp_

                    push!(uu,uj_+ej_)
                    push!(tt,tj_[1]+tj_[2])

                    tit+=j_iter_-j_iter
                    push!(iters, tit/j)
                    j=0
                    tit=0

                    tout=pop!(tstops)

                elseif tout==tj[1]+tj[2]
                    save_step=true
                    tout=pop!(tstops)

                end                
  
                if (dts[1]==0)
                    cont=false
                    break
                end
            end
  
            if save_step==true || save_everystep !=false || (cont==false)
                push!(uu,uj+ej)
                push!(tt,tj[1]+tj[2])
                push!(iters, tit/j)
                save_step=false
            end

        end # end while

        dts[1]=dt
        dts[2]=zero(tType2)

        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)

        if (myoutputs==true)
            return(sol,iters,step_number[])
        else
            return(sol)
        end

 
        
    end #else adaptive
  
  
  end


