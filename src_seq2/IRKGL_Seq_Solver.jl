
struct NireODESol{uType,tType,fType}
    t::Array{tType,1}
    u::Array{uType,1}
    iter::Array{Float64,1}
    retcode::Symbol
    f::fType
 end


struct IRKGL_Cache{uT,tT,fT,pT}
    odef::fT # function defining the ODE system
    p::pT # parameters and so
    b::Vector{tT}
    c::Vector{tT}
    mu::Array{tT,2}
    nu::Array{tT,2}
    U::Vector{uT}
    U_::Vector{uT}
    L::Vector{uT}
    dU::Vector{uT}
    F::Vector{uT}
    Dmin::Array{tT,1}
    maxiters::Int64
    step_number::Array{Int64,0}
    initial_interp::Array{Int64,0}
    length_u::Int64
    nrmdigits::Array{Int64, 0}
end


function Rdigits(x::Real,r::Real)

    mx=r*x
    mxx=mx+x
    return mxx-mx

end



function IRKGL_sek(u0::uType, dt::tType, t0::tType, tf::tType, f::fType, p::pType; 
                   s=8,m=1,initial_interp=1,myoutputs=false,nrmbits=0, save_everystep=true,
                   saveat=[],gamma=[], maxiters=100, adaptive=false) where {uType, tType, fType, pType}

    trace = false

    u0type=typeof(u0)
    uu = u0type[]
    tt = tType[]

    tType2=tType

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_interp

    uiType=eltype(uType)

    dtprev=zero(tType2) 
    signdt=sign(tf-t0)

    dts=[dt,dtprev,signdt]

    (b, c, mu, nu) = IRKGLCoefficients(s,dt)

    length_u = length(u0)
    dims = size(u0)
    indices=1:length_u

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

     Dmin=Array{uiType}(undef,length_u)
     Dmin.=zero(uiType)
 
     nrmdig = Array{Int64, 0}(undef)
     if (nrmbits > 0)
          nrmdig[] = 2^nrmbits
     else
         nrmdig[] = 0
     end

 
     irkgl_cache = IRKGL_Cache(f,p,b,c,mu,nu,
                               U, U_, L, dU, 
                               F,Dmin,maxiters,step_number,
                               init_interp,length_u,nrmdig)


     iters = Float64[]
                                
     push!(uu,copy(u0))
     push!(tt,t0)
     push!(iters,0.)
                                                           
     tj = [t0, zero(t0)]
     uj = copy(u0)
     ej= zero(u0)

     tj_=similar(tj)
     uj_=similar(uj)
     ej_=similar(ej)
     L_=deepcopy(L)
     dts_=similar(dts)

     tstops=[]

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
                               
 #      m=1
        save_everystep=false
        if tstops[end]==t0 pop!(tstops) end 
                              
     end
     
     tout= Array{uiType,0}(undef)
     tout[]=pop!(tstops)
     save_step=false    

    if adaptive == true

      println("Error: adaptive ==true")
    
    else

     cont=true
        
     while cont
          
            tit=0
            j=0
  
            for i in 1:m

                j+=1  
                step_number[]+= 1

                (status,j_iter) = IRKGLstep_fixed!(tj,tf,uj,ej,dts,irkgl_cache)
  
                if (status=="Failure")
                    println("Fail")
 #                   sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= ReturnCode.Failure)
                    sol = NireODESol(tt,uu,iters,:Successs,f)
                    return(sol)
                end
  
                tit+=j_iter

                
                if tout[]<tj[1]+tj[2]

                    println("emaitza itzuli ...")

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
    
                    dts_[1]=abs(tj[1]+tj[2]-tout[])
                    dts_[2]=zero(dts[2])
                    dts_[3]=-dts[3]
    
                    init_interp=0
                    (status_,j_iter_) = IRKGLstep_fixed!(tj_,tf, uj_,ej_,dts_,irkgl_cache)
                    
                    if (status=="Failure")
                        println("Fail- Computing tout")
#                        sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                        sol = NireODESol(tt,uu,iters,:Successs,f)
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
    
                    tout[]=pop!(tstops)
    
                elseif tout[]==tj[1]+tj[2]
                
                    save_step=true
                    tout[]=pop!(tstops)
    
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

    end


 #       sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)

    sol = NireODESol(tt,uu,iters,:Successs,f)
     
    return(sol)

  
  
end


