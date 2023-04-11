

include("VecArray_def.jl")

struct NireODESol{uType,tType,fType}
    t::Array{tType,1}
    u::Array{uType,1}
    iter::Array{Float64,1}
    retcode::Symbol
    f::fType
 end


struct IRKGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}
    odef::fType # function defining the ODE system
    p::pType # parameters and so
    b::Vec{s_,floatType}
    c::Vec{s_,floatType}
    mu::VecArray{s_,floatType,2}
    nu::VecArray{s_,floatType,2}
    U::VecArray{s_,floatType,dim_}
    U_::VecArray{s_,floatType,dim_}
    L::VecArray{s_,floatType,dim_}
    dU::VecArray{s_,floatType,dim_}
    F::VecArray{s_,floatType,dim_}
    Dmin::Array{floatType,dim}
    maxiters::Int64
    step_number::Array{Int64,0}
    initial_interp::Array{Int64,0}
    length_u::Int64
    nrmdigits::Array{floatType, 0}
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




function IRKGL_simd_2(u0::uType, dt::tType, t0::tType, tf::tType, f::fType, p::pType; 
                      s=8,m=1,dim=1, initial_interp=1,myoutputs=false,nrmbits=0, save_everystep=true,
#                      saveat=[], maxiters=100, adaptive=false) where {floatType<: Union{Float32,Float64},uType, tType, fType, pType}
                      saveat=[], maxiters=100, adaptive=false) where {uType, tType, fType, pType}
                      

    trace=false

    tspan=(t0,tf)
    dtprev=zero(tType)
    signdt=sign(tspan[2]-tspan[1]) 

    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    init_interp =  Array{Int64,0}(undef)
    init_interp[] = initial_interp


    uiType=eltype(u0)               
    if signdt==1 
        t0=tspan[1]
        tf=tspan[2]   # forward
    else
        t0=tspan[2]
        tf=tspan[1]   # backward
    end

    dts=[dt,dtprev,signdt]     

    (b_, c_, mu_, nu_) = IRKGLCoefficients(s,dt)
    length_u = length(u0)
    dims = size(u0)
    indices=1:length_u

    c = vload(Vec{s,uiType}, c_, 1)
    b = vload(Vec{s,uiType}, b_, 1)
    nu=VecArray{s,uiType,2}(nu_)
    mu=VecArray{s,uiType,2}(mu_)

    uu = uType[]
    tt = tType[]

    zz=zeros(uiType, s, dims...)
    U=VecArray{s,uiType,length(dims)+1}(zz)
    U_=deepcopy(U)
    L=deepcopy(U)
    dU=deepcopy(U)
    F=deepcopy(U)

#    Dmin = zero(u0)
    Dmin=Array{uiType}(undef,length_u)
    Dmin.=zero(uiType)


    nrmdig = Array{uiType, 0}(undef)
    if (nrmbits > 0)
        nrmdig[] = uiType(2^nrmbits)
    else
        nrmdig[] = zero(uiType)
    end


    irkgl_cache = IRKGL_SIMD_Cache(f,p,b,c,mu,nu,
                                   U, U_, L, dU,
                                   F,Dmin,maxiters,step_number,
                                   init_interp,length_u,nrmdig)


    iters = Float64[]
    push!(uu,copy(u0))
    push!(tt,t0)
    push!(iters,0.)
                                                                               
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej=zero(u0)
    
    
    tj_=similar(tj)
    uj_=similar(uj)
    ej_=similar(ej)
    L_=deepcopy(L)
    dts_=similar(dts)

    tstops=tType[]
                               
    if isempty(saveat)     
                               
        tstops=tType[Inf]
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
                                   
    if adaptive
     
        println("Error: adaptive ==true")

    else # adaptive=false

       cont=true

       while cont

          tit=0
          j=0
  
          for i in 1:m

              j+=1  
              step_number[] += 1

              (status,j_iter) = IRKGLstep_SIMD_fixed!(tj,tf, uj,ej,dts,irkgl_cache)
  
              if (status=="Failure")
                  println("Fail")
 #                 sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= ReturnCode.Failure)
                  sol = NireODESol(tt,uu,iters,:Successs,f)
                  if (myoutputs==true)
                    return(sol,iters, step_number[])
                  else
                    return(sol)
                  end
              end
  
              tit+= j_iter
  
              if tout[]<tj[1]+tj[2]

                tj_.=tj
                uj_.=uj
                ej_.=ej
                L_.data.=L.data
                init_interp_=init_interp

                dts_[1]=abs(tj[1]+tj[2]-tout[])
                dts_[2]=zero(dts[2])
                dts_[3]=-dts[3]

                init_interp=0
                (status_,j_iter_) = IRKGLstep_SIMD_fixed!(tj_,tf, uj_,ej_,dts_,irkngl_cache)
                
                if (status=="Failure")
                    println("Fail- Computing tout")
 #                   sol=DiffEqBase.build_solution(prob,alg,tt,uu, retcode= ReturnCode.Failure)
                    sol = NireODESol(tt,uu,iters,:Successs,f)
                    if (myoutputs==true)
                      return(sol,iters, step_number[])
                    else
                      return(sol)
                    end
                end

                L.data.=L_.data
                init_interp=init_interp_

                push!(uu,uj_+ej_)
                push!(tt,tj_[1]+tj_[2])

                tit+=j_iter_-j_iter
                push!(iters, tit/j)
                j=0
                tit=0

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

        end # while

#        sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= ReturnCode.Success)
        sol = NireODESol(tt,uu,iters,:Successs,f)

        if (myoutputs==true)
            return(sol,iters,step_number[])
        else
           return(sol)
        end
 
    end # adaptive-else

end

