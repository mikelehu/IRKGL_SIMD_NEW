

#
#  IRKGLstep_adap!
#  IRKGLstep_fixed!



function IRKGLstep_adap!(ttj,uj,ej,prob,dts,irkgl_cache::IRKGL_Cache{uType,tType2,fType,pType}) where {uType,tType2,fType,pType}
    trace = true
    f = irkgl_cache.odef
    p = irkgl_cache.p
    b = irkgl_cache.b
    c = irkgl_cache.c
    hi1 = irkgl_cache.hi1
    mu = irkgl_cache.mu
    nu = irkgl_cache.nu
    nu1 = irkgl_cache.nu1
    nu2 = irkgl_cache.nu2
    U = irkgl_cache.U
    Uz = irkgl_cache.Uz
    U_ = irkgl_cache.U_
    L = irkgl_cache.L
    Lz = irkgl_cache.Lz
    E = irkgl_cache.E
    D = irkgl_cache.D
#    Ws = irkgl_cache.Ws
#    Wz = irkgl_cache.Wz
    gamma=irkgl_cache.gamma
    theta=irkgl_cache.theta
    omega=irkgl_cache.omega
    beta=irkgl_cache.beta
    alpha=irkgl_cache.alpha
    g=irkgl_cache.g
    d=irkgl_cache.d
    F = irkgl_cache.F
    Dmin = irkgl_cache.Dmin
    step_number = irkgl_cache.step_number[]
    tau= irkgl_cache.tau[]
    first3iters = irkgl_cache.first3iters
    initial_interp = irkgl_cache.initial_interp[]
    len = irkgl_cache.length_u
    nrmdigits=irkgl_cache.nrmdigits[]
    s = length(b)
    maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
    tj = ttj[1]
    te = ttj[2]
    indices=1:len
    flzero = zero(tType2)


    uEtype=eltype(uj)
    R=uEtype(2)^(precision(uEtype)-20)

    dt=dts[1]
    dtz=dts[2]
    sdt=dts[3]
    tf=prob.tspan[2]


    if (initial_interp==1)
        for is in 1:s
            for k in indices
                dUik = muladd(nu[is,1], L[1][k], ej[k])
                for js in 2:s
                   dUik = muladd(nu[is,js], L[js][k], dUik)
                end
                U[is][k] =  uj[k]  + dUik
            end
        end
    else
        for is in 1:s
            for k in indices
               U[is][k] = uj[k] + ej[k]
            end
        end
    end

    trace ? println("Urratsa tj=",tj+te,", iter=", 1, ",dt=", Float32(dt)) : nothing

    j_iter = 0
    Dmin .= Inf
    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU =true
    dtj1=dt
    fase3=false
    Dtmin=Inf
    plusDt=2

   @inbounds while (iter && j_iter<maxiters) 
#   @inbounds while (j_iter<40)   #Float64 j_iter<22, BigFloat<40

        j_iter += 1
        dtj0=dtj1

        trace ? println("*** ", j_iter, ". iterazioa:") : nothing

        for is in 1:s
            f(F[is], U[is], p,  tj  + dt*c[is])
            for k in indices
                Uz[is][k] = U[is][k]
                L[is][k] = dt*(b[is]*F[is][k])
            end
        end

        for is in 1:s
            for k in indices
               dUik = muladd(mu[is,1], L[1][k], ej[k])
               for js in 2:s
                  dUik = muladd(mu[is,js], L[js][k], dUik)
               end
               U[is][k] =  uj[k] + dUik
            end
        end
  
        diffU = false
        iter= false 
  
        for k in indices
            Uk =[u[k] for u in U]
            Uzk =[u[k] for  u in Uz]
            DY = maximum(abs.(Uk-Uzk))
  
            if DY>0
                diffU = true
                if DY< Dmin[k]
                   Dmin[k]=DY
                   iter=true
                end
            end
        end
  
        if (!iter && diffU && plusIt)
            iter=true
            plusIt=false
        else
            plusIt=true
        end

        trace ? println("iter=", iter, ",diffU=", diffU, ",plusIt=", plusIt, ", norm(Dmin)=",norm(Dmin)) : nothing

#=
        for k in indices
               Ws[k] = gamma[s,1]*F[1][k]
              for js in 2:s
                   Ws[k] = muladd(gamma[s,js], F[js][k],  Ws[k])
              end
              Ws[k]=dt^(1-s)*Ws[k]
        end
=#

        for k in indices
            E[k] = g[1]*F[1][k]
            D[k] = d[1]*F[1][k]
            for js in 2:s
                E[k] = muladd(g[js], F[js][k],  E[k])
                D[k] = muladd(d[js], F[js][k],  D[k])
            end
            E[k]=dt*E[k]
        end

#        aux=abs(tau/norm(Ws)^(1/s))
        E2 = dot(E,E)
        DE = dot(D,E)
        #=
        B = E2^(1/(2s-2))
        Bder = DE*B/((s-1)*E2) 
        aux = dt-(B-tau)/Bder
        =#
        fz = log(E2) - (2s-2)*log(tau)
        fderz = 2*dt*DE/E2
        aux = dt*exp(-fz/fderz)
        
        dtj1=sdt*min(aux,abs(tf-(ttj[1]+ttj[2])))
        Dt=abs((dtj0-dtj1)/dtj1)
            
        if (0<Dt<0.0001) && !fase3

            Dtmin=min(Dt,Dtmin)

            if (Dt>Dtmin && plusDt==1) 
                fase3=true 
                plusDt-=1
            elseif (Dt>Dtmin && plusDt>1) 
                  plusDt-=1 
                else 
                  plusDt=2
            end  

            if trace
                println("eguneratu", " fase3=",fase3, ", dt=",Float32(dtj1) ,", dtj0=", Float32(dtj0),
            " ,dtj1=", Float32(dtj1), ", Dtmin=", Float32(Dtmin),", Dt=", Float32(Dt), ", plusDt=", plusDt)
                println()
            end
            for is in 1:s
                for k in indices
                    U_[is][k]=U[is][k]
                end
            end
          
            lambda=dtj1/dt-1

            for i in 1:s

                sumbetai=0 
                for j in 1:s+1
                    aux=muladd(lambda,theta[i,j],omega[i,j])
                    beta[i,j]=1/aux
                    sumbetai+=beta[i,j]
                end
                %
                for j in 1:s+1
                    alpha[i,j]=beta[i,j]/sumbetai
                end
            end

            for is in 1:s
                for k in indices
                    dUik = muladd(alpha[is,2],U_[1][k], ej[k])
                    for js in 3:s+1
                        dUik = muladd(alpha[is,js], U_[js-1][k], dUik)
                    end
                       U[is][k] = muladd(alpha[is,1],uj[k],dUik)
                end
            end

            dt=dtj1

        else  
            if !fase3 
                plusDt=2 
                Dtmin=min(Dt,Dtmin)
            if trace         
                println("baztertu", " fase3=",fase3, ", dt=",Float32(dt) ,", dtj0=", Float32(dtj0)," ,dtj1=", Float32(dtj1),
                ", Dtmin=", Float32(Dtmin),", Dt=", Float32(Dt), ", plusDt=", plusDt)
                println()
            end
            else
                if trace
                println("fase3=", fase3, ",dt=",Float32(dt), " finkoarekin")
                println()
                end
            end


        end



  end   # while

#  println("amaitu da:", "j_iter=", j_iter,",iter=", iter,",diffU=", diffU, ",plusIt=", plusIt)
  
 
  if (iter && j_iter==maxiters) 
      println("Failure !!!. Maximum number of iterations.  Step=",tj+te, " dt=", dt)
      return("Failure",0)
  end


  @inbounds if (j_iter<maxiters && diffU)    #s=8 j_iter>22
    j_iter += 1

    for is in 1:s
        f(F[is], U[is], p,  tj  + dt*c[is])
        for k in indices
            L[is][k] = dt*(b[is]*F[is][k])
        end
   end
  end 


  @inbounds for k in indices    #Batura konpentsatuaren parekoa

    L_sum = L[1][k]
    for is in 2:s
        L_sum+=L[is][k]
    end
    res = Base.TwicePrecision(uj[k], ej[k]) + L_sum

    uj[k] = res.hi
    ej[k] = res.lo
  end

  res = Base.TwicePrecision(tj, te) + sdt*dt
  ttj[1] = res.hi
  ttj[2] = res.lo

  dts[1]=sdt*min(abs(dt),abs(tf-(ttj[1]+ttj[2])))
  dts[2]=dt

#    println("Urratsa: dts[1]=",dts[1], " dts[2]=",dts[2])
#     println("urratsa tj=",tj+te)
#     println("")

    return  ("Success",j_iter)


end





function IRKGLstep_fixed!(ttj,uj,ej,prob,dts,irkgl_cache::IRKGL_Cache{uType,tType2,fType,pType}) where {uType,tType2,fType,pType}

 f = irkgl_cache.odef
 p = irkgl_cache.p
 b = irkgl_cache.b
 c = irkgl_cache.c
 mu = irkgl_cache.mu
 nu = irkgl_cache.nu
 nu1 = irkgl_cache.nu1
 nu2 = irkgl_cache.nu2
 U = irkgl_cache.U
 Uz = irkgl_cache.Uz
 L = irkgl_cache.L
 Lz = irkgl_cache.Lz
 F = irkgl_cache.F
 Dmin = irkgl_cache.Dmin
 step_number = irkgl_cache.step_number[]
 first3iters = irkgl_cache.first3iters
 initial_interp = irkgl_cache.initial_interp[]
 len = irkgl_cache.length_u
 s = length(b)
 maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
 tj = ttj[1]
 te = ttj[2]
 indices=1:len
 flzero = zero(tType2)

 dt=dts[1]
 dtprev=dts[2]
 sdt=dts[3]
 tf=prob.tspan[2]


 if (initial_interp==1)
     for is in 1:s
         for k in indices
             dUik = muladd(nu[is,1], L[1][k], ej[k])
             for js in 2:s
                dUik = muladd(nu[is,js], L[js][k], dUik)
             end
             U[is][k] =  uj[k]  + dUik
         end
     end
 else
     for is in 1:s
         for k in indices
            U[is][k] = uj[k] + ej[k]
         end
     end
 end


j_iter = 0  # counter of fixed_point iterations
Dmin .= Inf

iter = true # Initialize iter outside the for loop
plusIt=true
diffU = true

for is in 1:s
  for k in indices
     Lz[is][k] = L[is][k]
  end
end

#println("Urratsa tj=",tj+te,", iter=", 1, ",dt=", Float32(dt))

@inbounds while (j_iter<maxiters && iter)   #&&  j_iter<22)
#@inbounds while (j_iter<44)  #  #s=8 j_iter<22

      iter = false
      j_iter += 1

      for is in 1:s
          f(F[is], U[is], p,  tj  + dt*c[is])
          for k in indices
              Uz[is][k] = U[is][k]
              L[is][k] = dt*(b[is]*F[is][k])
          end
      end

      for is in 1:s
          for k in indices
             dUik = muladd(mu[is,1], L[1][k], ej[k])
             for js in 2:s
                dUik = muladd(mu[is,js], L[js][k], dUik)
             end
             U[is][k] =  uj[k] + dUik
          end
      end

      diffU = false

      for k in indices
          Uk =[u[k] for u in U]
          Uzk =[u[k] for  u in Uz]
          DY = maximum(abs.(Uk-Uzk))

          if DY>0
              diffU = true
              if DY< Dmin[k]
                 Dmin[k]=DY
                 iter=true
              end
          end
      end

      if (!iter && diffU && plusIt)
          iter=true
          plusIt=false
      else
          plusIt=true
      end

#  println("j_iter=",j_iter,",iter=", iter, ",diffU=", diffU, ",plusIt=", plusIt, ", norm(Dmin)=",norm(Dmin))

  end # while

  if (!iter && (j_iter==maxiters))
     println("Failure !!!  Step=",tj+te, " dt=", dts[1])
     return("Failure",0)
 end

  @inbounds if (j_iter<maxiters && diffU)    #s=8 j_iter>22
          j_iter += 1

          for is in 1:s
              f(F[is], U[is], p,  tj  + dt*c[is])
              for k in indices
                  L[is][k] = dt*(b[is]*F[is][k])
              end
         end
  end


  @inbounds for k in indices    #Batura konpentsatuaren parekoa

      L_sum = L[1][k]
      for is in 2:s
          L_sum+=L[is][k]
      end
      res = Base.TwicePrecision(uj[k], ej[k]) + L_sum

      uj[k] = res.hi
      ej[k] = res.lo
   end


   res = Base.TwicePrecision(tj, te) + sdt*dt
   ttj[1] = res.hi
   ttj[2] = res.lo

   dts[1]=sdt*min(abs(dt),abs(tf-(ttj[1]+ttj[2])))
   dts[2]=dt

   return  ("Success",j_iter)


end
