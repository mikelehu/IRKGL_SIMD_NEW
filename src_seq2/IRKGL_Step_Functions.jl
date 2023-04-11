#
#  IRKGLstep_fixed!


function IRKGLstep_fixed!(ttj,tj1, uj,ej,dts,irkgl_cache::IRKGL_Cache{uType,tType2,fType,pType}) where {uType,tType2,fType,pType}

    trace = false

    
     f = irkgl_cache.odef
     p = irkgl_cache.p
     b = irkgl_cache.b
     c = irkgl_cache.c
     mu = irkgl_cache.mu
     nu = irkgl_cache.nu
     U = irkgl_cache.U
     U_ = irkgl_cache.U_
     L = irkgl_cache.L
     dU = irkgl_cache.dU
     F = irkgl_cache.F
     Dmin = irkgl_cache.Dmin
     step_number = irkgl_cache.step_number[]
     initial_interp = irkgl_cache.initial_interp[]
     len = irkgl_cache.length_u
     s = length(b)
     maxiters = (step_number==1 ? 10+irkgl_cache.maxiters : irkgl_cache.maxiters )
     tj = ttj[1]
     te = ttj[2]
     indices=1:len
    
     dt=dts[1]
     dtprev=dts[2]
     signdt=dts[3]
     sdt=dt*signdt
    

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
    
    
    trace ? println("Urratsa tj=",tj+te,", iter=", 1, ",dt=", Float32(dt)) : nothing
    trace ? println("Norm L=", norm(L), ", Norm U=", norm(U), ", Norm ej=", norm(ej), ", Norm p =", p) : nothing
    trace ? println("tj=", tj, ",dts=", dts, "sdt=", sdt) : nothing
    
    @inbounds while (j_iter<maxiters && iter)  
    
          iter = false
          j_iter += 1
    
          for is in 1:s
              f(F[is], U[is], p,  tj  + sdt*c[is])
              for k in indices
                  U_[is][k] = U[is][k]
                  L[is][k] = sdt*(b[is]*F[is][k])
              end
          end

          trace ? println("Norm(F)=", norm(F)) : nothing
    
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
              DY = abs(U[1][k]-U_[1][k])
              for is in 2:s 
                  DY=max(abs(U[is][k]-U_[is][k]),DY)
              end 
    
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
    
          if trace == true
            Dminaux=replace(Dmin, Inf=>0)
            println("iter=", iter, ",diffU=", diffU, ",plusIt=", plusIt, ", norm(Dmin)=",norm(Dminaux)) 
          end
    
      end # while
    
      if (!iter && (j_iter==maxiters))
         println("Failure !!!  Step=",tj+te, " dt=", dts[1])
         return("Failure",0)
     end
    
      @inbounds if (j_iter<maxiters && diffU)    #s=8 j_iter>22
              j_iter += 1
    
              for is in 1:s
                  f(F[is], U[is], p,  tj  + sdt*c[is])
                  for k in indices
                      L[is][k] = sdt*(b[is]*F[is][k])
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
    
    
       res = Base.TwicePrecision(tj, te) + sdt
       ttj[1] = res.hi
       ttj[2] = res.lo
    
       dts[1]=min(abs(dt),abs(tj1-(ttj[1]+ttj[2])))
       dts[2]=dt
    
       return  ("Success",j_iter)
    
    
    end
    
