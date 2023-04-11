#   
#  IRKNGL Step functions

#      IRKNGLstep_SIMD_fixed! 


 function IRKNGLstep_SIMD_fixed!(ttj,tj1, uj,ej,dts, irknglcache::IRKNGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}) where {floatType,fType,pType,s_,dim,dim_}

    f = irknglcache.odef
    p = irknglcache.p
    b = irknglcache.b
    c = irknglcache.c
    mu = irknglcache.mu
    nu = irknglcache.nu
    U = irknglcache.U
    U_ = irknglcache.U_
    L = irknglcache.L
    dU = irknglcache.dU
    F = irknglcache.F
    Dmin = irknglcache.Dmin
    step_number = irknglcache.step_number[]
    initial_interp = irknglcache.initial_interp[]
    len = irknglcache.length_u
    lenq = irknglcache.length_q
    s = length(b)
    len = length(uj)
    maxiters = (step_number==1 ? 10+irknglcache.maxiters : irknglcache.maxiters )
    tj = ttj[1]
    te = ttj[2]


    indices=1:len
    indices1 = 1:lenq
    indices2 = (lenq+1):len

    dt=dts[1]
    dtprev=dts[2]
    signdt=dts[3]
    sdt=dt*signdt


    if (initial_interp==1)

      for k in indices2
         Lk = getindex_(L,k)
         dUk = muladd(nu[1], Lk[1], ej[k])
         for is in 2:s
             dUk = muladd(nu[is], Lk[is], dUk)
         end
         setindex_!(U, uj[k]+dUk, k)
      end

    else
    
       for k in indices2
         uej = uj[k] + ej[k]
         setindex_!(U, uej, k)
       end
       
    end

    f(F, U, p, tj + sdt*c, 1)
    
    for k in indices1
        Fk = getindex_(F,k)
        Lk =sdt*(b*Fk)
        setindex_!(L, Lk, k)
        dUk = muladd(mu[1], Lk[1], ej[k])
        for is in 2:s
            dUk = muladd(mu[is], Lk[is], dUk)
        end
        setindex_!(U, uj[k]+dUk, k)
    end

    j_iter = 0  # counter of fixed_point iterations
    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true

    @inbounds while (j_iter<maxiters && iter)

         iter = false
         j_iter += 1

         U_.data .= U.data

         f(F, U, p, tj + sdt*c, 2)

         for k in indices2
                 Fk = getindex_(F,k)
                 Lk =sdt*(b*Fk)
                 setindex_!(L, Lk, k)
                 dUk = muladd(mu[1], Lk[1], ej[k])
                 for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
                 end
                 setindex_!(U, uj[k] + dUk, k)
         end

         f(F, U, p, tj + sdt*c, 1)

         for k in indices1
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 setindex_!(L, Lk, k)
                 dUk = muladd(mu[1], Lk[1], ej[k])
                 for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
                 end
                 setindex_!(U, uj[k]+dUk, k)
         end


        diffU = false

        for k in indices1   # Hemen indices1 jarri liteke, q'=v, v'=f(q,t) moduko ED-a dela suposatuz

             Uk = getindex_(U,k)
             Uk_ = getindex_(U_,k)
             DY = maximum(abs(Uk-Uk_))

             if DY>0
                 diffU = true
                 if DY< Dmin[k]
                    Dmin[k]=DY
                    iter=true
                 end
             end
         end

         if (!iter && diffU && plusIt)  #
             iter=true
             plusIt=false
         else
             plusIt=true
         end

    end # while

    if (!iter && (j_iter==maxiters))
        println("Fail !!!  Step=",tj+te, " dt=", dts[1])
        return("Failure",0)
    end

    @inbounds if (j_iter<maxiters && diffU)

         j_iter += 1

         f(F, U, p, tj + sdt*c, 2)

         for k in indices2
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 dUk = muladd(mu[1], Lk[1], ej[k])
                 for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
                 end
                 setindex_!(U, uj[k]+dUk, k)
                 setindex_!(L, Lk, k)
         end

         f(F, U, p, tj + sdt*c, 1)

         for k in indices1
                 Fk = getindex_(F,k)
                 Lk = sdt*(b*Fk)
                 setindex_!(L, Lk, k)
         end

    end


     @inbounds for k in indices    #Batura konpentsatuaren parekoa
         Lk = getindex_(L,k)
         L_sum = sum(Lk)
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


