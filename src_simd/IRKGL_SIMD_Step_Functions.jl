#   
#  IRKGL Step functions

#      IRKGLstep_SIMD_fixed!    



function IRKGLstep_SIMD_fixed!(ttj,tj1, uj,ej,prob,dts,irkglcache::IRKGL_SIMD_Cache{floatType,fType,pType,s_,dim,dim_}) where {floatType,fType,pType,s_,dim,dim_}

    f = irkglcache.odef
    p = irkglcache.p
    b = irkglcache.b
    c = irkglcache.c
    mu = irkglcache.mu
    nu = irkglcache.nu
    U = irkglcache.U
    U_ = irkglcache.U_
    L = irkglcache.L
    dU = irkglcache.dU
    F = irkglcache.F
    Dmin = irkglcache.Dmin
    step_number = irkglcache.step_number[]
    initial_interp = irkglcache.initial_interp[]
    len = irkglcache.length_u
    s = length(b)
    maxiters = (step_number==1 ? 10+irkglcache.maxiters : irkglcache.maxiters )
    tj = ttj[1]
    te = ttj[2]
    indices=1:len


    dt=dts[1]
    dtprev=dts[2]
    signdt=dts[3]
    sdt=dt*signdt


    if (initial_interp==1)

      for k in indices
         Lk = getindex_(L,k)
         dUk = muladd(nu[1], Lk[1], ej[k])
         for js in 2:s
             dUk = muladd(nu[js], Lk[js], dUk)
         end
         setindex_!(U, uj[k]+dUk, k)
      end

    else

       for k in indices
         uej = uj[k] + ej[k]
         setindex_!(U, uej, k)
       end

    end


    j_iter = 0  # counter of fixed_point iterations
    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt=true
    diffU = true


    @inbounds while (j_iter<maxiters && iter ) 

         iter = false
         j_iter += 1

         U_.data .= U.data
         f(F, U, p, tj + sdt*c)

         diffU = false

         for k in indices

             Fk = getindex_(F,k)
             Lk = sdt*(b*Fk)
             dUk = muladd(mu[1], Lk[1], ej[k])
             for is in 2:s
                     dUk = muladd(mu[is], Lk[is], dUk)
             end
             Uk = uj[k]+dUk
             setindex_!(U, Uk, k)
             setindex_!(L, Lk, k)
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

         if (!iter && diffU && plusIt)
             iter=true
             plusIt=false
         else
             plusIt=true
         end

     end # while

     if (!iter && (j_iter==maxiters))
        println("Failure !!!  Step=",tj+te, " dt=", dts[1])
        return("Failure",0)
     end

     @inbounds if (j_iter<maxiters && diffU) #  && j_iter>22)   #s=8 j_iter>22
             j_iter += 1
             f(F, U, p, tj + sdt*c)

             for k in indices
                 Fk = getindex_(F, k)
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
