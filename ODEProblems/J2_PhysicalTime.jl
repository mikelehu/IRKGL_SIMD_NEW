
function Energy(u,p)
    μ = p[1]
    R = p[2]
    ϵ = p[3]
    q = u[1:3]
    v = u[4:6]
    r = norm(q)
    v2 = dot(v,v)
    sinth2 = (u[3]/r)^2
    aux = (R/r)^2
    return 0.5*v2 - μ/r * (1 + 0.5*ϵ*aux * (1 - 3*sinth2))
end


function J2ODE!(du,u,p,t)
    μ = p[1]
    R = p[2]
    ϵ = p[3]
    x = u[1]
    y = u[2]
    z = u[3]
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    r2 = x^2+y^2+z^2
    r = sqrt(r2)
    aux1 = -μ/(r*r2)
    F = 1.5 - 7.5*(z/r)^2
    G =  3 + F
    aux2 = ϵ*R^2/r2
    aux3 = aux1*(1+aux2*F)
    aux4 = aux1*(1+aux2*G)
    du[4]=aux3*x
    du[5]=aux3*y
    du[6]=aux4*z
    return nothing
end


# mikel 2023-04-11

function J2ODE!(du,u,p,t,part)
    
    if part==1    # Evaluate dq/dt   
       du[1] = u[4]
       du[2] = u[5]
       du[3] = u[6]
    else         # Evaluate dv/dt
      μ = p[1]
      R = p[2]
      ϵ = p[3]
      x = u[1]
      y = u[2]
      z = u[3]
    	r2 = x^2+y^2+z^2
    	r = sqrt(r2)
    	aux1 = -μ/(r*r2)
    	F = 1.5 - 7.5*(z/r)^2
    	G =  3 + F
    	aux2 = ϵ*R^2/r2
    	aux3 = aux1*(1+aux2*F)
    	aux4 = aux1*(1+aux2*G)
    	du[4]=aux3*x
    	du[5]=aux3*y
    	du[6]=aux4*z
    end
    return nothing
end


#
# functions for SecondOrderODEProblem 
#

function J2ODE2nd!(ddu,du,u,p,t)

     μ = p[1]
     R = p[2]
     ϵ = p[3]
     x = u[1]
     y = u[2]
     z = u[3]
     r2 = x^2+y^2+z^2
     r = sqrt(r2)
     aux1 = -μ/(r*r2)
     F = 1.5 - 7.5*(z/r)^2
     G =  3 + F
     aux2 = ϵ*R^2/r2
     aux3 = aux1*(1+aux2*F)
     aux4 = aux1*(1+aux2*G)
     ddu[1]=aux3*x
     ddu[2]=aux3*y
     ddu[3]=aux4*z
 
    return nothing

end



#
# functions for DynamicalODEProblem 
#


function J2ODEv!(dv,q,v,p,t)

     μ = p[1]
     R = p[2]
       ϵ = p[3]
     x = u[1]
     y = u[2]
     z = u[3]
     r2 = x^2+y^2+z^2
     r = sqrt(r2)
     aux1 = -μ/(r*r2)
     F = 1.5 - 7.5*(z/r)^2
     G =  3 + F
     aux2 = ϵ*(R/r)^2
     aux3 = aux1*(1+aux2*F)
     aux4 = aux1*(1+aux2*G)
     dv[1]=aux3*x
     dv[2]=aux3*y
     dv[3]=aux4*z
 
    return nothing

end


function J2ODEq!(dq,q,v,p,t)

     dq[1] = v[1]
     dq[2] = v[2]
     dq[3] = v[3]
     
     return nothing

end



