



function NbodyEnergy3REL(u,Gm)

     N = length(Gm)

     norm12=sqrt(u[1,1,1]^2+u[2,1,1]^2)
     norm23=sqrt(u[1,2,1]^2+u[2,2,1]^2)
     norm31=sqrt(u[1,3,1]^2+u[2,3,1]^2)

     T=Gm[1]*(u[1,1,2]^2+u[2,1,2]^2)+Gm[2]*(u[1,2,2]^2+u[2,2,2]^2)+Gm[3]*(u[1,3,2]^2+u[2,3,2]^2)
     U=Gm[1]*Gm[2]/norm12+Gm[2]*Gm[3]/norm23+Gm[3]*Gm[1]/norm31

     return T/2-U


end


function NbodyODE3REL!(du,u,Gm,t)
     
#     N = length(Gm)

#     q12 = u[:,1,1]
#     q23 = u[:,2,1]
#     q31 = u[:,3,1]

#     aux12 = q12*dot(q12, q12)^-1.5
#     aux23 = q23*dot(q23, q23)^-1.5
#     aux31 = q31*dot(q31, q31)^-1.5

#     aux12 = q12
#     aux23 = q23
#     aux31 = q31

#     dv1 = Gm[2]*aux12 - Gm[3]*aux31
#     dv2 = Gm[3]*aux23 - Gm[1]*aux12 
#     dv3 = Gm[1]*aux31 - Gm[2]*aux23  

#     @. du[:,1,1]=u[:,2,2]-u[:,1,2]
#     @. du[:,2,1]=u[:,3,2]-u[:,2,2]
#     @. du[:,3,1]=u[:,1,2]-u[:,3,2]

#     du[:,1,2].= dv1
#     du[:,2,2].= dv2
#     du[:,3,2].= dv3 

    return nothing

end




function NbodyODE3REL!(du,u,Gm,t,part)
    
   
     N = length(Gm)

    if part==1    # Evaluate dq/dt

     @. du[:,1,1]=u[:,2,2]-u[:,1,2]
     @. du[:,2,1]=u[:,3,2]-u[:,2,2]
     @. du[:,3,1]=u[:,1,2]-u[:,3,2]

    else         # Evaluate dv/dt
         
     q12 = u[:,1,1]
     q23 = u[:,2,1]
     q31 = u[:,3,1]

     aux12 = q12*dot(q12, q12)^-1.5
     aux23 = q23*dot(q23, q23)^-1.5
     aux31 = q31*dot(q31, q31)^-1.5

     dv1 = Gm[2]*aux12 - Gm[3]*aux31
     dv2 = Gm[3]*aux23 - Gm[1]*aux12 
     dv3 = Gm[1]*aux31 - Gm[2]*aux23   

     du[:,1,2].= dv1
     du[:,2,2].= dv2
     du[:,3,2].= dv3 
      
    end # if

    return nothing

end






