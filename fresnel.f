      subroutine fresnel(cost,n1,n2,sflag,weight,Rr,xp,yp,zcur,
     +                   sint,zmax,rbins,Tr,ddr)
      
      implicit none
      
      integer sflag,ir,rbins
      real cost,n1,n2,crit,rf1,rf2,rft,weight,transmitted
      real r,xp,yp,ddr,Rr(0:rbins-1),zcur,cost2,sint,zmax
      real Tr(0:rbins-1),specref,costt
      
      
      if(sflag.lt.1) then
      
!     do code for surface reflection
      rft=((n2-n1)/(n1+n2))**2
      specref=(1.-rft)*rft
      weight=weight*rft
      sflag=1
      else
            crit = n2/n1
            if(zcur.lt.1.0E-7) then
!            print *, 'bot',cost
            costt=cost
            cost=abs(cost)
            end if
            if(sint.gt.crit) then
            
!                 total internal reflection

            else
            
                  cost2=sqrt(1-(n1/n2*sint)**2)
                  rf1=abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
                  rf2=abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
                  rft=0.5*(rf1+rf2)
                  
                  if(zcur.gt..99999*2.*zmax) then
                  
!                  top reflection
!                        print *, cost
                        transmitted=(1.-rft)*weight
                        weight=weight*rft
                        r=sqrt(xp*xp+yp*yp)
                        ir=floor(r/ddr)
                        Rr(ir)=Rr(ir) + transmitted
                        cost=-cost
                  
                  else
                  
!                 bottom reflection                                          
                  transmitted=(1.-rft)*weight
                  weight=weight*rft
                  r=sqrt(xp*xp+yp*yp)
                  ir=floor(r/ddr)
                  Tr(ir)=Tr(ir) + transmitted
                  cost=-costt
                  
                  

                  end if
                  
            end if
            
      end if
                    
      return
      end subroutine fresnel
