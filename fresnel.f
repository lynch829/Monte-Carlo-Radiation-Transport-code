      subroutine fresnel(sflag,sint,cost,n1,n2,zcur,
     & ddr,ycur,weight,xmax,ymax,zmax,xp,yp,zp,xcur,
     & Rr,rbins,Tr,twopi,tflag,sinp,cosp)
      
      implicit none
      

      integer tflag,sflag,ir,rbins
      real cost2,f1,f2,rf,crit,n2,n1,transmitted,weight,specref
      real sint,cost,xp,yp,zp,xmax,ymax,zmax,twopi,sinp,cosp
      real sint2,theta1,theta2
      real zcur,xcur,ycur,ddr,Rr(0:rbins-1),rt,Tr(0:rbins-1)

      !need to adjust code so that can take light from outside
      if(sflag.eq.0)then
      rf=abs((n2-n1)**2/(n2+n1)**2)
      specref=rf
      weight=1.-rf
      sflag=1
      else
      
c***** crtitcal angle      
            crit=n2/n1
            ! correct angle for bottom side
       if(zcur.lt.1E-7) then
       cost=90.-(180.*acos(cost)/3.14-90.)
       cost=cos(3.14*cost/180.)
       end if
c*****      checks if angle is less than critical angle if so then TIR 
            if(cost.gt.crit)then

                  !tir rf=1.
                  !update direction properly
                 
              else
                  !no TIR
!                  sint2=n1*sint
!         theta1=asin(sint)
!         theta2=asin(sint2)
!          rf=0.5*(  (sin(theta1-theta2))**2/(sin(theta1+theta2))**2 + 
!     +              (tan(theta1-theta2))**2/(tan(theta1+theta2))**2  )
!     
                  cost2= sqrt(1-((n2/n1)*sint)**2)
            f1=(abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2)))**2
            f2=(abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost)))**2
!                   reflection coefficent
                  rf=.5*(f1+f2)
!                  print *, 'rf', f1,f2,cost2
                  !diffuse reflectance
                  if(zcur.gt..9999999*(2*zmax))then

                  !adjust weight of photon as some of photon transmitted and bin.
                  
                        transmitted=(1-rf)*weight
                        weight=weight*rf
                        
                        rt=sqrt(xp**2+yp**2)
                        ir=floor(rt/(ddr))
!                        if (ir.eq.0) print *, transmitted,Rr(0)
                        if (ir.gt.rbins-1) ir=rbins-1
                        Rr(ir)=Rr(ir) + transmitted
                        
                  ! transmission on bottom side      
                  else
                  !adjust weight of photon as some of photon transmitted and bin.
                  
                        transmitted=(1.-rf)*weight
                        weight=rf*weight
!                        print *, transmitted,rf
                        rt=sqrt(xp**2+yp**2)
                        ir=floor(rt/(ddr)-1)
                        if(ir.lt.0) ir=0
                        if (ir.gt.rbins-1) ir=rbins-1
                        Tr(ir)=Tr(ir) + transmitted
                  end if
            endif
      endif
      return
      end
