      subroutine fresnel(sint,cost,sinp,cosp,nxp,nyp,nzp,rflag
     +            ,iseed,n1,n2,xp,yp,zp,xcur,ycur,Rr,ddr,weight
     +            ,xmax,ymax,zmax,zcur,rbins,Tr)
      
      implicit none

      integer rflag,iseed,ir,rbins
      real sint,cost,nxp,nyp,nzp,tir,theta1,theta2,sint2
      real sinp,cosp,n1,n2,ran,ran2,xp,yp,xcur,ycur,ddr,weight
      real zmax,ymax,xmax,zp,zcur,rp,Rr(0:rbins-1),Tr(0:rbins-1)


      rflag=0

c**** Internal reflection and refraction at upper boundary
      if (sint.gt.n1/n2) then
         rflag=1   ! critical angle, photon fresnel reflected

      else
         sint2=n2*sint
         theta1=asin(sint)
         theta2=asin(sint2)
         tir=0.5*(  (sin(theta1-theta2))**2/(sin(theta1+theta2))**2 + 
     +              (tan(theta1-theta2))**2/(tan(theta1+theta2))**2  )
            ran=ran2(iseed)
         if(ran.gt.tir ) then
          xp=xcur-xmax
              yp=ycur-ymax
              zp=zcur-zmax
              if(zcur.lt.1E-7) then
              rp=sqrt(xp*xp+yp*yp)
            ir=floor(rp/ddr)
            if(ir.gt.rbins-1) ir=rbins-1
            Tr(ir)=Tr(ir)+weight
            weight=0.
              else
            rp=sqrt(xp*xp+yp*yp)
            ir=floor(rp/ddr)
            if(ir.gt.rbins-1) ir=rbins-1
            Rr(ir)=Rr(ir)+weight
            weight=0.
            end if
            sint=sint2                 !  photon refracted
            cost=(1.-sint*sint)
            if(cost.le.0.)then
              cost=1.
            else
              cost=-sqrt(cost)
            endif

            nxp=sint*cosp  
            nyp=sint*sinp
            nzp=cost
            
         elseif(ran.le.tir) then
            rflag=1                    !  fresnel reflection

         endif
      endif

                    
      return
      end
