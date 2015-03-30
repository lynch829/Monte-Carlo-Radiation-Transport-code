      subroutine fresnel(sint,cost,sinp,cosp,nxp,nyp,nzp,tflag
     +            ,iseed,n1,n2,xp,yp,zp,xcur,ycur,Rr,ddr,weight
     +       ,xmax,ymax,zmax,zcur,rbins,Tr,sflag,rflag,noise,cnt)
      
      implicit none

      include 'grid.txt'

      integer tflag,iseed,ir,rbins,sflag,rflag,cnt,yn,xn
      real sint,cost,nxp,nyp,nzp,tir,theta1,theta2,sint2,cost2
      real sinp,cosp,n1,n2,ran,ran2,xp,yp,xcur,ycur,ddr,weight
      real f1,f2,noise(1:cnt,1:cnt)
      real zmax,ymax,xmax,zp,zcur,rp,Rr(0:rbins-1),Tr(0:rbins-1)



c**** Internal reflection and refraction at upper boundary

      if(sflag.eq.0) then
!            xn=int(nxg*(xp+xmax)/(2.*xmax))+1
!            yn=int(nyg*(yp+ymax)/(2.*ymax))+1
!            cost=noise(xn,yn)
!            sint=sqrt(1-cost*cost)
!            nxp=sint*cosp
!            nyp=sinp*sint
!            nzp=cost
      !lookup pos and get pertubed normal and calculate specular ref using fresnel
            cost2=sqrt(1-(((n1/n2)*sint)**2))
            f1=((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
            f2=((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
      tir=((f1+f2)*.5)/100.
      if(ran2(iseed).lt.tir)then
            tflag=1
            weight = 0.
            sflag=1
      else
            sflag=1
      end if

      else
      
            xp=xcur-xmax
            yp=ycur-ymax
            zp=zcur-zmax
            if(sint.gt.n2/n1) then
            !reflect photon
            cost=-cost
            sint=1.-cost*cost
            if(sint.lt.0.)then
                  sint=0.
            else
                  sint=sqrt(sint)
            end if
            
            else

!            xn=int(nxg*(xp+xmax)/(2.*xmax))+1
!            yn=int(nyg*(yp+ymax)/(2.*ymax))+1
!            cost=cos(acos(noise(xn,yn))-(3.141592654*.5))
!            sint=sqrt(1-cost*cost)
                  cost2=cos(asin((n2/n1)*sint))
                  f1=abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
                  f2=abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
                  tir=(0.5*(f1+f2))/100.


                  ran=ran2(iseed)
                  if(ran.gt.tir) then
                        if(zcur.gt..9999999*(2.*zmax)) then
                        ! top escape
                        rp=sqrt(xp*xp+yp*yp)
                        ir=floor(rp/ddr)
                        if(ir.gt.rbins-1)ir=rbins-1
                        Rr(ir)=Rr(ir)+weight
                        weight=0.
                        tflag=1
                        else

                        ! bottom escape
                        rp=sqrt(xp*xp+yp*yp)
                        ir=floor(rp/ddr)
                        if(ir.gt.rbins-1)ir=rbins-1
                        Tr(ir)=Tr(ir)+weight
                        weight=0.
                        tflag=1
                        end if
                        
                        !update position
                        sint=n2*sint
                        cost=(1.-sint*sint)
                        if(cost.lt.1.)then
                              cost=0.
                        else
                              cost=sqrt(cost)
                        end if
                        
                  else
                  ! reflect photon
                        cost=-cost
                        sint=1.-cost*cost
                        if(sint.lt.0.)then
                            sint=0.
                        else
                           sint=sqrt(sint)
                        end if
                        
                  end if
            end if
            nxp=sint*cosp
            nyp=sint*sinp
            nzp=cost
      end if
                  


























!      if (sint.gt.n1/n2) then
!         rflag=1   ! critical angle, photon fresnel reflected

!      else
!         sint2=n2*sint
!         theta1=asin(sint2)
!         theta2=asin(sint)
!         tir=0.5*(  (sin(theta1-theta2))**2/(sin(theta1+theta2))**2 + 
!     +              (tan(theta1-theta2))**2/(tan(theta1+theta2))**2  )
!!       cost2=sqrt(1-((n1/n2)*sint)**2)
!!           tir=.5*(abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
!!     +       + abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2)
!           
!           
!            ran=ran2(iseed)
!         if(ran.gt.tir ) then
!          xp=xcur-xmax
!              yp=ycur-ymax
!              zp=zcur-zmax
!              if(zcur.lt.1E-7) then
!              rp=sqrt(xp*xp+yp*yp)
!            ir=floor(rp/ddr)
!            if(ir.gt.rbins-1) ir=rbins-1
!            Tr(ir)=Tr(ir)+weight
!            weight=0.
!              else
!            rp=sqrt(xp*xp+yp*yp)
!            ir=floor(rp/ddr)
!            if(ir.gt.rbins-1) ir=rbins-1
!            Rr(ir)=Rr(ir)+weight
!            weight=0.
!            end if
!            sint=sint2                 !  photon refracted
!            cost=(1.-sint*sint)
!            if(cost.le.0.)then
!              cost=1.
!            else
!              cost=-sqrt(cost)
!            endif

!            nxp=sint*cosp  
!            nyp=sint*sinp
!            nzp=cost
!            
!         elseif(ran.le.tir) then
!            rflag=1                    !  fresnel reflection

!         endif
!      endif

                    
      return
      end
