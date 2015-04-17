      subroutine fresnel(sint,cost,sinp,cosp,nxp,nyp,nzp,tflag,
     + iseed,n1,n2,xp,yp,zp,xcur,ycur,Rr,ddr,weight,dda,noise,
     + rd_ang,xmax,ymax,zmax,zcur,rbins,Tr,sflag,td_ang,abins
     + ,cnt)
      
      implicit none

      include 'grid.txt'

      integer tflag,iseed,ir,rbins,sflag,abins,ia,cnt,xn,yn
      real*8 sint,cost,nxp,nyp,nzp,tir,cost2,f1,f2
      real*8 sinp,cosp,n1,n2,ran,xp,yp,xcur,ycur,ddr
      real*8 weight,theta2,sint2,theta1,Tr(0:rbins-1)
      real*8 zmax,ymax,xmax,zp,zcur,rp,Rr(0:rbins-1)
      real*8 dda,td_ang(0:rbins-1,0:abins),rando,angle
      real*8 rd_ang(0:rbins-1,0:abins),noise(1:cnt,1:cnt)
      real ran2
     


c**** Internal reflection and refraction at upper boundary

      if(sflag.lt.1.) then
            !noisy surface
            xn=int(nxg*(xp+xmax)/(2.*xmax))+1
            yn=int(nyg*(yp+ymax)/(2.*ymax))+1

            cost=noise(xn,yn)
            sint=sqrt(1.-cost*cost)
            
            if(sint.gt.n1/n2) then
            !photon reflects
            tflag=1
            weight = 0.
            sflag=1
            else
!            photon refracts            
            !fresnel coefficents
            cost=abs(cost)
                  !fresnel coefficients
                  if(n1.eq.n2)then!equal refractive indices
                  tir=0.
                  elseif(cost.gt.1.-1.0E-8)then!cost is straight down
                  tir=(n2-n1)**2/(n2+n1**2)
                  elseif(cost.lt.1.0E-6)then!oblique angle
                  tir=1.0
                  else
            cost2=sqrt(1.-(n2*sint/n1)**2)
            f1=abs((n2*cost-n1*cost2)/(n2*cost+n1*cost2))**2
            f2=abs((n2*cost2-n1*cost)/(n2*cost2+n1*cost))**2
            tir=(0.5*(f1+f2))
                  end if
            end if
      if(ran2(iseed).lt.tir)then
                 
            tflag=1
            weight = 0.
            sflag=1
      else
            sflag=1
            ! update position
            sint=(1./n2)*sint
            cost=(1.-sint*sint)
            if(cost.lt.1.)then
                  cost=0.
            else
                  cost=-sqrt(cost)
            end if
      end if
       

       
!       tir=(n2-n1)**2/(n2+n1**2)
!       weight=weight-tir
       sflag=1


      else
            !update postion
            xp=xcur-xmax
            yp=ycur-ymax
            zp=zcur-zmax

            
!            if(cost.lt.0.)sint=cost

            if(sint.gt.n2/n1) then
            !reflect photon

            cost=-cost
            sint=1.-cost*cost
            if(cost.lt.0.)then
                  sint=0.
            else
                  sint=sqrt(sint)
            end if
            else
                  cost=abs(cost)
                  angle=acos(cost)
                  if(zcur.gt.zmax)then
                              !noisy surface
                  xn=int(nxg*(xp+xmax)/(2.*xmax))+1
                  yn=int(nyg*(yp+ymax)/(2.*ymax))+1

                  cost=cos(angle+acos(noise(xn,yn)))
                  sint=sqrt(1.-cost*cost)
                  end if
            
                  !fresnel coefficients
                  if(n1.eq.n2)then!equal refractive indices
                  tir=0.
                  elseif(cost.gt.1.-1.0E-8)then!cost is straight down
                  tir=(n2-n1)**2/(n2+n1**2)
                  elseif(cost.lt.1.0E-6)then!oblique angle
                  tir=1.0
                  else
                  
                  cost2=sqrt(1.-(n1*sint/n2)**2)
                  f1=abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
                  f2=abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
                  tir=(0.5*(f1+f2))
                  end if
                  

                  ran=ran2(iseed)
                 
                        if(ran.gt.tir) then
                              if(zcur.lt.1.0E-8) then
                        ! bin bottom escape
                                    rp=sqrt(xp*xp+yp*yp)
                                    ir=floor(rp/ddr)
                                    if(ir.gt.rbins-1)ir=rbins-1
                                    Tr(ir)=Tr(ir)+weight
!                                    weight=0.
                                     sint=n1*sint
                                    ia=floor(asin(sint)/dda)
                                    if(ia.lt.0.)ia=abins
                                    if(ia.gt.abins-1)ia=abins
                                    td_ang(ir,ia)=td_ang(ir,ia)+weight 
                                    tflag=1
                              else

                                    ! bin top escape
                                     rp=sqrt(xp*xp+yp*yp)
                                     ir=floor(rp/ddr)
                                     if(ir.gt.rbins-1)ir=rbins-1
                                     Rr(ir)=Rr(ir)+weight
!                                     weight=0.
                                     sint=(1./n1)*sint
                                    ia=floor(asin(sint)/dda)
                                    if(ia.lt.0.)ia=abins
                                    if(ia.gt.abins-1)ia=abins
                                    rd_ang(ir,ia)=rd_ang(ir,ia)+weight
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
      !update direction vectors
      nxp=sint*cosp
      nyp=sint*sinp
      nzp=cost
      end if
             
      return
      end
