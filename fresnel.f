      subroutine fresnel(sint,cost,sinp,cosp,nxp,nyp,nzp,tflag,
     + Amat,ybins,iseed,n1,n2,xp,yp,zp,xcur,ycur,Rr,ddr,weight,
     + xbins,xmax,ymax,zmax,zcur,rbins,Tr,sflag,rflag,noise,cnt
     + ,ddx,ddy,Amat1,o,Amat2)
      
      implicit none

      include 'grid.txt'

      integer tflag,iseed,ir,rbins,sflag,rflag,cnt,yn,xn,xbins
      integer ix,iy,ybins,o
      real*8 sint,cost,nxp,nyp,nzp,tir,theta1,theta2,sint2,cost2
      real*8 sinp,cosp,n1,n2,ran,xp,yp,xcur,ycur,ddr,weight
      real*8 f1,f2,noise(1:cnt,1:cnt),ddx,ddy
      real*8 zmax,ymax,xmax,zp,zcur,rp,Rr(0:rbins-1),Tr(0:rbins-1)
      real ran2
      real*8 Amat(0:xbins,0:ybins),Amat1(0:xbins,0:ybins)
      real*8 Amat2(0:xbins,0:ybins)


c**** Internal reflection and refraction at upper boundary

      if(sflag.lt.1) then
            xn=int(nxg*(xp+xmax)/(2.*xmax))+1
            yn=int(nyg*(yp+ymax)/(2.*ymax))+1
            cost=-1.+noise(xn,yn)
            sint=sqrt(1-cost*cost)
      !lookup pos and get pertubed normal and calculate specular ref using fresnel
            cost2=sqrt(1-(((n2/n1)*sint)**2))
            f1=((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
            f2=((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
      tir=((f1+f2)*.5)/100.
            nxp=sint*cosp
            nyp=sinp*sint
            nzp=cost
      if(ran2(iseed).lt.tir)then

!                 ix=floor(xcur/ddx)
!                 iy=floor(ycur/ddy)
!                 if(ix.gt.xbins) ix=xbins 
!                 if(iy.gt.ybins) iy=ybins
!                 if(ix.lt.0) ix=0 
!                 if(iy.lt.0) iy=0
!                 Amat1(ix,iy)=Amat1(ix,iy)+1.
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
!            cost=1.+noise(xn,yn)
!            sint=sqrt(1-cost*cost)
!            print *, cost,sint

                  cost2=cos(asin((n2/n1)*sint))
                  f1=abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
                  f2=abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
                  tir=(0.5*(f1+f2))/100.


                  ran=ran2(iseed)
                  if(ran.gt.tir) then
                        if(zcur.gt..9999999*(2.*zmax)) then
                        ! top escape
!                        rp=sqrt(xp*xp+yp*yp)
!                        ir=floor(rp/ddr)
!                        if(ir.gt.rbins-1)ir=rbins-1
!                        Rr(ir)=Rr(ir)+weight



                 ix=floor(xcur/ddx)
                 iy=floor(ycur/ddy)
                 if(ix.gt.xbins) ix=xbins 
                 if(iy.gt.ybins) iy=ybins
                 if(ix.lt.0) ix=0 
                 if(iy.lt.0) iy=0
                 if(o.eq.1)then
                  Amat(ix,iy)=Amat(ix,iy)+1.
                 elseif(o.eq.2) then
                  Amat1(ix,iy)=Amat1(ix,iy)+1.
                 elseif(o.eq.3) then
                  Amat2(ix,iy)=Amat2(ix,iy)+1.
                 end if
                        weight=0.
                        tflag=1
                        else

                        ! bottom escape
!                        rp=sqrt(xp*xp+yp*yp)
!                        ir=floor(rp/ddr)
!                        if(ir.gt.rbins-1)ir=rbins-1
!                        Tr(ir)=Tr(ir)+weight
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
                 
      return
      end
