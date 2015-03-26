      subroutine sourceph(xp,yp,zp,nxp,nyp,nzp,sint,
     +          cost,sinp,cosp,phi,xmax,ymax,zmax,twopi,
     +            xcell,ycell,zcell,nxg,nyg,nzg,iseed)

      implicit none

      include 'photon.txt'

      integer xcell,ycell,zcell,nxg,nyg,nzg,iseed,i,cnt,j,nlow
      real xmax,ymax,zmax,twopi,w,lam,phigauss,r1,flu
      real ran
      real ran2


!      zp=zmax-1E-7
!      w=0.2
c***** emit photon from a circle on the surface of skin
c      xp=xmax+1.
c      yp=ymax+1.
c      zp=zmax-(1E-3)
c      w=0.2             !radius of light illumination
c      do while((xp**2+yp**2).gt.w**2)
c            xp=(2*ran2(iseed)-1.)*w
c            yp=(2*ran2(iseed)-1.)*w
c      end do

c      phi=0.
c      cosp=1.
c      sinp=0.
c      cost=-1.
     
c**** emit photon from a gaussian beam on surface of skin

!      r1=w*sqrt(-log(1-ran2(iseed)))
!      phigauss=twopi*ran2(iseed) 
!      xp=r1*cos(phigauss)
!      yp=r1*sin(phigauss)    

c**** emit photons in a pencil beam 

      xp=0.
      yp=0.
      zp=zmax-0.001

c***** Set photon direction cosines for direction of travel(into skin) *********

      phi=0.
      cosp=1.
      sinp=0.
      cost=-1.

      nxp=0.  
      nyp=0.
      nzp=-1.
      

c*************** Linear Grid *************************
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
c*****************************************************

      return
      end

