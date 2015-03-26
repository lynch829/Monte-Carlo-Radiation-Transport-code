      subroutine sourceph(xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,fi,fq,fu,fv,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed)

      implicit none

      include 'photon.txt'

      integer xcell,ycell,zcell,nxg,nyg,nzg,iseed
      real xmax,ymax,zmax,twopi
      real ran2

c***** emit photon isotropically from origin
      xp=0.
      yp=0.
      zp=0.

      cost=2.*ran2(iseed)-1.
      sint=(1.-cost*cost)
      if(sint.le.0.)then
        sint=0.
      else
        sint=sqrt(sint)
      endif

      phi=twopi*ran2(iseed)
      cosp=cos(phi)
      sinp=sin(phi)

c***** Set photon direction cosines for direction of travel *********
      nxp=sint*cosp  
      nyp=sint*sinp
      nzp=cost

c***** Set Stokes fluxes ********************************************
      fi=1.
      fq=0.
      fu=0.
      fv=0.

c*************** Linear Grid *************************
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
c*****************************************************

      return
      end

