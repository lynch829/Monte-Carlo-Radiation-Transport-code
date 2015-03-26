      program mcpolar

      implicit none

      include 'grid.txt'
      include 'photon.txt'

c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,tflag
      real*8 nscatt
      real kappa,albedo,hgg,pl,pc,sc,xmax,ymax,zmax
      real pi,twopi,fourpi,g2,delta

      real ran2

c**** Read in parameters from the file input.params
      open(10,file='input.params',status='old')
          read(10,*) nphotons
          read(10,*) iseed
          read(10,*) kappa
          read(10,*) albedo
          read(10,*) hgg
          read(10,*) pl
          read(10,*) pc
          read(10,*) sc
          read(10,*) xmax
          read(10,*) ymax
          read(10,*) zmax
          close(10)

c***** Set up constants, pi and 2*pi  ********************************
      pi=4.*atan(1.)
      twopi=2.*pi
      fourpi=4.*pi

      iseed=-abs(iseed)  ! Random number seed must be negative for ran2

      g2=hgg*hgg  ! Henyey-Greenstein parameter, hgg^2

c**** Initialize arrays to zero *************************************
      call iarray(xface,yface,zface,rhokap)

c***** Set up density grid *******************************************
      call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,kappa)

c***** Set small distance for use in optical depth integration routines 
c***** for roundoff effects when crossing cell walls
      delta=1.e-3*(2.*xmax/nxg)


c**** Loop over nph photons from each source *************************
        nscatt=0
        do j=1,nphotons

          if(mod(j,10000).eq.0)then
             print *, j,' scattered photons completed'
          end if

c***** Release photon from point source *******************************
          call sourceph(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  fi,fq,fu,fv,xmax,ymax,zmax,twopi,
     +                  xcell,ycell,zcell,nxg,nyg,nzg,iseed)

c****** Find scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,tflag,iseed,delta)
     
c******** Photon scatters in grid until it exits (tflag=1) 
          dowhile(tflag.eq.0)
                 if( (ran2(iseed).lt.albedo) ) then
c************ Scatter photon into new direction and update Stokes parameters
                   call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi,iseed)
                   nscatt=nscatt+1
                else
                   goto 100
                endif

c************ Find next scattering location
              call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,tflag,iseed,delta)
          end do

100      continue

        end do      ! end loop over nph photons

      print*,'Avereage number of scatterings = ',sngl(nscatt/nphotons)

      stop
      end
