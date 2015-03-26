      program mcpolar

      implicit none

      include 'grid.txt'
      include 'photon.txt'

c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,tflag,i,ph1,ph2,ph3
      integer cnt,q,qz,numberrun,sflag,pflag
      real*8 nscatt
      real kappa,albedo,hgg,xmax,ymax,zmax,zpos,th,d
      real pi,twopi,fourpi,g2,delta,ydect,xdect,rdect,dectang
      real chance,terminate,weight,absorb,n1,n2,mua,mus,xcur,ycur,zcur
      character(len=70) :: fn
      integer ix,iy,iz,xbins,ybins,zbins,rbins,ir
      real ddz,ddx,ddy,ddr,r
      ! bin array for cartesian bins
c      real, allocatable :: Amat(:,:,:)
      ! bin arrays for cylindrical bins
      real, allocatable :: Amat(:,:),Az(:),Rr(:),Tr(:)
      real ran2,start,finish
      
c**** Read in parameters from the file input.params
      open(10,file='input.params',status='old')
          read(10,*) nphotons
          read(10,*) iseed
          read(10,*) mua
          read(10,*) mus
          read(10,*) hgg
          read(10,*) n1
          read(10,*) n2
          read(10,*) xmax
          read(10,*) ymax
          read(10,*) zmax
          read(10,*) ydect
          read(10,*) xdect
          read(10,*) rdect
          read(10,*) zpos
          read(10,*) dectang
          close(10)
          
          kappa=mua+mus
          albedo=mus/kappa

c***** Set up bins ***************************************************
  
!      abins=100
!      dda=int((pi/(2.*abins)))


!      **** Cartesian Bins *****

       ! set # of bins for x,y, and z directions
!      ybins=500
!      xbins=500
!      zbins=500
       ! set size of bins
!      ddy=((2*ymax)/ybins)
!      ddx=((2*xmax)/xbins)           
!      ddz=((2*zmax)/zbins) 
       ! allocate bin array based on # of bins
!      allocate(Amat(-xbins/2:xbins/2,-ybins/2:ybins/2,-zbins/2:zbins/2)) 
    
!      **** Cylindrical bins ****

      ! set # of bins for radial and vert direction
      zbins=500
      rbins=500
      ! set size of bins
      ddz=((2.)/zbins)
      ddr=((2.)/rbins) 
      ! allocate bin array based on # of bins     
      allocate(Amat(0:rbins-1,0:zbins))
      allocate(Az(0:zbins))
      allocate(Rr(0:rbins-1))
      allocate(Tr(0:rbins-1))
      
      
      ! clear arrays
      Amat=0.
      Az=0.
      Rr=0.
      Tr=0.
      
c***** Set up constants, pi and 2*pi  ********************************
      pi=4.*atan(1.)
      twopi=2.*pi
      fourpi=4.*pi
      ph1 = 0
      ph2 = 0
      ph3 = 0
      ! set terminate and chance values for roulette sub-routine
      terminate=0.00001
      chance=0.1
      ! set number of photons run to 0
      numberrun=0
      iseed=-abs(iseed)  ! Random number seed must be negative for ran2
      ! set detector angle into radians
      dectang=(dectang*pi)/180.
      g2=hgg*hgg  ! Henyey-Greenstein parameter, hgg^2

c**** Initialize arrays to zero *************************************
      call iarray(xface,yface,zface,rhokap,jmean)

c***** Set up density grid *******************************************
      call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,kappa)

c***** Set small distance for use in optical depth integration routines 
c***** for roundoff effects when crossing cell walls
      delta=1.e-5*(2.*xmax/nxg)


c**** Loop over nph photons from each source *************************
        nscatt=0
        call cpu_time(start)
        do while(numberrun.lt.nphotons)
            ! set weight for photon launch and sflag to 0 for surface reflections
            weight=1.0
            sflag=0
          if(mod(numberrun,10000).eq.0)then
             print *, numberrun,' scattered photons completed'
          end if

c***** Release photon from point source and drop weight if mismatched boundry *******************************
          call sourceph(xp,yp,zp,nxp,nyp,nzp,sint,
     +          cost,sinp,cosp,phi,xmax,ymax,zmax,twopi,
     +            xcell,ycell,zcell,nxg,nyg,nzg,iseed)

            pflag=1

      call fresnel(cost,n1,n2,sflag,weight,Rr,xp,yp,zcur,
     +                   sint,zmax,rbins,Tr,ddr)

c****** Find scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,Rr,Tr
     +                  ,xface,yface,zface,rhokap,zpos,phi,dectang,pi,
     +               twopi,xcell,ycell,zcell,tflag,iseed,delta,th,rbins,
     +                jmean,ph1,ph2,ph3,ydect,xdect,rdect,ddr,numberrun,
     +            weight,sint,n2,n1,cost,sflag,xcur,ycur,zcur,cosp,sinp)
     
c******** Photon scatters in grid until it exits (tflag=1) 

          do while(tflag.eq.0)
                 if( (ran2(iseed).lt.albedo) ) then
                 
c************ drop weight

                  absorb=weight*(1-albedo)
                  weight=weight-absorb

!                  ***** Cartesian bins *****
!                 ix=floor(xcur/ddx)+1
!                 iy=floor(ycur/ddy)+1
!                 if(iz.gt.zbins/2) iz=zbins/2 
!                 if(iy.gt.ybins/2) iy=ybins/2
!                 Amat(ix,iy,iz)=Amat(ix,iy,iz)+absorb

!                 ***** Cylindrical bins *****
                  r=sqrt(xp**2+yp**2)
                  ir=floor(r/ddr)
                  iz=floor(zcur/ddz)
                  if(iz.gt.zbins) iz=zbins
                  if(iz.lt.0) iz=0
                  if(ir.gt.rbins-1) ir=rbins-1
                  Amat(ir,iz)=Amat(ir,iz)+absorb

c************ Scatter photon into new direction and update Stokes parameters

                   call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  hgg,g2,pi,twopi,iseed)
     
c************ uses russian roulette to kill off photons

                   if(weight.lt.terminate)then
                        if(ran2(iseed).le.chance)then
                              weight=weight/chance
                        else
                        
                  absorb=weight
                  weight=weight-absorb

!                  ***** Cartesian bins *****
!                 ix=floor(xcur/ddx)+1
!                 iy=floor(ycur/ddy)+1
!                 if(iz.gt.zbins/2) iz=zbins/2 
!                 if(iy.gt.ybins/2) iy=ybins/2
!                 Amat(ix,iy,iz)=Amat(ix,iy,iz)+absorb

!                 ***** Cylindrical bins *****
                  r=sqrt(xp**2+yp**2)
                  ir=floor(r/ddr)
                  iz=floor(zcur/ddz)
                  if(iz.gt.zbins) iz=zbins
                  if(iz.lt.0) iz=0
                  if(ir.gt.rbins-1) ir=rbins-1
                  Amat(ir,iz)=Amat(ir,iz)+absorb
                              goto 100
                        endif
                   endif
                   nscatt=nscatt+1
                else
               
c************ drop weight 

                  absorb=weight*(1-albedo)
                  weight=weight-absorb

!                  ***** Cartesian bins *****
!                 ix=floor(xcur/ddx)+1
!                 iy=floor(ycur/ddy)+1
!                 if(iz.gt.zbins/2) iz=zbins/2 
!                 if(iy.gt.ybins/2) iy=ybins/2
!                 Amat(ix,iy,iz)=Amat(ix,iy,iz)+absorb

!                 ***** Cylindrical bins *****
                  r=sqrt(xp**2+yp**2)
                  ir=floor(r/ddr)
                  iz=floor(zcur/ddz)
                  if(iz.gt.zbins) iz=zbins
                  if(iz.lt.0) iz=0
                  if(ir.gt.rbins-1) ir=rbins-1
                  Amat(ir,iz)=Amat(ir,iz)+absorb
                  
                   goto 100
                endif

c************ Find next scattering location
             call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,Rr,Tr
     +                  ,xface,yface,zface,rhokap,zpos,phi,dectang,pi,
     +               twopi,xcell,ycell,zcell,tflag,iseed,delta,th,rbins,
     +                jmean,ph1,ph2,ph3,ydect,xdect,rdect,ddr,numberrun,
     +            weight,sint,n2,n1,cost,sflag,xcur,ycur,zcur,cosp,sinp)

      
          end do
100      continue
            numberrun=numberrun+1
        end do      ! end loop over nph photons
        print *, numberrun,' scattered photons completed'
            call cpu_time(finish)
            print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
      

      jmean=jmean*((pi*0.2**2)/(nphotons*((2*xmax)**3/(201**2))))
      print*, ' '  
      print*,'For detector at',xdect,ydect,zpos
      print *,'with radius',rdect
      print *,'Detects',ph1,'photons'
      print *, '',ph2
      print*,'Avereage number of scatterings = ',sngl(nscatt/nphotons)
      
      qz=201
c      open(13,file='numphotons.dat',status='unknown')
c    write(13,*) ph1
c      close(13)
c******** writes out 20 slices so that a gif(using GIMP) can be made of fluence through media
!      do q = 1,20
!            write(fn,"(i0,a)") q, '.dat'
!            open(12,file=fn)
!c            open(12,file='density.dat',status='unknown')
!            do i = 1,201
!            write(12,*) (jmean(i,j,qz),j=1,201)
!            end do
!            qz=qz-5
!           
!            close(12)
!      end do

!       open(25,file="absorb.dat")
!       do i=0,rbins-1
!       write(25,*) (Amat(i,j),j=0,zbins)
!       end do
       open(26,file='fres-same-sur1.dat')
       do i=0,zbins
            do j=0,rbins-1
                  Az(i)=Az(i)+Amat(j,i)
             end do
       end do
       Az=Az/(ddz*numberrun*mua)
       do i =0,zbins
            write(26,*) Az(i)
       end do
            close(26)
!            close(25)
            
            open(18,file="rr.dat")
            open(19,file="Tr1.dat")
       do i=0,rbins-1
            write(18,*) i,Rr(i)/(numberrun*twopi*ddr*ddr*(i-.5))
            write(19,*) i,Tr(i)/(numberrun*twopi*ddr*ddr*(i+.5))
            if(i.eq.rbins-1) then
                  j=rbins-1
            do while(j.ge.0)
            write(18,*) -j,Rr(j)/(numberrun*twopi*ddr*ddr*(j-.5))
            write(19,*) -j,Tr(j)/(numberrun*twopi*ddr*ddr*(j+.5))
                  j=j-1
            end do
            end if
       end do
            close(18)
            close(19)

            deallocate(Amat,Az,Rr,Tr)
      stop
      end
