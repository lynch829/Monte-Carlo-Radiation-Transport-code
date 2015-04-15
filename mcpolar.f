      program mcpolar

      implicit none

      include 'grid.txt'
      include 'photon.txt'

c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,tflag,i
      integer cnt,q,qz,numberrun,sflag,pflag,abins,io
      real*8 kappa,albedo,hgg,xmax,ymax,zmax,pi,twopi,fourpi,g2,delta
      real*8 chance,terminate,weight,absorb,n1,n2,mua,mus,xcur,ycur,zcur
      character(len=70) :: fn
      integer ix,iy,iz,xbins,ybins,zbins,rbins,ir
      real*8 ddz,ddx,ddy,ddr,r,dda,dalpha,domega,nscatt,rd,td
      ! bin array for cartesian bins
c      real, allocatable :: Amat(:,:,:)
      ! bin arrays for cylindrical bins
      real*8, allocatable :: Amat(:,:),Az(:),Rr(:),Tr(:),td_ang(:,:)
      real*8, allocatable ::rrang(:),rd_ang(:,:),angular(:),noise(:,:)
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
          close(10)
          
          kappa=mua+mus
          albedo=mus/kappa
          
c***** read in noise file and get size then allocate noise array based on size of file
          open(32,file='noisesmooth.dat')
          do 
            read(32,*,IOSTAT=io)

          if (io < 0) then
            close(32)
            ! allocates the arrays and inits the variables
            allocate(noise(1:cnt,1:cnt))
            noise = 0.0
            exit
          else 
            cnt = cnt + 1
          end if
          end do
          open(33,file='noisesmooth.dat')
          do i=1,cnt
            read(33,*) (noise(i,j),j=1,cnt)
          end do
          close(33)

c***** Set up bins ***************************************************

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
      zbins=204
      rbins=500
      ! set size of bins
      ddz=((2.*zmax)/zbins)
      ddr=((2.*xmax)/rbins) 
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
      ! set terminate and chance values for roulette sub-routine
      terminate=0.0001
      chance=0.1
      ! set number of photons run to 0
      numberrun=0
      iseed=-abs(iseed)  ! Random number seed must be negative for ran2
      g2=hgg*hgg  ! Henyey-Greenstein parameter, hgg^2

!      **** Angle bins ****    
   
      abins=30
!      dda=(pi/2.)/abins
      dda=((pi/2.)/abins)
      allocate(td_ang(0:rbins-1,0:abins),rrang(0:rbins-1))
      allocate(rd_ang(0:rbins-1,0:abins),angular(0:abins-1))
      td_ang=0.
      rd_ang=0.
      rrang=0.
      angular=0.
    
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
            xcur=xp+xmax
            ycur=yp+ymax
            zcur=zp+zmax
            pflag=0

      call fresnel(sint,cost,sinp,cosp,nxp,nyp,nzp,tflag,
     + iseed,n1,n2,xp,yp,zp,xcur,ycur,Rr,ddr,weight,dda,noise,
     + rd_ang,xmax,ymax,zmax,zcur,rbins,Tr,sflag,td_ang,abins
     + ,cnt)

c****** Find scattering location
        call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,Rr,Tr,kappa
     + ,pflag,xface,yface,zface,rhokap,phi,pi,zbins,twopi,
     + xcell,ycell,zcell,tflag,iseed,delta,rbins,Amat,jmean,numberrun
     + ,weight,sint,n2,n1,cost,sflag,xcur,ycur,cosp,sinp,abins,td_ang
     + ,dda,rd_ang,ddr,zcur,noise,cnt)

c******** Photon scatters in grid until it exits (tflag=1) 

          do while(tflag.eq.0)
                                  
c************ drop weight

                  absorb=weight*(mua/kappa)
                  weight=weight*albedo
!!                  ***** Cartesian bins *****
!!                 ix=floor(xcur/ddx)+1
!!                 iy=floor(ycur/ddy)+1
!!                 if(iz.gt.zbins/2) iz=zbins/2 
!!                 if(iy.gt.ybins/2) iy=ybins/2
!!                 Amat(ix,iy,iz)=Amat(ix,iy,iz)+absorb

!!                 ***** Cylindrical bins *****
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
                   if(weight.le.terminate)then
                        if(ran2(iseed).le.chance)then
                              weight=weight/chance
                        else
                  absorb=weight
                  weight=0.
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
               
c************ Find next scattering location
           call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,Rr,Tr,kappa
     + ,pflag,xface,yface,zface,rhokap,phi,pi,zbins,twopi,
     + xcell,ycell,zcell,tflag,iseed,delta,rbins,Amat,jmean,numberrun
     + ,weight,sint,n2,n1,cost,sflag,xcur,ycur,cosp,sinp,abins,td_ang
     + ,dda,rd_ang,ddr,zcur,noise,cnt)

      
          end do
100      continue
        if(xcur.lt.1E-8.or.xcur.gt.2.*xmax) print *, 'Error xcur:',xcur
        if(ycur.lt.1E-8.or.ycur.gt.2.*ymax) print *, 'Error ycur:',ycur
!         if(zcur.le.1E-8)print *, numberrun,zcur,cost
            numberrun=numberrun+1
        end do      ! end loop over nph photons
        print *, numberrun,' scattered photons completed'
            call cpu_time(finish)
            if(finish-start.ge.60.)then
            print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
            else
            print*, 'time taken ~',floor(finish-start/60.),'s'
            end if
      jmean=jmean*((nxg**3)/(xmax*ymax*zmax*numberrun))
      print*, ' '  

      print*,'Avereage number of scatterings = ',sngl(nscatt/nphotons)
c******* Fluence normalisation and output
      qz=130
c******** writes out 20 slices so that a gif(using GIMP) can be made of fluence through media
      do q = 1,20
            write(fn,"(i0,a)") q, '.dat'
            open(12,file=fn)
c            open(12,file='density.dat',status='unknown')
            do i = 1,201
            write(12,*) (jmean(i,qz,j),j=1,201)
            end do
            qz=qz-3
           
            close(12)
      end do
! path length fluence normaiastion and outputs it
!       open(25,file="jmean.dat")
!       do i=1,204
!            do j=1,zbins
!                  Az(i)=Az(i)+jmean(j,i)
!            end do
!       end do
!       az=2.*az*(0.1*0.01*mua)
!       do i=1,204	
!       write(25,*) Az(i)
!       end do
!       close(25)


! normalises fluence and outputs it.
       open(26,file='flunoise.dat')
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


c******** loops to normalise reflectance and transmittance arrays *****
! for Td_a
       open(49,file='td_a.dat')
       do i=0,abins-1
            do j=0,rbins-1
                 dalpha=twopi*(j+0.5)*ddr*ddr
                 domega=fourpi*sin(((i+0.5)*dalpha))*sin(dalpha/2.)
                 if(j.eq.0)domega=twopi*sin(pi/120.)*(pi/60.)
                  angular(i)=angular(i)+(td_ang(j,i)/abs(domega))
            end do
       end do
       do i=0,abins-1
!            dalpha=twopi*(i+0.5)*ddr*ddr
!            domega=twopi*sin(((i+0.5)*dalpha))*sin(dalpha/2.) 
!            domega=1.
            write(49,*) angular(i)/(numberrun)
       end do
       close(49)
       angular=0.
       
       ! for Rd_a

       open(50,file='rd_a.dat')
       do i=0,abins-1
            do j=0,rbins-1
             dalpha=twopi*(j+0.5)*ddr*ddr
            domega=fourpi*sin(((i+0.5)*dalpha))*sin(dalpha/2.)
            if(j.eq.0)domega=twopi*sin(pi/120.)*(pi/60.)
                  angular(i)=angular(i)+(rd_ang(j,i)/abs(domega))
            end do
       end do
       do i=0,abins-1
!            dalpha=twopi*(i+0.5)*ddr*ddr
!            domega=fourpi*sin(((i+0.5)*dalpha))*sin(dalpha/2.)
!            domega=1. 
            write(50,*) angular(i)/(numberrun)
       end do
       close(50)
       rrang=0.
       
c********** total Rd and Td *******************
       do i=0,rbins-1
            do j=0,abins-1
                  rrang(i)=rrang(i)+rd_ang(i,j)
            end do
       end do
       do i=0,rbins-1
            rd=rd+rrang(i)
       end do
       rrang=0.
       
       do i=0,rbins-1
            do j=0,abins-1
                  rrang(i)=rrang(i)+td_ang(i,j)
            end do
       end do
       do i=0,rbins-1
            td=td+rrang(i)
       end do
       print *, 'Total diffuse reflectance:',rd/numberrun
       print *, 'Total transmitance:',td/numberrun
       
            !for transmittance and reflectance radially dept.
            open(18,file="rr.dat")
            open(19,file="Tr3.dat")
       do i=0,rbins-1
!       print *, ((i+.5)+(1/(12*(i+0.5))))
            write(18,*) i,(Rr(i)/(numberrun*twopi*ddr*ddr*(
     +           (i+.5))))
            write(19,*) i,(Tr(i)/(numberrun*twopi*ddr*ddr*(
     +           (i+.5))))
            if(i.eq.rbins-1) then
                  j=rbins-1
            do while(j.ge.0)
            write(18,*) -j,(Rr(j)/(numberrun*twopi*ddr*ddr*(
     +           (j+.5))))
            write(19,*) -j,(Tr(j)/(numberrun*twopi*ddr*ddr*(
     +           (j+.5))))
                  j=j-1
            end do
            end if
       end do
            close(18)
            close(19)

            deallocate(Amat,Az,Rr,Tr,rrang,td_ang,rd_ang,angular)
      stop
      end
