      subroutine tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,Rr,Tr
     +                  ,xface,yface,zface,rhokap,zpos,phi,dectang,pi,
     +               twopi,xcell,ycell,zcell,tflag,iseed,delta,th,rbins,
     +                jmean,ph1,ph2,ph3,ydect,xdect,rdect,ddr,numberrun,
     +            weight,sint,n2,n1,cost,sflag,xcur,ycur,zcur,cosp,sinp)

      implicit none

      include 'grid.txt'

      integer tflag,iseed,xcell,ycell,zcell,ph1,ph2,ph3,sflag,numberrun
      real xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,phi,ddr,Rr(0:rbins-1)
      real xpos,ypos,zpos,tpos,flu,n1,n2,sint,cost,weight
      real ran2,pi,twopi

      integer celli,cellj,cellk,rbins
      real tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
      real dx,dy,dz,smax,delta,ydect,xdect,rdect,dectang,th,cosp,sinp
      real Tr(0:rbins-1)
c      real jmean(nxg,nyg,nzg)



c***** tflag=0 means photon is in envelope
      tflag=0
!      dsx=0.
!      dsy=0.
!      dsz=0.
!      dx=0.
!      dy=0.
!      dz=0.
c**** generate random optical depth tau
      tau=-alog(ran2(iseed))

c***** set the cumulative distance and optical depth (d and taurun) 
c***** along the photon path to zero.  set the current photon coordinates.
c***** note that the origin of the (xcur,ycur,zcur) system is at the 
c***** bottom corner of the grid.
      taurun=0.    
      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax
345   continue      
      celli=xcell
      cellj=ycell
      cellk=zcell 
        

      if (cellk.lt.1 )cellk=1 
      d=0.
c***** calculate smax -- maximum distance photon can travel
      if(nxp.gt.0.) then
         dsx=(2.*xmax-xcur)/nxp
      elseif(nxp.lt.0.) then
         dsx=-xcur/nxp
      elseif(nxp.eq.0.) then
         dsx=1.e2*xmax
      endif

      if(nyp.gt.0.) then
         dsy=(2.*ymax-ycur)/nyp
      elseif(nyp.lt.0.) then
         dsy=-ycur/nyp
      elseif(nyp.eq.0.) then
         dsy=1.e2*ymax
      endif

      if(nzp.gt.0.) then
         dsz=(2.*zmax-zcur)/nzp
      elseif(nzp.lt.0.) then
         dsz=-zcur/nzp
      elseif(nzp.eq.0.) then
         dsz=1.e2*zmax
      endif

      smax=amin1(dsx,dsy,dsz)
      if(smax.lt.delta) then
         tflag=1
         return
      endif
       
c***** integrate through grid
      do while((taurun.lt.tau).and.(d.lt.(.999*smax)))

c***** find distance to next x, y, and z cell walls.  
c***** note that dx is not the x-distance, but the actual distance along 
c*****the direction of travel to the next x-face, and likewise for dy and dz.
         if(nxp.gt.0.) then
            dx=(xface(celli+1)-xcur)/nxp
            if(dx.lt.delta) then
               xcur=xface(celli+1)
               celli=celli+1
               dx=(xface(celli+1)-xcur)/nxp
            endif
         elseif(nxp.lt.0.) then
            dx=(xface(celli)-xcur)/nxp
            if(dx.lt.delta) then
               xcur=xface(celli)
               dx=(xface(celli-1)-xcur)/nxp
               celli=celli-1
            endif
         elseif(nxp.eq.0.) then
            dx=1.e2*xmax
         endif

         if(nyp.gt.0.) then
            dy=(yface(cellj+1)-ycur)/nyp
            if(dy.lt.delta) then
               ycur=yface(cellj+1)
               cellj=cellj+1
               dy=(yface(cellj+1)-ycur)/nyp
            endif
         elseif(nyp.lt.0.) then
            dy=(yface(cellj)-ycur)/nyp
            if(dy.lt.delta) then
               ycur=yface(cellj)
               if(cellj-1.eq.0) cellj=2
               dy=(yface(cellj-1)-ycur)/nyp
               cellj=cellj-1
            endif
         elseif(nyp.eq.0.) then
            dy=1.e2*ymax
         endif

         if(nzp.gt.0.) then
            dz=(zface(cellk+1)-zcur)/nzp
            if(dz.lt.delta) then
               zcur=zface(cellk+1)
               cellk=cellk+1
               dz=(zface(cellk+1)-zcur)/nzp
            endif
         elseif(nzp.lt.0.) then
            if (cellk.gt.nzg+3) then
            cellk=nzg+3
            print *, 'cellk > 204'
            end if
            dz=(zface(cellk)-zcur)/nzp
            if(dz.lt.delta) then
               zcur=zface(cellk)
               if(cellk.lt.2) then
               cellk=2
               print *, 'cellk < 2'
               end if
               dz=(zface(cellk-1)-zcur)/nzp
               cellk=cellk-1
               if(cellk.lt.1.) then
                  print *, 'cellk<1'
                  cellk=1
               end if
            endif
         elseif(nzp.eq.0.) then
            dz=1.e2*zmax
         endif

c***** distances are only zero if photon is on cell wall.  if it is 
c***** on cell wall then set to arbitrary large distance, since we will
c***** in fact hit another wall
         if( (dx.eq.0.) .or. ((abs(dx)).lt.(delta)) ) dx=1.e2*xmax
         if( (dy.eq.0.) .or. ((abs(dy)).lt.(delta)) ) dy=1.e2*ymax
         if( (dz.eq.0.) .or. ((abs(dz)).lt.(delta)) ) dz=1.e2*zmax
            
c***** find distance to next cell wall -- minimum of dx, dy, and dz
         dcell=amin1(dx,dy,dz)
         if(dcell.le.0.) then
            print *,'tauint2: dcell < 0'
         endif
         if(dx.lt.0.) dcell=amin1(dy,dz)
         if(dy.lt.0.) dcell=amin1(dx,dz)
         if(dz.lt.0.) dcell=amin1(dx,dy)

c***** optical depth to next cell wall is 
c***** taucell= (distance to cell)*(opacity of current cell)
         taucell=dcell*rhokap(celli,cellj,cellk)

c***** if taurun+taucell>tau then scatter at distance d+d1.  
c***** update photon position and cell.  
c***** if taurun+taucell<tau then photon moves distance dcell 
c***** (i.e. ends up on next cell wall) and update photon position
c***** and cell.
         if((taurun+taucell).ge.tau) then
            d1=(tau-taurun)/rhokap(celli,cellj,cellk)
            d=d+d1
            taurun=taurun+taucell
            xcur=xcur+d1*nxp
            ycur=ycur+d1*nyp
            zcur=zcur+d1*nzp
            jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + d1

c*************** Linear Grid ************************
            celli=int(nxg*xcur/(2.*xmax))+1
            cellj=int(nyg*ycur/(2.*ymax))+1
            cellk=int(nzg*zcur/(2.*zmax))+1
c****************************************************

         else

            d=d+dcell
            taurun=taurun+taucell
            xcur=xcur+dcell*nxp
            ycur=ycur+dcell*nyp
            zcur=zcur+dcell*nzp
            jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + dcell

c*************** Linear Grid ************************
            celli=int(nxg*xcur/(2.*xmax))+1
            cellj=int(nyg*ycur/(2.*ymax))+1
            cellk=int(nzg*zcur/(2.*zmax))+1
c****************************************************

          endif
!      if(numberrun.gt.217703) print *, xcur,ycur,zcur,d
      end do
      
c***** calculate photon final position.  if it escapes envelope then
c***** set tflag=1.  if photon doesn't escape leave tflag=0 and update 
c***** photon position.
      if(ycur.gt..99999*2.*ymax) tflag=1
      if(ycur.lt.1E-6) tflag=1
      if(xcur.gt..99999*2.*xmax) tflag=1
      if(xcur.lt.1E-6) tflag=1
      if((d.ge.(.99999*smax))) then
            d=smax
            if(zcur.ge..9999999*(2*zmax).or.zcur.le.1e-8)then
            if(zcur.lt.1E-8) then
!            if(cost.lt.0.) print *, 90.-abs(180.*cost/3.14)
            end if
          call fresnel(sflag,sint,cost,n1,n2,zcur,
     & ddr,ycur,weight,xmax,ymax,zmax,xp,yp,zp,xcur,
     & Rr,rbins,Tr,twopi,tflag,sinp,cosp)
     
      cost=cos(pi-asin(sint))
!      sint=(1.-cost*cost)
!      if(sint.le.0) then
!            sint=1.
!      else
!            sint=sqrt(sint)
!      end if
            xp=xcur-xmax
            yp=ycur-ymax
            zp=zcur-zmax
     
!            nxp=sint*cosp  
!            nyp=sint*sinp   
            nzp=cost
            
            goto 345
            else
                  ! exits out a side, so ignore
                 tflag=1
            endif
      else
         xp=xp+d*nxp
         yp=yp+d*nyp
         zp=zp+d*nzp
         xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
         ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
         zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
         
      endif

      return
      end
