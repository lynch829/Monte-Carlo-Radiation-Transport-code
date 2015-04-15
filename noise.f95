program perlin

implicit none

real, allocatable :: noise(:,:)
real :: ran2,smooth,inter,r,p,b
integer :: i,j,iseed,ngrid

ngrid=204

allocate(noise(1:ngrid,1:ngrid))
iseed=15628628
noise=0.
!**** circle code
!do i=1,ngrid
!      do j=1,ngrid
!            r=sqrt((real(j)-103.)**2+((real(i)-103.)**2))
!            p=-r
!            if((r.gt.(90.)).and.(r.lt.(100.))) then
!            noise(i,j)=noise(i,j) + .5
!            elseif((p.gt.-(100.+100.)).and.(p.lt.-(100.+100.)))then
!                        noise(i,j)=noise(i,j) + .5
!                        else
!            noise(i,j)=noise(i,j)
!            end if
!      end do
!end do

! adds 3 dots
!do i=1,ngrid
!      do j=1,ngrid
!            if(sqrt(real((i-100.)**2+(j-100.)**2)).lt.20.) then
!            noise(i,j)=0.
!            elseif(sqrt(real((i-150.)**2+(j-150.)**2)).lt.20.) then
!            noise(i,j)=0.
!            elseif(sqrt(real((i-40.)**2+(j-30.)**2)).lt.15.) then
!            noise(i,j)=0.            
!            else
!            noise(i,j)=noise(i,j)+ran2(iseed)
!            end if
!      end do
!end do

!just noise
do i=1,ngrid
      do j=1,ngrid
            noise(i,j)=noise(i,j)+ran2(iseed)
      end do
end do



do j=1,ngrid
      do i=1,ngrid
            noise(i,j)=smooth(i,j,noise,ngrid) 
      end do
end do


open(1,file='noisesmooth.dat')
do i=1,ngrid
           write(1,*) ((noise(i,j)),j=1,ngrid)
end do
close(1)
deallocate(noise)

end program perlin

function smooth(x,y,noise,ngrid)

integer, intent(in) :: x,y,ngrid
integer :: ix,iy
real :: corner,side,center
real, intent(in) :: noise(1:ngrid,1:ngrid)

ix=x
iy=y

if(ix-1.lt.1) ix=2
if(iy-1.lt.1) iy=2
if(ix+1.gt.ngrid) ix=ngrid-1
if(iy+1.gt.ngrid) iy=ngrid-1
corner=(noise(ix-1,iy-1)+noise(ix+1,iy-1)+noise(ix-1,iy+1)+noise(ix+1,iy+1))/16.
side=(noise(ix-1,iy)+noise(ix+1,iy)+noise(ix,iy-1)+noise(ix,iy+1))/8.
 center=noise(ix,iy)/4.
 smooth=corner+side+center
return
end
