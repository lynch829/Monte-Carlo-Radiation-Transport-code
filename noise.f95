program perlin

implicit none

real, allocatable :: noise(:,:)
real :: ran2,smooth,inter
integer :: i,j,iseed,ngrid

ngrid=204

allocate(noise(1:ngrid,1:ngrid))
iseed=15628628
do i=1,ngrid
      do j=1,ngrid
            noise(i,j)=noise(i,j) + (2.*ran2(iseed)-1.)
      end do
end do



do j=1,ngrid
      do i=1,ngrid
            noise(i,j)=smooth(i,j,noise,ngrid) 
      end do
end do

!do j=1,ngrid
!do i=1,ngrid
!      noise(i,j)=inter(j,j,noise,ngrid)
!end do
!end do

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
side=(noise(ix-1,iy)+noise(ix+1,1)+noise(ix,iy-1)+noise(ix,iy+1))/6.
 center=noise(ix,y)/4.
 smooth=corner+side+center
return
end
