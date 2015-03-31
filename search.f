       subroutine search(cnt,cdf,nlow,ran,iseed)
            implicit none

            integer nup,nlow,cnt,middle,iseed
            real*8  ran,ran2
            real*8 cdf(cnt-1)
            

c search by bisection algorithm. finds bracketing indice nlow and nlow+1 that bracket a random number in the cdf
            nup = cnt
            nlow = 1
            middle = int((nup+nlow)/2.)
            ran = ran2(iseed)

            do while((nup - nlow) .gt. 1)

                  middle = int((nup+nlow)/2.)

                  if (ran .gt. cdf(middle)) then
                        nlow = middle

                  else
                        nup = middle

                  end if
            end do
      return
      end 
