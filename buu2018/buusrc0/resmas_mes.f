      subroutine resmas_mes(vrel,idm,density,mres,mmes)
      implicit none

      integer idm
      real*8 mres,mmes,density

      real*8 ssdown,ssup,max_dist,fu,aa
      real*8 dist_mes,rn,vrel
      integer i

      include "common"
      include "cominput"

c      write(*,*) 'in resmas_mes start'
c      call f77flush()
c  --------------  nolimit   hw
      ssdown = 2.*pmass + 0.001
      if (icbro .gt. 1 .and. idm .eq. 3) ssdown = 0.09
c  --------------  nolimit   hw
      ssup = mres - rmass
      max_dist = 0.
      do i = 1,50
        mmes = ssdown + dble(i)/50. * (ssup-ssdown)
        aa = dist_mes(vrel,idm,density,mmes,mres)
        max_dist = max(max_dist,aa)
c       if (idm .eq. 5)  write(*,*)
c    1                  ssdown, ssup, mmes, aa
      end do
 10   continue
      mmes = ssdown + rn(iseed)*(ssup-ssdown)
      fu = dist_mes(vrel,idm,density,mmes,mres)
c     write(46,*)  '  in resmas_mes  ',  mmes, fu, max_dist
      if (fu.lt.max_dist*rn(iseed)) goto 10
      return
      end
