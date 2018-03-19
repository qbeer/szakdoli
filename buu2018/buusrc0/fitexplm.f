
************************************************************************
*                                                                      *
      subroutine fitexplms(pit,epic,ikin,iti ,e0,b,imin,imax)
*                                                                      *
************************************************************************
      implicit none
      integer maxbin
      real*8 x,x2,y,xy,sum,b,e0,xi,yi,a
      integer ii,ikin,iti,imin,imax
      parameter     (maxbin =   50)
      include"common"
*----------------------------------------------------------------------*
      real*8  pit(11,maxbin),epic(maxbin)
*----------------------------------------------------------------------*
      x   = 0.0
      x2  = 0.0
      y   = 0.0
      xy  = 0.0
      sum = 0.0
      b   = 0.0
      e0  = 0.0
      do 50  ii = 1,ikin
        pit(11,ii) = 0.0
  50  continue
      if(pit(iti ,imin) .lt. 1.-5) return
      do 100  ii = imin,imax
        if(pit(iti ,ii) .lt. .1)                                goto 100
        yi  = log(pit(iti ,ii))
        xi  = epic(ii)
        x   = x  + xi
        x2  = x2 + xi**2
        y   = y  + yi
        xy  = xy + xi * yi
        sum = sum + 1.0
 100  continue
      if(sum .lt. 0.5) return
      x   = x / sum
      x2  = x2/ sum
      y   = y / sum
      xy  = xy/ sum
      a   = (x * y - xy) / (x * x - x2)
      b   = y - a * x
      do 200  ii = 1,ikin
        pit(11,ii) = exp( a * epic(ii) + b)
 200  continue
      e0 = -1000./a
      return
      end
