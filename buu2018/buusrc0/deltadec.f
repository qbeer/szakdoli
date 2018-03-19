************************************************************************
*                                                                      *
      subroutine deltadec(dt,num,igamma)
*                                                                      *
************************************************************************
      implicit none
      real*8 dt,ww,dendel
      integer num,igamma,k,iendel,ixx,iyy,izz
      include 'common'
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*   loop over all testparticles                                        *
      do 1000 k = 1,maxpar
        if(abs(id(1,k)) .le. 1)                                goto 1000
        if(id(1,k) .gt. nres+1)                                goto 1000
c        write(*,*) "test 1 for deldil passed"
        if(time-dtim(k) .le. dt)                               goto 1000
c        write(*,*) "test 2 for deldil passed"
        if((id(2,k)+2)/2 .ne. 1)                               goto 1000
c        write(*,*) "test 3 for deldil passed"
        ixx = nint(r(1,k))
        iyy = nint(r(2,k))
        izz = nint(r(3,k))
        dendel = 0.0
       if(iabs(ixx).le.maxx.and.iabs(iyy).le.maxx.and.iabs(izz).le.maxz)
     &      dendel = rhb(ixx,iyy,izz)/rho0
        iendel = nint(10.0*dendel)
        ww = 1.0
cc        if(ideldil.eq.1) call deldil(dt,k,id(1,k),e(k),
cc     &                 r(1,k),r(2,k),r(3,k),p(1,k),p(2,k),p(3,k),ww,0,
cc     &                                            ndedel(k),iendel)
cc        if(igamma.eq.1) call delgam(dt,k,id(1,k),e(k),
cc     &                 r(1,k),r(2,k),r(3,k),p(1,k),p(2,k),p(3,k),ww,0)
 1000 continue
      return
      end
