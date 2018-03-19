
************************************************************************
*                                                                      *
      subroutine pdenp(minnum,maxnum,num,out,fac)
*                                                                      *
*       purpose:   determine momentumspace-density from momenta of     *
*                  pseudoparticles   pl-p//  plane                     *
*       variables:                                                     *
*         minnum - number of first pseudoparticle      (integer,input) *
*         maxnum - number of last pseudoparticle       (integer,input) *
*         num    - number of pseudoparticles/nucleon   (integer,input) *
*         out    - number of pseudopart. out of range (integer,output) *
*         fac    - factor                                 (real,input) *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      integer out,iz,ix,minnum,maxnum,num,i
      real*8 fac
*
      do 200 iz = -24,24
        do 100 ix = -20,20
          prho(ix,iz) = 0.0
  100   continue
  200 continue
*
      out = 0
*
      do 300 i = minnum,maxnum
        ix =  int( fac  * sqrt(p(1,i)**2+p(2,i)**2) )
        iz = nint( fac  * p(3,i) )
        if( ix .le. -20 .or. ix .ge. 20 .or.
     &      iz .le. -24 .or. iz .ge. 24   )    then
          out = out + 1
        else
          if(ix.le.1) then
          prho( 0,iz) = prho( 0,iz) + 1.0/2.25/pi
          prho( 1,iz) = prho( 1,iz) + 1.0/2.25/pi
          prho(-1,iz) = prho(-1,iz) + 1.0/2.25/pi
          else
          prho(ix,iz) = prho(ix,iz) + 1.0/2./pi/float(ix)
          prho(-ix,iz) = prho(-ix,iz) + 1.0/2./pi/float(ix)
          end if
        end if
  300 continue
*
      do 500 iz = -24,24
        do 400 ix = -20,20
          prho(ix,iz) = prho(ix,iz) / dble(num)
  400   continue
  500 continue
*
      return
      end
