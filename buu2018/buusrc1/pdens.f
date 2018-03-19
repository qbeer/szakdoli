
************************************************************************
*                                                                      *
      subroutine pdens(minnum,maxnum,num,out,fac)
*                                                                      *
*       purpose:   determine momentumspace-density from momenta of     *
*                  pseudoparticles   px-pz plane                       *
*       variables:                                                     *
*         minnum - number of first pseudoparticle      (integer,input) *
*         maxnum - number of last pseudoparticle       (integer,input) *
*         num    - number of pseudoparticles/nucleon   (integer,input) *
*         out    - number of pseudopart. out of range (integer,output) *
*         fac    - facor                                  (real,input) *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      integer out,iz,ix,iy,minnum,maxnum,num,i
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
        ix = nint( fac * p(1,i) )
        iy = nint( fac * p(2,i) )
        iz = nint( fac * p(3,i) )
        if( ix .le. -20 .or. ix .ge. 20 .or.
     &      iz .le. -24 .or. iz .ge. 24 .or.
     &      iy .le. -4  .or. iy .ge. 4    )    then
          out = out + 1
        else
          prho(ix,iz) = prho(ix,iz) + 1.0
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
