
************************************************************************
*                                                                      *
      subroutine inir(minnum,maxnum,num,radi,x0,z0,
     &                iseed,mass,mspr,msto,surf)
*                                                                      *
*       purpose:     providing initial conditions for coordinate       *
*                    distribution of testparticles                     *
*                    according to woods-saxon with saa=0.15 and cut=1.0*
*       variables:   (all input)                                       *
*         minnum  - first testparticle treated in one run    (integer) *
*         maxnum  - last testparticle treated in one run     (integer) *
*         num     - number of testparticles per nucleon      (integer) *
*         radi    - radius of nucleus "fm"                      (real) *
*         x0,z0   - displacement of center of nucleus in x,z-          *
*                   direction "fm"                              (real) *
*         iseed   - seed for random-number generator         (integer) *
*         mass    - total mass of the system                 (integer) *
*         mspr    - proton number in the one nucleus         (integer) *
*         msto    - total number in the one nucleus          (integer) *
*                                                                      *
************************************************************************
      implicit none
      include"common"
*
*-----------------------------------------------------------------------
****  surface parameter for test particle distribution
*      next mass dependence is obtained by considering the ground
*      state stability, only for t0+t3+yukawa+coulomb with rpot=0.3333
*      and smu=2.175
*
      integer minnum,maxnum,num,iseed,mass,mspr,msto
      real*8 cut , testfu, rdis, x,y,z , radi
      logical flag
      real*8    saa
      integer idnum, icount,irun,init,i, ii
      real*8    rtest, rn, rram, rwod, rsqr, x0,z0
      real*8    surf
      real*8  aa, rdeut, xmax, ymax, rho, zeta, seta
      write(*,*)'in inir surf = ', surf

      saa = surf
      radi = radi
!      write(20,*)'radi = ', radi
!      write(20,*)'saa = ', saa

*-----------------------------------------------------------------------
*     target-id = 1 and projectile-id = -1
*     proton-id = 1 and neutron-id = -1
*     nucleon-id= 1(and delta-id = 2, ...)
*
*
*  id(1,i) = particle type (nucleon,delta, nstar...)
*  id(2,i) = particle charge
*  id(3,i) = last colliding partner of particle i
*  id(4,i) = for resonances: how many times the meson was created
*  id(5,i) = number of pion absorption
*  id(6,i) = abs: number of barion collision; sign: target or projectile
*
*
      if (minnum .eq. 1) then
        idnum = 1
      else
        idnum = -1
      end if
      flag = .true.
      icount = 0
      do while(flag)
        icount = icount + 1
        rtest = radi + 0.001*float(icount)
        testfu= 1.0/(1.0 + exp((rtest-radi)/saa))

**************** gauss density distr for C******************
       if(mass.eq.13) then
         testfu=(1.+1.082*(rtest/1.692)**2)*exp(-(rtest/1.692)**2)
       endif
************************************************************

        if(testfu .lt. 1.0e-06) then
          flag = .false.
          cut  = rtest
        end if
      end do

*-----------------------------------------------------------------------
*     identification of testparticles and assigment of restmass
*     occupation of coordinate-space
*     put particles in their position in coordinate-space
*     (shift and relativistic contraction)
*
*     loop over all parallel runs:
c      write(*,*)'inir0 before loop',maxb,minnum,maxnum
      do 400 irun = 1,num
        init     = (irun-1)*maxb
        do 300 i = minnum+init,maxnum+init
          p(1,i) = 0.0
          p(2,i) = 0.0
          p(3,i) = 0.0
          id(1,i) = 1
          if(mspr.eq.-1) id(1,i) = -1
          id(3,i) = 0
          id(4,i) = 0
          id(5,i) = 0
          id(6,i) = idnum
          e(i)  = rmass
          if(mspr.eq.-1) then
            id(2,i) = -1
          else if(i.lt.minnum+init+mspr) then
            id(2,i) = 1
          else
            id(2,i) = 0
          end if
*
         if (maxnum-minnum .eq. 1) then         !    deuteron
         if (init + minnum .eq. i) then         !    deuteron
           aa = 6.16 - 1.0                        !   Hulthen wave fct
           rdeut = 0.5 / 0.231                     !    fm
           xmax = (1./(2.*aa+1.))**(1./aa)
           ymax = xmax*(1.-xmax**aa)**2
   50      continue
           x = rn(iseed)
           y = x*(1. - x**aa)**2  / ymax
           if (y  .lt. rn(iseed) ) goto 50
           rho = -log(x)*rdeut
           zeta = (2.*rn(iseed) - 1.)
           seta = sqrt(abs(1.-zeta**2))
           z  = rho * zeta
           x =  (6.283*rn(iseed))
           y =  rho * seta * sin(x)
           x =  rho * seta * cos(x)
           r(1,i) = x  + x0
           r(2,i) = y
           r(3,i) = z  + z0
           r(1,i+1) = -x  + x0
           r(2,i+1) = -y
           r(3,i+1) = -z  + z0
         endif
         goto 300
         endif
          if(msto.gt.1) then
  200     continue
          x =(1.0 - 2.0 * rn(iseed))*cut
          y =(1.0 - 2.0 * rn(iseed))*cut
          z =(1.0 - 2.0 * rn(iseed))*cut
          rsqr = x*x+y*y+z*z
          rdis = sqrt(rsqr)
          if(rdis.gt.cut) goto 200

          rram = rn(iseed)
          rwod = 1.0/(1.0 + exp((rdis-radi)/saa))

**************** gauss density distr for C ******************
          if(mass.eq.13) then
           rwod=(1.+1.082*(rdis/1.692)**2)*exp(-(rdis/1.692)**2)
          endif
************************************************************

          if(rram.gt.rwod) goto 200
*
          r(1,i) = x  + x0
          r(2,i) = y
          r(3,i) = z  + z0
*
          else
          r(1,i) = x0
          r(2,i) = 0.0
          r(3,i) = z0
          end if

          if(i.eq.6424) then
            write(*,*)'ini inir ', r(1,i),r(2,i),r(3,i)
            write(*,*) x0, z0
          end if
*
  300   continue
  400 continue
*

      return
      end



