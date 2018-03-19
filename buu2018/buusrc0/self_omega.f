************************************************************************
      subroutine self_omega(vrel, dens,m_omega,rself,iself,sgamma)
      implicit none

      include "common"
      include "cominput"

      real*8 totom(7,2000),totro(7,2000)
      COMMON/self_stored/ totom,totro
      integer idm
      real*8 vrel,dens,m_omega,rself,iself, coll_bro
      real*8 exmass, fgamma,sgamma

      idm = 5
      exmass = m_omega - 3.*pmass
      if (exmass .le. .01) then
          fgamma = .0
      else
          fgamma = (1. - exp(-10.*exmass))/0.97
      endif
cccccccccccc matter calculation usually 50 Mev Shift
      coll_bro  =  .0
      if(icbro.ge.1) call broadening(vrel,dens,idm, m_omega,coll_bro)
      rself = .0
      if(ivecmatt.eq.1)
     &  rself = dens / rho0 * 2. * omass * omemshift  ! omemshift shift
c     write(*,*)  '  self omega ', dens,0.00849*fgamma, coll_bro
      iself = - omass * (owidth*fgamma + coll_bro)
      sgamma = owidth
      return
      end

      subroutine broadening(vrel,dens,idm, m_meson,coll_bro)
      implicit none

      include "common"
      include "cominput"

      real*8 totom(7,2000),totro(7,2000)
      COMMON/self_stored/ totom,totro

      integer idm
      real*8 vrel,dens,m_meson, coll_bro, cross1 ,eps, pf, vrel2, ef
      real*8 help,cross,qq2,dummyf(1:9),ss,bwdist_selfen,isofac,vrel1
      real*8 p00, e00
      integer jres

      coll_bro = 0.0
      if(dens .le. 0.001) return
         
      pf = (1.5* pi**2 * dens)**0.33333
      ef = sqrt(rmass**2+pf**2)

      vrel1 = abs(vrel)
      if (vrel1.ge.1.0) then
        write(*,*) 'vrel>1 in self_meson!!',vrel
        vrel1 = 0.99
c        stop
      end if

      if (vrel1 .lt. 1.e-3)  then
        vrel2 = (ef/pf)**3-(rmass/pf)**3 -3.0*rmass**2*(ef-rmass)/pf**3
      else if (vrel1.gt.pf/ef) then
        vrel2 = vrel1 + 1.0/vrel1*(0.3333-rmass**2/pf**2+(rmass/pf)**3*
     &           atan(pf/rmass))
      else
         p00 = rmass*vrel1/sqrt(1.0-vrel1**2)
         e00 = rmass/sqrt(1.0-vrel1**2)
         
         vrel2 = (p00/pf)**3 * vrel1+ 
     &     1.0/vrel1*(0.3333-rmass**2*p00/pf**3
     &    +(rmass/pf)**3*atan(p00/rmass))
     &    + 1.0/pf**3 * ( ef**3-e00**3 - 
     &     3.0*rmass**2 * (ef-e00)+0.3333*vrel1**2 *(ef**3-e00**3))
      end if

      eps = .001
      if (vrel1 .lt. 1.e-3)  vrel1 = 1.e-3
      if (vrel1 .gt. 0.95)   vrel1 = 0.95

      ss = m_meson**2 + rmass**2 + 2.*rmass*m_meson/sqrt(1.-vrel1**2)
      cross = 0.0
      do jres=1,nres
        help = bwdist_selfen(ss,idm,jres,m_meson,dummyf)
        qq2= 0.25*(ss-rmass**2+m_meson**2)**2/ss-m_meson**2
        isofac = (resprop2(jres,1) + 1.)/2.
        cross1 = 4.0*pi/(qq2+eps) * isofac * (resprop2(jres,3)+1.)/6.
     &     * help * hbc**2      ! fm^2
        cross = cross + cross1
c        write(*,*) 'dummyf in broadening:',m_meson,jres,help,dummyf
c    1   , '  jres, cross ', jres, cross1, cross, isofac,qq2,vrel1,
c    2   ss
c         write(mspfpri,*) '  jres, cross ', jres, cross1, cross
      end do

      coll_bro =  dens*abs(vrel1)/sqrt(1-vrel1**2)*cross * hbc
c      write(*,*) 'cross in broaden',m_meson,cross,ss,vrel,coll_bro
      return
      end

