      function dist_mes(vrel,idm,density,mmes,mres)
      implicit none

      integer idm
      real*8 vrel,density,mmes,mres,dist_mes
      real*8 bb, pp, ppp, spect_mes

      include "common"

c     write(*,*) ' in dist_mes ', idm, mres
      pp = (mres**2 - (rmass+mmes)**2) * (mres**2 - (rmass-mmes)**2)
      ppp = sqrt(pp/4.)
c     bb = mres**2 + rmass**2 - mmes**2

      dist_mes = ppp * spect_mes(vrel,idm,density,mmes)

      return
      end



      function spect_mes(vrel,idm,density,mmes)
      implicit none

      include "common"

      integer idm
      real*8 vrel,density,mmes,spect_mes
      real*8 rself,iself,mpole, sgam

      mpole = 0.0
      if (idm.eq.3) then
        mpole = romas
        call self_rho(vrel,density,mmes,rself,iself,sgam)
      else if (idm.eq.5) then
        mpole = omass
        call self_omega(vrel,density,mmes,rself,iself,sgam)
      end if

      spect_mes = - 2.*mmes * iself/pi /
     &   ((mmes**2 - mpole**2 - rself)**2 + iself**2)

      return
      end


      subroutine self_rho(vrel,dens,mrho,rself,iself,sgamma)
      implicit none

      include "common"
      include "cominput"

      real*8 totom(7,2000),totro(7,2000)
      COMMON/self_stored/ totom,totro

      real*8 density,dens,mrho,rself,iself, vrel, coll_bro
      real*8 sgamma, totro6_lim, bg, cg, dgdm
      real*8 delta,aaa, totro6, eps, mrho1, mrho2, corr
      integer ii, idm

      eps = .0005
c     mrho_lim  = 0.38
      mrho_lim  = 0.427
      g_col_null = .02
      gam_null  = 0.01
c     rself = dens * 0.5 * atan((mrho-1.3)/0.2)

      rself = .0
      if(ivecmatt.eq.1) then
        rself = dens / rho0 * rhomshift ! rhomshift masshift
        rself = 2.*romas*rself   ! + rself**2
      endif

c-----------------------------------
c                       that's the original version
      IF (icbro .lt. 2) THEN
        if (mrho .gt. 2.*pmass+ eps) then
          mrho1 = mrho
          totro6 = romas * rowidth *
     &     ((mrho1**2-4.*pmass**2)/(romas**2-4.*pmass**2))**1.5 *
     &     romas/mrho1*
     &((.09+(.25*romas**2-pmass**2))/(.09+(.25*mrho1**2-pmass**2)))**2
        else
          mrho1 = 2.*pmass+ eps
          totro6 = .001
        endif
        sgamma = 2.*totro6 / (romas+mrho)
c
      ELSEIF (icbro .eq. 2) THEN
c-----------------------------------
c              now the replacement
c----------------------
        mrho1 = mrho
        if (mrho .gt. .380 ) then
          totro6 = romas * rowidth *
     &     ((mrho1**2-4.*pmass**2)/(romas**2-4.*pmass**2))**1.5 *
     &     romas/mrho1*
     &((.09+(.25*romas**2-pmass**2))/(.09+(.25*mrho1**2-pmass**2)))**2
        else
          mrho2 = 0.38
          totro6 = romas * rowidth *
     &      ((mrho2**2-4.*pmass**2)/(romas**2-4.*pmass**2))**1.5 *
     &      romas/mrho2
          totro6 = totro6 * (mrho/mrho2)**1.5
          if (mrho .lt.  2.*pmass+ eps) mrho1 = 2.*pmass+ eps
        endif
        sgamma = 2.*totro6 / (romas+mrho)
c
      ELSE
c
        mrho1 = mrho
        if (mrho .gt. mrho_lim ) then
          totro6 = romas * rowidth *
     &     ((mrho1**2-4.*pmass**2)/(romas**2-4.*pmass**2))**1.5 *
     &     romas/mrho1
          sgamma = 2.*totro6 / (romas+mrho1)
        else
          mrho2 = mrho_lim
          totro6_lim = romas * rowidth *
     &      ((mrho2**2-4.*pmass**2)/(romas**2-4.*pmass**2))**1.5 *
     &      romas/mrho2
          dgdm = totro6_lim*(2.0*mrho2**2+4.*pmass**2) /
     1           (mrho2*(mrho2**2-4.*pmass**2))
          cg   = dgdm / (totro6_lim - gam_null)
          bg   = dgdm / (cg * exp(cg*mrho2))
          totro6 = gam_null + bg*exp(cg*mrho)
          sgamma = 2.*totro6 / (romas+mrho2)
          if (mrho .lt.  2.*pmass+ eps) mrho1 = 2.*pmass+ eps
        endif
      ENDIF
c-------------------------------------
c
      idm = 3
c     write(*,*)  ' in self_rho-  before broadening ', mrho
      coll_bro = .0
      if(icbro.ge.1) call broadening(vrel,dens,idm, mrho1, coll_bro)
c      write(*,*)   ' in rho_wid:  after broadening ',
c     1             vrel,dens,idm, mrho1, coll_bro
c      write(*,*)  'rho_wid: icbro ', icbro
      corr = 1.
      if (icbro .eq. 3) corr = dens/(rho0+dens)+0.001
      iself = - totro6*corr -  romas * coll_bro
      if (icbro.eq.3) iself = iself - dens/rho0 * romas *
     1                                 g_col_null
c


c to be used with Gyuri's tables (amplread must be called!)
c      subroutine self_rho(0.,dens,mrho,rself,iself)
c      implicit none
c
c      include "common"
c
c      real*8 totom(7,2000),totro(7,2000)
c      COMMON/self_stored/ totom,totro
c
c      real*8 density,dens,mrho,rself,iself
c      real*8 delta
c      integer ii
c
c      density = 0.   !!!!!!!!!!!
c      rself = 0.
c      iself = .117  ! !!!
c
c      ii = int(1000.0*mrho-141.0) + 1
c      delta = 1000.0*mrho-141.0 + 1.0 - ii
c
c      if (ii.gt.999 .or. ii.lt.1) ii = 999
c
c      rself = - (density*((1.-delta)*totro(3,ii) + delta*totro(3,ii+1))
c     &   + (1.-delta)*totro(5,ii) + delta*totro(5,ii+1))
c      iself = - (density*((1.-delta)*totro(4,ii) + delta*totro(4,ii+1))
c     &   + (1.-delta)*totro(6,ii) + delta*totro(6,ii+1))
c      return
c      end
      return
      end


      subroutine print_spect
      implicit none

      include "common"

      real*8 mmes,spect_rho,spect_omega,spect_mes,vrel,density
      real*8 rself, rself2, iself, iself2, sgam
      real*8 vrel0, a1,a2,a3,a4,a5
      vrel = 0.596
      vrel0 = .0
      density = 0.32
      write(mspfpri,*)'# spectral function of rho, omega at dens.=0.16'
      write(mspfpri,*) '# mass       rho(re, im, A) ',
     1           '         omega(re, im, A)'
c     do mmes = 0.3, 2.0, 0.025
c       call self_rho(vrel,density,mmes,rself,iself,sgam)
c       spect_rho = spect_mes(vrel,3,density,mmes)
c       call self_omega(vrel,density,mmes,rself2,iself2, sgam)
c       spect_omega = spect_mes(vrel,5,density,mmes)
c       write(mspfpri,'(f8.3,3x,2f8.4,f10.5,2f8.4,f10.5)')
c    1   mmes, rself,iself,spect_rho, rself2,iself2,spect_omega
c     end do
c     do mmes = 0.6, 0.85, 0.002    !     0.002
c       spect_rho = spect_mes(vrel,3,density,mmes)
c       call self_omega(vrel,density,mmes,rself2,iself2, sgam)
c       a5 = 3.14/mmes
c       a1   = a5 * spect_mes(vrel,5,.600,mmes)
c       a2   = a5 * spect_mes(vrel,5,.300,mmes)
c       a3   = a5 * spect_mes(vrel0,5,.0,mmes)
c       a4   = spect_mes(vrel,5,0.3,mmes)
c       write(mspfpri,'(f8.3,3x,5f10.5)')
c    1  mmes, a1,a2,a3
c     end do
c     do mmes = 0.02, 1.00, 0.002    !     0.02
c       call self_rho(vrel,density,mmes,rself,iself,sgam)
c       a4 = sgam*0.5*(mmes+romas)
c       a5 = 3.14/2./mmes
c       a1   = a5 * spect_mes(vrel,3,.599,mmes)
c       a2   = a5 * spect_mes(vrel,3,.320,mmes)
c       a3   = a5 * spect_mes(vrel0,3,.0,mmes)
c       write(mspfpri,'(f8.3,3x,3f10.5,5x,3f10.6)')
c    1  mmes, a1,a2,a3, rself,iself, a4
c     end do
c     write(*,*) '  start self_rho  '
c     call self_rho(vrel0,0.0,1.4193182,rself,iself,sgam)
c
      return
      end
