      function dens_bin(xxx,yyy,zzz)
*     returns the bin corresponding to the local density
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      include 'com_cont_epair'
      integer ix,iy,iz,dens_bin
      real*8 xxx,yyy,zzz, density
      ix = nint(xxx)
      iy = nint(yyy)
      iz = nint(zzz)
      density = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and.
     &       iabs(iz).le. maxz)
     &  density = sqrt(rhob_4(0,ix,iy,iz)**2-rhob_4(1,ix,iy,iz)**2
     &                -rhob_4(2,ix,iy,iz)**2-rhob_4(3,ix,iy,iz)**2)/rho0
      dens_bin = nint(density/densste+1.0)
      dens_bin = min(dens_bin,maxde)
      return
      end

************************************************************************
      subroutine piNcrossmodfact(xxx,yyy,zzz,px,py,pz,vmass,
     &     factrho,factome,factrhome)
*     returns the factors for changing the piNdilep crossections imedium
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      include 'cominput'
      real*8 factrho,factome,factrhome,xxx,yyy,zzz,px,py,pz,vmass
      real*8 deriv(0:3),j0,j1,j2,j3,dens,ener,vrel,ppx,ppy,ppz
      real*8 rselfrho,iselfrho,sgammarho,rselfome,iselfome,sgammaome
      real*8 gammarho0,factrhome0,factrhome1,betlrfx,betlrfy,betlrfz

      if(isplipi .eq. 1) then
        call splinint1(xxx,yyy,zzz,deriv)
        j0    = deriv(0)
        j1    = deriv(1)
        j2    = deriv(2)
        j3    = deriv(3)
      else if(isplipi .eq. 0) then
        call linint1(xxx,yyy,zzz,deriv)
        j0    = deriv(0)
        j1    = deriv(1)
        j2    = deriv(2)
        j3    = deriv(3)
      end if
      dens=sqrt(j0**2-j1**2-j2**2-j3**3)
      if(j0 .gt. 1.0e-6) then
        betlrfx = j1/j0
        betlrfy = j2/j0
        betlrfz = j3/j0
      else
        betlrfx = 0.0
        betlrfy = 0.0
        betlrfz = 0.0
      end if

c      write(*,*) 'piNdilepfactors',betlrfx,betlrfy,betlrfz,px,py,pz,dens
      ppx=px
      ppy=py
      ppz=pz
      ener=sqrt(vmass**2+ppx**2+ppy**2+ppz**2)
      call lorentz(betlrfx,betlrfy,betlrfz,ppx,ppy,ppz,ener)
      vrel = sqrt(ppx**2+ppy**2+ppz**2)/ener   ! vrel (vel. rel. to medium)

      if(vrel.ge.1.0)
     & write(*,*)'dilsupp',vrel,betlrfx,betlrfy,betlrfz,ppx,ppy,ppz,ener
      call self_omega(vrel, dens,vmass,rselfome,iselfome,sgammaome)
      call self_rho(vrel, dens,vmass,rselfrho,iselfrho,sgammarho)

c      write(*,*) 'facts',vmass,vrel,rselfome,iselfome,rselfrho,iselfrho

      factome = ((vmass**2-omass**2)**2+vmass**2*owidth**2)/
     &          ((vmass**2-omass**2-rselfome)**2+iselfome**2)
c      factome = ((vmass**2-omass**2)**2+vmass**2*owidth**2)
      gammarho0 = 0.0
      if(vmass.gt.2.0*pmass+0.001) 
     & gammarho0 =  rowidth * romas/vmass *
     &      ((vmass**2-4.*pmass**2)/(romas**2-4.*pmass**2))**1.5*
     &((.09+(0.25*romas**2-pmass**2))/(.09+(0.25*vmass**2-pmass**2)))**2
      factrho  = ((vmass**2-romas**2)**2+vmass**2*gammarho0**2)/
     &          ((vmass**2-romas**2-rselfrho)**2+iselfrho**2)
c      factrho  = ((vmass**2-romas**2)**2+vmass**2*gammarho0**2)

      factrhome0 = ((vmass**2-omass**2)*(vmass**2-romas**2)+
     &             vmass*owidth*vmass*gammarho0)/
     &           (((vmass**2-omass**2)**2+vmass**2*owidth**2)*
     &            ((vmass**2-romas**2)**2+vmass**2*gammarho0**2))
      factrhome1 = ((vmass**2-omass**2-rselfome)*
     &               (vmass**2-romas**2-rselfrho)+iselfome*iselfrho)/
     &            (((vmass**2-omass**2-rselfome)**2+iselfome**2)*
     &             ((vmass**2-romas**2-rselfrho)**2+iselfrho**2))
      factrhome = factrhome1/factrhome0
c      factrhome = 1.0/factrhome0
c      write(*,*) 'factors',gammarho0,factrho,factrhome0,factrhome1

      return
      end

************************************************************************
      function dgdm_bary_dalitz(idres,charge,srt,dilmass)
c-----------------------------------------------------------------------
      implicit none
      include 'common'

      integer idres,charge
      real*8 srt,dilmass,dgdm_bary_dalitz

      real*8 matrix_t,matrix_l
      integer spin,parity
      real*8 lambda

      real*8 g_em(1:nres,0:1)
      data ((g_em(idres,charge),charge=0,1),idres=1,nres)
c    &   /0.599, 0.599,         ! D1232
c    &   0.098, 0.139,          ! N1440
c    &   0.238, 0.238,          ! N1520
     &   /1.98, 1.980,         ! D1232
     &   0.098, 0.139,          ! N1440
     &   0.719, 0.793,          ! N1520
     &   0.549, 0.634,          ! N1535
     &   0.285, 0.315,          ! N1650
     &   1.52,  0.678,          ! N1675
     &   2.74,  0.971,          ! N1680
     &   0.193, 0.126,          ! N1700
     &   0.019, 0.030,          ! N1710
     &   0.386, 0.193,          ! N1720
     &   0.0,   0.0,            ! N2000
     &   0.0,   0.0,            ! N2080
     &   0.0,   0.0,            ! N2190
     &   0.0,   0.0,            ! N2220
     &   0.0,   0.0,            ! N2250
     &   0.162, 0.162,          ! D1600
     &   0.162, 0.162,          ! D1620
     &   0.549, 0.549,          ! D1700
     &   0.0,   0.0,            ! D1900
     &   0.713, 0.713,          ! D1905
     &   0.066, 0.066,          ! D1910
     &   0.0,   0.0,            ! D1920
     &   0.479, 0.479,          ! D1930
     &   0.0,   0.0/            ! D1950

      integer parity_res(1:nres)
      data (parity_res(idres),idres=1,nres)
     &   /+1,                   ! D1232
     &   +1,                    ! N1440
     &   -1,                    ! N1520
     &   -1,                    ! N1535
     &   -1,                    ! N1650
     &   -1,                    ! N1675
     &   +1,                    ! N1680
     &   -1,                    ! N1700
     &   +1,                    ! N1710
     &   +1,                    ! N1720
     &   +1,                    ! N2000
     &   -1,                    ! N2080
     &   -1,                    ! N2190
     &   +1,                    ! N2220
     &   -1,                    ! N2250
     &   +1,                    ! D1600
     &   -1,                    ! D1620
     &   -1,                    ! D1700
     &   -1,                    ! D1900
     &   +1,                    ! D1905
     &   +1,                    ! D1910
     &   +1,                    ! D1920
     &   -1,                    ! D1930
     &   +1/                    ! D1950

      if (charge.lt.0 .or. charge.gt.1) then
        dgdm_bary_dalitz = 0.
      else
        spin = resprop2(idres,3) ! 2*spin!
        parity = parity_res(idres)

        matrix_t = 0.
        matrix_l = 0.

        if (spin.eq.1) then
          if (parity.eq.+1) then
            matrix_t = 1./(2.*rmass**4) * (srt**2-rmass**2)**2 *
     &         ( (srt+rmass)**2 - dilmass**2 )
            matrix_l = dilmass**2/(2.*rmass**4) * (srt+rmass)**2 *
     &         ( (srt+rmass)**2 - dilmass**2 )
          else
            matrix_t = 1./(2.*rmass**4) * (srt**2-rmass**2)**2 *
     &         ( (srt-rmass)**2 - dilmass**2 )
            matrix_l = dilmass**2/(2.*rmass**4) * (srt-rmass)**2 *
     &         ( (srt-rmass)**2 - dilmass**2 )
          end if
        else if (spin.eq.3) then
          if (parity.eq.+1) then
            matrix_t = 1./(12.*srt**2*rmass**2) *
     &         ( (srt-rmass)**2 - dilmass**2 ) *
     &         (3.*srt**4 + 6.*srt**3*rmass + 4.*srt**2*rmass**2
     &         + 2.*srt*rmass**3 + rmass**4 - 2.*srt*rmass*dilmass**2
     &         - 2.*rmass**2*dilmass**2 + dilmass**4)
            matrix_l = dilmass**2/(3.*rmass**2) *
     &         ( (srt-rmass)**2 - dilmass**2 )
          else
            matrix_t = 1./(12.*srt**2*rmass**2) *
     &         ( (srt+rmass)**2 - dilmass**2 ) *
     &         (3.*srt**4 - 6.*srt**3*rmass + 4.*srt**2*rmass**2
     &         - 2.*srt*rmass**3 + rmass**4 + 2.*srt*rmass*dilmass**2
     &         - 2.*rmass**2*dilmass**2 + dilmass**4)
            matrix_l = dilmass**2/(3.*rmass**2) *
     &         ( (srt+rmass)**2 - dilmass**2 )
          end if
        else if (spin.eq.5) then
          if (parity.eq.+1) then
            matrix_t = 1./(480.*srt**4*rmass**4) *
     &         ( (srt-rmass)**2 - dilmass**2 ) *
     &         ( (srt+rmass)**2 - dilmass**2 )**2 *
     &         (2.*srt**4 - 4.*srt**3*rmass + 3.*srt**2*rmass**2
     &         - 2.*srt*rmass**3 + rmass**4 + 2.*srt*rmass*dilmass**2
     &         - 2.*rmass**2*dilmass**2 + dilmass**4)
            matrix_l = dilmass**2/(120.*srt**2*rmass**4) *
     &         ( (srt-rmass)**2 - dilmass**2 ) *
     &         ( (srt+rmass)**2 - dilmass**2 )**2
          else
            matrix_t = 1./(480.*srt**4*rmass**4) *
     &         ( (srt-rmass)**2 - dilmass**2 )**2 *
     &         ( (srt+rmass)**2 - dilmass**2 ) *
     &         (2.*srt**4 + 4.*srt**3*rmass + 3.*srt**2*rmass**2
     &         + 2.*srt*rmass**3 + rmass**4 - 2.*srt*rmass*dilmass**2
     &         - 2.*rmass**2*dilmass**2 + dilmass**4)
            matrix_l = dilmass**2/(120.*srt**2*rmass**4) *
     &         ( (srt-rmass)**2 - dilmass**2 )**2 *
     &         ( (srt+rmass)**2 - dilmass**2 )
          end if
        end if

        dgdm_bary_dalitz = alfa**2 * g_em(idres,charge)**2 *
     &     sqrt(lambda(srt**2,rmass**2,dilmass**2)) /
     &     (6.*pi*srt**3*dilmass) *
     &     (2.*matrix_t + matrix_l)
c        write(*,*) 'dgdm_bary_dalitz',idres,charge,parity,srt,dilmass,
c     &     dgdm_bary_dalitz
c        write(*,*) '  .....> ',alfa,g_em(idres,charge),
c     &     sqrt(lambda(srt**2,rmass**2,dilmass**2)),
c     &     1./(6.*pi*srt**3*dilmass),
c     &     matrix_t,matrix_l
      end if
      return
      end


************************************************************************
      function dgdm_omega_dalitz(m_omega,dilmass)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      real*8 dgdm_omega_dalitz,m_omega,dilmass

      real*8 cutoff
      parameter (cutoff=0.65)   ! GeV
      real*8 Gamma
      parameter (Gamma=0.075)  ! GeV
      real*8 Gamma_phot
      parameter (Gamma_phot=0.000757)
      real*8 formfac,lambda

c      formfac = 1./( (1.-dilmass**2/cutoff**2)**2 + Gamma**2/cutoff**2 )
      formfac = 1.
      dgdm_omega_dalitz = 2.*alfa/(3.*pi) * Gamma_phot/dilmass *
     &   ( lambda(m_omega**2,pmass**2,dilmass**2)
     &   / (m_omega**2-pmass**2)**2 )**1.5 * formfac
c      write(*,*) 'dgdm_omega_dalitz',m_omega,dilmass,dgdm_omega_dalitz
      return
      end

************************************************************************
      function dgdm_eta_dalitz(m_eta,dilmass)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      real*8 dgdm_eta_dalitz,dilmass,m_eta

      real*8 cutoff
      parameter (cutoff=0.77)   ! GeV
      real*8 Gamma_phot
      parameter (Gamma_phot=5.086e-7)
      real*8 formfac

      formfac = 1./(1.-dilmass**2/cutoff**2)**2
      dgdm_eta_dalitz = 4.*alfa/(3.*pi) * Gamma_phot/dilmass *
     &   (1. - dilmass**2/m_eta**2)**3 * formfac
c      write(*,*) 'dgdm_eta_dalitz',dilmass,formfac,dgdm_eta_dalitz
      return
      end

************************************************************************
      function dgdm_pion_dalitz(m_pion,dilmass)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      real*8 dgdm_pion_dalitz,dilmass,m_pion

      real*8 cutoff
      parameter (cutoff=5.5)   ! GeV^-2
      real*8 Gamma_phot
      parameter (Gamma_phot=8.e-9)
      real*8 formfac

      formfac = (1.+dilmass**2*cutoff)**2
      dgdm_pion_dalitz = 4.*alfa/(3.*pi) * Gamma_phot/dilmass *
     &   (1. - dilmass**2/m_pion**2)**3 * formfac
c      write(*,*) 'dgdm_pion_dalitz',dilmass,formfac,dgdm_pion_dalitz
      return
      end

************************************************************************
      real*8 function lambda(a,b,c)
c-----------------------------------------------------------------------
      implicit none
      real*8 a,b,c
      lambda = a**2 + b**2 + c**2 - 2.*(a*b + b*c + c*a)
      return
      end

************************************************************************
      subroutine test_acc
c     test the hades acceptance filter
      implicit none
      real*8 dilmass,pt_dil,rap_dil,pz_dil,e_dil,acc
      real*8 getHadesPairAcceptance
      integer ii,jj,kk

      do ii=1,6
        dilmass = float(ii)*0.2
        do jj=-4,4
          pz_dil = float(jj)*0.2
          do kk=-4,4
            pt_dil = float(kk)*0.2
            e_dil = sqrt(dilmass**2 + pz_dil**2 + pt_dil**2)
            rap_dil = 0.5 * log( (e_dil+pz_dil)/(e_dil-pz_dil) )
            acc = getHadesPairAcceptance(dilmass,pt_dil,rap_dil,0)
            write(*,*) '->',dilmass,pt_dil,rap_dil,acc
          end do
        end do
      end do
      return
      end

************************************************************************
      function mass_bin_num(vmass,i_reg)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 vmass
      integer mass_bin_num,ibin,i_reg

      mass_bin_num = 1
      if(ndlmas.le.1) return
c if vmass smaller than the lower value of the smallest mass bin
      if(vmass.lt.1.5*qy(0,1)-0.5*qy(0,nqt*ny*nf+1)) then
           mass_bin_num = 0
c if vmass larger than the upper value of the largest mass bin
      else if(vmass.gt.1.5*qy(0,(ndlmas-1)*nqt*ny*nf+1)-
     &           0.5*qy(0,(ndlmas-2)*nqt*ny*nf+1)) then
        mass_bin_num=ndlmas+1
      else
        do ibin = 2,ndlmas
          if(vmass.lt.qy(0,(ibin-1)*nqt*ny*nf+1)) goto 11
        end do
        ibin = ndlmas
 11     continue
        if(vmass.gt.0.5*qy(0,(ibin-1)*nqt*ny*nf+1)+
     &            0.5*qy(0,(ibin-2)*nqt*ny*nf+1)) then
          mass_bin_num=ibin
        else
          mass_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function JP_mass_bin_num(vmass)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 vmass
      integer JP_mass_bin_num,ibin,nqtyf

      nqtyf=JP_qt*JP_y*JP_f
      JP_mass_bin_num = 1
      if(JP_dlmas.le.1) return

c if vmass smaller than the lower value of the smallest mass bin
      if(vmass.lt.1.5*JP_qy(0,1)-0.5*JP_qy(0,nqtyf+1)) then
            JP_mass_bin_num = 0
c if vmass larger than the upper value of the largest mass bin
      else if(vmass.gt.1.5*JP_qy(0,(JP_dlmas-1)*nqtyf+1)-
     &           0.5*JP_qy(0,(JP_dlmas-2)*nqtyf+1)) then
         JP_mass_bin_num=JP_dlmas+1
      else
        do ibin = 2,JP_dlmas
          if(vmass.lt.JP_qy(0,(ibin-1)*nqtyf+1)) goto 11
        end do
        ibin = JP_dlmas
 11     continue
        if(vmass.gt.0.5*JP_qy(0,(ibin-1)*nqtyf+1)+
     &            0.5*JP_qy(0,(ibin-2)*nqtyf+1)) then
           JP_mass_bin_num=ibin
        else
           JP_mass_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function IR_mass_bin_num(vmass)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 vmass
      integer IR_mass_bin_num,ibin,nqtyf

      nqtyf=IR_qt*IR_y*IR_f
      IR_mass_bin_num = 1
      if(IR_dlmas.le.1) return

c if vmass smaller than the lower value of the smallest mass bin
      if(vmass.lt.1.5*IR_qy(0,1)-0.5*IR_qy(0,nqtyf+1)) then
            IR_mass_bin_num = 0
c if vmass larger than the upper value of the largest mass bin
      else if(vmass.gt.1.5*IR_qy(0,(IR_dlmas-1)*nqtyf+1)-
     &           0.5*IR_qy(0,(IR_dlmas-2)*nqtyf+1)) then
         IR_mass_bin_num=IR_dlmas+1
      else
        do ibin = 2,IR_dlmas
          if(vmass.lt.IR_qy(0,(ibin-1)*nqtyf+1)) goto 11
        end do
        ibin = IR_dlmas
 11     continue
        if(vmass.gt.0.5*IR_qy(0,(ibin-1)*nqtyf+1)+
     &            0.5*IR_qy(0,(ibin-2)*nqtyf+1)) then
           IR_mass_bin_num=ibin
        else
           IR_mass_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function mass_bin_size(imass)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 mass_bin_size
      integer imass

      if(imass.eq.1) then
        mass_bin_size = qy(0,nqt*ny*nf+1)-qy(0,1)
      else if(imass.eq.ndlmas) then
        mass_bin_size=
     &      qy(0,(ndlmas-1)*nqt*ny*nf+1)-qy(0,(ndlmas-2)*nqt*ny*nf+1)
      else
        mass_bin_size=
     &    0.5*(qy(0,imass*nqt*ny*nf+1)-qy(0,(imass-2)*nqt*ny*nf+1))
      end if
      return
      end

************************************************************************
      function qt_bin_num(qt)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 qt
      integer qt_bin_num,ibin

c if qt smaller than the lower value of the smallest qt bin
      if(qt.lt.1.5*qy(1,1)-0.5*qy(1,ny*nf+1)) then
             qt_bin_num = 0
c if qt larger than the upper value of the largest qt bin
      else if(qt.gt.1.5*qy(1,(nqt-1)*ny*nf+1)-
     &         0.5*qy(1,(nqt-2)*ny*nf+1)) then
              qt_bin_num=nqt+1
      else
        do ibin = 2,nqt
          if(qt.lt.qy(1,(ibin-1)*ny*nf+1)) goto 11
        end do
        ibin = nqt
 11     continue
        if(qt.gt.0.5*qy(1,(ibin-1)*ny*nf+1)+
     &         0.5*qy(1,(ibin-2)*ny*nf+1)) then
          qt_bin_num=ibin
        else
          qt_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function JP_qt_bin_num(qt)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 qt
      integer JP_qt_bin_num,ibin,nyf

      nyf = JP_y*JP_f
c if qt smaller than the lower value of the smallest qt bin
      if(qt.lt.1.5*JP_qy(1,1)-0.5*JP_qy(1,nyf+1)) then
             JP_qt_bin_num = 0
c if qt larger than the upper value of the largest qt bin
      else if(qt.gt.1.5*JP_qy(1,(JP_qt-1)*nyf+1)-
     &         0.5*JP_qy(1,(JP_qt-2)*nyf+1)) then
              JP_qt_bin_num=JP_qt+1
      else
        do ibin = 2,JP_qt
          if(qt.lt.JP_qy(1,(ibin-1)*nyf+1)) goto 11
        end do
        ibin = JP_qt
 11     continue
        if(qt.gt.0.5*JP_qy(1,(ibin-1)*nyf+1)+
     &         0.5*JP_qy(1,(ibin-2)*nyf+1)) then
          JP_qt_bin_num=ibin
        else
          JP_qt_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function IR_qt_bin_num(qt)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 qt
      integer IR_qt_bin_num,ibin,nyf

      nyf = IR_y*IR_f
c if qt smaller than the lower value of the smallest qt bin
      if(qt.lt.1.5*IR_qy(1,1)-0.5*IR_qy(1,nyf+1)) then
             IR_qt_bin_num = 0
c if qt larger than the upper value of the largest qt bin
      else if(qt.gt.1.5*IR_qy(1,(IR_qt-1)*nyf+1)-
     &         0.5*IR_qy(1,(IR_qt-2)*nyf+1)) then
              IR_qt_bin_num=IR_qt+1
      else
        do ibin = 2,IR_qt
          if(qt.lt.IR_qy(1,(ibin-1)*nyf+1)) goto 11
        end do
        ibin = IR_qt
 11     continue
        if(qt.gt.0.5*IR_qy(1,(ibin-1)*nyf+1)+
     &         0.5*IR_qy(1,(ibin-2)*nyf+1)) then
          IR_qt_bin_num=ibin
        else
          IR_qt_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function qt_bin_size(iqt)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 qt_bin_size
      integer iqt

      if(iqt.eq.1) then
        qt_bin_size = qy(1,ny*nf+1)-qy(1,1)
      else if(iqt.eq.nqt) then
        qt_bin_size=
     &      qy(1,(nqt-1)*ny*nf+1)-qy(1,(nqt-2)*ny*nf+1)
      else
        qt_bin_size=
     &      0.5*(qy(1,iqt*ny*nf+1)-qy(1,(iqt-2)*ny*nf+1))
      end if
      return
      end

************************************************************************
      function y_bin_num(rap)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 rap
      integer y_bin_num,ibin

c if rap smaller than the lower value of the smallest rapidity bin
      if(rap.lt.1.5*qy(2,1)-0.5*qy(2,nf+1)) then
         y_bin_num = 0
c if rap larger than the upper value of the largest rapidity bin
      else if(rap.gt.1.5*qy(2,(ny-1)*nf+1)-
     &          0.5*qy(2,(ny-2)*nf+1)) then 
          y_bin_num=ny+1
      else
        do ibin = 2,ny
          if(rap.lt.qy(2,(ibin-1)*nf+1)) goto 11
        end do
        ibin = ny
 11     continue
        if(rap.gt.0.5*qy(2,(ibin-1)*nf+1)+
     &          0.5*qy(2,(ibin-2)*nf+1)) then
          y_bin_num = ibin
        else
          y_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function JP_y_bin_num(rap)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 rap
      integer JP_y_bin_num,ibin

c if rap smaller than the lower value of the smallest rapidity bin
      if(rap.lt.1.5*JP_qy(2,1)-0.5*JP_qy(2,JP_f+1)) then
         JP_y_bin_num = 0
c if rap larger than the upper value of the largest rapidity bin
      else if(rap.gt.1.5*JP_qy(2,(JP_y-1)*JP_f+1)-
     &          0.5*JP_qy(2,(JP_y-2)*JP_f+1)) then 
          JP_y_bin_num=JP_y+1
      else
        do ibin = 2,JP_y
          if(rap.lt.JP_qy(2,(ibin-1)*JP_f+1)) goto 11
        end do
        ibin = JP_y
 11     continue
        if(rap.gt.0.5*JP_qy(2,(ibin-1)*JP_f+1)+
     &          0.5*JP_qy(2,(ibin-2)*JP_f+1)) then
          JP_y_bin_num = ibin
        else
          JP_y_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function IR_y_bin_num(rap)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 rap
      integer IR_y_bin_num,ibin

c if rap smaller than the lower value of the smallest rapidity bin
      if(rap.lt.1.5*IR_qy(2,1)-0.5*IR_qy(2,JP_f+1)) then
         IR_y_bin_num = 0
c if rap larger than the upper value of the largest rapidity bin
      else if(rap.gt.1.5*IR_qy(2,(IR_y-1)*IR_f+1)-
     &               0.5*IR_qy(2,(IR_y-2)*IR_f+1)) then 
          IR_y_bin_num=IR_y+1
      else
        do ibin = 2,IR_y
          if(rap.lt.IR_qy(2,(ibin-1)*IR_f+1)) goto 11
        end do
        ibin = IR_y
 11     continue
        if(rap.gt.0.5*IR_qy(2,(ibin-1)*IR_f+1)+
     &            0.5*IR_qy(2,(ibin-2)*IR_f+1)) then
          IR_y_bin_num = ibin
        else
          IR_y_bin_num = ibin -1 
        end if
      end if
      return
      end

************************************************************************
      function y_bin_size(iy)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      real*8 y_bin_size
      integer iy

      if(iy.eq.1) then
        y_bin_size = qy(2,nf+1)-qy(2,1)
      else if(iy.eq.ny) then
        y_bin_size=
     &      qy(2,(ny-1)*nf+1)-qy(2,(ny-2)*nf+1)
      else
        y_bin_size=
     &      0.5*(qy(2,iy*nf+1)-qy(2,(iy-2)*nf+1))
      end if
      return
      end

