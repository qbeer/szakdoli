***********************************************************************
*                                                                      *
*                                                                      *
*       purpose:    calculating K_minus  production from:              *
*                                n+n  collision      (entry kminu_dbb) *
*                                n+pi collision      (entry kminu_dpi) *
*                                n+pi collision      (entry kminu_hpi) *
*                                perturbative, only the K- is stored   *
*--------------------------------------------------------------------  *
*                   initialization by subroutine phi_in                *
*--------------------------------------------------------------------  *
*                   final cross section by entry kaonout               *
*--------------------------------------------------------------------  *
*                                                                      *
*         yref    - -1 * target rapidity in the frame                  *
*         nyk     - number of rapidities               (integer,input) *
*         nqtk    - number of transverse momenta       (integer,input) *
*         nfk     - number of angles                   (integer,input) *
*                                                                      *
*         ireac: 1-> nn      ->  nn K+K-                               *
*                2-> nD      ->  nn K+K-                               *
*                3-> DD      ->  nn K+K-                               *
*                4-> RR      ->  nn K+K-                               *
*                5-> n pi    ->  n  K+K-                               *
*                6-> D pi    ->  n  K+K-                               *
*                7-> R pi    ->  n  K+K-                               *
*                8-> L pi    ->  n  K-                                 *
*                9-> S pi    ->  n  K-                                 *
*                                                                      *
************************************************************************
      subroutine kminu_dbb(id1,id2,beta,gamma,srt,xxx,yyy,zzz,sig0,irun
     &  ,i1,i2)
*       variables:                                                     *
*         ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*----------------------------------------------------------------------*
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*----------------------------------------------------------------------*
      integer irun, ireac, inkrun, ink
      integer i1, i2,id2,id1, ntag
      integer inkmin, kl, ich
      integer  bin_dens,bin_time,ix,iy,iz
      real*8 srt, xxx,yyy,zzz, s,s0
      real*8 szig, szigma, sig0
      real*8 xx,yy,zz
      real*8 gamma, phase, density, denst
      real*8 pxx,pxy,pxz, vx,vy,vz, vxx,vyy,vzz, tmas_m, tmas_p
      real*8 mmass, srtmin, pmax, pmax2, amass, traf, sprim
      real*8 valkapi, ranp, pka, rr, pnupr2, pnupr
      real*8 us(0:3), qqk(0:3),pn3(0:3),ps3(0:3),ps2(0:3),pn2(0:3)
      real*8 pn23(0:3), pn2c(0:3), pn3c(0:3),pn2s(0:3),pn3s(0:3)
      real*8 qqks(0:3)
      real*8 pnup1, pnup2, s20, s30
      real*8 rn
*----------------------------------------------------------------------*
      real*8  beta(3)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
c
c      write(*,*)  ' start kminu_dbb ',id1,id2,beta,gamma,
c     +            srt,xxx,yyy,zzz,sig0,irun,i1,i2
      if(id1.gt.nres+1 .or. id2.gt.nres+1)                      return
      ireac = 4
      if(id1.eq.1 .and. id2.eq.1) ireac = 1
      if(id1*id2.eq.2)            ireac = 2
      if(id1.eq.2 .and. id2.eq.2) ireac = 3
cc
      tmas_m = xkmas
      tmas_p = xkmas
      us(0) = gamma
      us(1) = gamma * beta(1)
      us(2) = gamma * beta(2)
      us(3) = gamma * beta(3)
c      write(*,*)  ' us - beta ',us,beta
      pxx = us(1) * xkmas
      pxy = us(2) * xkmas
      pxz = us(3) * xkmas
      density  = .0
      if (i_kminu_pot  .gt. 0) then
         ich = -1
c        if (nthw  .ge. 15)
c    1   write(*,*)  ' in kminu_dbb vor  kao2 ', pxx, pxy, pxz
         call gradukao2(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  density,ich, vx,vy,vz,vxx,vyy,vzz)
         tmas_m = mmass
      endif
      if (ikaonpot .gt. 0) then
         ich =  1
!          call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
!      1                  density, ich, vx,vy,vz,vxx,vyy,vzz)	HS!
         call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  ich, vx,vy,vz,vxx,vyy,vzz)
         tmas_p = mmass
      endif
c      write(*,*)  ' kminu_dbb - srt ',srt,
c     1             density, rmass , tmas_m , tmas_p
      srtmin = 2.*rmass + tmas_m + tmas_p
      s0  = srtmin**2
      if (srt .lt. srtmin+.0001)                   return
      s     = srt**2
      amass = 2.*rmass + tmas_p
      pmax2=.25*(s-(amass+tmas_m)**2)*(s-(amass-tmas_m)**2)/s
c     pmax  = the maximal kminus  momentum
      pmax  = sqrt(pmax2)
      traf  = gamma / (gamma+1.0)
***                                                                  ***
      sprim = srt - srtmin
      if (i_kminu_cr .eq. 0)
     1 szigma = exp( 0.54+2.28*log(sprim)-3.03*log(3.0+sqrt(sprim))) ! hw
      if (i_kminu_cr .eq. 1)
     1 szigma = 0.3*(1.-s0/s)**3*(s/s0)**0.9 ! sibirtsev  nucl-th 9805021
      szig = szigma / sig0
c      write(*,*) '   kminu_bb szig ', srt, sprim, szig, pmax
c***********************************************************************
      inkrun = (max_kminu/num)
      ink = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 10 kl=1,inkrun                    !        iphi_bb       hw
        ink = ink + 1
        if (nx_kminu(0,ink) .eq. 0) goto 12
        if (valkapi .gt. p_kminu(4,ink)) then
           inkmin = ink
           valkapi  = p_kminu(4,ink)
        endif
  10  continue
      ink = inkmin
      if (szig .lt. valkapi) goto 99
  12  continue
c      write(*,*)   '  kminus found in kminu_dbb',ink
   15 continue
          ranp     = rn(iseed)
          xx = ranp**2 * (1.00-ranp**2)**2 / 0.1482
          if(xx .lt. rn(iseed)) go to 15
          pka       = pmax*ranp
c         if (ink .eq. 3197) write(*,*) ' in BB kminus', ink,
c    1                    s, pmax, ranp, amass, tmas_m, tmas_p
c
  20      xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 20
          qqk(0)  = sqrt(pka**2 + tmas_m**2)
          qqk(1)  = pka * xx / rr
          qqk(2)  = pka * yy / rr
          qqk(3)  = pka * zz / rr                   !    K- momentum
c            nucleon
          ps3(0) = sqrt(amass**2+pka**2)
          ps3(1) = -qqk(1)
          ps3(2) = -qqk(2)
          ps3(3) = -qqk(3)
          s30 = ps3(0)**2 - ps3(1)**2 - ps3(2)**2 - ps3(3)**2
c          write(*,*)  '  ps3, s30 ',ps3,s30,pka,amass
          amass = tmas_p + rmass
          pnupr2=.25*(s30-(amass+rmass)**2)*(s30-(amass-rmass)**2)/s30
          if (pnupr2 .lt. -.001)
     1      write(*,*) 'kinemat1 error in kminu_dbb' , pnupr2
          pnupr2 = max(.0, pnupr2)
          pnupr    = sqrt(pnupr2)
  115  continue
       ranp     = rn(iseed)
       xx = ranp**2 * sqrt(1.00-ranp**2) / 0.385
       if(xx .lt. rn(iseed)) go to 115
       pnup1 = ranp * pnupr
  30      xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 30
*   p3:                    nucleon momentum in colliding  system i1,i2
          pn3(1)    =  pnup1*xx/rr
          pn3(2)    =  pnup1*yy/rr
          pn3(3)    =  pnup1*zz/rr
          pn3(0)    =  sqrt(rmass**2 + pnup1**2)
          ps2(0)    =  sqrt(amass**2 + pnup1**2)
          ps2(1)    = - pn3(1)
          ps2(2)    = - pn3(2)
          ps2(3)    = - pn3(3)
          s20 = ps2(0)**2 - ps2(1)**2 - ps2(2)**2 - ps2(3)**2
          pnup2=.25*(s20-(tmas_p+rmass)**2)*(s20-(tmas_p-rmass)**2)/s20
          if (pnup2 .lt. -.001)
     1     write(*,*) 'kinemat2 error in kminu_dbb', pnup2
          pnup2 = max(.0, pnup2)
c         write(*,*)  ' pnup2 ',pnup2, amass, rmass
          pnup2 = sqrt(pnup2)
  31      xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 31
          pn2(1)    =  pnup2*xx/rr   ! the second nucleon
          pn2(2)    =  pnup2*yy/rr
          pn2(3)    =  pnup2*zz/rr
          pn2(0)    = sqrt(rmass**2+pn2(1)**2+pn2(2)**2+pn2(3)**2)
          call lorentz_hw(ps2, pn2, pn23)
          call lorentz_hw(ps3, pn23, pn2c)
          call lorentz_hw(ps3, pn3 , pn3c)
          call lorentz_hw(us , pn2c, pn2s)
          call lorentz_hw(us , pn3c, pn3s)
          call lorentz_hw(us , qqk , qqks)
c          if (ink .eq. 3197)  then
c           write(*,*)  '  us   ',us
c           write(*,*)  '  ps2  ',ps2
c           write(*,*)  '  ps3  ',ps3
c           write(*,*)  '  pn2  ',pn2
c           write(*,*)  '  pn3  ',pn3
c           write(*,*)  '  pn23 ',pn23
c           write(*,*)  '  pn2c ',pn2c
c           write(*,*)  '  pn3c ',pn3c
c           write(*,*)  '  pn2s ',pn2s
c           write(*,*)  '  pn3s ',pn3s
c           write(*,*)  '  qqks ',qqks
c           write(*,*)  '  qqk  ',qqk
c         endif
          if(ipauli.eq.1) call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,
     1                        pn2s(1),pn2s(2),pn2s(3))
          szig = szig * (1.0- phase)
          if(ipauli.eq.1) call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,
     1                        pn3s(1),pn3s(2),pn3s(3))
          szig = szig * (1.0- phase)
c
c          write(*,*)'kminu_dbb e',ink,qqks(1),qqks(2),qqks(3),szig,ireac
          p_kminu(0,ink)= density
          p_kminu(1,ink)= qqks(1)
          p_kminu(2,ink)= qqks(2)
          p_kminu(3,ink)= qqks(3)
          p_kminu(4,ink)= szig
          nx_kminu(5,ink)= 1001 * nint(200+zzz)
c         nx_kminu(5,ink)= 1001 * nthw
c         write(*,*) ' nx_kminu5 dbb',ink,zzz, nx_kminu(5,ink)
          r_kminu(1,ink)= xxx
          r_kminu(2,ink)= yyy
          r_kminu(3,ink)= zzz
          nx_kminu(0,ink) = 1
          nx_kminu(1,ink) = ireac
          nx_kminu(2,ink) = i1
          nx_kminu(3,ink) = i2
c         if (szig .lt. .0) then
c            write(*,*)  '  kminus negativ dbb', ink, ireac
c            stop
c         endif
c-----------------------------------------------------------------------
c           density dependence
          ix = nint(xxx)
          iy = nint(yyy)
          iz = nint(zzz)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
          bin_dens = int(denst*10.)
          bin_time = int(time*1.)
c          kmi_dbb(bin_dens,bin_time)=kmi_dbb(bin_dens,bin_time)+
c     1   p_kminu(4,ink)/real(num*isubs)
c      if(ink .eq. 20039) write(*,*) 'kminus found in BB coll',
c     1     irun,ink,ireac,szig,phase, qqks
   99 continue
      return
      end
************************************************************************
      subroutine kminu_dpi(ede,pdx,pdy,pdz,srt,xxx,yyy,zzz,szig0,i1,id1,
     +     i2,id2,iz1,iz2,irun)
*       variables:    1 = pion       2 = baryon                        *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*----------------------------------------------------------------------*
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*----------------------------------------------------------------------*
      integer irun, ireac, inkrun, ink
      integer i1, i2,id2,iz1,iz2,id1
      integer inkmin, kl, ich
      integer  bin_dens,bin_time,ix,iy,iz
      real*8 ede,pdx,pdy,pdz,srt, xxx,yyy,zzz, s,s0, srt0
      real*8 szig0, szig, szigr
      real*8 xx,yy,zz, qqabs, qqx,qqy,qqz, betax,betay,betaz
      real*8 density, denst
      real*8 p00,pxx,pxy,pxz, vx,vy,vz, vxx,vyy,vzz, tmas_m, tmas_p
      real*8 mmass, traf
      real*8 valkapi, ranp, pka, rr
      real*8 pbeta, transf, q0, qq2, gammx, rnaa, xrest2
      real*8 rn
*----------------------------------------------------------------------*
c      write(*,*) 'in kminu_dpi ', ede,pdx,pdy,pdz,srt,xxx,yyy,zzz,
c     +      szig0,i1,id1, i2, iz1,iz2,irun
      if(id1.ne.1)                      return
c      if(id2.ne.1)                      return
! 	ireac = 10
c      if(id2.gt.2)                      return
      ireac = 7
      if(id2.eq.1) ireac = 5
      if(id2.eq.1) ireac = 6

      tmas_m = xkmas
      tmas_p = xkmas
      p00 = sqrt(ede**2 - pdx**2 -pdy**2 -pdz**2)
      pxx = pdx / p00 * xkmas
      pxy = pdy / p00 * xkmas
      pxz = pdz / p00 * xkmas
      density = .0
      if (i_kminu_pot .gt. 0) then
         ich = -1
         call gradukao2(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  density, ich, vx,vy,vz,vxx,vyy,vzz)
         tmas_m = mmass
      endif
      if (ikaonpot  .gt. 0) then
         ich =  1
!          call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
!      1                  density, ich, vx,vy,vz,vxx,vyy,vzz)	HS!
         call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  ich, vx,vy,vz,vxx,vyy,vzz)
         tmas_p = mmass
      endif
c      if (abs (nthw -31) .lt. 1)
c     1 write(*,*) ' kminu_dpi - mass', density, tmas_m,tmas_p
      srt0  = rmass + tmas_m + tmas_p
      s0    = srt0**2
      if(srt.le.srt0)                                      return
      xrest2=  (rmass + tmas_p)**2
      s     = srt**2
c       q0  : energy of the K-
      q0      = 0.5 * (s + tmas_m**2 - xrest2) / srt
      if(q0 .le. tmas_m)                                    return
      qq2     = q0**2 - tmas_m**2
      qqabs   = sqrt(qq2)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gammx = ede / srt
      traf  = gammx / (gammx+1.0)
*----------------------------------------------------------------------*
c     write(*,*)   '  in  kminu_dpi   gammx ', gammx,traf
c...  cross sections:
      if (id2.eq.1) szig = 1.45 * (1. - s0/s)**2 * (s/s0)**.26 ! sibirtsev, nucl-th 9805021
      if (id2.eq.2) szig = 0.75*1.45 * (1. - s0/s)**2 * (s/s0)**.26
      if (id2.gt.2) szig = 1.45 * (1. - s0/s)**2 * (s/s0)**.26 ! piN másolva
      szigr  = szig/szig0
c***********************************************************************
c      write(*,*)' in kminu_dpi  szig', szig, ede, betax, betay, betaz
      inkrun = (max_kminu/num)
      ink = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 810 kl=1,inkrun                    !        iphi_bb       hw
        ink = ink + 1
        if (nx_kminu(0,ink) .eq. 0) goto 812
        if (valkapi .gt. p_kminu(4,ink)) then
           inkmin = ink
           valkapi  = p_kminu(4,ink)
        endif
  810 continue
      ink = inkmin
      if (szigr .lt. valkapi) goto 898
  812 continue
c     write(*,*)  '  in kminu_dpi  szig 2 ', ink, valkapi
  817  continue
          ranp     = rn(iseed)
          rr      = ranp*ranp
          rnaa = 2.598 * rr * sqrt(1.-rr)
          if(rnaa .lt. rn(iseed)) go to  817
          pka       = qqabs*ranp
  860 continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 860
c         write(*,*) 'kminus_dpi',szig,traf,q0, ink, irun, kl
          qqx= xx*pka/rr
          qqy= yy*pka/rr
          qqz= zz*pka/rr
*   kminu_ momentum in observable system      PAULI
          pbeta  = betax*qqx + betay*qqy + betaz*qqz
          transf = gammx * (traf * pbeta + q0)
          p_kminu(0,ink)= density
          p_kminu(1,ink)= qqx + betax * transf
          p_kminu(2,ink)= qqy + betay * transf
          p_kminu(3,ink)= qqz + betaz * transf
          p_kminu(4,ink)= szigr
          nx_kminu(5,ink)= 1001 * nint(200+zzz)
c         nx_kminu(5,ink)= 1001 * nthw
          r_kminu(1,ink)= xxx   !HS war auskommentiert??!
          r_kminu(2,ink)= yyy
          r_kminu(3,ink)= zzz
          nx_kminu(0,ink) = 1
          nx_kminu(1,ink) = ireac
          nx_kminu(2,ink) = i2
          nx_kminu(3,ink) = 0
c-----------------------------------------------------------------------
c           density dependence
          ix = nint(xxx)
          iy = nint(yyy)
          iz = nint(zzz)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
        bin_dens = int(denst*10.)
        bin_time = int(time*1.)
c        kmi_dpi(bin_dens,bin_time)=kmi_dpi(bin_dens,bin_time)+
c     &   p_kminu(4,ink)/real(num*isubs)
*
c      if(ink .eq. 20039) write(*,*) 'kminus found in piB coll',
c     1     irun,ink,ireac
c      write(*,*) 'kminus found in kminu_dpi', irun, ink, ireac, szigr,
c    1                qqabs,p_kminu(1,ink),p_kminu(2,ink),p_kminu(3,ink)
  898 continue
c************************
c       write(*,*) 'end of kminu_dpi'
c         if (szigr .lt. .0) then
c            write(*,*)  '  kminus negativ in dpi ', ink,ireac
c            stop
c         endif
      return
      end
************************************************************************
      subroutine kminu_hpi(ede,pdx,pdy,pdz,srt_0,xxx,yyy,zzz,szig0,
     1                i1,id1, i2, ihyp,iz1,iz2,irun)
*       variables:    1 = pion       2 = hyperon                       *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*     cross-sections are taken from Chung et al., Phys. lett. B401(97)1*
*----------------------------------------------------------------------*
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*----------------------------------------------------------------------*
      integer irun, ihyp, ireac, inkrun, ink
      integer i1, i2,iz1,iz2,id1, ntag
      integer inkmin, kl, ich
      integer  bin_dens,bin_time,ix,iy,iz
      real*8 h1
      real*8 ede,pdx,pdy,pdz,srt_0, xxx,yyy,zzz, srt0,srtp
      real*8 szig0, szig, pin2, szigr
      real*8 xx,yy,zz, qqabs, qqx,qqy,qqz, betax,betay,betaz
      real*8 px,py,pz, e3, fac, phase, density, denst
      real*8 p00,pxx,pxy,pxz, vx,vy,vz, vxx,vyy,vzz, tmas_m
      real*8 mmass, amass, traf
      real*8 valkapi, rr
      real*8 srt_in, srt_fin, s_in, s_fin, dsrtf
      real*8 rmass2
      real*8 pbeta, transf, q0, qq2, gammx
      real*8 rn
*----------------------------------------------------------------------*
      write(*,*) 'in kminu_hpi ', ede,pdx,pdy,pdz,srt_0,xxx,yyy,zzz,
     +     szig0,i1,id1, i2, ihyp,iz1,iz2,irun
      if(iz1+iz2 .gt. 0)                                          return
      if(iz1+iz2 .lt. -1)                                         return
      if(id1.ne.1)                                                return
      mmass = xkmas
      p00 = sqrt(ede**2 - pdx**2 -pdy**2 -pdz**2)
      pxx = pdx / p00 * mmass
      pxy = pdy / p00 * mmass
      pxz = pdz / p00 * mmass
      density = .0
      if (i_kminu_pot  .gt. 0) then
         ich = -1
         call gradukao2(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  density, ich, vx,vy,vz,vxx,vyy,vzz)
      endif
      write(*,*) 'kminu_hpi after gradukao2',ihyp
      tmas_m  = mmass
      if(ihyp .eq. 1) amass = xlmas
      if(ihyp .eq. 2) amass = xsmas
      dsrtf  = xkmas - tmas_m
      srt_in  = srt_0
      s_in     = srt_in**2
      pin2= ((s_in-pmass**2-amass**2)**2-(2.*pmass*amass)**2)/
     1        (4.*s_in)
      if (pin2 .le. .0001)  return
      srt_fin = srt_0 + dsrtf
      s_fin  = srt_fin**2
      srt0 = rmass+xkmas
      if(srt_fin .le. srt0)                   return
      rmass2= rmass**2
c       q0  : energy of the K-
      q0   = 0.5 * (s_fin + xkmas**2 - rmass2) / srt_fin  ! in medium wavenum
      qq2     = q0**2 - xkmas**2
      if(qq2 .le. .0)                           return
      qqabs   = sqrt(qq2)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gammx = ede / srt_0
      traf  = gammx / (gammx+1.0)
*----------------------------------------------------------------------*
c     write(*,*)   '  in  kminu_hpi   gammx ', gammx,traf
c...  cross sections:
      szig = 0.0
      srtp  = srt_fin
      h1   = srt_fin - srt0
      if(ihyp .eq. 1) then         ! lambda+pi
        ireac = 5
        if (iz1 .eq. 0) then
          if (srtp .lt. 1.5) then
             szig = 24. * sqrt(h1) / (1. +h1*h1)
          else
             szig = 354.* exp(-3.*srtp)
     1              + 2.07 / (1. + ((srtp-1.775)/0.065)**2)
          endif
        else
          szig = 170.*exp(-2.444*srtp) +
     1             3.23*exp(-((srtp-1.7)/0.1)**2)
     2           + 0.91*exp(-((srtp-2.04)/0.172)**2)
        endif
        szig = qq2 / pin2 * szig
        szig = min(szig, 100.)
      elseif (ihyp .eq. 2) then
        ireac = 6
        if (iz1 .eq. 0) then
           szig = 42.9*exp(-2.207*srtp) +
     1            .6 * exp(-((srtp-1.65)/.169)**2)
        elseif (iz1 .eq. 1) then
           if (srtp .lt. 1.5) then
             szig = 0.106 * sqrt(h1) / (h1+.005)**2
           else
             szig = 55.4 * exp(-2.26*srtp +11.275*srtp**(-3.26)) +
     1            0.335 * exp(-((srtp-2.06)/0.15)**2)
           endif
c          szig = min(szig, 10.)
        else
           if (srtp .lt. 1.5) then
             szig = 0.106 * sqrt(h1) / (h1+.005)**2
           else
             szig = 3572.*exp(-4.5*srtp)
     +           + 1.4 / (((srtp-1.690)/0.031)**2+1.)
     +           + 0.4 / (((srtp-1.830)/0.050)**2+1.)
           endif
        endif
        szig = qq2 / pin2 * szig
        szig = min(szig, 100.)
      endif
      szigr  = szig/szig0
c      szigr  = szig/szig0  * p_hyp(4,i2)
c***********************************************************************
      write(*,*)  '  in kminu_hpi  szig  ', ihyp, szig
      inkrun = (max_kminu/num)
      ink = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 710 kl=1,inkrun                    !        iphi_bb       hw
        ink = ink + 1
        if (nx_kminu(0,ink) .eq. 0) goto 712
        if (valkapi .gt. p_kminu(4,ink)) then
           inkmin = ink
           valkapi  = p_kminu(4,ink)
        endif
  710 continue
      ink = inkmin
      if (szigr .lt. valkapi) goto 98
  712 continue
c     write(*,*)  '  in kminu_hpi  szig 2 ', ink, valkapi
   60 continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 60
c         write(*,*) 'kminus_hpi',szig,traf,q0, ink, irun, kl
          qqx= xx*qqabs/rr
          qqy= yy*qqabs/rr
          qqz= zz*qqabs/rr
*   kminu_ momentum in observable system      PAULI
          pbeta  = betax*qqx + betay*qqy + betaz*qqz
          transf = gammx * (traf * pbeta + q0)
          p_kminu(0,ink)= density
          p_kminu(1,ink)= qqx + betax * transf
          p_kminu(2,ink)= qqy + betay * transf
          p_kminu(3,ink)= qqz + betaz * transf
          nx_kminu(5,ink)= 1001 * nint(200+zzz)
c         nx_kminu(5,ink)= 1001 * nthw
c         write(*,*) ' nx_kminu5 hpi',ink,zzz, nx_kminu(5,ink)
          r_kminu(1,ink)= xxx
          r_kminu(2,ink)= yyy
          r_kminu(3,ink)= zzz
          nx_kminu(0,ink) = 1
          nx_kminu(1,ink) = ireac
          nx_kminu(2,ink) = 0
          nx_kminu(3,ink) = 0
c-----------------------------------------------------------------------
c           density dependence
          ix = nint(xxx)
          iy = nint(yyy)
          iz = nint(zzz)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0

*   p3:                    nucleon momentum in observable system
          e3 = sqrt(rmass**2 + qqabs**2)
          transf = gammx * (-traf * pbeta + e3)
          px= -qqx + betax * transf
          py= -qqy + betay * transf
          pz= -qqz + betaz * transf
          if(ipauli.eq.1) 
     &      call pauli(i2,ntag,iseed,phase,xxx,yyy,zzz,px,py,pz)
          fac = szigr*(1.0-phase)
          p_kminu(4,ink)= fac

          bin_dens = int(denst*10.)
          bin_time = int(time*1.)
c          kmi_hpi(bin_dens,bin_time)=kmi_hpi(bin_dens,bin_time)+
c     1   p_kminu(4,ink)/real(num*isubs)
*
      if(ink .eq. 20039) write(*,*) 'kminus found in piY coll',
     1     irun,ink,ireac,szig,phase, p_kminu(1,ink), transf,
     2     betax, gammx, traf, q0, qqabs, srt_fin, srt_0, dsrtf
c      write(33,*) ' kminus in hypi-pi',
c    1         irun, ink, ireac, srt, tmas_m, szig, phase
 98   continue
c******
c      write(*,*) 'end of kminu_hpi'
c         if (fac  .lt. .0) then
c            write(*,*)  '  kminus negativ in hpi ', ink, szigr, pin2,
c    1       srt_in, ireac
c            stop
c         endif
      return
************************************************************************
      end
