************************************************************************
*                                                                      *
      subroutine kaondin
*                                                                      *
*       purpose:    calculating kaon production from:                  *
*                               n+n  collision      (entry kaondbb)    *
*                               n+pi collision      (entry kaondpi)    *
*                   -------------------------------------------------  *
*                   initialization by subroutine kaonin                *
*                   -------------------------------------------------  *
*                   final cross section by entry kaonout               *
*                   -------------------------------------------------  *
*                                                                      *
*         yref    - -1 * target rapidity in the frame                  *
*         nyk     - number of rapidities               (integer,input) *
*         nqtk    - number of transverse momenta       (integer,input) *
*         nfk     - number of angles                   (integer,input) *
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*-----------------------------------------
      common /nthwhw/  nthw, isu_hw
      integer nthw, isu_hw
      common /counthw/ ihw
      integer ihw
      common /correlkaon/ ikairun(0:maxkaon)
      integer ikairun

      integer  inkrun, inkmin
      integer  bin_dens,bin_time
      integer  kpl_dbb(0:999,0:999),kpl_dpi(0:999,0:999)
      integer  lams_dbb(0:999,0:999),lams_dpi(0:999,0:999)
      integer  i, ii,iy,i1,i2,ireac,ix,iz,maxde,id2,jj
      integer  i61,i62,id6,id62,iz1,iz2,iz12,id1, idd1, ihyp
      integer  irun,ink,ntag,kl,ikaonbb,ikaonpi
      real*8   gamma,srt,srt1,xxx,yyy,zzz, s0,s1, dsrtk, srtp, sprim
      real*8   sig0,facreac,em12,em22,tmass2,emm2,s,pmax2,pmax,denst,
     1         pmaxl,xsmas2, ddmass, fac2, p00,pxx,pxy,pxz
      real*8   szigzw,szigrk,szig,szigx, traf,spr, sig_tsu, sig_ts1
      real*8   xx,yy,zz,rr,transf,qqx,qqy,qqz
      real*8   e3,pbeta,ede,pdx,pdy,pdz,szig0,srt0,phase,szigww
      real*8   xlmas2,gamm,betax,betay,betaz,q0,qq2,qqabs, q1, q01
      real*8   rn,fackaondi2,pkaonbb, pkaonpi
      real*8   pmaxs,pmaxs2,szigko,srtmin,etotal,valkabb,valkapi
      real*8   ranp,pka,eka,enupr,pnupr,factn1,factn2, pkaon
      real*8   vx,vy,vz, vxx,vyy,vzz, tmass, density, mmass
      real*8   scal_hartnack
*----------------------------------------------------------------------*
      integer  ihyp0, n_hyp, ihypx, ic_hyp
      integer ic_rnd, ic_rnd0, ic_hyp0, ic_hypx, izd
      real*8   szig00, pp, sig_tsu0, xmas
      real*8   dmassx, hmassx, rnaa
c----------------
      real*8  s001, s002, xr, sigtp, siglp, sigsp, sigt1, sigl1, sigs1
      real*8  spp_plp, spp_0sp, spp_psp, spp_psn
*----------------------------------------------------------------------*
      parameter (maxde=20)
*----------------------------------------------------------------------*
      real*8  p3(3),beta(3), us(0:3)
*----------------------------------------------------------------------*
      save ikaonbb,ikaonpi,pkaon
*----------------------------------------------------------------------*
      pkaon = 1./real(pkaonnum)
      pkaonbb = 1.
      pkaonpi = 1.
      valkabb = 0.0
      valkapi = 0.0
      do ii=1,maxkaon
        ika(1,ii)=0
        ika(5,ii)=0
        pkao(4,ii)  = .0
      enddo
      return
*----------------------------------------------------------------------*
      entry kaondbb(id1,id2,beta,gamma,srt,xxx,yyy,zzz,sig0,i61,i62,irun
     &  ,i1,i2,etotal,kpl_dbb,lams_dbb)
*       variables:                                                     *
*         ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*     cross-sections are taken from cassing et al., phys. lett.        *
*                                 = cugnon et al.,nucl. phys. a422 635 *
*                                 = lang et al., nucl. phys. a541 507  *
*     for channels including delta we use a factor facreac multiplying *
*                  the n+n cross sections (taken from ko and randrup)  *
*              delta + n:     facreac=0.75                             *
*              delta + delta: facreac=0.5                              *
*----------------------------------------------------------------------*
*                                                                      *
*     hw:  ikaoncr = 4
*          for N+N collision tsushima for lambda and sigma production  *
*          for higher states isospin average of tsushima  are used     *
*----------------------------------------------------------------------*
c     return            !   TEST
      write(*,*)  '  kaon mass kaondbb2', pkaon
      if(rn(iseed) .ge. pkaon)                    return
      scal_hartnack = 1.0     !    1.34
c
      write(*,*)  '  kaon mass kaondbb1'
c
      us(0) = gamma
      us(1) = gamma * beta(1)
      us(2) = gamma * beta(2)
      us(3) = gamma * beta(3)
      pxx = us(1) * xkmas
      pxy = us(2) * xkmas
      pxz = us(3) * xkmas
      if (ikaonpot .gt. 0) then
!          call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
!      1                  density,1, vx,vy,vz,vxx,vyy,vzz)		!HS
         call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  1, vx,vy,vz,vxx,vyy,vzz)
      else
         density = .0
	 mmass = xkmas
      endif
      tmass  = mmass
      dsrtk = mmass - xkmas
      srtmin = rmass + tmass + xlmas
cc      write(*,*)  '  kaon mass kaondbb ',tmass, srt, srtmin
      s     = srt**2
      srtp  = srt - dsrtk + .005
      sprim = srtp**2
      s001 = (rmass + xkmas + xlmas)**2
      s002 = (rmass + xkmas + xsmas)**2
      ddmass =  1.232
      pmax2 = .25*(s-(rmass+xlmas+tmass)**2)*
     1            (s-(rmass+xlmas-tmass)**2)/s
      if(pmax2 .le. 0.0)                                  return
      pmax  = sqrt(pmax2)
      dmassx = rmass
***                                                                  ***
      traf  = gamma / (gamma+1.0)
c
      ireac = 0
      if(id1+id2 .eq. 2) ireac = 1
      if(id1+id2 .eq. 3) ireac = 2
      if(id1+id2 .ge. 4) ireac = 4
      if(id1.eq.2 .and. id2.eq.2) ireac = 3
      facreac = 1.0
      if(id1+id2 .eq. 3)          facreac = 0.75
      if(id1.eq.2 .and. id2.eq.2) facreac = 0.5
      em12  = rmass**2
      em22  = xlmas**2
      tmass2= tmass**2
      emm2  = (rmass+xlmas)**2
      if (ireac .eq. 0  .or. ireac .eq. 4)  return
c     IF (ireac .eq. 1) RETURN           !            TEST
      ihyp = 0
      ihyp0 = 0
      hmassx = xlmas
      dmassx = rmass
      if (ikaoncr .eq. 4 .or. ikaoncr .le. 6) then
c--------------------------------------------------------
c            =4         !   Tsushima  nucl-th/9801063
c            =5         !   NN  Sibirtsev
c--------------------------------------------------------
         iz12 = id(2,i1) + id(2,i2)
         rnaa = rn(iseed)
         ic_rnd = nint(rnaa)
         rnaa = rn(iseed)
         ic_rnd0 = nint(rnaa)
c
         sig_tsu =  .0
         sig_tsu0 = .0
         if (ireac .eq. 1) then
          if (ikaoncr .eq. 5) then
           if (sprim .lt. s001) return
           xr = s001 / sprim
           spp_plp = 0.732 * (1.-xr)**1.8 * xr**1.5
c  means : s(pp->K(p,0),(L/S),(n,p)
           if (sprim .le. s002) then
              spp_0sp = .0
              spp_psp = .0
           else
              xr = s002 / sprim
              spp_0sp =  0.338 * (1.-xr)**2.25 * xr**1.35
              spp_psp =  0.275 * (1.-xr)**1.98 * xr
           endif
              spp_psn =  2.*spp_psp
           if (iz12 .eq. 2) then
              siglp = spp_plp
              sigsp = spp_psp + spp_psn
              sigtp = sigsp + siglp
              ihyp  =  1
              if (rn(iseed)*sigtp .gt. siglp) ihyp = 2
              ic_hyp = iz12 - 1
              if (ihyp .gt. 1 ) ic_hyp = iz12-1-ic_rnd
              sigs1 = spp_0sp         !   K0  production
              sigt1 = sigs1
              ihyp0  =  2
              ic_hyp0 = 1
           elseif (iz12 .eq. 1) then
              siglp = 2.5 * spp_plp
              sigsp = 2.5 * (spp_0sp + spp_psp)
              sigtp = sigsp + siglp
              ihyp  =  1
              ic_hyp = iz12 - 1
              if (rn(iseed)*sigtp .gt. siglp) ihyp = 2
              if (ihyp .gt. 1 ) ic_hyp = iz12-1-ic_rnd
              sigt1 = sigtp
              ihyp0  =  1
              if (rn(iseed)*sigt1 .gt. sigl1) ihyp0 = 2
              ic_hyp0 = iz12
              if (ihyp0.gt. 1 ) ic_hyp0 = iz12 - ic_rnd0
           elseif  (iz12 .eq. 0) then
              sigsp = spp_0sp
              sigtp = sigsp
              ihyp = 2
              ic_hyp = -1
              sigl1 = spp_plp
              sigs1 = spp_psp + spp_psn
              sigt1 = sigl1 + sigs1
              ihyp0  =  1
              if (rn(iseed)*sigt1 .gt. sigl1) ihyp0 = 2
              ic_hyp0 = iz12
              if (ihyp0.gt. 1 ) ic_hyp0 = iz12 - ic_rnd0
           endif
           sig_tsu = sigtp
           sig_tsu0 = sigt1
          elseif (ikaoncr .eq. 4)  then
           call kaon_nu_nu (sprim, iz12, sig_tsu, ihyp, ic_hyp)
           if (ihyp .gt. 1 ) ic_hyp = iz12-1-ic_rnd
c          write(*,*) ' tsushima1 ', sprim, iz12, sig_tsu, ihyp
           call kaon_nu_nu (sprim, 2-iz12, sig_tsu0, ihyp0, ic_hyp0)!K0 product.
            ic_hyp0 = -ic_hyp0
           if (ihyp0.gt. 1 ) ic_hyp0 = iz12 - ic_rnd0
          elseif (ikaoncr .eq. 6)  then
           call kaon_nu_sib (sprim, iz12, sig_tsu, ihyp, ic_hyp)
           if (ihyp .gt. 1 ) ic_hyp = iz12-1-ic_rnd
c          write(*,*) ' tsushima1 ', sprim, iz12, sig_tsu, ihyp
           call kaon_nu_sib (sprim, 2-iz12, sig_tsu0, ihyp0, ic_hyp0)!K0 product
            ic_hyp0 = -ic_hyp0
           if (ihyp0.gt. 1 ) ic_hyp0 = iz12 - ic_rnd0
          endif
         elseif(ireac .ge. 2) then
           fac2 = 1.
           if (ireac .gt. 2) fac2 = 2./3.
           ihyp  =  0
           ihyp0 =  0
           ic_hyp = 0
           ic_hyp0= 0
	   s0 = 6.504
	   s1 = 6.904
           if (sprim .le. s0) return
           sig_tsu = 2.085*(sprim/s0 -1.)**2.227 * (s0/sprim)**2.511
           sig_ts1 = .0
	   if (sprim .gt. s1) sig_ts1 =
     1          19.77*(sprim/s1 -1.)**2.799 * (s1/sprim)**6.303
           if (rn(iseed) .gt. sig_tsu/(sig_ts1+sig_tsu))  ihyp = 1
           if (rn(iseed) .gt. sig_tsu/(sig_ts1+sig_tsu)) ihyp0= 1
              rnaa = iz12 -1.5
              izd = 1
              if (rnaa .gt. .0)  izd = 0
           if (ihyp .eq. 1) ic_hyp = ic_rnd  - izd
              rnaa = iz12 -0.5
              izd = 1
              if (rnaa .gt. .0)  izd = 0
           if (ihyp0.eq. 1) ic_hyp0= ic_rnd0  - izd
           sig_tsu =  fac2 * (sig_ts1+sig_tsu)
           sig_tsu =  sig_tsu   ! That's a reduction for Delta's
           sig_tsu0= sig_tsu
	 endif
        szig = sig_tsu / (sig0*pkaon)
        szig00 = sig_tsu0 / (sig0*pkaon)
        write(*,*) ' tsushima2 ', s, i1, i2, sig_tsu, sig_tsu0,
     1               ireac, tmass, pmax, sig_tsu, sig_ts1, szig,
     2               iz12, ihyp, ic_hyp, ihyp0, ic_hyp0
        goto 1200
      endif
c
      szigzw = 0.8 * pmax2**2
      szigrk = 0.036 * (pmax+pmaxs)/tmass
      szigww = 1.1 * (srt-srtmin)**1.71/
     & (0.8*(srt-srtmin)**1.7+4.5*sqrt(srt-srtmin)+1.0)
      szigko = 0.8 * (srt-srtmin)**1.6
      if(srt.ge.2.75) szigko=
     & 0.8*(2.75-srtmin)**1.6+0.25*(srt-2.75)**1.1/(0.5+(srt-2.75)**1.5)
      if(ikaoncr.eq.1) szig=szigzw
      if(ikaoncr.eq.2) szig=szigko
      if(ikaoncr.eq.0 .or. ikaoncr.eq.3) szig=szigww
      szig   = szig/sig0 * facreac / pkaon
      ihyp0  = 0
      ic_hyp0  = 0
      szig00 = szig
 1200 continue
c***********************************************************************
c     if (szig00 + szig  .lt. 1.e-7) return
      inkrun = (maxkaon/num)
      szigx = szig
      ihypx = ihyp
c-----------------------------------------
      do 1250 n_hyp = 0,1
        if (n_hyp .eq. 1) then
                szigx = szig00
                ihypx = ihyp0
        endif
        if (szigx  .lt. 1.e-7)  goto  1250
        inkmin = 0
        valkabb  = 1000.0
        ink = (irun-1) * inkrun
        do 10 kl=1,inkrun                    !        iphi_bb       hw
          ink = ink + 1
          if (ika(1,ink) .eq. 0) goto 12
          if (valkabb .gt. pkao(4,ink)) then
            inkmin = ink
            valkabb  = pkao(4,ink)
          endif
   10  continue
   13  ink = inkmin
       if (szigx .lt. valkabb) goto 99
   12  continue
c     write(*,*)  ' kaons dbb found  at ',ink, ireac, szigx
   15  continue
          ranp     = rn(iseed)
c------------         this is not good  for kaons
c         if((ranp**2-ranp**3)*6.75 .lt. rn(iseed)) go to 15
c-------------
          if(ranp**2*sqrt(1.-ranp)*3.94 .lt. rn(iseed)) go to 15
c       pmax  = the maximal kaon momentum
        pmax = 0.0
        hmassx = xlmas
        dmassx = rmass
        if (ihypx .eq. 1) hmassx = xsmas
        if (ihypx .eq. 2) dmassx = ddmass
        if (ihypx .eq. 3) then
          dmassx = ddmass
          hmassx = xsmas
        endif
        pmax2=(s-(dmassx+hmassx+tmass)**2)*
     1      (s-(dmassx+hmassx-tmass)**2) / s
        if(pmax2 .gt. 0.0) pmax = 0.5 * sqrt(pmax2)
          pka       = pmax*ranp
  20      xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 20
          eka      = sqrt(pka**2 + tmass**2)
          spr      = s - 2.0*srt*eka + tmass**2
          enupr    = 0.5*(spr + rmass**2 - xlmas**2)/sqrt(spr)
          pnupr    = sqrt(enupr**2-rmass**2)
          qqx      = pka * xx / rr
          qqy      = pka * yy / rr
          qqz      = pka * zz / rr
  30      xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 30
          factn1   = spr + sqrt(spr)*(srt-eka)
          factn2=pnupr*(qqx*xx+qqy*yy+qqz*zz)/factn1/rr-enupr/sqrt(spr)
*   p3:                    nucleon momentum in i1-i2-c.m. system
          p3(1)    = pnupr*xx/rr + factn2*qqx
          p3(2)    = pnupr*yy/rr + factn2*qqy
          p3(3)    = pnupr*zz/rr + factn2*qqz
          e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p3:                    nucleon momentum in observable system
          pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
          transf = gamma * (traf * pbeta + e3)
          p3(1)  = p3(1) + beta(1) * transf
          p3(2)  = p3(2) + beta(2) * transf
          p3(3)  = p3(3) + beta(3) * transf
c         write(*,*) '  start pauli  in  kaondi ',
c    1    xxx,yyy,zzz,p3(1),p3(2),p3(3)
          if(ipauli.eq.1 .and.id1.eq.1)
     &     call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
*   kaon momentum in observable system
          write(*,*) '  end  pauli  in  kaondi ', ntag, phase
cc          fackaondi2 = szigx*pkaonbb*(1.0-phase)
          fackaondi2 = szigx*(1.0-phase)
          if (n_hyp .eq. 0) szig   = fackaondi2
          if (n_hyp .eq. 1) szig00 = fackaondi2
          if(ikaondi.eq.2.and.fackaondi2.lt.rn(iseed))         goto 99
          if(ikaondi.eq.2.and.fackaondi2.gt.1.1) write(*,*)
     &         'hiba in kaonbb, fackaondi2> 1.1',fackaondi2,szigx
          pbeta  = beta(1)*qqx + beta(2)*qqy + beta(3)*qqz
          transf = gamma * (traf * pbeta + eka)
          pkao(1,ink)= qqx + beta(1) * transf
          pkao(2,ink)= qqy + beta(2) * transf
          pkao(3,ink)= qqz + beta(3) * transf
          pkao(4,ink)= szigx*(1.0-phase) * scal_hartnack  !  1.34 = Hartnack
c         write(*,*)  ' kaon momentum ', ink,(pkao(i,ink),i=1,4),
c    1    beta, transf, pka, pmax
c         write(*,*) 'kaondbb',i,irun, kl, ink, i1,i2,id1,id2,
c    1                szigx, phase, valkabb
          rkao(1,ink)= xxx
          rkao(2,ink)= yyy
          rkao(3,ink)= zzz
          ika(1,ink) = 1 + n_hyp
          ika(2,ink) = ireac
          ika(3,ink) = i1
          ika(4,ink) = i2
c         ika(6,ink)= 1001 * nint(200+zzz)
c         ika(6,ink)= 1001 * nthw
          ikairun(ink) = irun
          pkao(5,ink)= qqx + beta(1) * transf
          pkao(6,ink)= qqy + beta(2) * transf
          pkao(7,ink)= qqz + beta(3) * transf
          if(ikaondi.eq.2) then
            p(1,i1)  = p3(1)
            p(2,i1)  = p3(2)
            p(3,i1)  = p3(3)
            p(1,i2)  = etotal*beta(1) - pkao(1,ink) - p3(1)
            p(2,i2)  = etotal*beta(2) - pkao(2,ink) - p3(2)
            p(3,i2)  = etotal*beta(3) - pkao(3,ink) - p3(3)
          end if
          if (n_hyp .eq. 1) goto 99
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
	 kpl_dbb(bin_dens,bin_time)=kpl_dbb(bin_dens,bin_time)+1
c-------------  creation time ---------
c          ika(6,ink)= 1001 * nthw 				!HS end
          goto  99
c-----------------------------------------------------------------------
 99   continue
 1250 continue      !  end of n_hyp
      if (i_kminu .eq. 1) then                  !  storage of hyperons
          n_hyp = 0
 3000 continue !                         do loop for K+, K0
          n_hyp = n_hyp + 1
          if (n_hyp .eq. 1) then
               ihypx = ihyp
               pp    = szig
               ic_hypx = ic_hyp
          else
               ihypx = ihyp0
               pp    = szig00
               ic_hypx = ic_hyp0
          endif
          if (szigx .lt. 1.e-7)  goto 3030
          hmassx = xlmas
          if (ihypx .eq. 1) hmassx = xsmas
          if (ihypx .eq. 2) dmassx = rmass
          if (ihypx .eq. 3) then
             dmassx = ddmass
             hmassx = xsmas
          endif
          pmax2=(s-(dmassx+hmassx+tmass)**2)*
     1            (s-(dmassx+hmassx-tmass)**2) / s
          if(pmax2 .le. 0.0) goto 3030
          pmax = 0.5 * sqrt(pmax2)
c
      inkrun = (max_kminu/num)
      ink = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 3020 kl=1,inkrun                    !        look for place     hw
        ink = ink + 1
        if (nx_hyp(0,ink) .eq. 0) goto 3024
        if (valkapi .gt. p_hyp(4,ink)) then
           inkmin = ink
           valkapi  = p_hyp(4,ink)
        endif
 3020 continue
 3022 ink = inkmin
      if (pp .lt. valkapi) goto 3030
 3024 continue
c     write(*,*)' sigma storage in dbb',ink,n_hyp,ihypx,ic_hypx,pp,xxx
c-------------------------------------
 3015  continue
          ranp     = rn(iseed)
          rr      = ranp*ranp
          rnaa = 2.598 * rr * sqrt(1.-rr)
          if(rnaa .lt. rn(iseed)) go to 3015
          pka       = pmax*ranp
 3102  continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 3102
          qqx= xx*pka/rr
          qqy= yy*pka/rr
          qqz= zz*pka/rr
*   kaon momentum in observable system
          traf  = gamma / (gamma+1.0)
          q01 = sqrt(pka**2 + hmassx**2)
          pbeta  = beta(1)*qqx + beta(2)*qqy + beta(3)*qqz
          transf = gamma * (traf * pbeta + q01)
          p_hyp(0,ink)= hmassx
          p_hyp(1,ink)= qqx + beta(1) * transf
          p_hyp(2,ink)= qqy + beta(2) * transf
          p_hyp(3,ink)= qqz + beta(3) * transf
          p_hyp(4,ink)= pp * scal_hartnack
          r_hyp(1,ink)= xxx
          r_hyp(2,ink)= yyy
          r_hyp(3,ink)= zzz
          ii = ihypx - 2*(ihypx/2)
          nx_hyp(0,ink) = ii + 1
          nx_hyp(1,ink) = ic_hypx
          nx_hyp(2,ink) = i1
          nx_hyp(3,ink) = i2
          nx_hyp(4,ink) = 1	!ireac for prod in BB

c-----------------------------------------------------------------------
	if(nx_hyp(1,ink).eq.0) then
c           density dependence
          ix = nint(xxx)
          iy = nint(yyy)
          iz = nint(zzz)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
	  bin_dens = int(denst*10.)
	  bin_time = int(time*1.)
	 lams_dbb(bin_dens,bin_time)=lams_dbb(bin_dens,bin_time)+1
	endif
c-----------------------------------------------------------------------

      write(*,*) ' produced - hyperon in bb', nthw, nx_hyp(1,ink),
     1    ink, p_hyp(4,ink), nx_hyp(0,ink), phase, szig00, szig
 3030 continue
          if (n_hyp .le. 1) goto 3000
      endif                             !    end storage of hyperons
c-----------------------------------------------------------------------
      write(*,*) 'end of kaondbb '
      return
************************************************************************
      entry kaondpi(ede,pdx,pdy,pdz,srt,xxx,yyy,zzz,szig0,i2,id2,id6,
     &        id62,iz1,iz2,i1,irun,kpl_dpi,lams_dpi)
*       variables:                                                     *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*     cross-sections are taken from cassing et al., phys. lett.        *
*                                 = cugnon et al.,nucl. phys. a422 635 *
*      id2 type of baryon, id6, id62 collision partner, i1 = meson
*----------------------------------------------------------------------*
c     return     !  TEST
c     write(*,*) '  start kaondpi ',nthw
c     return
      scal_hartnack = 1.00    !    1.34
      idd1 = ipi(1,i1)
      if (idd1 .ne. 1) return
      p00 = sqrt(ede**2 - pdx**2 -pdy**2 -pdz**2)
      pxx = pdx / p00 * xkmas
      pxy = pdy / p00 * xkmas
      pxz = pdz / p00 * xkmas
      if (ikaonpot .gt. 0) then
!          call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
!      1                  density,1, vx,vy,vz,vxx,vyy,vzz)		!HS
         call gradukaon(-1, xxx,yyy,zzz, pxx,pxy,pxz,mmass,
     1                  1, vx,vy,vz,vxx,vyy,vzz)
       else
            density = .0
       	    mmass = xkmas
      endif
      tmass  = mmass
      dsrtk = mmass - xkmas
      srtp = srt - dsrtk + .005
      srt0  = xlmas+tmass
      srt1  = xsmas+tmass
      if(srt.le.srt0)                                return
      tmass2= tmass**2
      xlmas2= xlmas**2
      xsmas2= xsmas**2
      s     = srt**2
c       q0  : energy of the kaon
      q0      = 0.5 * (s + tmass2 - xlmas2) / srt
      if(q0 .le. tmass)     return
      pmaxl   = .0
      if(q0 .gt. tmass) pmaxl = sqrt(q0**2 - tmass2)
      q1      = 0.5 * (s + tmass2 - xsmas2) / srt
      pmaxs   = .0
      if(q1 .gt. tmass) pmaxs = sqrt(q1**2 - tmass2)
      q01    =  q0
      qqabs   = pmaxl
*----------------------------------------------------------------------*
      facreac = 1.0
      sig_tsu  = .0
      sig_tsu0 = .0
      ireac = 5
      if(id2.eq.2)     ireac = 6
      ihyp = 0
      ihyp0 = 0
      if (ikaoncr .eq. 4) then      !   Tsushima  nucl-th/9602005
         if(ireac .eq. 5) then
           call pi_N_kplu(srtp, iz1, iz2, ihyp,  sig_tsu)
           call pi_N_kplu(srtp,-iz1,1-iz2, ihyp0,  sig_tsu0)
         elseif(ireac .eq. 6) then
c          return                !  test
           call pi_D_kplu(srtp, iz1, iz2, ihyp,  sig_tsu)
           call pi_D_kplu(srtp,-iz1,1-iz2, ihyp0,  sig_tsu0)
         else
            return
         endif
           if (sig_tsu + sig_tsu0 .lt. 1.e-7) return
           if (ihyp .eq. 1) qqabs = pmaxs
           if (ihyp .eq. 1) q01  = q1
c     if (nthw .eq. 37) write(*,*) ' density ', density, mmass,
c    1 tmass, srt, srt0, s, q0, q1, pmaxl, pmaxs, ihyp, ireac,
c    2 qqabs, q01
         szig = sig_tsu/szig0
         szig00 = sig_tsu0/szig0
      else
c
c***********************************************************************
        facreac = 1.0
        if(id2.eq.2) facreac = 0.75
        ireac = 5
        if(id2.gt.1) ireac = 6
        if(idd1.eq.2) ireac = 7
cc        if(srt.le.1.70) szig = 2.47*(srt-srt0)/szig0
cc        if(srt.gt.1.70) szig = 0.0225/(srt-1.60)/szig0
cc        if(iz1.eq.1 .and. iz2.eq.0) szig = 4.0 * szig
cc        if(iz1.eq.0 .and. iz2.eq.1) szig = 2.0 * szig
        if(srt.le.1.70) szig = .9*(srt-srt0)/.091
        if(srt.gt.1.70) szig = 0.09/(srt-1.60)
        if(iz1.eq.0 .and. iz2.eq.1) szig = 0.5 * szig
        if(id2.gt.1) szig=szig*pidelka
        szig   = szig/szig0 * facreac
        ihyp0   = 0
        szig00 = szig
      endif
c       write(*,*) ' tsushima7 ', s, i1, i2, sig_tsu, sig_tsu0,
c    1               ireac, iz1, iz2, ihyp, ihyp0


c***********************************************************************
 1300 continue
c     write(*,*) '  kaon  search  begins ', szig
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gamm  = ede / srt
      traf  = gamm / (gamm+1.0)
      inkrun = (maxkaon/num)
c-------------------------------------------------------
      szigx  =  szig
      qqabs   = pmaxl
      do 1350  n_hyp = 0,1
        if (n_hyp .eq. 1) then
           szigx  =  szig00
           qqabs   = pmaxs
        endif
        if (szigx .lt. 1.e-7) goto 1350
        ink = (irun-1) * inkrun
        inkmin = 0
        valkapi  = 1000.0
        do 710 kl=1,inkrun                    !        iphi_bb       hw
          ink = ink + 1
          if (ika(1,ink) .eq. 0) goto 712
          if (valkapi .gt. pkao(4,ink)) then
             inkmin = ink
             valkapi  = pkao(4,ink)
          endif
  710  continue
  713  ink = inkmin
      if (szigx .lt. valkapi) goto 96
  712  continue
c     write(*,*) '  kaon found at ',ink, ireac, szigx
   60  continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 60
          qqx= xx*qqabs/rr
          qqy= yy*qqabs/rr
          qqz= zz*qqabs/rr
*   kaon momentum in observable system
          pbeta  = betax*qqx + betay*qqy + betaz*qqz
          transf = gamm * (traf * pbeta + q01)
          pkao(1,ink)= qqx + betax * transf
          pkao(2,ink)= qqy + betay * transf
          pkao(3,ink)= qqz + betaz * transf
          pkao(4,ink)= szigx   * scal_hartnack     !  1.34 = Hartnack
          write(*,*) 'kaondpi',irun, kl, ink, n_hyp, szigx, valkapi
          rkao(1,ink)= xxx
          rkao(2,ink)= yyy
          rkao(3,ink)= zzz
          ika(1,ink) = 1 + n_hyp
          ika(2,ink) = ireac
          ika(3,ink) = i2
c         ika(6,ink)= 1001 * nint(200+zzz)
c         ika(6,ink)= 1001 * nthw
          ikairun(ink) = irun
          pkao(5,ink)= qqx + betax * transf
          pkao(6,ink)= qqy + betay * transf
          pkao(7,ink)= qqz + betaz * transf
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
	 kpl_dpi(bin_dens,bin_time)=kpl_dpi(bin_dens,bin_time)+1
   96 continue
 1350 continue       !    END n_hyp
c======================================
      if (i_kminu .eq. 1) then                  !  storage of hyperons
          n_hyp = 0
 2000 continue !                         do loop for K+, K0
          n_hyp = n_hyp + 1
          if (n_hyp .eq. 1) then
               ihypx = ihyp
               pp    = szig
               ic_hyp = iz1 + iz2 - 1
          else
               ihypx = ihyp0
               pp    = szig00
               ic_hyp =  iz1 + iz2
          endif
          if (abs(ic_hyp) .gt. 1) goto 2030
          if (pp .lt. 1.e-7)      goto 2030
      inkrun = (max_kminu/num)
      ink = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 2020 kl=1,inkrun                    !        look for place     hw
        ink = ink + 1
        if (nx_hyp(0,ink) .eq. 0) goto 2024
        if (valkapi .gt. p_hyp(4,ink)) then
           inkmin = ink
           valkapi  = p_hyp(4,ink)
        endif
 2020 continue
 2022 ink = inkmin
      if (pp .lt. valkapi) goto 2030
 2024 continue
c     write(*,*)  ' sigma storage ', ink, n_hyp, pp, xxx
c-------------------------------------
 2102  continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 2102
          qqx= xx*qqabs/rr
          qqy= yy*qqabs/rr
          qqz= zz*qqabs/rr
*   kaon momentum in observable system
          if (ihypx .eq. 0) xmas = xlmas
          if (ihypx .eq. 1) xmas = xsmas
          q01 = sqrt(qqabs**2 + xmas**2)
          pbeta  = betax*qqx + betay*qqy + betaz*qqz
          transf = gamm * (traf * pbeta + q01)
          p_hyp(0,ink)= xmas
          p_hyp(1,ink)= qqx + betax * transf
          p_hyp(2,ink)= qqy + betay * transf
          p_hyp(3,ink)= qqz + betaz * transf
          p_hyp(4,ink)= pp * scal_hartnack           !!  hartnack like
          r_hyp(1,ink)= xxx
          r_hyp(2,ink)= yyy
          r_hyp(3,ink)= zzz
          nx_hyp(0,ink) = ihypx + 1
          nx_hyp(1,ink) = ic_hyp
          nx_hyp(2,ink) = i2
          nx_hyp(3,ink) = -i1
	  nx_hyp(4,ink) = 2   !ireac for prod in piN


c-----------------------------------------------------------------------
	if(nx_hyp(1,ink).eq.0) then
c           density dependence
          ix = nint(xxx)
          iy = nint(yyy)
          iz = nint(zzz)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
	  bin_dens = int(denst*10.)
	  bin_time = int(time*1.)
	 lams_dpi(bin_dens,bin_time)=lams_dpi(bin_dens,bin_time)+1
	endif
c-----------------------------------------------------------------------

      write(*,*) ' produced - hyperon in piN',
     1    nthw, ink, nx_hyp(0,ink), p_hyp(1,ink),p_hyp(4,ink)
 2030 continue
          if (n_hyp .le. 1) goto 2000
      endif                             !    end storage of hyperons
          goto  98
c-----------------------------------------------------------------------
cc          if(iabs(ix).gt.maxx) ix=ix/iabs(ix)*maxx
cc          if(iabs(iy).gt.maxx) iy=iy/iabs(iy)*maxx
cc          if(iabs(iz).gt.maxz) iz=iz/iabs(iz)*maxz
cc          ithermvar = ithermvar + 1
cc          thermok(1,ink) = ithermvar
cc          thermok(2,ink) = ithermvar
cc          if(ithermo.eq.1) call tmunu(ix,iy,iz)
c          if(ithermo.eq.1) write(mterpri,'(i8,f6.2,3i4,14e11.3)')
c     &      ithermvar,time,ix,iy,iz,(avp(jj,ix,iy,iz),jj=1,14)
c-----------------------------------------------------------------------
 98   continue
c......................
      write(*,*) 'end of kaondpi'
      return
      end
c**********************************************************************
      subroutine kaon_nu_nu(ss, izz, sig, ihyp, ic_sig)
c
c          probability       for charge sigma(ic_sig)
c                            for hyperon (ihyp = 0,1)
c          sig  =  sig K+
c
      implicit none
      include 'cominput'
      real*8 ss, sig
      integer izz, ihyp, ic_sig
      real*8 s0,ap,bp,cp,an,bn,cn
      real*8 s1,appp,bppp,cppp,appn,bppn,cppn
      real*8 apnp,bpnp,cpnp,apnn,bpnn,cpnn, annn,bnnn,cnnn
      real*8 s2, annl, bnnl, cnnl, s3, apps, bpps, cpps
      real*8 sig1,sig1a,sig2,fac,sig3,rn
      integer ic_sig0
      data s0,ap,bp,cp,an,bn,cn /6.504,1.879,2.176,5.264,
     1                           2.812, 2.121, 4.893/
      data s1,appp,bppp,cppp,appn,bppn,cppn
     1         /6.904,5.321,2.735,8.510,1.466, 2.743, 3.271/
       data apnp,bpnp,cpnp,apnn,bpnn,cpnn, annn,bnnn,cnnn
     1   /11.02,2.782,7.674, 6.31,2.773,7.820,7.079,2.76,8.164/
       data s2, annl, bnnl, cnnl, s3, apps, bpps, cpps
     1   /8.085, 6.166, 2.842, 1.960, 8.531, 10.00, 2.874, 2.543/
c     data from Tsushima et al., nucl-th/9801063
	     sig = .0
             ihyp = 0
             ic_sig = 0
             ic_sig0= 0
      if (ss .lt. s0) return
c
      if (izz .eq. 0)  then
	     sig = .0
      elseif (izz .eq. 1) then
             sig =  an * (ss/s0 -1.)**bn * (s0/ss)**cn
      elseif (izz .eq. 2)  then
             sig = ap * (ss/s0 -1.)**bp * (s0/ss)**cp
      endif
c
      if (ss .lt. s1) return
             sig1 = .0
      if (izz .eq. 0)  then
         ic_sig0 = -1
	 sig1 = annn * (ss/s1 -1.)**bnnn * (s1/ss)**cnnn
      elseif (izz .eq. 1) then
         sig1 = apnp * (ss/s1 -1.)**bpnp * (s1/ss)**cpnp
         sig1a =apnn * (ss/s1 -1.)**bpnn * (s1/ss)**cpnn
         sig1  = sig1+sig1a
         if (rn(iseed)  .lt.  sig1a / sig1) ic_sig0 = -1
      elseif (izz .eq. 2) then
	 sig1 = appp * (ss/s1 -1.)**bppp * (s1/ss)**cppp
	 sig1a= appn * (ss/s1 -1.)**bppn * (s1/ss)**cppn
         sig1  = sig1+sig1a
         if (rn(iseed)  .gt.  sig1a / sig1) ic_sig0 =  1
      endif
c        sig1  = sig1 / 4.      !   reduction of Tsushima's estimate
         sig = sig + sig1
         if(rn(iseed) .lt. sig1/ sig) then
                 ihyp = 1
                 ic_sig = ic_sig0
         endif
c
             sig2 = .0
      if (ss .lt. s2) return
      if (izz .eq. 0)  then
         fac = 1.0
      elseif (izz .eq. 1) then
         fac = .3333
      elseif (izz .eq. 2) then
         fac = 0.3333
      endif
         sig2 = fac * annl * (ss/s2 -1.)**bnnl * (s2/ss)**cnnl
         if(rn(iseed) .gt. sig / (sig2 + sig)) ihyp = 2
         sig = sig + sig2
c
      if (ss .lt. s3) return
             sig3 = .0
      if (izz .eq. 0)  then
         fac = 0.5
      elseif (izz .eq. 1) then
         fac = .5
      elseif (izz .eq. 2) then
         fac = 1.1666
      endif
         sig3 = fac * apps * (ss/s3 -1.)**bpps * (s3/ss)**cpps
         if(rn(iseed) .gt. sig / (sig3 + sig)) ihyp = 3
         sig = sig + sig3
      return
      end
c
      subroutine kaon_nu_sib(ss, izz, sig, ihyp, ic_sig)
c
c          probability       for charge sigma(ic_sig)
c                            for hyperon (ihyp = 0,1)
c          sig  =  sig K+
c
      implicit none
      integer izz, ihyp, ic_sig
      real*8 ss, sig
      real*8 s0,ap,bp,cp,an,bn,cn
      real*8 s1,appp,bppp,cppp,appn,bppn,cppn
      real*8 apnp,bpnp,cpnp,apnn,bpnn,cpnn, annn,bnnn,cnnn
      real*8 s2, annl, bnnl, cnnl, s3, apps, bpps, cpps
      real*8 sig1,sig1a,sig2,sig3, fac,rn
      integer ic_sig0
      include 'cominput'
      data s0,ap,bp,cp,an,bn,cn /6.504, 0.732, 1.8, 3.3,
     1                           1.830, 1.800, 3.300/
      data s1,appp,bppp,cppp,appn,bppn,cppn
     1         /6.904,5.321,2.735,8.510,1.466, 2.743, 3.271/
       data apnp,bpnp,cpnp,apnn,bpnn,cpnn, annn,bnnn,cnnn
     1   /11.02,2.782,7.674, 6.31,2.773,7.820,7.079,2.76,8.164/
       data s2, annl, bnnl, cnnl, s3, apps, bpps, cpps
     1   /8.085, 6.166, 2.842, 1.960, 8.531, 10.00, 2.874, 2.543/
c     data from Tsushima et al., nucl-th/9801063
	     sig = .0
             ihyp = 0
             ic_sig = 0
             ic_sig0= 0
      if (ss .lt. s0) return
c
      if (izz .eq. 0)  then
	     sig = .0
      elseif (izz .eq. 1) then
             sig =  an * (ss/s0 -1.)**bn * (s0/ss)**cn
      elseif (izz .eq. 2)  then
             sig = ap * (ss/s0 -1.)**bp * (s0/ss)**cp
      endif
c
      if (ss .lt. s1) return
             sig1 = .0
      if (izz .eq. 0)  then
         ic_sig0 = -1
	 sig1 = annn * (ss/s1 -1.)**bnnn * (s1/ss)**cnnn
      elseif (izz .eq. 1) then
         sig1 = apnp * (ss/s1 -1.)**bpnp * (s1/ss)**cpnp
         sig1a =apnn * (ss/s1 -1.)**bpnn * (s1/ss)**cpnn
         sig1  = sig1+sig1a
         if (rn(iseed)  .lt.  sig1a / sig1) ic_sig0 = -1
      elseif (izz .eq. 2) then
	 sig1 = appp * (ss/s1 -1.)**bppp * (s1/ss)**cppp
	 sig1a= appn * (ss/s1 -1.)**bppn * (s1/ss)**cppn
         sig1  = sig1+sig1a
         if (rn(iseed)  .gt.  sig1a / sig1) ic_sig0 =  1
      endif
c        sig1  = sig1 / 4.      !   reduction of Tsushima's estimate
         sig = sig + sig1
         if(rn(iseed) .lt. sig1/ sig) then
                 ihyp = 1
                 ic_sig = ic_sig0
         endif
c
             sig2 = .0
      if (ss .lt. s2) return
      if (izz .eq. 0)  then
         fac = 1.0
      elseif (izz .eq. 1) then
         fac = .3333
      elseif (izz .eq. 2) then
         fac = 0.3333
      endif
         sig2 = fac * annl * (ss/s2 -1.)**bnnl * (s2/ss)**cnnl
         if(rn(iseed) .gt. sig / (sig2 + sig)) ihyp = 2
         sig = sig + sig2
c
      if (ss .lt. s3) return
             sig3 = .0
      if (izz .eq. 0)  then
         fac = 0.5
      elseif (izz .eq. 1) then
         fac = .5
      elseif (izz .eq. 2) then
         fac = 1.1666
      endif
         sig3 = fac * apps * (ss/s3 -1.)**bpps * (s3/ss)**cpps
         if(rn(iseed) .gt. sig / (sig3 + sig)) ihyp = 3
         sig = sig + sig3
      return
      end
c
      subroutine pi_N_kplu(sss, izpi, izba, ihyp,  sig)
      implicit none
      integer izpi, ihyp, izba
      real*8 sss, sig
      real*8 xmkaon,xmsigma,xmlambda,thr,thr1,sig1,fac,dtr,rn
      integer ic_sig0
      include 'cominput'
      ihyp = 0
      sig = .0
      xmkaon = 0.494
      xmsigma = 1.193
      xmlambda = 1.116
      thr  = xmkaon + xmlambda
      thr1 = xmkaon + xmsigma
      if (sss .le. thr) return
           if (izpi + izba .eq. 1) then
	      fac = 1.0
	      if(izpi .eq. 0) fac = .5
c             sig =  fac* .02279  * (sss-thr)**0.3893 /
              sig =  fac* .02000  * (sss-thr)**0.3893 /
     1                 ((sss-1.700)**2 + .01031 )
c             sig =  .007665 * (sss-thr)**0.1341 /
c    1                 ((sss-1.720)**2 + .007826)
           endif
      if (sss .le. thr1) return
            sig1 = .0
            dtr  = sss - thr1
      if (izba .eq. 1) then
        if(izpi .eq. 1) sig1 =
     1     .03591 * dtr**.9541  /((sss-1.890)**2+.01548)
     2   + .1594  * dtr**.01056 /((sss-3.000)**2+.9412)
        if(izpi .eq. 0) sig1 =
     1     .003978* dtr**0.5848 /((sss-1.740)**2+.006670)
     2   + .047090* dtr**2.1650 /((sss-1.905)**2+.006358)
        if(izpi .eq.-1) sig1 =
     1      .009803* dtr**.6021  /((sss-1.742)**2+.006583)
     2    + .00651 * dtr**1.4728 /((sss-1.940)**2+.006248)
      endif
      if (izba .eq. 0) then
        if(izpi .ge. 0) sig1 =
     1     .05014* dtr**1.2878 /((sss-1.730)**2+.006455)
      endif
         if(rn(iseed) .gt. sig / (sig1 + sig)) ihyp = 1
         sig = sig + sig1
      return
      end
c
      subroutine pi_D_kplu(sss, izpi, izba, ihyp,  sig)
      implicit none
      integer izpi, ihyp, izba
      real*8 sss, sig
      real*8 xmkaon,xmsigma,xmlambda,thr,thr1,sig1,fac,dtr,rn
      include 'cominput'
      ihyp = 0
      sig = .0
      xmkaon = 0.494
      xmsigma = 1.193
      xmlambda = 1.116
      thr  = xmkaon + xmlambda
      thr1 = xmkaon + xmsigma
      if (sss .le. thr) return
           if (izpi + izba .eq. 1) then
	      fac = 1.0
	      if(izpi .eq. 0) fac = .6666
	      if(izpi .eq. 1) fac = .3333
              sig =  fac* .00988  * (sss-thr)**0.7866 /
     1                 ((sss-1.720)**2 + .004852)
           endif
      if (sss .le. thr1) return
             sig1 = .0
            dtr  = sss - thr1
      if (izpi + izba .gt. 2) return
      if (izpi + izba .lt. 0) return
        sig1= .01052* dtr**0.8140 /((sss-1.725)**2+.007713) !  pi0+D0
        fac = 1.0
      if (izpi .eq. -1)  then
        if(izba .eq. 2) fac = .75
        if(izba .eq. 1) fac = .5
      endif
      if (izpi .eq.  1)  then
        if(izba .eq. 0) fac = .25
        if(izba .eq.-1) fac = 1.5
      endif
         sig1 = fac * sig1
         if(rn(iseed) .gt. sig / (sig1 + sig)) ihyp = 1
         sig = sig + sig1
      return
      end

