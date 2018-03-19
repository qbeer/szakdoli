************************************************************************
*                                                                      *
      subroutine ksi_kminu(mass,ksi_KN)
*                                                                      *
*     purpose:  ksi_minus production in  K- + N -->  KSI + K+          *
*           perturbative (do not decrease the K- prop)
*           what is 3.355 ? Clebsh?                                    *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      integer ksi_KN(0:999,0:999),bin_dens,bin_time
      real*8 radmax, radcol, e_thr, ksimas, sigpr, prob, valkapi
      real*8 vx,vy,vz, vxx,vyy,vzz, pxx1, pxy1,pxz1,p00,mmass1
      real*8 pxx2, pxy2,pxz2,mmass2,denst,crossfact
      real*8 xx, yy, zz, s, srt, srtp, pcm, rr, r2
      real*8 x1, y1, z1, px1, py1, pz1, em1, em12, e1, e1cm
      real*8 dx, dy, dz, rsqare, px2, py2, pz2, em2, em22, e2, x2, y2,z2
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt, b21,dxm
      real*8 rn,phase,p2beta,e2cm,qx1,qy1,qz1,qx2,qy2,qz2,dxp
      real*8 pn2(0:3), ps2(0:3), pn23(0:3)
      integer mass, id1, id2, iz2, inksi, maxloop, inkmin, kl
      integer maxk, irun, ink, ini, i1, i2, ii, ix, iy, iz, jj
      integer ireac
      parameter (radmax = 0.07 * 10.0) ! -> 0.15 mb (max-cross section)
      parameter (e_thr  = 1.816) ! ->  massksi + xkmas
      parameter (ksimas = 1.321)
      parameter (crossfact = 100.) !increase the crossX, then devide prob
*----------------------------------------------------------------------*
      write(*,*)  '  start  ksi_production kminus'
      maxk = max_kminu/num
      numk_minu = 0
      ireac     = 1
*     loop over all parallel runs
      do 1000 irun = 1,num
        ink = (irun-1) * maxk
        ini = (irun-1) * mass
        do 800 ii  = 1,maxk
          i1  = ii + ink
          if(nx_kminu(0,i1) .eq. 0 )                 goto 800
          numk_minu = numk_minu + 1
          x1  = r_kminu(1,i1)
          y1  = r_kminu(2,i1)
          z1  = r_kminu(3,i1)
          px1 = p_kminu(1,i1)
          py1 = p_kminu(2,i1)
          pz1 = p_kminu(3,i1)
          em1 = xkmas
          em12= em1**2
          e1  = sqrt( em12 + px1**2 + py1**2 + pz1**2 )
*     look for a scattering pseudonucleon in the same run
          i2  = ini
  600     i2  = i2 + 1
          if(i2.eq.nx_kminu(2,i1) .or. i2.eq.nx_kminu(3,i1))
     1                                          goto 600
          if(i2 .gt. ini+mass)                            goto 800
          id2 = id(1,i2)
          if (id2 .ne. 1) goto 600        !only protons
          x2 = r(1,i2)
          dx     = x1 - x2
            if(nbound.eq.1) then
              dxp = dx+2.0*boxx
              dxm = dx-2.0*boxx
              if(abs(dx) .gt. abs(dxp)) dx=dxp
              if(abs(dx) .gt. abs(dxm)) dx=dxm
            end if
          if (abs(dx) .gt. radmax)                              goto 600
          y2 = r(2,i2)
          dy     = y1 - y2
            if(nbound.eq.1) then
              dxp = dy+2.0*boxx
              dxm = dy-2.0*boxx
              if(abs(dy) .gt. abs(dxp)) dy=dxp
              if(abs(dy) .gt. abs(dxm)) dy=dxm
            end if
          if (abs(dy) .gt. radmax)                              goto 600
          z2 = r(3,i2)
          dz     = z1 - z2
            if(nbound.eq.1) then
              dxp = dz+2.0*boxz
              dxm = dz-2.0*boxz
              if(abs(dz) .gt. abs(dxp)) dz=dxp
              if(abs(dz) .gt. abs(dxm)) dz=dxm
            end if
          if (abs(dz) .gt. radmax)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. radmax**2)                            goto 600
*         now particles are close enough to each other !
          px2    = p(1,i2)
          py2    = p(2,i2)
          pz2    = p(3,i2)
          em2    = e(i2)
          em22   = e(i2)**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                           - (pz1+pz2)**2
          if (s .lt. e_thr**2)  goto 600
*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if (brel .gt. radmax)                                goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
*
*   now  the kminus  will collide and produce ksi_s in  this time step
************************************************************************
          iz2 = id(2,i2)
          srt = sqrt(s)
**  modified medium masses mmass1,2 of the reaction products Xi- and K+/0 **
      p00 = sqrt((e1+e2)**2 - (px1+px2)**2 -(py1+py2)**2 -(pz1+pz2)**2)
      pxx1 = (px1+px2) / p00 * ksimas
      pxy1 = (py1+py2) / p00 * ksimas
      pxz1 = (pz1+pz2) / p00 * ksimas

      if (i_ksi_pot .ge. 1) then
         call graduxi(-1, x1,y1,z1,pxx1,pxy1,pxz1,mmass1,-1,
     1                   vx,vy,vz,vxx,vyy,vzz)
      else
          mmass1 = ksimas
      endif
      pxx2 = (px1+px2) / p00 * xkmas
      pxy2 = (py1+py2) / p00 * xkmas
      pxz2 = (pz1+pz2) / p00 * xkmas

      if (ikaonpot .gt. 0) then
         call gradukaon(-1,x1,y1,z1,pxx2,pxy2,pxz2,mmass2,
     1                  1, vx,vy,vz,vxx,vyy,vzz)
       else
       	    mmass2 = xkmas
      endif
******************************
      srtp = sqrt(s) - (mmass1 + mmass2)
      xx = srtp + .00001
      sigpr = 2.24 * xx**2.16 * exp(-8.23*xx)  ! sigma in fm^2
c------------------------------------------------------------
      radcol = sqrt(sigpr*crossfact/pi)
*     write(*,*)   '  radcol, brel ', radcol, brel, srt, sigpr
      if (brel .gt. radcol)                        goto 600
!       prob  = p_kminu(4,i1)
!       prob  = p_kminu(4,i1)*sigpr/(pi*brel*2)
      prob  = p_kminu(4,i1)
!          write(54,*) "xi+K+" ,srtp, sigpr, id2, iz2
c    ----
      maxloop = (max_ksi/num)
      inksi = (irun-1) * maxloop
      inkmin = 0
      valkapi  = 1000.0
      do 10 kl=1, maxloop
        inksi = inksi + 1
        if (nx_ksi(0,inksi) .eq. 0) goto 12
        if (valkapi .gt. p_ksi(4,inksi)) then
           inkmin = inksi
           valkapi  = p_ksi(4,inksi)
        endif
  10  continue
      inksi = inkmin
       if (prob .lt. valkapi) goto 600
  12  continue
      write(*,*)   '  ksi found in ksi_kminu',inksi
c
  100         xx = 0.5 - rn(iseed)
              yy = 0.5 - rn(iseed)
              zz = 0.5 - rn(iseed)
              pcm = sqrt(0.25*(s-mmass2**2+mmass1**2)**2/s -mmass1**2)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 0.25 .or. r2 .lt. 0.000001) goto 100
              rr = sqrt(r2)
              ps2(0) = (e1 + e2)
              ps2(1) = (px1+px2)
              ps2(2) = (py1+py2)
              ps2(3) = (pz1+pz2)
c
              pn2(1) = pcm * xx/rr
              pn2(2)   = pcm * yy/rr
              pn2(3)   = pcm * zz/rr
              pn2(0)   = sqrt(mmass1**2 + pcm**2)
*             lorentz-transformation into lab frame
              call lorentz_hw(ps2, pn2, pn23)

c-----------------------------------------------------------------------
c           density dependence
          ix = nint(x1)
          iy = nint(y1)
          iz = nint(z1)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
         bin_dens = int(denst*10.)
         bin_time = int(time*1.)
         ksi_KN(bin_dens,bin_time)=ksi_KN(bin_dens,bin_time)+1
c-----------------------------------------------------------------------
c
              nx_ksi(0,inksi) = 1
              nx_ksi(1,inksi) = ireac
              p_ksi(1,inksi) =  pn23(1)
              p_ksi(2,inksi) =  pn23(2)
              p_ksi(3,inksi) =  pn23(3)
              p_ksi(4,inksi) =  (prob/crossfact) * 3.355
*
************************************************************************
*                                                                      *
  800   continue
 1000 continue
      write(*,*) '  end ksi_kminu prod '
      return
      end

************************************************************************
*                                                                      *
      subroutine ksi_km_hyp_eta(ksi_hyp_eL,ksi_hyp_eS)
*                                                                      *
*     purpose: ksi_minus production in  K- + Y -->  KSI + eta            *
*              since K0  are absent
*              it is assumed they are of the same number
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      real*8 radmax, e_thr, ksimas, sig_0 ,mmass2,crossfact
      real*8 vx,vy,vz, vxx,vyy,vzz,p00,pxx1,pxy1,pxz1,mmass1

      parameter (radmax = 0.287*10.) !2.6
      parameter (e_thr  = 1.869) ! ->  massksi + etamass
      parameter (ksimas = 1.321)
      parameter (crossfact = 100.)

      integer bin_dens,bin_time,ix,iy,iz
      integer ksi_hyp_eL(0:999,0:999),ksi_hyp_eS(0:999,0:999)

      real*8 tmas, srt, srtp, sigpr, rsqare, rn, denst
      real*8 radcol, prob, valkapi
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt, b21
      real*8 x1,x2,y1,y2,z1,z2, px1,py1,pz1, px2,py2,pz2, pp, e1,e2
      real*8 em1, em2, em12, em22, dx,dy,dz, s, pcm, xx,yy,zz, r2,rr
      real*8 pn2(0:3), ps2(0:3), pn23(0:3)
      integer irun, ii, i1, iz1, iz2
      integer maxhyp, inhyp, i2, id2, ireac
      integer  inksi, maxloop, inkmin, kl, inkma
c
      maxhyp = max_kminu / num
c     loop over all parallel runs
      do 1000 irun = 1,num
        inkma = (irun-1) * maxhyp
        inhyp = (irun-1) * maxhyp
        do 800 ii  = 1,maxhyp
          i1  = ii + inkma
          if(nx_kminu(0,i1).eq. 0 )                 goto 800
c
          x1         = r_kminu(1,i1)
          y1         = r_kminu(2,i1)
          z1         = r_kminu(3,i1)
          px1        = p_kminu(1,i1)
          py1        = p_kminu(2,i1)
          pz1        = p_kminu(3,i1)
          em1        = xkmas
          pp     = sqrt(px1**2+py1**2+pz1**2)
          em12       = em1**2
          e1         = sqrt(em1**2+px1**2+py1**2+pz1**2)
c
          i2  = inhyp
  600     i2  = i2 + 1
          if (i2 .gt. inhyp + maxhyp)                      goto 800
          if (nx_hyp(0,i2) .eq. 0)     goto 600
          iz2 = nx_hyp(1,i2)
          if (iz2  .ne. 0) goto  600    !only charge zero is allowed for the hyperon when producing an eta
          x2 = r_hyp(1,i2)
          dx     = x1 - x2
          if (abs(dx) .gt. radmax)                              goto 600
          y2 = r_hyp(2,i2)
          dy     = y1 - y2
          if (abs(dy) .gt. radmax)                              goto 600
          z2 = r_hyp(3,i2)
          dz     = z1 - z2
          if (abs(dz) .gt. radmax)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. radmax**2)                            goto 600
*         now particles are close enough to each other !
          id2  = nx_hyp(0,i2) - 1
          if (id2 .eq. 0) tmas = xlmas
          if (id2 .eq. 1) tmas = xsmas
          px2    = p_hyp(1,i2)
          py2    = p_hyp(2,i2)
          pz2    = p_hyp(3,i2)
          em2    = tmas
          em22   = em2**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                       - (pz1+pz2)**2
          if (s .lt. e_thr**2)  goto 600
*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if (brel .gt. radmax)                                goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
c------------------------------------------------------
*            OK, there is a collision
c------------------------------------------------------
          srt = sqrt(s)
** modified medium masses mmass1,2 of the reeaction products Xi- and eta ***
      p00 = sqrt((e1+e2)**2 - (px1+px2)**2 -(py1+py2)**2 -(pz1+pz2)**2)
      pxx1 = (px1+px2) / p00 * ksimas
      pxy1 = (py1+py2) / p00 * ksimas
      pxz1 = (pz1+pz2) / p00 * ksimas

      if (i_ksi_pot .ge. 1) then
         call graduxi(-1, x1,y1,z1,pxx1,pxy1,pxz1,mmass1,-1,
     1                   vx,vy,vz,vxx,vyy,vzz)
      else
          mmass1 = ksimas
      endif
      mmass2 = emass
******************************
        srtp = sqrt(s) - (mmass1 + mmass2)+ .00001
        if (id2 .eq. 0) then
          sigpr = 0.24 - (0.20*exp(-50.*srtp)) ! sigma in fm^2
	else
          sigpr = 0.26 - (0.25*exp(-8.0*srtp)) ! sigma in fm^2
	endif
c------------------------------------------------------------
      radcol = sqrt(sigpr*crossfact/pi)
*     write(*,*)   '  radcol, brel ', radcol, brel, srt, sigpr
      if (brel .gt. radcol)                        goto 600
!       prob  = p_hyp(4, i2) * p_kminu(4,i1)
      prob  = p_hyp(4, i2)*p_kminu(4,i1)
      ireac = 19 + id2
!          write(54,*) "xi+eta" ,srtp, sigpr, id2, iz2
c    ----
      maxloop = (max_ksi/num)
      inksi = (irun-1) * maxloop
      inkmin = 0
      valkapi  = 1000.0
      do 10 kl=1, maxloop
        inksi = inksi + 1
        if (nx_ksi(0,inksi) .eq. 0) goto 12
        if (valkapi .gt. p_ksi(4,inksi)) then
           inkmin = inksi
           valkapi  = p_ksi(4,inksi)
        endif
  10  continue
      inksi = inkmin
       if (prob .lt. valkapi) goto 600
  12  continue
c     write(*,*)   '  ksi found in ksi_km_hyp',inksi
c
  100         xx = 0.5 - rn(iseed)
              yy = 0.5 - rn(iseed)
              zz = 0.5 - rn(iseed)
              pcm = sqrt(0.25*(s-mmass2**2+mmass1**2)**2/s -mmass1**2)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 0.25 .or. r2 .lt. 0.000001) goto 100
              rr = sqrt(r2)
              ps2(0) = (e1 + e2)
              ps2(1) = (px1+px2)
              ps2(2) = (py1+py2)
              ps2(3) = (pz1+pz2)
c
              pn2(1) = pcm * xx/rr
              pn2(2)   = pcm * yy/rr
              pn2(3)   = pcm * zz/rr
              pn2(0)   = sqrt(mmass1**2 + pcm**2)
*             lorentz-transformation into lab frame
              call lorentz_hw(ps2, pn2, pn23)

c-----------------------------------------------------------------------
c           density dependence
          ix = nint(x1)
          iy = nint(y1)
          iz = nint(z1)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
         bin_dens = int(denst*10.)
         bin_time = int(time*1.)
          if (id2 .eq. 0) then
           ksi_hyp_eL(bin_dens,bin_time)=ksi_hyp_eL(bin_dens,bin_time)+1
       endif
       if (id2 .eq. 1) then
         ksi_hyp_eS(bin_dens,bin_time)=ksi_hyp_eS(bin_dens,bin_time)+1
       endif
c-----------------------------------------------------------------------
c
              nx_ksi(0,inksi) = 1
              nx_ksi(1,inksi) = ireac
              p_ksi(1,inksi) =  pn23(1)
              p_ksi(2,inksi) =  pn23(2)
              p_ksi(3,inksi) =  pn23(3)
              p_ksi(4,inksi) =  (prob/crossfact)  * 3.355
*
************************************************************************
*                                                                      *
  800   continue
 1000 continue
      write(*,*) '  end ksi_km_hyp prod '
      return
      end

************************************************************************
*                                                                      *
      subroutine ksi_km_hyp(ksi_hyp_piL,ksi_hyp_piS)
*                                                                      *
*     purpose: ksi_minus production in  K- + Y -->  KSI + pi            *
*              since K0  are absent
*              it is assumed they are of the same number
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      real*8 radmax, e_thr, ksimas, sig_0 ,mmass2, denst,crossfact
      real*8 vx,vy,vz, vxx,vyy,vzz,p00,pxx1,pxy1,pxz1,mmass1

      parameter (radmax = 0.564*10.)  !10mb ! -> 5mb (max-cross section)  :radmax = 0.4
      parameter (e_thr  = 1.459) ! ->  massksi + pionmass
      parameter (ksimas = 1.321)
      parameter (sig_0 = 1.20)
      parameter (crossfact = 100.)

      integer bin_dens,bin_time,ix,iy,iz
      integer ksi_hyp_piL(0:999,0:999),ksi_hyp_piS(0:999,0:999)

      real*8 tmas, srt, srtp, sigpr, rsqare, rn
      real*8 radcol, prob, valkapi
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt, b21
      real*8 x1,x2,y1,y2,z1,z2, px1,py1,pz1, px2,py2,pz2, pp, e1,e2
      real*8 em1, em2, em12, em22, dx,dy,dz, s, pcm, xx,yy,zz, r2,rr
      real*8 pn2(0:3), ps2(0:3), pn23(0:3)
      integer irun, ii, i1, iz2
      integer maxhyp, inhyp, i2, id2, ireac
      integer  inksi, maxloop, inkmin, kl, inkma
c
      maxhyp = max_kminu / num
c     loop over all parallel runs
      do 1000 irun = 1,num
        inkma = (irun-1) * maxhyp
        inhyp = (irun-1) * maxhyp
        do 800 ii  = 1,maxhyp
          i1  = ii + inkma
          if(nx_kminu(0,i1) .eq. 0 )                        goto 800
c
          x1         = r_kminu(1,i1)
          y1         = r_kminu(2,i1)
          z1         = r_kminu(3,i1)
          px1        = p_kminu(1,i1)
          py1        = p_kminu(2,i1)
          pz1        = p_kminu(3,i1)
          em1        = xkmas
          pp     = sqrt(px1**2+py1**2+pz1**2)
          em12       = em1**2
          e1         = sqrt(em1**2+px1**2+py1**2+pz1**2)
c
          i2  = inhyp
  600     i2  = i2 + 1
          if (i2 .gt. inhyp + maxhyp)                           goto 800
          if (nx_hyp(0,i2) .eq. 0)                              goto 600
          iz2 = nx_hyp(1,i2)
!          if (iz2  .eq. 1)                  goto 600
          x2 = r_hyp(1,i2)
          dx     = x1 - x2
          if (abs(dx) .gt. radmax)                              goto 600
          y2 = r_hyp(2,i2)
          dy     = y1 - y2
          if (abs(dy) .gt. radmax)                              goto 600
          z2 = r_hyp(3,i2)
          dz     = z1 - z2
          if (abs(dz) .gt. radmax)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. radmax**2)                            goto 600
*         now particles are close enough to each other !
          id2  = nx_hyp(0,i2) - 1
          if (id2 .eq. 0) tmas = xlmas
          if (id2 .eq. 1) tmas = xsmas
          px2    = p_hyp(1,i2)
          py2    = p_hyp(2,i2)
          pz2    = p_hyp(3,i2)
          em2    = tmas
          em22   = em2**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                       - (pz1+pz2)**2
          if (s .lt. e_thr**2)                                goto 600
!       always possible cause of exothermic reaction!
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if (brel .gt. radmax)                                	goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            	goto 600
c------------------------------------------------------
*            OK, there is a collision
c------------------------------------------------------
c     write(*,*)' ksi production by K-', irun,i1,i2,nx_kminu(1,i1),
c    1           nx_kminu(2,i1),nx_kminu(3,i1),nx_kminu(2,i1),
c    2           x1,x2,y1,y2,z1,z2
          srt = sqrt(s)
*** modified medium masses mmass1,2 of the reeaction products Xi- and pi ***
      p00 = sqrt((e1+e2)**2 - (px1+px2)**2 -(py1+py2)**2 -(pz1+pz2)**2)
      pxx1 = (px1+px2) / p00 * ksimas
      pxy1 = (py1+py2) / p00 * ksimas
      pxz1 = (pz1+pz2) / p00 * ksimas

      if (i_ksi_pot .ge. 1) then
         call graduxi(-1, x1,y1,z1,pxx1,pxy1,pxz1,mmass1,-1,
     1                   vx,vy,vz,vxx,vyy,vzz)
      else
          mmass1 = ksimas
      endif
      mmass2 = pmass
******************************
          srtp = sqrt(s) - (mmass1 + mmass2) + .00001
!         sigpr = sig_0 * srtp**(-0.308)       ! sigma in fm^2
        if (id2 .eq. 0) then
!         if (id2 .eq. 1) then
!            if (iz2 .eq. -1)  sigpr = 2./3. * sigpr
!            if (iz2 .eq.  0)  sigpr = 1./3. * sigpr
!         endif
        sigpr =(3.47/4.)*(0.9/srtp**(0.2))*
     &                  ((srtp/1.611)+1.)**(-2.)       ! sigma in fm^2
      else
        sigpr =(31.8/12.)*(exp(-140.*(srtp-0.03))+1.)*
     &  (1.-((srtp/1.688)+1.)**(-2.))**0.6 *((srtp/1.688)+1.)**(-3.4)! sigma in fm^2
      endif
c------------------------------------------------------------
      radcol = sqrt(sigpr*crossfact/pi)
*     write(*,*)   '  radcol, brel ', radcol, brel, srt, sigpr
      if (brel .gt. radcol)                        goto 600
      prob  = p_hyp(4, i2) * p_kminu(4,i1)
      ireac = 5 + id2
!          write(54,*) "xi+pi" ,srtp, sigpr, id2, iz2
c    ----
      maxloop = (max_ksi/num)
      inksi = (irun-1) * maxloop
      inkmin = 0
      valkapi  = 1000.0
      do 10 kl=1, maxloop
        inksi = inksi + 1
        if (nx_ksi(0,inksi) .eq. 0) goto 12
        if (valkapi .gt. p_ksi(4,inksi)) then
           inkmin = inksi
           valkapi  = p_ksi(4,inksi)
        endif
  10  continue
      inksi = inkmin
       if (prob .lt. valkapi) goto 600
  12  continue
c     write(*,*)   '  ksi found in ksi_km_hyp',inksi
c

  100         xx = 0.5 - rn(iseed)
              yy = 0.5 - rn(iseed)
              zz = 0.5 - rn(iseed)
              pcm = sqrt(0.25*(s-mmass2**2+mmass1**2)**2/s -mmass1**2)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 0.25 .or. r2 .lt. 0.000001) goto 100
              rr = sqrt(r2)
              ps2(0) = (e1 + e2)
              ps2(1) = (px1+px2)
              ps2(2) = (py1+py2)
              ps2(3) = (pz1+pz2)
c
              pn2(1) = pcm * xx/rr
              pn2(2)   = pcm * yy/rr
              pn2(3)   = pcm * zz/rr
              pn2(0)   = sqrt(mmass1**2 + pcm**2)
*             lorentz-transformation into lab frame
              call lorentz_hw(ps2, pn2, pn23)

c-----------------------------------------------------------------------
c           density dependence
          ix = nint(x1)
          iy = nint(y1)
          iz = nint(z1)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      if (id2 .eq. 0) then
         ksi_hyp_piL(bin_dens,bin_time)=ksi_hyp_piL(bin_dens,bin_time)+1
      endif
      if (id2 .eq. 1) then
    	 ksi_hyp_piS(bin_dens,bin_time)=ksi_hyp_piS(bin_dens,bin_time)+1
      endif
c-----------------------------------------------------------------------
c
              nx_ksi(0,inksi) = 1
              nx_ksi(1,inksi) = ireac
              p_ksi(1,inksi) =  pn23(1)
              p_ksi(2,inksi) =  pn23(2)
              p_ksi(3,inksi) =  pn23(3)
              p_ksi(4,inksi) =  (prob/crossfact) * 3.355
c
*
************************************************************************
*                                                                      *
  800   continue
 1000 continue
      write(*,*) '  end ksi_km_hyp prod '
      return
      end

************************************************************************
*                                                                      *
      subroutine ksi_pion(ksi_pi_L,ksi_pi_S)
*                                                                      *
*     purpose: ksi_minus production in  pi + Y -->  KSI + K+           *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      real*8 radmax, e_thr, ksimas, sig_lam, sig_sig,crossfact
      real*8 vx,vy,vz, vxx,vyy,vzz, pxx1, pxy1,pxz1,p00,mmass1
      real*8 pxx2, pxy2,pxz2,mmass2, denst, ka,ma,mb,mc

      parameter (radmax = 0.3989*10.)  !  5mb
      parameter (e_thr  = 1.816) ! ->  massksi + xkmas
      parameter (ksimas = 1.321)
      parameter (sig_lam = 1.3, sig_sig = 2.0)
      parameter (crossfact = 100.)

      integer bin_dens,bin_time,ix,iy,iz
      integer ksi_pi_L(0:999,0:999),ksi_pi_S(0:999,0:999)

      real*8 tmas, srt, srtp, sigpr, rsqare, rn
      real*8 radcol, prob, valkapi
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt, b21
      real*8 x1,x2,y1,y2,z1,z2, px1,py1,pz1, px2,py2,pz2, pp, e1,e2
      real*8 em1, em2, em12, em22, dx,dy,dz, s, pcm, xx,yy,zz, r2,rr
      real*8 pn2(0:3), ps2(0:3), pn23(0:3)
      integer irun, inp, ii, i1, iz1, iz12
      integer maxhyp, inhyp, i2, id2, ireac
      integer inksi, maxloop, inkmin, kl
c
      maxhyp = max_kminu / num
c     loop over all parallel runs
      do 1000 irun = 1,num
        inp   = (irun-1) * maxp
        inhyp = (irun-1) * maxhyp
        do 800 ii  = 1,maxp
          i1  = ii + inp
          if(ipi(1,i1) .ne. 1 )                         goto 800
c
          iz1  =  ipi(2,i1)
          x1         = rpi(1,i1)
          y1         = rpi(2,i1)
          z1         = rpi(3,i1)
          px1        = ppi(1,i1)
          py1        = ppi(2,i1)
          pz1        = ppi(3,i1)
          em1        = epi(i1)
          pp     = sqrt(px1**2+py1**2+pz1**2)
          em12       = epi(i1)**2
          e1         = sqrt(em1**2+px1**2+py1**2+pz1**2)
c
          i2  = inhyp
  600     i2  = i2 + 1
          if (i2 .gt. inhyp + maxhyp)                      goto 800
          if (nx_hyp(0,i2) .eq. 0)     goto 600
          iz12 = nx_hyp(1,i2) + iz1
          if (iz12 .ne. 0)             goto 600      !only pi+Y-, pi- Y+ or pi0 Y0 are allowed
!           if (iz12 .lt.-1)             goto 600
          x2 = r_hyp(1,i2)
          dx     = x1 - x2
          if (abs(dx) .gt. radmax)                              goto 600
          y2 = r_hyp(2,i2)
          dy     = y1 - y2
          if (abs(dy) .gt. radmax)                              goto 600
          z2 = r_hyp(3,i2)
          dz     = z1 - z2
          if (abs(dz) .gt. radmax)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. radmax**2)                            goto 600
*         now particles are close enough to each other !
          id2  = nx_hyp(0,i2) - 1
          if (id2 .eq. 0) tmas = xlmas
          if (id2 .eq. 1) tmas = xsmas
          px2    = p_hyp(1,i2)
          py2    = p_hyp(2,i2)
          pz2    = p_hyp(3,i2)
          em2    = tmas
          em22   = em2**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                       - (pz1+pz2)**2
          if (s .lt. e_thr**2)  goto 600
*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if (brel .gt. radmax)                                goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
c------------------------------------------------------
*            OK, there is a collision
c------------------------------------------------------
c     write(*,*)' ksi production by pion', irun,i1,i2,nx_kminu(1,i1),
c    1           nx_kminu(2,i1),nx_kminu(3,i1),nx_kminu(2,i1),
c    2           x1,x2,y1,y2,z1,z2
          srt = sqrt(s)
** modified medium masses mmass1,2 of the reeaction products Xi- and K+/0 **
      p00 = sqrt((e1+e2)**2 - (px1+px2)**2 -(py1+py2)**2 -(pz1+pz2)**2)
      pxx1 = (px1+px2) / p00 * ksimas
      pxy1 = (py1+py2) / p00 * ksimas
      pxz1 = (pz1+pz2) / p00 * ksimas

      if (i_ksi_pot .ge. 1) then
         call graduxi(-1, x1,y1,z1,pxx1,pxy1,pxz1,mmass1,-1,
     1                   vx,vy,vz,vxx,vyy,vzz)
      else
          mmass1 = ksimas
      endif
      pxx2 = (px1+px2) / p00 * xkmas
      pxy2 = (py1+py2) / p00 * xkmas
      pxz2 = (pz1+pz2) / p00 * xkmas

      if (ikaonpot .gt. 0) then
         call gradukaon(-1,x1,y1,z1,pxx2,pxy2,pxz2,mmass2,
     1                  1, vx,vy,vz,vxx,vyy,vzz)
       else
       	    mmass2 = xkmas
      endif
******************************
      srtp = sqrt(s) -  (mmass1 + mmass2) + .00001

      ma = 0.494
      mb = 1.116
      mc = 0.135

      ka = sqrt( ((s-(ma+mb)**2)*(s-(ma-mb)**2))/
     1               ((s-(mc+mb)**2)*(s-(mc-mb)**2)))

      if (id2 .eq. 0) then    !only Lambda
!            sigpr = (1./3.)*sig_lam * srtp**0.465 !  sigma in fm**2
        sigpr=ka*(3.47/4.)*(0.9/srtp**(0.2))*
     1             ((srtp/1.611)+1.)**(-2.)       ! sigma in fm^2
!           sigpr = sigpr * real(1-iz1) / 3.
         elseif (id2 .eq. 1) then      !only all Sigmas +-0
!            sigpr = sig_sig * srtp**0.85    !  sigma in fm**2
!            sigpr  = (1. - 0.2*iz1) * sigpr
        sigpr=ka*(31.8/12.)*(exp(-140.*(srtp-0.03))+1.)*
     &  (1.-((srtp/1.688)+1.)**(-2.))**0.6 *((srtp/1.688)+1.)**(-3.4) ! sigma in fm^2
      elseif (id2 .gt. 1 .or. id2 .lt. 0) then
          sigpr  = 0.0
      endif
c------------------------------------------------------------
      radcol = sqrt(sigpr*crossfact/pi)
*     write(*,*)   '  radcol, brel ', radcol, brel, srt, sigpr
      if (brel .gt. radcol)                        goto 600
!       prob  = p_hyp(4, i2) * sigpr / (pi * radmax**2)
      prob  = p_hyp(4, i2)
!       prob  = p_hyp(4, i2)
      ireac = 10 + id2
!          write(54,*) "xiK+p" ,srtp, sigpr, id2
c    ----
      maxloop = (max_ksi/num)
      inksi = (irun-1) * maxloop
      inkmin = 0
      valkapi  = 1000.0
      do 10 kl=1, maxloop
        inksi = inksi + 1
        if (nx_ksi(0,inksi) .eq. 0) goto 12
        if (valkapi .gt. p_ksi(4,inksi)) then
           inkmin = inksi
           valkapi  = p_ksi(4,inksi)
        endif
  10  continue
      inksi = inkmin
      if (prob .lt. valkapi) goto 600
  12  continue
      write(*,*)   '  ksi found in ksi_pion',inksi
c
  100         xx = 0.5 - rn(iseed)
              yy = 0.5 - rn(iseed)
              zz = 0.5 - rn(iseed)
              pcm = sqrt(0.25*(s-mmass2**2+mmass1**2)**2/s -mmass1**2)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 0.25 .or. r2 .lt. 0.000001) goto 100
              rr = sqrt(r2)
              ps2(0) = (e1 + e2)
              ps2(1) = (px1+px2)
              ps2(2) = (py1+py2)
              ps2(3) = (pz1+pz2)
c
              pn2(1) = pcm * xx/rr
              pn2(2)   = pcm * yy/rr
              pn2(3)   = pcm * zz/rr
              pn2(0)   = sqrt(mmass1**2 + pcm**2)
*             lorentz-transformation into lab frame
              call lorentz_hw(ps2, pn2, pn23)
c
c-----------------------------------------------------------------------
c           density dependence
          ix = nint(x1)
          iy = nint(y1)
          iz = nint(z1)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      if (id2 .eq. 0) then
         ksi_pi_L(bin_dens,bin_time)=ksi_pi_L(bin_dens,bin_time)+1
      endif
      if (id2 .eq. 1) then
    	 ksi_pi_S(bin_dens,bin_time)=ksi_pi_S(bin_dens,bin_time)+1
      endif
c-----------------------------------------------------------------------

              nx_ksi(0,inksi) = 1
              nx_ksi(1,inksi) = ireac
              p_ksi(1,inksi) =  pn23(1)
              p_ksi(2,inksi) =  pn23(2)
              p_ksi(3,inksi) =  pn23(3)
              p_ksi(4,inksi) =  (prob/crossfact)  * 3.355
*
************************************************************************
*                                                                      *
  800   continue
 1000 continue
      write(*,*) '  end ksi_pion prod '
      return
      end


