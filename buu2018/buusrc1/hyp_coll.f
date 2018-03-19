************************************************************************
*                                                                      *
      subroutine hyp_coll
*              collisions of baryons and  hyperons                     *
*              includes also K- production                             *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      common /nthwhw/  nthw
      integer nthw
      real*8 pirk,s,pcm,xx,yy,zz,r2,rr,betax,betay,betaz,gamma
      real*8 x1, y1, z1, px1, py1, pz1, em1, em12, e1
      real*8 dx, dy, dz, rsqare, px2, py2, pz2, em2, em22,e2, x2,y2,z2
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2,ddlt,b21,dxm
      real*8 rn,phase,dxp
      real*8 q1(0:3), q2(0:3), us(0:3), pn2s(0:3), pn1s(0:3), pn2(0:3)
      real*8 sig00, radka, s0, s1, ratio, ax1r
      real*8 sqs, em3, pcm3, pcm32, em1r, delhyp
      integer nx0, nx1, nx0r, nx1r, ntag, nx2,nx3
      integer maxk, irun, ink, inkk, ini, i1, i2, ii
      real*8  thr_kmin, valkapi, pp, pkmax2, pkmax, pka, ranp, sigk
      integer  inkmin, kl, inkrun, ich
      real*8  vx,vy,vz,vxx,vyy,vzz, pxx, pxy, pxz, mmass, dkmass
      parameter (pirk   = 0.798)  !  30 mb
c     parameter (pirk   = 0.977)  !  30 mb
      parameter (delhyp = 2.5  )  ! 200mb
*----------------------------------------------------------------------*
*    pirk = 0.565 fm                corresponds  to 10.2 mb            *
*----------------------------------------------------------------------*
      write(*,*) 'call hyp_coll'
c     return
      maxk = max_kminu/num
      kanum = 0
*     loop over all parallel runs
      do 1000 irun = 1,num
        ink = (irun-1) * maxk
        ini = (irun-1) * maxb
        do 800 ii  = 1,maxk
          i1  = ii + ink
          nx0 = nx_hyp(0,i1)
          nx1 = nx_hyp(1,i1)
          if(nx0 .eq. 0 )                            goto 800
          if(nx0 .gt. 2 ) write(*,*) 'hiba hyp_coll',nx0
          kanum = kanum + 1
          x1  = r_hyp(1,i1)
          y1  = r_hyp(2,i1)
          z1  = r_hyp(3,i1)
          px1 = p_hyp(1,i1)
          py1 = p_hyp(2,i1)
          pz1 = p_hyp(3,i1)
          if (nx0 .eq. 1) em1 = xlmas
          if (nx0 .eq. 2) em1 = xsmas
          em12= em1**2
          e1  = sqrt( em12 + px1**2 + py1**2 + pz1**2 )
*     look for an scattering pseudonucleon in the same run
          i2  = ini
  600     i2  = i2 + 1
          if(i2 .gt. ini+maxb)                                 goto 800
          if(i2.eq.nx_hyp(2,i1) .or. i2.eq.nx_hyp(3,i1) .or.
     &      id(1,i2).eq.0)       goto 600
          x2 = r(1,i2)
          dx     = x1 - x2
            if(nbound.eq.1) then
              dxp = dx+2.0*boxx
              dxm = dx-2.0*boxx
              if(abs(dx) .gt. abs(dxp)) dx=dxp
              if(abs(dx) .gt. abs(dxm)) dx=dxm
            end if
          if (abs(dx) .gt. delhyp)                              goto 600
          y2 = r(2,i2)
          dy     = y1 - y2
            if(nbound.eq.1) then
              dxp = dy+2.0*boxx
              dxm = dy-2.0*boxx
              if(abs(dy) .gt. abs(dxp)) dy=dxp
              if(abs(dy) .gt. abs(dxm)) dy=dxm
            end if
          if (abs(dy) .gt. delhyp)                              goto 600
          z2 = r(3,i2)
          dz     = z1 - z2
            if(nbound.eq.1) then
              dxp = dz+2.0*boxz
              dxm = dz-2.0*boxz
              if(abs(dz) .gt. abs(dxp)) dz=dxp
              if(abs(dz) .gt. abs(dxm)) dz=dxm
            end if
          if (abs(dz) .gt. delhyp)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. delhyp**2)                            goto 600
*         now particles are close enough to each other !
          px2    = p(1,i2)
          py2    = p(2,i2)
          pz2    = p(3,i2)
          em2    = e(i2)
          em22   = e(i2)**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                           - (pz1+pz2)**2
*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          pcm = sqrt(0.25*(s-em1**2+em2**2)**2/s -em2**2)
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          sqs = sqrt(s)
          s0  = sqs-em1-em2
          sig00 = .099/s0**1.39 + 8.8 + 105.*s0 -93.*s0*s0
          sig00 = min (200.,sig00)
          sig00 = max(sig00,30.)
          radka =  sqrt(0.1*sig00 / pi)
          if (brel .gt. radka)                                  goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
*
*   now  there  will a collision in this time step
************************************************************************
          if (nx0 .eq. 1) then
             nx0r = 2
             ax1r = 2.98 * rn(iseed)-1.49
             nx1r = nint(ax1r)
             em1r = xsmas
          else
             nx0r = 1
             nx1r = 0
             em1r = xlmas
          endif
          pcm3  = .0
          pcm32 =  0.25*(s-em1r**2+em2**2)**2/s -em2**2
          if (pcm32 .gt. .0) pcm3 = sqrt(pcm32)
          ratio =  pcm3/ pcm * 0.5
         if (rn(iseed) .lt. ratio) then
              nx2 = nx0r
              nx3 = nx1r
              em3 = em1r
              pcm = pcm3
         else
              nx2 = nx0
              nx3 = nx1
              em3 = em1
         endif
c         if (i1 .eq. 10004)  then
c         write(*,*) ' sigma+ coll1 ', i2, nx0, em3,px1,py1,pz1,
c    1     em2, px2,py2,pz2,  sqs, ax1r, nx1r,nx0r, pcm,pcm3,em1,em3,
c    2     '  pcm ', s, e1, e2 , p12
c         endif
c       --------------     end inelastic
          betax = (px1+px2) / (e1+e2)
          betay = (py1+py2) / (e1+e2)
          betaz = (pz1+pz2) / (e1+e2)
          gamma  = 1.0 / sqrt(1.0-betax**2-betay**2-betaz**2)
          us(0) = gamma
          us(1) = gamma * betax
          us(2) = gamma * betay
          us(3) = gamma * betaz
          if (i_kminu_pot  .gt. 0) then
             ich = -1
             pxx = us(1) * xkmas
             pxy = us(2) * xkmas
             pxz = us(3) * xkmas
             call gradukao2(-1, x1,y1,z1, pxx,pxy,pxz,mmass,ich,
     1                          vx,vy,vz,vxx,vyy,vzz)
          else
             mmass = xkmas
          endif
          thr_kmin = 2.*rmass  + mmass
          dkmass   = xkmas - mmass
c------------------------------------------
c-------------------   kminus production -------
c-----------------------------------------
c     goto 2000
c-----------------------------------------
      if (sqs .lt. thr_kmin)   goto  2000
      s1 = (sqs+dkmass)**2
      s0 = (thr_kmin+dkmass)**2
      if (nx0 .eq. 1) sigk = 18.*(s1/s0-1.)**1.90*(s0/s1)**4.93
      if (nx0 .eq. 2) sigk = 19.*(s1/s0-1.)**1.63*(s0/s1)**5.86
      pp   = sigk / (sigk+sig00) * p_hyp(4,i1)
      inkrun = (max_kminu/num)
      inkk = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 10 kl=1,inkrun                    !        iphi_bb       hw
        inkk = inkk + 1
        if (nx_kminu(0,inkk) .eq. 0) goto 12
        if (valkapi .gt. p_kminu(4,inkk)) then
           inkmin = inkk
           valkapi  = p_kminu(4,inkk)
        endif
  10  continue
      inkk = inkmin
c     write(*,*)  ' sqs - thr_kmin ', sqs ,thr_kmin, pp, valkapi
       if (pp .lt. valkapi) goto 2000
  12  continue
        pkmax2=((s1-4.*rmass**2-xkmas**2)**2 - (4.*rmass*xkmas)**2)
     1         / (4.*s1)
c     write(*,*)   '  kminus found in kminu_coll',inkk, pp, thr_kmin,
c    1             xkmas,mmass,rmass,s,s1, pkmax2
        if (pkmax2 .le. .0) goto  2000
        pkmax = sqrt(pkmax2)
   15  continue
          ranp     = rn(iseed)
          if((ranp**2-ranp**3)*6.75 .lt. rn(iseed)) go to 15
          pka = pkmax * ranp
  31      xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 31
          pn2(1)    =  pka *xx/rr   ! the second nucleon
          pn2(2)    =  pka *yy/rr
          pn2(3)    =  pka *zz/rr
          pn2(0)    = sqrt(xkmas**2+pn2(1)**2+pn2(2)**2+pn2(3)**2)
          call lorentz_hw(us , pn2, pn2s)
c
          p_kminu(1,inkk)= pn2s(1)
          p_kminu(2,inkk)= pn2s(2)
          p_kminu(3,inkk)= pn2s(3)
          p_kminu(4,inkk)= pp
          r_kminu(1,inkk)= x1
          r_kminu(2,inkk)= y1
          r_kminu(3,inkk)= z1
          nx_kminu(0,inkk) = 1
          nx_kminu(1,inkk) = 7
          nx_kminu(2,inkk) = i2
          nx_kminu(3,inkk) = i1
c            write(*,*)  '  kminus found in NY', inkk, pp
c
 2000     continue
  100         xx = 1.0 - 2.*rn(iseed)
              yy = 1.0 - 2.*rn(iseed)
              zz = 1.0 - 2.*rn(iseed)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 1.0 .or. r2 .lt. 0.000001) goto 100
              rr = sqrt(r2)
c
              q1(1)  = pcm * xx/rr
              q1(2)  = pcm * yy/rr
              q1(3)  = pcm * zz/rr
c
              q2(1)   = -pcm * xx/rr
              q2(2)   = -pcm * yy/rr
              q2(3)   = -pcm * zz/rr
*             lorentz-transformation into lab frame
              q2(0)  = sqrt (em2**2 + q2(1)**2 + q2(2)**2 + q2(3)**2)
              call lorentz_hw(us , q2, pn2s)
              ntag = 0
              if(id(1,i2).eq.1 .and. ipauli.eq.1)
     &           call pauli(i2,ntag,iseed,phase,r(1,i2),r(2,i2),r(3,i2),
     &                         pn2s(1),pn2s(2),pn2s(3))
*
c             write(*,*) '  pauli in kaoncoll ', i1, i2, ntag
              if(ntag .eq. -1) go to 600
c
              q1(0)  = sqrt (em3**2 + q1(1)**2 + q1(2)**2 + q1(3)**2)
              call lorentz_hw(us , q1, pn1s)
              p_hyp(1,i1) = pn1s(1)
              p_hyp(2,i1) = pn1s(2)
              p_hyp(3,i1) = pn1s(3)
              nx_hyp(0,i1) = nx2
              nx_hyp(1,i1) = nx3
              nx_hyp(2,i1) = i2
              nx_hyp(3,i1) = i2
              icollh = icollh + 1
*
************************************************************************
*                                                                      *
              if (i1 .eq. 10004) write(*,*) ' sigma+ coll ',p_hyp(1,i1)
  800   continue
 1000 continue
      return
      end
