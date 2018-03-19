************************************************************************
*                                                                      *
      subroutine kaoncoll(collkplus,collkplusd,colkpl)
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
      common /nthwhw/  nthw, isu_hw
      integer nthw, isu_hw

      real*8 s,pcm,xx,yy,zz,r2,rr,betax,betay,betaz,gamma
      real*8 x1,y1,z1,px1, py1, pz1, em1, em12, e1, e1cm, p1beta,transf
      real*8 dx, dy, dz, rsqare, px2, py2, pz2, em2, em22, e2, x2, y2,z2
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt, b21,dxm
      real*8 rn,phase,p2beta,e2cm,qx1,qy1,qz1,qx2,qy2,qz2,dxp
      real*8 pklab, prel, radka, sig_rn, sig_exc
      real*8 q0, xmpion2,xmnucl2,xmkaon2, pmax, pk, ek ,denst
      real*8 rx, x_s, xm2, sqs, sig_in, sig_el, sig_t, ff, fmax, delkao
      integer ntag, chk, chn, ch2, ch_fin, scatt,countcolk
      integer collkplus(2*maxkaon)
      integer collkplusd(2*maxkaon)
      integer ixn,iyn,izn
      integer maxk, irun, ink, ini, i1, i2, ii, ix, iy, iz, jj
      integer ikpel,ikpinel
      parameter (delkao = 1.5  )
! 	save icollk,ikpel,ikpinel
*----------------------------------------------------------------------*
*    pirk = 0.565 fm                corresponds  to 10.2 mb            *
*----------------------------------------------------------------------*
***      write(*,*) 'call kaoncoll'
c     return
      maxk = maxkaon/num
      kanum = 0

c  222 format(5e12.4)

*     loop over all parallel runs
      do 1000 irun = 1,num
        ink = (irun-1) * maxk
        ini = (irun-1) * maxb
        do 800 ii  = 1,maxk
          i1  = ii + ink
          if(ika(1,i1) .eq. 0)          goto 800
          chk  = 2 - ika(1,i1)
          kanum = kanum + 1
          x1  = rkao(1,i1)
          y1  = rkao(2,i1)
          z1  = rkao(3,i1)
          px1 = pkao(1,i1)
          py1 = pkao(2,i1)
          pz1 = pkao(3,i1)
          em1 = xkmas
          em12= em1**2
          e1  = sqrt( em12 + px1**2 + py1**2 + pz1**2 )
*     look for an scattering pseudonucleon in the same run
          i2  = ini
  600     i2  = i2 + 1
          if(i2 .gt. ini+maxb)     goto 800
         if(i2.eq.ika(3,i1).or.i2.eq.ika(4,i1).or.id(1,i2).eq.0)goto 600
          chn   =   id(2,i2)
          if (chn .lt. 0)  goto 600
          if (chn .gt. 1)  goto 600
          x2 = r(1,i2)
          dx     = x1 - x2
            if(nbound.eq.1) then
              dxp = dx+2.0*boxx
              dxm = dx-2.0*boxx
              if(abs(dx) .gt. abs(dxp)) dx=dxp
              if(abs(dx) .gt. abs(dxm)) dx=dxm
            end if
          if (abs(dx) .gt. delkao)                              goto 600
          y2 = r(2,i2)
          dy     = y1 - y2
            if(nbound.eq.1) then
              dxp = dy+2.0*boxx
              dxm = dy-2.0*boxx
              if(abs(dy) .gt. abs(dxp)) dy=dxp
              if(abs(dy) .gt. abs(dxm)) dy=dxm
            end if
          if (abs(dy) .gt. delkao)                              goto 600
          z2 = r(3,i2)
          dz     = z1 - z2
            if(nbound.eq.1) then
              dxp = dz+2.0*boxz
              dxm = dz-2.0*boxz
              if(abs(dz) .gt. abs(dxp)) dz=dxp
              if(abs(dz) .gt. abs(dxm)) dz=dxm
            end if
          if (abs(dz) .gt. delkao)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. delkao**2)                            goto 600
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
          sqs = sqrt(s)
c --------------------
          ch2  = abs(chk + chn - 1)
          if (ch2 .eq. 1) then
            if (sqs .le. 1.8)  then
               x_s = sqs - 1.432
               sig_el =  11.35 + x_s * (9.49 - 68.7*x_s*x_s)  ! cross sections in mb
            else
               x_s    = log(sqs/1.8)
               sig_el = 11.41 + x_s * (-24.4 + 18.6*x_s)
            endif
            sig_in = .0
            x_s   = sqs - 1.57
            if (x_s .gt. 0)
     1      sig_in  = 15.04 - 2.04 / ((log(sqs))**2-.15)
            sig_in  = max(.0, sig_in)
            sig_exc = .0
            sig_t   = sig_el + sig_in
          else
            if (sqs .le. 1.8)  then
               sig_t =  13.23 + 0.167 /(0.124**2+(sqs-1.92)**2)
            else
               sig_t = 21.8 - 3.81*log(sqs)
            endif
            x_s = sqs - 1.57
            sig_in = 620.* x_s**2.937 * exp(-5.122*x_s)
            x_s = sqs - 1.432
            sig_exc = 66.5 * x_s**0.979 * exp(-3.718*x_s)
            sig_el  = max (.0, sig_t-2.*sig_in-sig_exc)
          endif

!           write(51,222) sqs,sig_el,sig_in,sig_exc,sig_t

c---------------------------------
c         sig_t = .5 *sig_t
c         sig_exc = .5 * sig_exc
c         sig_in = .5 *sig_in
c
c         sig_in = .0
c         sig_el = .0                !     increase - reduction
c---------------------------------
c-hw           end cross section
          radka = sqrt(0.1 * sig_t / pi)
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if (brel .gt. radka)                            goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                        goto 600
*
*   now  the kaon will collide in this time step
          colkpl(i1) = colkpl(i1) + 1
c           density dependence
          ixn = nint(x1)
          iyn = nint(y1)
          izn = nint(z1)
      if(iabs(ixn).le.maxx.and.iabs(iyn).le.maxx.and.iabs(izn).le.maxz)
     &      denst = rhb(ixn,iyn,izn) / rho0
           collkplus(i1) = int(time*2.)
           collkplusd(i1)= int(denst*20.)
************************************************************************
          em2 = rmass     !     deexcitation  of resonances
          pcm = sqrt(0.25*(s-em1**2+em2**2)**2/s -em2**2)
          scatt = 0
*
c-hw           inelastic (1pion production)
*
         sig_rn = rn(iseed)*sig_t
         ch_fin = chk
         if (sig_in .gt. sig_rn) then
            xmpion2 = pmass**2
            xmnucl2 = rmass**2
            xmkaon2 = xkmas**2
            q0    = (sqs-xkmas)**2
            fmax = sqrt((q0-xmpion2-xmnucl2)**2-4.*xmpion2*xmnucl2) /
     1             (xkmas * q0)
            xm2 = (rmass + pmass)**2
            pmax = sqrt((s-xmkaon2-xm2)**2-4.*xm2*xmkaon2)
     1          /(2.*sqs)
  210       continue
            rx = rn(iseed)
            pk = rx**(1./3.) * pmax
            ek = sqrt(xmkaon2+pk*pk)
            q0    = s - 2.*sqs*ek + xmkaon2
            ff = sqrt((q0-xmpion2-xmnucl2)**2-4.*xmpion2*xmnucl2) /
     1                (q0 * ek)
            if (ff .lt. fmax*rn(iseed)) goto 210
           scatt = 1
! 	        if(ika(2,i1).eq.77) then
           ikpinel = ikpinel + 1
!           	write(54,*)  'kaon inelastic scattering',i1,i2,ikpinel
! 		endif

c          write(54,*) ' kplus inelastic scattering ',ikpinel,time
           pcm = pk
            ch_fin = nint(rn(iseed))
         elseif(sig_in+sig_exc .gt. sig_rn)  then
            ch_fin = 1 - chk
         endif
c-hw
  100         xx = 1.0 - 2.*rn(iseed)
              yy = 1.0 - 2.*rn(iseed)
              zz = 1.0 - 2.*rn(iseed)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 1.0 .or. r2 .lt. 0.000001) goto 100
              rr = sqrt(r2)
              betax = (px1+px2) / (e1+e2)
              betay = (py1+py2) / (e1+e2)
              betaz = (pz1+pz2) / (e1+e2)
              gamma  = 1.0 / sqrt(1.0-betax**2-betay**2-betaz**2)
              qx1   = pcm * xx/rr
              qy1   = pcm * yy/rr
              qz1   = pcm * zz/rr
              qx2   = -pcm * xx/rr
              qy2   = -pcm * yy/rr
              qz2   = -pcm * zz/rr
*             lorentz-transformation into lab frame
              e2cm  = sqrt (em2**2 + qx2**2 + qy2**2 + qz2**2)
              p2beta = qx2*betax + qy2*betay + qz2*betaz
              transf = gamma * ( gamma * p2beta / (gamma + 1.0) + e2cm )
              qx2 = betax * transf + qx2
              qy2 = betay * transf + qy2
              qz2 = betaz * transf + qz2
              ntag = 0
              if(ipauli.eq.1 .and.id(1,i2).eq.1)
     &           call pauli(i2,ntag,iseed,phase,r(1,i2),r(2,i2),r(3,i2),
     &                         qx2,qy2,qz2)
*
!               write(54,*) '  pauli in kaoncoll ', i1, i2, ntag
              if(ntag .eq. -1) go to 600

              e1cm  = sqrt (em1**2 + qx1**2 + qy1**2 + qz1**2)
              p1beta = qx1*betax + qy1*betay + qz1*betaz
              transf = gamma * ( gamma * p1beta / (gamma + 1.0) + e1cm )
              pkao(1,i1) = betax * transf + qx1
              pkao(2,i1) = betay * transf + qy1
              pkao(3,i1) = betaz * transf + qz1
              if(scatt.eq.1) then
                if(ika(5,i1).ne.0) then
                nx_kminu(4,ika(5,i1)) = 6666
                ika(5,i1) = 6666   ! zm
                endif
              else

!             if(i1.eq.16) write(54,*)'vor k',i1,nx_kminu(6,i1),
!      &         nx_kminu(4,ika(5,i1)),ika(5,ika(5,i1))

               if(ika(5,i1).ne.0 .and.
     &           ika(5,ika(5,i1)).ne.0 .and.
     &           ika(5,ika(5,i1)).ne.9999 .and.
     &           nx_kminu(6,ika(5,i1)).ne.1 .and.
     &           nx_kminu(6,i1).ne.1 .and.
     &           ika(6,ika(5,i1)).ne.1) then
                nx_kminu(4,ika(5,i1)) = 7777
                ika(5,i1) = 7777   ! zm
                end if
              endif

!              if(i1.eq.16) write(54,*)'nach k',i1,nx_kminu(6,i1),
!      &         nx_kminu(4,ika(5,i1)),ika(5,ika(5,i1))

!              if(ika(2,i1).eq.77) then	
                icollk = icollk + 1
                ikpel = ikpel + 1
!             write(54,*) 'coll',nx_kminu(4,ika(5,i1)), i1
!             endif
!            collkplus(int(time*2.)) = collkplus(int(time*2.)) + 1
c          write(54,*) ' kplus elastic scattering ',ikpel,time

c              ika(6,i1) = 1000 * (ika(6,i1)/1000) + nint(200+z1)         !HS
c             ika(6,i1) = 1000 * (ika(6,i1)/1000) + nthw
              ika(1,i1)  = 2 - ch_fin
c             if(abs(pkao(1,i1)+qx2-px1-px2).gt.1.e-3 .or.
c    &           abs(pkao(2,i1)+qy2-py1-py2).gt.1.e-3 .or.
c    &           abs(pkao(3,i1)+qz2-pz1-pz2).gt.1.e-3)
c    &        write(*,*) 'hiba in kaoscoll, momentum noncoservation',
c    &        pkao(1,i1),qx2,px1,px2,pkao(2,i1),qy2,py1,py2,
c    &        pkao(3,i1),qz2,pz1,pz2
! 	  write(54,*) "absorbed kplus?",ika(1,i1),i1,time

                if(pkao(1,i1)**2+pkao(2,i1)**2+pkao(3,i1)**2.gt.2.5)
     &        write(*,*) 'warning in kaoscoll, momentum is large',
     &        pkao(1,i1),qx2,px1,px2,pkao(2,i1),qy2,py1,py2,
     &        pkao(3,i1),qz2,pz1,pz2

              if(ikaondi.eq.2) then
                p(1,i2) = qx2
                p(2,i2) = qy2
                p(3,i2) = qz2
              end if
c       write(*,*)  '  kaoncoll - charge ', ch2, chk, 2-ika(1,i1)
c             ix=nint(rkao(1,i1))
c             iy=nint(rkao(2,i1))
c             iz=nint(rkao(3,i1))
c             if(iabs(ix).gt.maxx) ix=ix/iabs(ix)*maxx
c             if(iabs(iy).gt.maxx) iy=iy/iabs(iy)*maxx
c             if(iabs(iz).gt.maxz) iz=iz/iabs(iz)*maxz
cc              ithermvar = ithermvar + 1
cc              thermok(2,i1) = ithermvar
cc              if(ithermo.eq.1) call tmunu(ix,iy,iz)
c              if(ithermo.eq.1) write(mterpri,'(i8,f6.2,3i4,14e11.3)')
c     &         ithermvar,time,ix,iy,iz,(avp(jj,ix,iy,iz),jj=1,14)
c-----------------------------------------------------------------------
*
************************************************************************
*                                                                      *
  800   continue
 1000 continue
!       if(kanum.ne.0)
!      & write(54,*) 'kaon elastic coll',icollk
!      & write(54,*) 'kaon elastic coll,# kaons',icollk,kanum
      if(kanum.eq.0 .and. icollk.gt.0)
     & write(*,*) 'hiba in kaoncoll, kanum=0 and icollk>0',icollk,kanum
      return
      end
