************************************************************************
*                                                                      *
      subroutine pionco(lppan,lpiro,lpisi,lprom)
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"com_kminu"
      include"cominput"
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      integer lppan, irun, inp, i1, ii, i2, jj, izzz, ix, iy, iz
      integer lpiro, lpisi, lprom, npat, iendel, idj, id1, id2,iz1,iz2
      real*8  x1, y1, z1, px1, py1, pz1, em1, e1, px2, py2, pz2, em2, e2
      real*8  dx, dy, dz, rsqare, px, py, pz, ee, s, srt, pt, yrap, phi
      real*8  p12, p1dr, p2dr, b12, a12, c12, b21, brel, t1, t2, ddlt
      real*8  xxx, yyy, zzz, densi, qpi2, bwmes, sigrho, sigsig, sigmax
      real*8  rnx, faclebs, pcm, xx, yy, zz, r2, rr, gamma, e1cm, transf
      real*8  rn, p1beta, betax, betay, betaz, path, sigome

      real*8  crx, cry, crz, cpx, cpy, cpz, cpot
      real*8  cfox, cfoy, cfoz, etot
      integer cid2


      lppan= 0
      lpiro= 0
      lpisi= 0
      lprom= 0
      write(*,*)'in pionco'
*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        do 800 ii  = 1,maxp-1
          i1  = ii + inp
          if(ipi(1,i1) .eq. 0)                                 goto 800
          if(ipi(1,i1) .eq.6 .and. i_phi.eq.1) call phi_KK(i1,irun)
          if(ipi(1,i1) .ne. 1 .and. mescol.eq.0)               goto 800
          id1 = ipi(1,i1)
          iz1 = ipi(2,i1)
          x1  = rpi(1,i1)
          y1  = rpi(2,i1)
          z1  = rpi(3,i1)
          px1 = ppi(1,i1)
          py1 = ppi(2,i1)
          pz1 = ppi(3,i1)
          em1 = epi(i1)
          e1  = sqrt( em1**2 + px1**2 + py1**2 + pz1**2 )
          if((id1.eq.2.or.id1.eq.4.or.id1.eq.5).and.ipi(2,i1).ne.0) then
            write(*,*) 'hiba pionco charge',id1,ipi(2,i1)
            stop
          else if((id1.eq.1.or.id1.eq.3).and.iabs(ipi(2,i1)).gt.1) then
            write(*,*) 'hiba pionco charge',id1,ipi(2,i1)
            stop
          else if((id1.eq.6).and.ipi(2,i1)/2 .ne.0) then
            write(*,*) 'hiba pionco charge',id1,ipi(2,i1)
            stop
          end if
*     look for a second pion in the same run
          do 600 jj  = ii+1,maxp
            i2  = jj + inp
            if(ipi(1,i2) .eq. 0)                               goto 600
            if(ipi(1,i2) .ne. 1 .and. mescol.eq.0)             goto 600
            id2 = ipi(1,i2)
            iz2 = ipi(2,i2)
           if((id2.eq.2.or.id2.eq.4.or.id2.eq.5).and.iz2.ne.0)then
              write(*,*) 'hiba pionco charge',id2,ipi(2,i2)
              stop
            else if((id2.eq.1.or.id2.eq.3).and.iabs(iz2).gt.1)then
              write(*,*) 'hiba pionco charge',id2,ipi(2,i2)
              stop
            else if((id2.eq.6).and.iz2/2 .ne.0) then
              write(*,*) 'hiba pionco charge',id2,ipi(2,i2),iz2
              stop
            end if
            dx  = x1 - rpi(1,i2)
            if (abs(dx) .gt. delpi)                            goto 600
            dy  = y1 - rpi(2,i2)
            if (abs(dy) .gt. delpi)                            goto 600
            dz     = z1 - rpi(3,i2)
            if (abs(dz) .gt. delpi)                            goto 600
            rsqare = dx**2 + dy**2 + dz**2
            if (rsqare .gt. delpi**2)                          goto 600
*         now particles are close enough to each other !
            px2 = ppi(1,i2)
            py2 = ppi(2,i2)
            pz2 = ppi(3,i2)
            em2 = epi(i2)
            e2  = sqrt ( em2**2 + px2**2 + py2**2 + pz2**2 )
*   is their impact parameter small enough?
*            write(*,*)'is imp par small enough 1'
            p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
            p1dr   = px1 * dx + py1 * dy + pz1 * dz
            p2dr   = px2 * dx + py2 * dy + pz2 * dz
            a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
            b12    = p1dr / em1 - p2dr * em1 / p12
            c12    = rsqare + ( p1dr / em1 )**2
            brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
*            write(*,*)'small enough 2'
            if (brel .gt. dispi)                               goto 600
            b21    = - p2dr / em2 + p1dr * em2 / p12
            t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
            t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
            ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
*            write(*,*)'closest approach '
            if ( abs(t1+t2) .gt. dt )                          goto 600
*   now  the pions may annihilate
            px     = px1 + px2
            py     = py1 + py2
            pz     = pz1 + pz2
            ee     = e1  + e2
            s      = ee**2 - px**2 - py**2 - pz**2
            izzz   = ipi(2,i1)+ipi(2,i2)
            srt    = sqrt(s)
            pt     = sqrt(px**2 + py**2)
            yrap   = 0.5 * log( (ee+pz)/(ee-pz) )
            phi    = atan2(py,px)
            xxx    = 0.5 * (x1+rpi(1,i2))
            yyy    = 0.5 * (y1+rpi(2,i2))
            zzz    = 0.5 * (z1+rpi(3,i2))
            ix = nint(xxx)
            iy = nint(yyy)
            iz = nint(zzz)
            densi = 0.0
            if(iabs(ix).le.maxx.and.iabs(iy).le.maxx.and.
     &         iabs(iz).le. maxz) densi = rhb(ix,iy,iz)
*
*************************************************************************
*        pion annihilation
            if(id1.eq.1 .and. id2.eq.1) then
              lppan  = lppan + 1
c              write(*,*)'vor imeson'
              if(imeson.eq.1)
     &          call pipiann(srt,pt,yrap,phi,densi)
            endif
c            write(*,*)'vor ipico'
************************************************************************
*           simple isotropic collision
            if(ipico.eq.1.and. idec2pi.eq.0) then
c              write(*,*) 'pionco, ipico1',ipico, idec2pi
              pcm = sqrt(0.25*(s-em1**2+em2**2)**2/s -em2**2)
  100         xx = 1.0 - 2.*rn(iseed)
              yy = 1.0 - 2.*rn(iseed)
              zz = 1.0 - 2.*rn(iseed)
              r2 = xx**2 + yy**2 + zz**2
              if(r2 .gt. 1.0 .or. r2 .lt. 0.000001) go to 100
              rr = sqrt(r2)

              px1   = pcm * xx/rr
              py1   = pcm * yy/rr
              pz1   = pcm * zz/rr
              betax = px / ee
              betay = py / ee
              betaz = pz / ee
              gamma  = 1.0 / sqrt(1.0-betax**2-betay**2-betaz**2)
*             lorentz-transformation into lab frame
              e1cm  = sqrt (em1**2 + px1**2 + py1**2 + pz1**2)
              p1beta  = px1*betax + py1*betay + pz1*betaz
              transf  = gamma * ( gamma * p1beta / (gamma + 1) + e1cm )
              ppi(1,i1) = betax * transf + px1
              ppi(2,i1) = betay * transf + py1
              ppi(3,i1) = betaz * transf + pz1
*
              ppi(1,i2) = px - ppi(1,i1)
              ppi(2,i2) = py - ppi(2,i1)
              ppi(3,i2) = pz - ppi(3,i1)
            end if
************************************************************************
*           phi production
            sigmax = 10.*pi*dispi**2
            if(i_phi.eq.1) then
c-------  pi+pi -> pi+phi
              if(id1.eq.1 .and. id2.eq.1) then
                call phi_pi_pipi(ee,px,py,pz, srt, xxx,yyy,zzz, sigmax,
     &           i1,i2, iz1,iz2, irun)
c                write(*,*) 'pionco call phi_pi_pipi',i_phi
              end if
C c-------  pi+rho -> phi
              if((id1.eq.1.and.id2.eq.3).or.(id2.eq.1.and.id1.eq.3))then
                call phi_pirho(ee,px,py,pz, srt, xxx,yyy,zzz, sigmax,
     &           i1,i2, iz1,iz2, irun)
c                write(*,*) 'pionco call phi_pirho',i_phi
              end if
            end if
************************************************************************
*            rho and sigma creation
*         rho   cross section at peak : 113.5 mb bmax:1.9 fm
*         sigma cross section at peak : 34.67 mb bmax:1.05 fm; m=0.8 GeV
            if((ipico.eq.2.or.idec2pi.gt.0).and.id1.eq.1.and.id2.eq.1)
     &         then
*              write(*,*) 'pionco, ipico2',ipico, idec2pi
              if(izzz.gt.1 .or. izzz.lt.-1)                     goto 400
              qpi2 = 0.25 * s - pmass**2
              sigrho=120.*pi/qpi2*hbc**2*
     &               bwmes(s,1,idec2pi,iresmode,1,0,iwidth,1)
              if(ipi(2,i1).eq.0 .and. ipi(2,i2).eq.0) sigrho=0.0
              if(izzz.ne.0) then
                sigsig= 0.0
              else
                faclebs = 0.333333
                if(ipi(2,i1).ne.0) faclebs = 0.6666667
                sigsig=40.*pi/qpi2*faclebs*hbc**2*
     &                    bwmes(s,2,idec2pi,iresmode,1,0,iwidth,1)
              endif
              rnx = rn(iseed)
              sigmax= 10.*pi*dispi**2
              if(rnx.gt.(sigrho+sigsig)/sigmax)                 goto 400
c-----------------------------------------------------------------------
*        sigma or rho is now created
c-----------------------------------------------------------------------
c           meson collision number

              idj = min(ipi(4,i1) , 50)
              mlife(1,idj) = mlife(1,idj) + 1
              idj = min(ipi(4,i2) , 50)
              mlife(1,idj) = mlife(1,idj) + 1
              pitomes = pitomes + 1
c-----------------------------------------------------------------------
              ppi(1,i1) = px
              ppi(2,i1) = py
              ppi(3,i1) = pz
              epi(i1)   = srt
              rpi(1,i1) = xxx
              rpi(2,i1) = yyy
              rpi(3,i1) = zzz
              if(rnx.lt.sigrho/sigmax) then
                ipi(1,i1) = 3
                lpiro     = lpiro + 1
*               write(*,*) 'pionco rho is created', i1, rpi(1,i1)
************************************************************************
*       meson dileptonic decay
                if(ivmesdil.eq.1 .and. izzz.eq.0) then
                  write(*,*) 'pionco call vmesdil'
                  gamma=ee/srt
c                  call vmesdil(srt,pt,yrap,phi,densi,gamma,50.,3,-9)
                end if
************************************************************************
              else
                ipi(1,i1) = 4
                lpisi     = lpisi + 1
              endif
              ipi(2,i1) = izzz
              ipi(4,i1) = 1
              ipi(5,i1) = -1
              ipi(1,i2) = 0
              if((ipi(1,i1).eq.2.or.ipi(1,i1).eq.4).and.ipi(2,i1).ne.0)
     &           then
                write(*,*)'hiba:pionco veg',ipi(1,i1), ipi(2,i1),sigsig
     &                ,sigrho/sigmax,rnx
                stop
              end if
*-----------------------------------------------------------------------
*                 path length distribution
              path=sqrt((rpie(1,i1)-x1)**2+(rpie(2,i1)-y1)**2+
     &                      (rpie(3,i1)-z1)**2)
              npat=min(30,nint(2.0*path))
              mpath(1,npat)=mpath(1,npat) + 1
              path=sqrt((rpie(1,i2)-rpi(1,i2))**2+
     &             (rpie(2,i2)-rpi(2,i2))**2+(rpie(3,i2)-rpi(3,i2))**2)
              npat=min(30,nint(2.0*path))
              mpath(1,npat)=mpath(1,npat) + 1
*-----------------------------------------------------------------------
*               lifetime distribution
              path=max( 0.0, time - rpie(6,i1))
              npat=min(50,nint(2.0*path))
              melet(1,npat)=melet(1,npat) + 1
              path=max( 0.0, time - rpie(6,i2))
              npat=min(50,nint(2.0*path))
              melet(1,npat)=melet(1,npat) + 1
*-----------------------------------------------------------------------
              iendel = nint(5.0*densi/rho0)
              if(iendel .gt. 50) iendel = 50
              mdens(ipi(1,i1),iendel) = mdens(ipi(1,i1),iendel) + 1
c-----------------------------------------------------------------------
              rpie(1,i1)= xxx
              rpie(2,i1)= yyy
              rpie(3,i1)= zzz
              rpie(4,i1)= densi
c-hw             rpie(5,i1)= 2.0 * pmass
              rpie(6,i1)= time
************** store coulomb pot
                crx = rpi(1,i1)
                cry = rpi(2,i1)
                crz = rpi(3,i1)

                cpx = ppi(1,i1)
                cpy = ppi(2,i1)
                cpz = ppi(3,i1)

                cid2 = ipi(2,i1)

                etot = sqrt(pmass**2+cpx**2+cpy**2+cpz**2)
                call emfoca(crx,cry,crz,cid2,
     &                       cfox,cfoy,cfoz,ncont,cpot)

c-hw               rpie(7,i1) = cpot


*********************************************
            end if
  400       continue
c-----------------------------------------------------------------------
*         omega cross section at peak : mb bmax: fm;
            if(mescol.eq.1.and.(id1*id2.eq.3).and.(id1+id2.eq.4)) then
c              write(*,*) 'pionco, pi+rho->omega',ipico, idec2pi
              if(izzz.ne.0)                                     goto 500
              qpi2 = 0.25*(s-em1**2+em2**2)**2/s -em2**2
*        gamin = gamtot/3. since 3 decay channels(+-,00,-+) are possible
              sigome=40./3.*pi/qpi2*hbc**2*
     &               bwmes(s,3,idec2pi,iresmode,2,0,iwidth,1)
              write(*,*) 'pionco,omega',sigome,sqrt(s),qpi2
              rnx = rn(iseed)
              sigmax= 10.*pi*dispi**2
              if(rnx.gt.sigome/sigmax)                          goto 500
c-----------------------------------------------------------------------
*        omega is now created
c-----------------------------------------------------------------------
c           meson collision number

              idj = min(ipi(4,i1) , 50)
              mlife(id1,idj) = mlife(id1,idj) + 1
              idj = min(ipi(4,i2) , 50)
              mlife(id2,idj) = mlife(id2,idj) + 1
              pitomes = pitomes + 1
c-----------------------------------------------------------------------
              ppi(1,i1) = px
              ppi(2,i1) = py
              ppi(3,i1) = pz
              epi(i1)   = srt
              rpi(1,i1) = xxx
              rpi(2,i1) = yyy
              rpi(3,i1) = zzz
              ipi(1,i1) = 5
              lprom     = lprom + 1
              ipi(2,i1) = izzz
              ipi(4,i1) = 1
              ipi(5,i1) = -1
              ipi(1,i2) = 0
c             write(*,*) 'pionco - omega created: ',i1,px,py,pz,
c    &           xxx,yyy,zzz,srt
              if((ipi(1,i1).eq.2.or.ipi(1,i1).ge.4).and.ipi(2,i1).ne.0)
     &           write(*,*) 'hiba:pionco',ipi(1,i1), ipi(2,i1), sigsig
              if((ipi(1,i2).eq.2.or.ipi(1,i2).ge.4).and.ipi(2,i2).ne.0)
     &           write(*,*) 'hiba:pionco',ipi(1,i2), ipi(2,i2), sigsig
*-----------------------------------------------------------------------
*                 path length distribution
              path=sqrt((rpie(1,i1)-x1)**2+(rpie(2,i1)-y1)**2+
     &                      (rpie(3,i1)-z1)**2)
              npat=min(30,nint(2.0*path))
              mpath(id1,npat)=mpath(id1,npat) + 1
              path=sqrt((rpie(1,i2)-rpi(1,i2))**2+
     &             (rpie(2,i2)-rpi(2,i2))**2+(rpie(3,i2)-rpi(3,i2))**2)
              npat=min(30,nint(2.0*path))
              mpath(id2,npat)=mpath(id2,npat) + 1
*-----------------------------------------------------------------------
*               lifetime distribution
              path=max( 0.0, time - rpie(6,i1))
              npat=min(50,nint(2.0*path))
              melet(id1,npat)=melet(id1,npat) + 1
              path=max( 0.0, time - rpie(6,i2))
              npat=min(50,nint(2.0*path))
              melet(id2,npat)=melet(id2,npat) + 1
*-----------------------------------------------------------------------
              iendel = nint(5.0*densi/rho0)
              if(iendel .gt. 50) iendel = 50
              mdens(ipi(1,i1),iendel) = mdens(ipi(1,i1),iendel) + 1
c-----------------------------------------------------------------------
              rpie(1,i1)= xxx
              rpie(2,i1)= yyy
              rpie(3,i1)= zzz
              rpie(4,i1)= densi
c-hw             rpie(5,i1)= 2.0 * pmass
              rpie(6,i1)= time
            end if
  500     continue
*
  600     continue
  800   continue
 1000 continue
*      write(*,*)'end of pionco'
      return
      end
