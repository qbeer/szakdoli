************************************************************************
*                                                                      *
      subroutine pionab(lmesa)
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
c      include"resdata"
      real*8 pot_hw
      real*8 pirp, pire, pirr, crospi, croset, crosro, pirkaon
      real*8 em12, em22, srt,  srtp, s, srt_p,gyoks
      real*8 srt_0, vx, vy, vz, pot1, pot2
      real*8 x1, y1, z1, px1, py1, pz1, em1, e1, e10, e20
      real*8 dx, dy, dz, rsqare, px2, py2, pz2, em2, e2, sig0, qq2, ecm
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt,pirb,b21
      real*8 ede, pdx, pdy, pdz, xx, yy,   rn, rdx, rdy, rdz
      real*8 bwdist,   path
      real*8 gamma, p1beta
c      real*8 sigkaon
      integer dimlmesc,id1,iz1,iz2,ibar
      integer irun, inp, ini,ii,i1,i2, idres, id6, id62,id2
      parameter     (pirp   = 2.52,  pire  = 1.38, pirkaon= 1.2)
      parameter(pirr = 1.3)
      parameter(dimlmesc = nres+9)
      integer lmesa(dimlmesc), npat, nresmo
*----------------------------------------------------------------------*
*    pirr = 1.3  fm                 corresponds  to 53 mb              *
*    pir0 = 2.07 fm                 corresponds  to 135 mb             *
*    pirp = 2.52 fm and pir2 = 1.49fm correspond to 200 mb and 70 mb   *
*    pire = 1.38 fm                 corresponds  to 60 mb              *
*----------------------------------------------------------------------*
      real*8 beta(3)

*----------------------------------------------------------------------*
*     variables needed for mom-dep stuff                               *
      real*8    betlrfx,betlrfy,betlrfz
      real*8    j02, pabs
      real*8    vecpot, scapot(1:2), potanal, finmass(1:2)
      real*8    meff1, meff2 , minmass, test
      integer  itcount
      logical channel(1:2), decision, stopflag
c     real*8    t0, tin, tout
      real*8    consen, phelp1, phelp2
      real*8    deltait, epsit, itmass, phelp
      real*8    pdxb, pdyb, pdzb, newmass, oldmass , pabsold
      real*8    help
      real*8    dummyf(1:9), maxcross, testsig, pres
      real*8     sigt, isofac, sig(1:nres)
      integer ichannel, multipl, i, qmes, qtot, qn
      real*8    srts
      integer totabs, nloopind
      real*8    plrf(1:3), rhap(1:3)

c      include"resdata1"

c     t0 = 0.0
c      tin = secnds(t0)
      totabs = 0
*----------------------------------------------------------------------*
      crospi = 10.0 * pi * pirp**2
      croset = 10.0 * pi * pire**2
      crosro = 10.0 * pi * pirr**2

      write(*,*)  ' pionab 1b ', dimlmesc

      call dens_4
      call potcalc
      do 10 ii=1,dimlmesc
        lmesa(ii) = 0
  10  continue

*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        ini = (irun-1) * maxb
        do 800 ii  = 1,maxp
          i1  = ii + inp
c         write(*,*) ' pionab 2 ', irun, ii, i1, ipi(1,i1), iseed
          if(ipi(1,i1) .eq. 0 )                                goto 800
c
          if(ipi(1,i1).eq.4 .and.ipi(2,i1).ne.0 ) then
             write(*,*) 'sigma charge is not zero',ipi(1,i1),ipi(2,i1)
             stop
          end if

          x1         = rpi(1,i1)
          y1         = rpi(2,i1)
          z1         = rpi(3,i1)
          px1        = ppi(1,i1)
          py1        = ppi(2,i1)
          pz1        = ppi(3,i1)
          em1        = epi(i1)
          meff1      = em1 +  mpot(i1)
          phelp1     = sqrt(px1**2+py1**2+pz1**2)
          em12       = epi(i1)**2
          e1         = sqrt(meff1**2+px1**2+py1**2+pz1**2)
          e10        = sqrt(em1**2+px1**2+py1**2+pz1**2)
          id1        = ipi(1,i1)
          iz1        = ipi(2,i1)
c          if(phelp1.gt.1.d20) then
c            write(*,*)'pionab p1',px1,py1,pz1,id1,em1,meff1
c            stop
c          end if
c          if(ipi(1,i1).eq.6) write(*,*)'pionab, ipi',(ipi(jj,i1),jj=1,7)
*     look for an absorbent pseudonucleon in the same run
          i2  = ini
c          write(*,*)'t1 ', irun, ii
  600     i2  = i2 + 1
c          write(*,*) 'pionab0',i2, id(1,i2), maxb, ini, maxp, num
          if(i2 .gt. ini+maxb)                                 goto 800
          if(id(1,i2).eq.0)                                    goto 600
c          write(*,*) 'pionab1',maxb, ini, maxp, num
          if(id(1,i2).eq.1) then
            if(id(2,i2)/2.ne.0) then
              write(*,*) 'hiba in pionab 0',id(1,i2),id(2,i2)
              stop
            end if
          else if(id(1,i2).le.nres+1 .and. id(1,i2).gt.1) then
            if(resprop2(abs(id(1,i2))-1,1).eq.1) then
              if(id(2,i2).gt.1 .or. id(2,i2).lt.0) then
                write(*,*)'charge violation in pionab 0'
                write(*,*)i2, id(1,i2), id(2,i2)
                stop
              end if
            else if(resprop2(abs(id(1,i2))-1,1).eq.3) then
              if(id(2,i2).gt.2 .or. id(2,i2).lt.-1) then
                write(*,*)'charge violation in pionab 0a'
                write(*,*)i2, id(1,i2), id(2,i2)
                stop
              end if
            end if
          end if
c          write(*,*)'t19 ', id(2,816)
          if(i2 .eq. ipi(3,i1) .and. id(8,i2).eq.i1)           goto 600
c          if(i2 .eq. ipi(3,i1) .and. id(8,i2).eq.ipi(6,i1).and.
c     &       id(8,i2).ne.0 )                                   goto 600
c          if(ipi(7,i1).eq.i2.and.id(3,i2).eq.ipi(3,i1).and.
c     &       id(3,i2) .ne.0)                                   goto 600

c--- selection of baryons no meson absorption in the model for id1.gt.4:
c          if(id(1,i2).gt.4 .and. id(1,i2).lt.nres+2)           goto 600
          if(id(1,i2).ge.2.and.id(1,i2).lt.nres+2.and.irescol.lt.2)
     &                                                         goto 600

c--- strange collide only with strange:
C          if(ipi(1,i1).eq.6 .and. (id(1,i2).ne.nres+2.and.
C     &                             id(1,i2).ne.nres+3))        goto 600
c          if(id(1,i2).eq.nres+2 .and. ipi(1,i1).ne.6)          goto 600
c          if(id(1,i2).eq.nres+3 .and. ipi(1,i1).ne.6)          goto 600

c          if(ipi(1,i1).eq.6) write(*,*) 'pionab, id',(id(jj,i1),jj=1,8),
c     &          (id(jj,i2),jj=1,7)
          dx     = x1 - r(1,i2)
          if (abs(dx) .gt. delpi)                              goto 600
          dy     = y1 - r(2,i2)
          if (abs(dy) .gt. delpi)                              goto 600
          dz     = z1 - r(3,i2)
          if (abs(dz) .gt. delpi)                              goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. delpi**2)                            goto 600
*         now particles are close enough to each other !
c          if(ipi(1,i1).eq.3) write(*,*) 'zm: rho+bar. are close enough'

          rhap(1)    = r(1,i2)
          rhap(2)    = r(2,i2)
          rhap(3)    = r(3,i2)

          id2        = id(1,i2)
          iz2        = id(2,i2)
          px2        = p(1,i2)
          py2        = p(2,i2)
          pz2        = p(3,i2)
          em2        = e(i2)
          meff2      = em2 + upot(i2)
          phelp2     = sqrt(px2**2+py2**2+pz2**2)
          em22       = e(i2)**2
          e2         = sqrt (meff2**2+px2**2+py2**2+pz2**2 )
          e20        = sqrt (em2**2+px2**2+py2**2+pz2**2 )
          betlrfx    = betlrfboo(i2,1)
          betlrfy    = betlrfboo(i2,2)
          betlrfz    = betlrfboo(i2,3)
          j02        = betlrfboo(i2,4)
          s          = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                           - (pz1+pz2)**2
          srt        = sqrt(s)
          srt_p  = sqrt( (e10+e20)**2 - (px1+px2)**2 - (py1+py2)**2
     1                 - (pz1+pz2)**2) + upot(i2) + mpot(i1)
          if(phelp2.gt.1.d20) then
            write(*,*)'pionab p2',px2,py2,pz2,id2,em2,meff2
            stop
          end if
c          write(*,*)'t2 ',  irun, ii, i2

*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( meff1 * meff2 / p12 ) ** 2
          b12    = p1dr / meff1 - p2dr * meff1 / p12
          c12    = rsqare + ( p1dr / meff1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if(ipi(1,i1) .eq. 1 )  pirb = pirp
          if(ipi(1,i1) .ge. 2 )  pirb = pire
          if(ipi(1,i1) .eq. 3 )  pirb = pirr
          if (brel .gt. pirb)                                  goto 600
c          if(ipi(1,i1).eq.3) write(*,*) 'zm: impact par. is small en.'
          maxcross = 10.*pi*pirb**2
          b21    = - p2dr / meff2 + p1dr * meff2 / p12
          t1     = ( p1dr / meff1 - b12 / a12 ) * e1 / meff1
          t2     = ( - p2dr / meff2 - b21 / a12 ) * e2 / meff2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
c          if(ipi(1,i1).eq.3) write(*,*) 'zm: closest point now'
*   now  the pion may be absorbed in this time step
*
          pdx  = px1 + px2
          pdy  = py1 + py2
          pdz  = pz1 + pz2
          pabs = (pdx**2 + pdy**2 +pdz**2)
          ede = e1 + e2

          if(pabs.gt.1.d20) then
            write(*,*)'pionab pd',pdx,pdy,pdz,ede,id1,id2
            stop
          end if
c          if(i_phi.eq.1 .and. (id1.eq.5.or.id1.eq.3) .and. iz2/2.eq.0) 
c     &       call mesphi(irun,i1,i2,id1,id2,iz2,srt,pdx,pdy,pdz,ede,
c     &         rhap(1),rhap(2),rhap(3),maxcross)

c          write(*,*)'t3 ', id(2,816)
*---------------------------------------------------------------------*
*      in the new version s contains the full s (including pots)      *
*      one has to extract the free s by subtracting the potentials    *
*      this is easier than in the other parts of the collision term   *
*      because the momentum of the resonance in the final channel is  *
*      equal to the sum of the baryon and the meson in the initial    *
*      channel.                                                       *
*      so in order to obtain the free srt one has to take the full    *
*      srt = m + uscalar(resonance in final channel                   *
*       ==> srt(free) = srt - uscalar                                 *
*      the resonance is created at the place of the baryon in the     *
*      initial channel

*     abs value of p of the resonance                                 *
*     sqrt(s(incl. pots)) corresponds to the eff. mass of the res.    *
*---------------------------------------------------------------------*
          consen  = ede
          pabsold = pabs
          gyoks   = sqrt(consen**2 - pabsold)
          pabsold = sqrt(pabsold)
          epsit = 1.0d-06
*
*         nresmo = 1 : nucleon
*         nresmo = 2 : delta, no potential yet

          do nresmo = 1, 2
            oldmass = dmass
            deltait  = 100000.0
            itcount = 0
            do while(abs(deltait).gt.epsit)
              itcount = itcount + 1
              pdxb    = pdx
              pdyb    = pdy
              pdzb    = pdz
              itmass  = oldmass
              ede  = sqrt(pdx**2 + pdy**2 + pdz**2 + itmass**2)
c              write(*,*)'pionab itmass',itmass,itcount,pabsold,j02,ede
*
*         boost back to the LRF of the baryon in the initial state    *
             call lorentz(-betlrfx,-betlrfy,-betlrfz,pdxb,pdyb,pdzb,ede)
              pabs   = sqrt(pdxb**2+pdyb**2+pdzb**2)
              plrf(1) = pdxb
              plrf(2) = pdyb
              plrf(3) = pdzb
              if(pabs.gt.1.d20) then
c                write(*,*)'pionab lrf',betlrfx,betlrfy,betlrfz,
c     &    pdxb,pdyb,pdzb,ede,pabsold,itmass,scapot(nresmo),vecpot
                stop
              end if
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*
              vecpot        = potanal(rho0,j02,plrf,nresmo,rhap)
              scapot(nresmo)  =  -itmass + sqrt(itmass**2 +
     +           2.0*sqrt(pabs**2+itmass**2)*vecpot + vecpot**2)
c             write(*,*)'pionab scapot',vecpot,itmass,scapot(nresmo),pabs
c              newmass = consen**2 - pabsold**2
              newmass = gyoks - scapot(nresmo)
              newmass = max(newmass,0.001)
              deltait = newmass - oldmass
              if(abs(deltait).le.epsit)then
                if(newmass.le.0.001) 
     &                 write(*,*)'pionab vigyazz newmass small',newmass
                finmass(nresmo) = newmass
              end if
              oldmass = newmass
              if(itcount.gt.1000) then
                finmass(nresmo) = 0.0001
                deltait = 0.01*epsit
              end if
            end do
          end do

*     das finmass(1), finmass(2) stellt auch gleichzeitig das wurzel *
*     s im freien stoss dar.                                         *
*     und nur dieses frei s wird ab jetzt noch gebraucht.            *

c          write(*,*)'t4 ', id(2,816)
          srtp    = max(finmass(1),finmass(2))
          minmass = rmass + pmass
          if (srtp .le. minmass  )                          goto 600
c          if(ipi(1,i1).eq.3) write(*,*) 'zm: srt is above threshold'
          channel(1) = .false.
          channel(2) = .false.
          if(finmass(1).gt.minmass) channel(1) = .true.
          if(finmass(2).gt.minmass) channel(2) = .true.

          if(ipi(3,i1).gt.0) id62= id(6,ipi(3,i1))
          if(ipi(3,i1).le.0) id62= 1

c            if(time.eq.11.5)   write(*,*)'stop 1'
c          write(*,*)'t4 ', id(2,816)
************************************************************************
*               pion nucleon -> kaon + lambda                          *
c          write(*,*) 'in pionab, before partprod'
          if((ikaon.eq.1 .or. i_phi.eq.1 . or. imeson.eq.1) .and.
     &          ((ipi(1,i1).eq.1 .or. ipi(1,i1).eq.2))) then
            sig0 = crospi
            rdx = r(1,i2)
            rdy = r(2,i2)
            rdz = r(3,i2)
c            ede = sqrt(srtp**2 + pabsold**2) ! energy of the resonance
            ede = e1 + e2
            id6 = id(6,i2)
            id2 = id(1,i2)
            vx = px1 / e10
            vy = py1 / e10
            vz = pz1 / e10
            if((i_phi.eq.1 .or. ikaon.eq.1)  .and.
     1                      ipi(1,i1).eq.1        )  then
              ibar = 0
              pot1=0.0
              pot2=0.0
c              write(*,*) "pionab potcalci1", id1,i1, ibar, vx, vy, vz
              if(ipipot.eq.1) call potcalc_i(i1, ibar, vx, vy, vz, pot1)
              ibar = 1
c              write(*,*) "pionab potcalci2", id2,i2,ibar,vx, vy, vz,pot1
c     &   ,newmass,pdx,pdy,pdz,ede
              call potcalc_i(i2, ibar, vx, vy, vz, pot2)
              vx = px2/e2
              vy = py2/e2
              vz = pz2/e2
c              write(*,*) "pionab potcalci3", id2,i2,ibar,vx, vy, vz,pot2
              call potcalc_i(i2, ibar, vx, vy, vz, pot_hw)
c              write(*,*) "pionab potcalci4", pot1,pot2,pot_hw,srt_p
c              pot_hw = 0.0
c              if(pot2.gt.-0.1 .and. pot2.lt.0.5) pot_hw = pot2
              srt_0 = srt_p - pot1 - pot_hw
c         write(*,*)'before phi/ka_dpi ',mpot(i1),i2,upot(i2),pot2,pot_hw
c              if(ikaon.eq.1  .and. brel .le. pirkaon
c     1          .and.   iphi_dec .ne. 1)            then
c                sigkaon = 10.*pirkaon**2 *pi
c                srt_0 = srt_p  - 0.666*pot2
c                call kaondpi(ede,pdx,pdy,pdz,srt_0,rdx,rdy,rdz,sigkaon,
c     &            i2,id2,id6,id62,ipi(2,i1),id(2,i2),i1,irun)
c                write(*,*) 'in pionab, after kaondpi', iseed
c              endif
c              write(*,*)'pionab, srt ',srt_0,srt_p,srt,srtp,pot1,pot_hw
              if(i_phi .eq. 1 .and. id1.eq.1) then
                 call phi_dpi(ede,pdx,pdy,pdz,srt,rdx,rdy,rdz,sig0,
     &                i2,id2,id6,ipi(2,i1),id(2,i2),ipi(1,i1),irun)
c     write(*,*) 'in pionab, after phi_dpi', iseed
              endif
              if(i_phi .eq. 1 .and. (id1.eq.2.or.id1.eq.5)) then
                 call phi_omegaeta(ede,pdx,pdy,pdz,srt,meff1,meff2,
     &                rdx,rdy,rdz,sig0,i2,id2,id6,iz1,iz2,id1,irun)
c     write(*,*) 'in pionab, after phi_omegaeta', iseed
              endif
              if(i_phi .eq. 1 .and. id1.eq.3) then
c      write(*,*) ' now calling phi_rho'
                 call phi_rho(ede,pdx,pdy,pdz,srt,meff1,meff2,
     &                rdx,rdy,rdz,sig0,i2,id2,id6,iz1,iz2,id1,irun)
              end if
              if(i_kminu .eq. 1) then
                 call kminu_dpi(ede,pdx,pdy,pdz, srt,rdx,rdy,rdz,sig0,
     &                i1,ipi(1,i1), i2, id2,iz1,iz2,irun)
c      write(*,*) 'in pionab, after kminu_dpi'
              endif
              if(i_kminu.eq.1 .and.  ipi(1,i1).eq.1 .and. 
     &             (id2.eq.nres+2 .or. id2.eq.nres+3)) then
c      write(*,*) 'in pionab, before kminu_hpi'
      call kminu_hpi(ede,pdx,pdy,pdz, srt,rdx,rdy,rdz,sig0,
     &            i1,id1, i2, id2-nres-1,iz1,iz2,irun)
      write(*,*) 'in pionab, after kminu_hpi'
                 goto 600
              endif
              if(i_JPsi .ge. 1 .and. id1.eq.1) then
                 call JPsi_pert_Npi(i1,i2,id1,id2,iz1,
     $                iz2,ede,pdx,pdy,pdz,srt,rdx,rdy,rdz,
     $                sig0,irun)
c     write(*,*) 'in pionab, after JPsi_pert_Npi', iseed
              endif              
           endif
           if(abs(id(1,i2)).eq.1) then
              if(id(2,i2)/2.ne.0) then
                 write(*,*) 'hiba in pionab 0',id(1,i2),id(2,i2)
                 stop
              end if
           else if(id(1,i2).le.nres+1) then
              if(resprop2(abs(id(1,i2))-1,1).eq.1) then
                 if(id(2,i2).gt.1 .or. id(2,i2).lt.0) then
                    write(*,*)'charge violation in pionab 0b'
                    write(*,*)i2, id(1,i2), id(2,i2)
                    stop
                 end if
              else if(resprop2(abs(id(1,i2))-1,1).eq.3) then
                 if(id(2,i2).gt.2 .or. id(2,i2).lt.-1) then
                    write(*,*)'charge violation in pionab 0c'
                    write(*,*)i2, id(1,i2), id(2,i2)
                    stop
                 end if
              end if
           end if

c            if(imeson.eq.1)
c     &            call vecmespi(ede,pdx,pdy,pdz,srt,rdx,rdy,rdz,sig0
c     &             ,id2,ipi(2,i1),id(2,i2))
          endif
c          write(*,*) 'in pionab, after vecmespi'
          call f77flush()
            if(abs(id(1,i2)).eq.1) then
              if(id(2,i2)/2.ne.0) then
                write(*,*) 'hiba in pionab 0',id(1,i2),id(2,i2)
                stop
              end if
            else if(id(1,i2).le.nres+1 .and. id(1,i2).gt.1) then
              if(resprop2(abs(id(1,i2))-1,1).eq.1) then
                if(id(2,i2).gt.1 .or. id(2,i2).lt.0) then
                  write(*,*)'charge violation in pionab 0d'
                  write(*,*)i2, id(1,i2), id(2,i2)
                  stop
                end if
              else if(resprop2(abs(id(1,i2))-1,1).eq.3) then
                if(id(2,i2).gt.2 .or. id(2,i2).lt.-1) then
                  write(*,*)'charge violation in pionab 0e'
                  write(*,*)i2, id(1,i2), id(2,i2)
                  stop
                end if
              end if
            end if

************************************************************************
*               pion baryon -> omega + nucleon                         *
C           if(iomega.eq.1 .and. ipi(1,i1).eq.1) then
C             sig0 = crospi
C             ede = e1 + e2
C             pdx = px1+px2
C             pdy = py1+py2
C             pdz = pz1+pz2
C             rdx = r(1,i2)
C             rdy = r(2,i2)
C             rdz = r(3,i2)
C             id2 = id(1,i2)
C             id6 = id(6,i2)
C             call omegapi(ede,pdx,pdy,pdz,srt,rdx,rdy,rdz,sig0,
C      &            i2,id2,id6,id62,ipi(2,i1),id(2,i2),ipi(1,i1))
C           endif
c          write(*,*) 'in pionab, after omegapi'
************************************************************************
*               rho nucleon -> phi nucleon      zm
           if( (ikaon.eq.1 .or. imeson.eq.1) .and. i_phi.eq.1 .and.
     &          ipi(1,i1).eq.3 .and. id(1,i2).le.2 ) then
              pot1 = 0.0
              pot2 = 0.0
              sig0 = crosro
              ede = e1 + e2
              pdx = px1+px2
              pdy = py1+py2
              pdz = pz1+pz2
              rdx = r(1,i2)
              rdy = r(2,i2)
              rdz = r(3,i2)
              id2 = id(1,i2)
              id6 = id(6,i2)
              vx = pdx / ede
              vy = pdy / ede
              vz = pdz / ede
c              write(*,*) ' in pionab rho-phi ', i1,i2,vz
C             ibar = 0
C           call potcalc_i(i1, ibar, vx, vy, vz, pot1)
C           ibar = 1
C           call potcalc_i(i2, ibar, vx, vy, vz, pot2)
           srt_0 = srt - pot1 - pot2
c            write(*,*) ' before phi_rho ', mpot(i1),upot(i2) , pot1,pot2
           call phi_rho(ede,pdx,pdy,pdz,srt_0, meff1,meff2, rdx,rdy,
     1      rdz,sig0,i2,id2,id6,iz1,iz2,id1,irun)
           endif
c          write(*,*) 'in pionab, after phi_rho'
************************************************************************

*   now  the pion is absorbed by the nucleon.
************************************************************************
*               pion nucleon bremsstrahlung                            *
          if(ipiNbrems.eq.1 .and. ipi(1,i1).eq.1 .and. id(1,i2).eq.1
     &      .and. ipi(2,i1)+id(2,i2).ge.0 .and. ipi(2,i1)+id(2,i2).lt.2) 
     &              then
            sig0 = crospi
            ede = e1 + e2
            pdx = px1+px2
            pdy = py1+py2
            pdz = pz1+pz2
            rdx = 0.5*(r(1,i2)+rpi(1,i1))
            rdy = 0.5*(r(2,i2)+rpi(2,i1))
            rdz = 0.5*(r(3,i2)+rpi(3,i1))
            beta(1) = pdx / ede
            beta(2) = pdy / ede
            beta(3) = pdz / ede
            gamma  = 1.0 / sqrt(1.0-beta(1)**2-beta(2)**2-beta(3)**2)
*         transformation of pion energy
            p1beta = px1*beta(1) + py1*beta(2) + pz1*beta(3)
            ecm    = gamma * ( e1 - p1beta)
c           write(*,*) 'bremslep in pionab!'
            call piNdilepCreate(ipi(2,i1),id(2,i2),irun,srt,sig0,
     &       beta,gamma,rdx,rdy,rdz)
          end if
c          write(*,*) 'in pionab, after bremslep'
************************************************************************
***            new new new (st 13.9.95)                            *****
***********************************************************************

          if(id(1,i2).gt.3 .and. id(1,i2).lt.nres+2)           goto 600
c          write(*,*)'t5 ', id(2,816)
          if((ipi(1,i1).eq.1).and.(id(1,i2).eq.1))then
            ichannel = 1     ! N pi coll.
            multipl  = 2
          else if((ipi(1,i1).eq.2).and.(id(1,i2).eq.1))then
            ichannel = 2     ! N eta coll.
            multipl  = 2
          else if(ipi(1,i1).eq.3 .and. (id(1,i2).eq.1))then
            ichannel = 4     ! N rho coll.
            multipl  = 6
          else if(ipi(1,i1).eq.4 .and. (id(1,i2).eq.1))then
            ichannel = 3     ! N sigma coll.
            multipl  = 2
          else if(ipi(1,i1).eq.5 .and. (id(1,i2).eq.1))then
            ichannel = 9     ! N omega coll.
            multipl  = 6
          else if((ipi(1,i1).eq.1).and. id(1,i2).eq.2 )then
            ichannel = 5     ! Delta(1232) pi coll.
            multipl  = 4
          else if((ipi(1,i1).eq.1).and.id(1,i2).eq.3)then
            ichannel = 6     ! N(1440) pi coll.
            multipl  = 2
          else if((ipi(1,i1).eq.6).and.id(1,i2).eq.nres+2)then
            ichannel = 8     ! K + Lambda coll.
            multipl  = 2
          else if((ipi(1,i1).eq.6).and.id(1,i2).eq.nres+3)then
            ichannel = 7     ! K + Sigma coll.
            multipl  = 2
          else
            goto 600
          end if

          if(ichannel.eq.7 .or. ichannel.eq.8) write(*,*) 'pionab kaon'

          if(id(1,i2).gt.3 .and. ipi(1,i1).ne.6) then
            write(*,*)'error in pionab id = ', id(1,i2),ipi(1,i1)
          end if
c          write(*,*)'t6 ', id(2,816)
*------- mom. or resonance , needed for in-medium width(option for
*------- later build ins)
          pres = 0.0

          call f77flush()
          if(id(1,i2).eq.1) then
            if(id(2,i2)/2.ne.0) then
              write(*,*) 'hiba in pionab 0',id(1,i2),id(2,i2)
              stop
            end if
          else if(id(1,i2).le.nres+1) then
            if(resprop2(abs(id(1,i2))-1,1).eq.1) then
              if(id(2,i2).gt.1 .or. id(2,i2).lt.0) then
                write(*,*)'charge violation in pionab 0f'
                write(*,*)i2, id(1,i2), id(2,i2)
                stop
              end if
            else if(resprop2(abs(id(1,i2))-1,1).eq.3) then
              if(id(2,i2).gt.2 .or. id(2,i2).lt.-1) then
                write(*,*)'charge violation in pionab 0g'
                write(*,*)i2, id(1,i2), id(2,i2)
                stop
              end if
            end if
          end if

*-------- charges
          qmes = ipi(2,i1)
          qn   = id(2,i2)
          qtot = qmes + qn
          sigt = 0.0
c          write(*,*)'t7 ', id(2,816)
c         write(*,*)'pionab *******************************'
c         write(*,*)'gesamt ladg = ',qtot
c         write(*,*)'in pionab mesosn : ', i1
c         write(*,*)'id1, id2 ', ipi(1,i1), ipi(2,i1)
c         write(*,*)'in pionab nucl : ', i2
c         write(*,*)'id1, id2' , id(1,i2), id(2,i2)
c         write(*,*)' *******************************'


******------- LOOP OVER ALL RESONANCES
c          write(*,*) 'in pionab, before loop over all Res.'
          nloopind = nres
          if(isdeltaall.eq.1) nloopind = 1
          if(isdeltaonly.eq.1) nloopind = 1
          do i = 1, nloopind
******     the energy of the res. is dep. if one uses Delta or
******     nucl. potential
c          write(*,*)'t8 ', id(2,816)
            if(resprop2(i,1).eq.1) then
              decision = channel(1)
            else if(resprop2(i,1).eq.3) then
             decision = channel(2)
            end if

            if(decision) then !!!!!! reaction is energetically possib.
              if(resprop2(i,1).eq.1) then
                srts = finmass(1)
              else if(resprop2(i,1).eq.3) then
                srts = finmass(2)
              end if
*------- momentum of the pion in the rest-frame of the res
              qq2= 0.25*(srts**2-em22+em12)**2/srts**2-em12

*------- calculate the isospin factors for the reactions ------------*
              if(ichannel.eq.1 .or. ichannel.eq.4 .or. ichannel.eq. 6
     &                      .or. ichannel.eq. 7)then
**       isospin 1/2 x 1
                if(resprop2(i,1).eq.3) then
                  if((qtot.eq.2).or.(qtot.eq.-1)) then
                    isofac=1.
                  else if(qmes.eq.0) then
                    isofac=2./3.
                  else
                    isofac=1./3.
                  end if
                else if(resprop2(i,1).eq.1) then
                  if((qtot.eq.2).or.(qtot.eq.-1)) then
                    isofac=0.
                  else if(qmes.eq.0) then
                    isofac=1./3.
                  else
                    isofac=2./3.
                  end if
                end if
              else if((ichannel.eq.2).or.(ichannel.eq.3)
     &          .or.(ichannel.eq.8).or.(ichannel.eq.9)) then
** 1/2 x 0
                isofac = 1.0
              else if(ichannel.eq.5)  then
** 3/2 x 1
                if(resprop2(i,1).eq.3) then
                  if((qtot.eq.3).or.(qtot.eq.-2)) then
                    isofac=0.
                  else if((qtot.eq.2).or.(qtot.eq.-1)) then
                    if(qmes.eq.0) then
                      isofac=3./5.
                    else
                      isofac=2./5.
                    end if
                  else
                    if(qmes.eq.0) then
                      isofac=1./15.
                    else if(qmes.eq.(-2*qtot+1)) then
                      isofac=2./5.
                    else
                      isofac=8./15.
                      if(qmes.ne.(2*qtot-1)) write(*,*)'Fehler im
     &                   Isospinteil von pionab',qmes,qtot,qn,
     &                   ichannel
                    end if
                  end if
                else
                  if((qtot.gt.1).or.(qtot.lt.0)) then
                    isofac=0.
                  else
                    if(qmes.eq.0) then
                      isofac=1./3.
                    else if(qmes.eq.(-2*qtot+1)) then
                      isofac=1./2.
                    else
                      isofac=1./6.
                      if(qmes.ne.(2*qtot-1)) write(*,*)'Fehler im
     &                    Isospinteil von pionab',qmes,qtot,qn,
     &                    ichannel
                    end if
                  end if
                end if
              end if
c                        write(*,*)'t9 ', id(2,816)
*----------------- end of isospin factors -------------------------*
c          write(*,*) 'in pionab, after isospin',i

              help = bwdist(srts**2,i,idec2pi,iresmode,ichannel,
     &                     0,iwidth,1,dummyf)

              sig(i) = 40.0*pi/qq2 * isofac * (abs(resprop2(i,3))+1.)/
     &                real(multipl) * help * hbc**2
c             write(*,*) "in pionab sig(i) ",i,sig(i),isofac,multipl,
c     &          help,inotwopi,ichannel,resprop2(i,1),resprop2(i,3),
c     &          qtot,qmes,srts,iresmode,iwidth,ipi(1,i)
*----------- sum up sigmas
              sigt = sigt + sig(i)
*--------- end of decision if
            end if
********** END OF RES. LOOP
          end do
c          write(*,*) 'in pionab, after loop over all Res.'

c          write(*,*)'t10 ', id(2,816)
********** CHECK FOR ABSORPTION
c       if (ipi(1,i1) .eq. 2) write(*,*) ' t11z ',maxcross,sigt,iseed
          if(maxcross.lt.sigt) maxcross = sigt
          xx = rn(iseed)

c only Delta and n1520 exists:
c             sig(2) = 0.
c             do i=4,nres
c               sig(i) = 0.
c             end do
c             sigt = sig(1) + sig(3)

          if(xx .le. sigt/maxcross) then
**--------------- meson will be absorbed
c          write(*,*)'t12 ', id(2,816)
            testsig = 0.0
            yy = rn(iseed)
            stopflag = .true.

            i = 0
            do while(stopflag)
              i = i + 1
              testsig = testsig + sig(i)
c                 write(*,*) "inside pionab testsig",i,sig(i),testsig,
c     &              sigt,yy
              if(yy .le. testsig/sigt) then
***     got it
c          write(*,*)'t13 ', id(2,816),i,resprop2(i,1),i2
                stopflag = .false.
                totabs = totabs + 1

                if(resprop2(i,1).eq.1) then
                  upot(i2)= scapot(1)
                else if(resprop2(i,1).eq.3) then
                  upot(i2)= scapot(2)
                end if
c                write(*,*) 'pionab 3',scapot(1),scapot(2)
              end if
            end do

          else
c          write(*,*)'t14 ', id(2,816)
***-------------- meson will not be absorbed
            goto 600
          end if
***-------------- END Of ABSORPTION CHECK IF
c          write(*,*)'t15 ', id(2,816)

          e(i2)    = srts
          id(1,i2) = i + 1
          id(2,i2) = qtot
          id(4,i2) = ipi(4,i1)
          id(5,i2) = id(5,i2) + 1
          id(7,i2) = ipi(1,i1)
c         write(*,*)'t17 ', id(2,816)
          id(3,i2) = 0
          id(8,i2) = 0

c          write(*,*)'t18 ', id(2,816)
          if(resprop2(abs(id(1,i2))-1,1).eq.1) then
            if(id(2,i2).gt.1 .or. id(2,i2).lt.0) then
              write(*,*)'charge violation in pionab 1'
              write(*,*)i2,id(1,i2),id(2,i2),qmes,qn,ipi(1,i1)
     &                                         ,ichannel
              stop
            end if
          else if(resprop2(abs(id(1,i2))-1,1).eq.3) then
            if(id(2,i2).gt.2 .or. id(2,i2).lt.-1) then
              write(*,*)'charge violation in pionab 2'
              write(*,*)i2, id(1,i2), id(2,i2)
              stop
            end if
          end if
c          write(*,*)'t19 ', id(2,816)

           if(ichannel.eq.1) then
             if(i.eq.1) lmesa(1)=lmesa(1)+1
             if(i.eq.2) lmesa(3)=lmesa(3)+1
             if(i.eq.4) lmesa(2)=lmesa(2)+1
             if(i.eq.3 .or. i.gt.4) lmesa(3)=lmesa(3)+1
           else if(ichannel.eq.2) then
             lmesa(4)=lmesa(4)+1
c           if(i.ne.3) then
c             write(*,*)'hiba in pionab etaabs',ichannel,i
c           end if
           else if(ichannel.eq.5) then
             if(i.eq.2) lmesa(6)=lmesa(6)+1
             if(i.eq.4) lmesa(5)=lmesa(5)+1
             if(i.eq.3 .or. i.gt.4) lmesa(6)=lmesa(6)+1
           else if(ichannel.eq.6) then
             if(i.eq.4) lmesa(7)=lmesa(7)+1
             if(i.ne.4) lmesa(8)=lmesa(8)+1
           else if(ichannel.eq.4) then
            if(i.eq.2) lmesa(10)=lmesa(10)+1
            if(i.eq.4) lmesa(9)=lmesa(9)+1
            if(i.eq.3 .or. i.gt.4) lmesa(10)=lmesa(10)+1
          else if(ichannel.eq.3) then
            if(i.eq.2) lmesa(12)=lmesa(12)+1
            if(i.eq.4) lmesa(11)=lmesa(11)+1
            if(i.eq.3 .or. i.gt.4) lmesa(12)=lmesa(12)+1
          else if(ichannel.eq.9) then
            if(i.eq.2) lmesa(13)=lmesa(13)+1
            if(i.eq.4) lmesa(14)=lmesa(14)+1
            if(i.eq.3 .or. i.gt.4) lmesa(13)=lmesa(13)+1
          else if(ichannel.eq.8) then
            lmesa(16)=lmesa(16)+1
          else if(ichannel.eq.7) then
            lmesa(17)=lmesa(17)+1
          end if

*lmesa-Legende
* n + pi -> delta    : 1
* n + pi -> 1535     : 2
* n + pi -> N*       : 3
* n + eta -> 1535    : 4
* d + pi -> 1535     : 5
* d + pi -> N*       : 6
* 1440 + pi -> 1535  : 7
* 1440 + pi -> N*    : 8
* n + rho -> n*      : 10
* n + rho -> 1535    : 9
* n + sigma -> N*    : 12
* n + sigma -> 1535  : 11
* n + omega -> N*    : 13
* n + omega -> 1535  : 14
* K + Lambda -> res  : 16
* K + Sigma -> res   : 17
*-----------------------------------------------------------------------
c          write(*,*) 'in pionab, old statistic routines'

*----------- old statistic routines --------------------------------*
          idres = id(1,i2) - 1
          call resprod(i2,0.5*(t1+t2),id(6,i2),id62,idres,2)
*-----------------------------------------------------------------------
c          write(*,*)'t22 ', id(2,816)
*                 path length distribution
          path=sqrt((rpie(1,i1)-x1)**2+(rpie(2,i1)-y1)**2+
     &                      (rpie(3,i1)-z1)**2)
          npat=min(30,nint(2.0*path))
          if (npat .lt. 0) then
cc        write(mspfpri,*) ' i1, ipi(1,i1)', i1, ipi(1,i1)
c         write(mspfpri,*)   rpie(1,i1), x1, rpie(2,i1),y1,
c    2                      rpie(3,i1),z1
          endif
          mpath(ipi(1,i1),npat)=mpath(ipi(1,i1),npat) + 1
c          write(*,*)'t23 ', id(2,816)
*-----------------------------------------------------------------------
*               lifetime distribution
          path=max( 0.0, time - rpie(6,i1))
          npat=min(50,nint(2.0*path))
          melet(ipi(1,i1),npat)=melet(ipi(1,i1),npat) + 1
*-----------------------------------------------------------------------
c          write(*,*)'t24 ', id(2,816)
          p(1,i2)= px1 + px2
          p(2,i2)= py1 + py2
          p(3,i2)= pz1 + pz2
          if(srtp.le.rmass+pmass) write(*,'(''hiba:pionab,mass'',f9.3)')
     &               srt

          ipi(1,i1) = 0
c          write(*,*)'t25 ', id(2,816)
********* check for energy conservation
          phelp = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
          test = consen -sqrt((e(i2)+ upot(i2))**2+phelp**2)
c          write(*,*)'t26 ', id(2,816)
          if(abs(test).gt.1e-05) then
            write(122,*)'hiba energy conservation in pionab ', test
          end if
c          write(*,*) 'in pionab, before continue 800'


  800   continue
 1000 continue

      write(*,*)'total pi abs : ', totabs
c      tout = secnds(t0)
c     tin = tout - tin
c     write(*,*)'time ellapsed in pionab = ',tin,'  sec.'
      return
      end
