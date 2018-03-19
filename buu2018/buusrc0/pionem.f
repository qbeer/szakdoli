************************************************************************
*                                                                      *
      subroutine pionem(decti,lmesc,flag, cres)
*                                                                      *
*     purpose: calculating meson production from resonance decays      *
*     -----------------------------------------------------------      *
*     possible final states:                                           *
*        N pi, N eta, Delta pi, N(1440) pi, N rho, N sigma             *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"

      common /nthwhw/  nthw

      integer  nhw0, nhw1, nhw2
      real*8     hw1, hw2, hw3, hw4, j01, beta(3)
      integer dimlmesc, nthw, k_hw, dummy
      parameter(dimlmesc = nres+9)
      real*8 ede, ema2, rm2, pm2, gamma, dt0, w, rnxx
      real*8 xx
      real*8 bmass,   eba, decti, stot, dti,  bm2, wmass
*      real*8 gam1pi, gameta, gam2pi, gamdpi, gamnrp, gamnro, gamnsi, wmass
      real*8 gam1pi, gameta,  gamdpi, gamnrp, gamnro, gamnsi,gam,gamnom
      real*8 gamsik, gamlak
      integer irun, inp, ini, ii, jj, id1, ib, ichannel
      integer ire, ibin, ntag, idj, id7, idn, ixx, iyy, izz, iendel
      real*8  rnx,   emd2,  phase, dendel
      real*8 rn, bwdist, pdx, pdy, pdz,   x1, x, y, z
      integer lmesc(1:dimlmesc), nptest1
      integer cres(1:nres,1:9), resin
      logical gridte
      integer qtest, imomch,iii
*     stuff needed for the mom-dep. potentials

*      number of channels possible depending on idec2pi
      integer momnum(0:1), monuch
*      max 9 poss. channels, 3:where is b (1or 2),
*                            4:where is m (1 or 2)
*                         1,2: id of mes or b
      integer momid(1:9,1:4), idhelp
*      masses of the particles in the final state
      real*8    mommass(1:2)
*     beta to boost from calc.fram into the LRF of resonance
      real*8    betlrfx(1:2), betlrfy(1:2),betlrfz(1:2)
*     beta to boost from calc.frame into the res frame of res.
      real*8    betacm(1:3)
*     maximal potentials needed in between
      real*8    maxbarpot, maxmespot, vecpot
      real*8    phx, phy, phz, eh, pabsh, vecmax1, vecmax2
      real*8    srtit, ema, massmin, remass

      logical channel, flag
      real*8    momp(0:3), momscapot(1:2)
      real*8    gamst(1:9)
      logical nope, negp, massflag
      integer nptest, j, i, ihelp
      real*8    deriv(0:4), j0,j1, j2, j3, u3, u4
      real*8    pin, potanal, potmes, pxc, pyc, pzc
      integer npot
      real*8    bapot, mespot, vecpotb, vecpotm
c     real*8    t0, tin, tout
      real*8    consen ,  consenout  , pinabs
      real*8    massba, massme, ratio(9)
      real*8    stlrfprop(1:4), check
      integer totdec

      real*8    crx, cry, crz, cpx, cpy, cpz, cpot
      real*8    cfox, cfoy, cfoz, etot
      integer cid2

      real*8    plrf(1:3), rhap(1:3), rhap2(1:3)
      real*8    density, vrel

*     0, 1 : corresponds to idec2pi:: number of poss. final cahnnels   *
      data    momnum(0) /2/
      data    momnum(1) /9/
*---------- possible final channels------------------------------------*
*     1 :  N pi
*     2 :  N eta
*     3 :  n sigma
*     4 :  N rho
*     5 :  Delta pi
*     6 :  N14   pi
*     7 :  Sigma K
*     8 :  Lambda K
*     9 :  N omega
*     initialized :   =, b-id, mes-id
      data (momid(1,j),j= 1,4) / 1,  1, 0, 0/
      data (momid(2,j),j= 1,4) / 1,  2, 0, 0/
      data (momid(3,j),j= 1,4) / 1,  4, 0, 0/
      data (momid(4,j),j= 1,4) / 1,  3, 0, 0/
      data (momid(5,j),j= 1,4) / 2,  1, 0, 0/
      data (momid(6,j),j= 1,4) / 3,  1, 0, 0/
      data (momid(7,j),j= 1,4) / 27, 6, 0, 0/
      data (momid(8,j),j= 1,4) / 26, 6, 0, 0/
      data (momid(9,j),j= 1,4) / 1,  5, 0, 0/

c      write(*,*) '  start pionem ', dimlmesc
c      call f77flush()
      call potcalc
      totdec = 0

c     t0 = 0.0
c      tin = secnds(t0)

*----------------------------------------------------------------------*
      do 30 ii=1,dimlmesc
        lmesc(ii) = 0
  30  continue

      do ii = 1, nres
        do iii = 1,9
          cres(ii,iii) = 0
        end do
      end do

      pm2  = pmass**2
      rm2  = rmass**2

*   loop over all parallel runs*
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        ini = (irun-1) * maxb
        do 800 ii  = 1,maxb

          jj = ini + ii
          if(abs(id(1,jj)) .le. 1 .or. id(1,jj).gt.nres+1)      goto 800

*----------------------------------------------------------------------*
*      randomize
      monuch = momnum(0)
      if(idec2pi .ge. 1) monuch = momnum(1)

*---------------------------------------------------------------------*
*     store resonance properties                                      *
          id1      = id(1,jj)
          ema      = e(jj) + upot(jj)
          ema2     = ema**2
          ede      = sqrt(ema2+p(1,jj)**2+p(2,jj)**2+p(3,jj)**2)
          x        = r(1,jj)
          y        = r(2,jj)
          z        = r(3,jj)
          gamma    = ede / ema
          dt0      = decti / gamma

          rhap(1) = x
          rhap(2) = y
          rhap(3) = z

*              kinematics for all the final states is calculated     *
*--------------------------------------------------------------------*
*              determine the free sqrt(s) for each final channel     *
*              according to the above kinematics                     *

          resin = id1-1
c          write(*,*) 'resonance type in pionem: ',resin,jj,irun
c          call f77flush()

          idn = id1 - 1
          gam = 0.0
          emd2 = e(jj)**2
          gam  = bwdist(emd2,idn,idec2pi,iresmode,
     +                           0,0,iwidth,0,ratio)

c          write(*,*) 'call baryon_dalitz:',jj,id1,e(jj),upot(jj)
c
          beta(1) = p(1,jj) / ede
          beta(2) = p(2,jj) / ede
          beta(3) = p(3,jj) / ede
          if(iresdalitz.eq.1 .and.id(2,jj).le.1 .and.id(2,jj).ge.0) then
c          write(*,*) 'pionem resdalitz',id1-1,id(2,jj),gam,e(jj),dt0,
c     &                beta,gamma,x,y,z
          call baryon_dalitz(idn,id(2,jj),gam,e(jj),dt0,beta,gamma,
     &                x,y,z)
c          write(*,*) 'pionem resdalitz',id1-1,id(2,jj),gam,e(jj)
          end if
            do i = 1, 9
              gamst(i) = gam*ratio(i)
c  n1520 can only decay into rho:
c              if (id1.eq.4 .and. i.ne.4) gamst(i) = 0.
c              if (id1.eq.4 .and. i.eq.4) gamst(i) = gam
            end do
c            if (idn.eq.1) then
c              write(*,*) 'gam, ratio:',gam,ratio(1)
c            end if

            check = 0.0
            do i = 1, 9
              check = check + gamst(i)
            end do

            if(abs(check-gam).gt.1.0e-06) then
              write(*,*)'pionem 2'
              write(*,*)check, gam,e(jj), e(jj)-rmass-pmass
              stop
            end if

 19         continue
            gam = 0.0
            do i = 1,9
              gam = gam + gamst(i)
            end do

          w   = exp(-dt0 * gam / hbc / dlife)
          rnxx = rn(iseed)
c          if(resin.eq.3) write(*,*) "pionem, 1520",dlife,gam,w,
c     &       rnxx
          if(rnxx .gt. w .or. flag) then

            totdec = totdec + 1
*----------------------------------------------------------------------*
*               the resonance may decay                                *
*                                                                      *

*------------- store the partial widths

            gam1pi=gamst(1)
            gameta=gamst(2)
            gamnsi=gamst(3)
            gamnro=gamst(4)
            gamdpi=gamst(5)
            gamnrp=gamst(6)
            gamsik=gamst(7)
            gamlak=gamst(8)
            gamnom=gamst(9)

c            if (id1.eq.2) then
c              write(*,*) 'branching ratios:'
c              write(*,*) ' gam1pi=', gam1pi
c              write(*,*) ' gameta=', gameta
c              write(*,*) ' gamnsi=', gamnsi
c              write(*,*) ' gamnro=', gamnro
c              write(*,*) ' gamdpi=', gamdpi
c              write(*,*) ' gamnrp=', gamnrp
c              write(*,*) ' gamsik=', gamsik
c              write(*,*) ' gamlak=', gamlak
c              write(*,*) ' gamnom=', gamnom
c            end if

c            if(abs(gam-gam1pi-gameta-gamdpi-gamnrp-gamnro-gamnsi)
c     &       .gt.1.e-6) then
c              write(*,*) 'hiba pionem gam',gam,id(1,jj),e(jj)
c              write(*,*)'diff = ', gam-gam1pi-gameta-gamdpi-gamnrp-
c     &                             gamnro-gamnsi
c              write(*,*)'ch ', channel
c              write(*,*)'gams', gamst
c            end if

*------------------- pick the particular channel --------------------*

 42         rnx = rn(iseed)
            ichannel = 0
            imomch   = 0
c        write(*,*) ' in pionem  gammas=', gam1pi,gameta,gamdpi,
c    1                gamnrp,gamnro,gamnsi

            if(rnx .lt. gam1pi/gam .or.flag ) then
*             nucleon pion may be produced
              ichannel = 1
              imomch   = 1
              ire      = 1
            else if(rnx .lt. (gam1pi+gameta)/gam) then
*             nucleon eta  may be produced
              ichannel = 2
              imomch   = 2
              ire      = 2
            else if(rnx .lt. (gam1pi+gameta+gamdpi)/gam) then
*             delta pion is produced
              ichannel = 6
              imomch   = 5
              ire      = 1
            else if(rnx .lt. (gam1pi+gameta+gamdpi+gamnrp)/gam) then
*             N(1440) pion is produced
              ichannel = 7
              imomch   = 6
              ire      = 1
            else if(rnx .lt. (gam1pi+gameta+gamdpi+gamnrp+gamnro)/gam)
     &                      then
*             rho N may be produced
c-hw                               this is the important rho production
              k_hw = 1          !       for on-off switching
              if (k_hw .eq. 1) then
              ichannel = 5
              imomch   = 4
              ire      = 3
              endif
            else if(rnx .lt. (gam1pi+gameta+gamnro+gamnrp+gamdpi+
     &                        gamnsi)/gam) then
*             sigma N may be produced
              ichannel = 4
              imomch   = 3
              ire      = 4
            else if(rnx .lt. (gam1pi+gameta+gamdpi+gamnrp+gamnro+gamnsi
     &              +gamlak)/gam) then
*             kaon + Lambda may be produced
              ichannel = 9
              imomch   = 8
              ire = 6
            else if(rnx .lt. (gam1pi+gameta+gamdpi+gamnrp+gamnro+gamnsi
     &              +gamlak+gamsik)/gam) then
*             kaon + Sigma may be produced
              ichannel = 8
              imomch   = 7
              ire = 6
            else if(rnx .lt. (gam1pi+gameta+gamdpi+gamnrp+gamnro+gamnsi
     &              +gamlak+gamsik+gamnom)/gam) then
*             omega N may be produced
              ichannel = 10
              imomch   = 9
              ire = 5
            end if

c           if (jj.eq.3626) write(*,*) 'ichannel',ichannel
c           call f77flush()


            if(ichannel .eq. 0) goto 42 !!security on the 10e-08 level

*                hopefully we have found a channel                    *
*---------------------------------------------------------------------*
*              do the max-pots. i.o. to find a mass                   *

c   it has no sence here on my opinion, perhaps in BB collision
c     momid(i,3) is the baryon place in momid(i,momid(i,3))=id1 
c     here it can simply be momid(i,1)

            x1 = rn(iseed)
            if(x1 .le. 0.5) then
              momid(imomch,3) = 1
              momid(imomch,4) = 2
            else
              momid(imomch,3) = 2
              momid(imomch,4) = 1
              idhelp     = momid(imomch,1)
              momid(imomch,1) = momid(imomch,2)
              momid(imomch,2) = idhelp
            end if

            channel      = .false.
            mommass(1)   = 0.0
            mommass(2)   = 0.0
            momp(0)      = 0.0
            momp(1)      = 0.0
            momp(2)      = 0.0
            momp(3)      = 0.0

*---------------------------------------------------------------------*
*     beta that boosts from calc.frame into the resframe of resonance *
          betacm(1)=   p(1,jj) / ede
          betacm(2)=   p(2,jj) / ede
          betacm(3)=   p(3,jj) / ede

*---------------------------------------------------------------------*
*     determine the density in the LRF  and the trafo-variables for   *
*     the boost into the LRF                                          *
          if(isplipi .eq. 1) then

            call splinint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

          else if(isplipi .eq. 0) then
            call linint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)
          end if


           if(j0 .gt. 1.0e-6) then
             betlrfx(1) = j1/j0
             betlrfy(1) = j2/j0
             betlrfz(1) = j3/j0
             betlrfx(2) = j1/j0
             betlrfy(2) = j2/j0
             betlrfz(2) = j3/j0
           else
             betlrfx(1) = 0.0
             betlrfy(1) = 0.0
             betlrfz(1) = 0.0
             betlrfx(2) = 0.0
             betlrfy(2) = 0.0
             betlrfz(2) = 0.0
           end if

              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba pionem1 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
           call lorentz(betlrfx(1),betlrfy(1),betlrfz(1),j1,j2,j3,j0)

           stlrfprop(1) = betlrfx(1)
           stlrfprop(2) = betlrfy(1)
           stlrfprop(3) = betlrfz(1)
           stlrfprop(4) = j0

           density = j0



*---------------------- baryon --------------------------------------*
           remass    = rmass
           if(imomch.le.4) then
             remass    = rmass
           else if(imomch .eq.5) then
             remass    = resprop1(1,1)
           else if(imomch.eq.6) then
             remass    = resprop1(2,1)
           else if(imomch.eq.7) then
             remass    = xsmas
           else if(imomch.eq.8) then
             remass    = xlmas
           else if(imomch.eq.9) then
             remass    = rmass
           end if

           i  = imomch !!!!!!!!!!
           pin       = 0.0
           plrf(1) = 0.0
           plrf(2) = 0.0
           plrf(3) = 0.0

           vecpot    = potanal(rho0,j0,plrf,momid(i,momid(i,3)),
     &                         rhap)
c          write(*,*) "pionem vecpot1",rho0,j0,plrf,momid(i,3),
c    1             momid(i,momid(i,3)), vecpot

           maxbarpot = -remass + sqrt(remass**2 +
     +                   2.0*remass*vecpot + vecpot**2)


*     for the upper bound in the mass determination one has to express
*     the max. scalar pot. in terms of the vector pot.. therefor one
*     needs the mom. in this frame that corresponds to p=0 in the LRF

c    eh in the local rest frame
          phx     = 0.0
          phy     = 0.0
          phz     = 0.0
          eh      = 5.0*remass +  maxbarpot
          
          call lorentz(-betlrfx(1),-betlrfy(1),-betlrfz(1),
     +                 phx,phy,phz,eh)
          call lorentz(  betacm(1),  betacm(2),  betacm(3),
     +                   phx,phy,phz,eh)
          pabsh   = sqrt(phx**2+phy**2+phz**2)

          vecmax1 = -sqrt(pabsh**2+(5.0*remass)**2)+
     +               sqrt(pabsh**2+(5.0*remass)**2+
     +              2.0*5.0*remass*maxbarpot+maxbarpot**2)


*-------------------------- meson -----------------------------------*

           if(imomch.eq.1)then
             remass    = pmass
           else if(imomch.eq.2)then
             remass    = emass
           else if(imomch .eq.3) then
             remass    = 0.800
           else if(imomch.eq.4) then
             remass    = romas
           else if(imomch.eq.5) then
             remass    = pmass
           else if(imomch.eq.6) then
             remass    = pmass
           else if(imomch.eq.7) then
             remass    = xkmas
           else if(imomch.eq.8) then
             remass    = xkmas
           else if(imomch.eq.9) then
             remass    = omass
           end if

           i = imomch !!!!!!!!!!!!!
           pin       = 0.0
           vecpot    = potmes( j0,pin,momid(i,momid(i,4)))
           maxmespot = -pmass + sqrt(pmass**2 +
     +                  2.0*pmass*vecpot + vecpot**2)

          phx     = 0.0
          phy     = 0.0
          phz     = 0.0
          eh      = remass + maxmespot
          call lorentz(-betlrfx(2),-betlrfy(2),-betlrfz(2),
     +                  phx,phy,phz,eh)
          call lorentz(  betacm(1),  betacm(2),  betacm(3),
     +                  phx,phy,phz,eh)
          pabsh   = sqrt(phx**2+phy**2+phz**2)
          vecmax2 = -sqrt(pabsh**2+remass**2)+sqrt(pabsh**2+remass**2+
     +               2.0*remass*maxmespot+maxmespot**2)

*-------------------------- together ---------------------------------*
*  maximal possible srts(free) for the particles in the final channel *
          srtit   =  e(jj)+upot(jj) - vecmax1 - vecmax2   !!!
          consen  =  e(jj) +  upot(jj)
c          write(*,*) "pionem srtit ",e(jj),upot(jj),vecmax1,vecmax2,
c     &       srtit,consen,i,jj

          if(srtit .lt. rmass+pmass) then
            write(*,*)'something wrong with srtit in pionem'
            write(*,*)'srtit   = ', srtit, rmass,pmass
            write(*,*)'jj      = ', jj
            write(*,*)'id(1,jj)= ', id(1,jj) , e(jj), upot(jj)
            write(*,*)'ema     = ', ema
            write(*,*)'vecmax1 = ',vecmax1, j0,vecpot
            write(*,*)'vecmax2 = ',vecmax2, j0
            srtit = rmass+pmass + .001                  !           hw
          end if

*------------- done with the dumb stuff ------------------------------*

        nope     =.true.
        negp     =.true.
        nptest1  = 0
        massflag =.false.
        i        = imomch
        channel  = .false.

        do while(.not.channel .and. nptest1.lt.10)
          nptest1 = nptest1 + 1

*            do the masses for the different channels
            if(i .eq. 1) then
              massmin = rmass + pmass
              if (srtit .ge. massmin ) then
                mommass(momid(i,3)) = rmass
                mommass(momid(i,4)) = pmass
                massflag = .true.
              end if
            else if(i .eq. 2) then
              massmin = rmass + emass
              if (srtit .ge. massmin) then
                mommass(momid(i,3)) = rmass
                mommass(momid(i,4)) = emass
                massflag  = .true.
              end if
            else if(i .eq. 3) then
              massmin = 2.0*pmass + rmass
              if(srtit .ge. massmin+0.001) then
                mommass(momid(i,3)) = rmass
                call resmasdi(srtit, 2.0*pmass, rmass, -2,
     +                      mommass(momid(i,4)))
                massflag = .true.
              end if
            else if((i .eq. 4) .or. (i .eq. 9)) then
              massmin = 2.0*pmass+rmass
              if(srtit .ge. massmin+0.001) then
                mommass(momid(i,3)) = rmass
c----------------------------------
c                call resmasdi(srtit, 2.0*pmass, rmass, -1,
c     +                      mommass(momid(i,4)))
c               write(*,77543) srtit,pabsh,eh,
c    &             betlrfx(1),betlrfy(1),betlrfz(1),
c    &             betacm(1),betacm(2),betacm(3)
c77543           format('before resmas_mes srtit,pabsh,eh ',
c     &             3e14.5,/,'beta = ',3f8.4,/,'betar = ',3f8.4)
c----------------------------------
                vrel = pabsh/eh   ! vrel (vel. rel. to medium) - needed later
c----------------------------------
c stlrfprop(4): density
c srtit: resonance mass + potentials (??)
c                write(*,*) 'before call resmas_mes',stlrfprop(4),srtit,
c     &             mommass(momid(i,4))
c                call f77flush()
c----------------------------------

                call resmas_mes(vrel,ire,stlrfprop(4),srtit,
     &             mommass(momid(i,4)))
c         write(*,*) 'now a rho-omega could be here: ',jj,
c     &   i,mommass(momid(i,4)),ire,stlrfprop(4),e(jj),srtit,nptest1
                call f77flush()
                massflag = .true.
              end if
            else if(i .eq. 5) then
              massmin = rmass+pmass+pmass
              if(srtit .ge. massmin+0.001) then
                mommass(momid(i,4)) = pmass
                call resmasdi(srtit, rmass+pmass, pmass, 1,
     +                      mommass(momid(i,3)))
                massflag = .true.
              end if
            else if(i .eq. 6) then
              massmin = rmass+pmass+pmass
              if(srtit .ge. massmin+0.001) then
                mommass(momid(i,4)) = pmass
                call resmasdi(srtit, rmass+pmass, pmass, 2,
     +                      mommass(momid(i,3)))
                massflag = .true.
              end if
            else if(i .eq. 7) then
              massmin = xsmas + xkmas
              if (srtit .ge. massmin) then
                mommass(momid(i,3)) = xsmas
                mommass(momid(i,4)) = xkmas
                massflag  = .true.
              end if
            else if(i .eq. 8) then
              massmin = xlmas + xkmas
              if (srtit .ge. massmin) then
                mommass(momid(i,3)) = xlmas
                mommass(momid(i,4)) = xkmas
                massflag  = .true.
              end if
c            else if(i .eq. 9) then
c              massmin = rmass + omass
c              if (srtit .ge. massmin) then
c                mommass(momid(i,3)) = rmass
c                mommass(momid(i,4)) = omass
c                massflag  = .true.
c              end if
            end if
*          do the momenta
            if(massflag) then
            npot = 1
            negp =.true.
            nptest = 0
            hw1 = .0
            hw2 = .0
            hw3 = .0
            hw4 = .0
            nhw2 = 2
            nhw1 = 1
            nhw0 = 0
            j01 = j0
            rhap2(1) = rhap(1)
            rhap2(2) = rhap(2)
            rhap2(3) = rhap(3)
c            write(*,*)  '  pionem : vor while  ', nptest, negp,i
            do while(negp .and. nptest.lt.100)
              nptest = nptest + 1

c              write(*,*)  '  pionem : vor momiter ', nptest
               call f77flush()
             call momiter( consen,npot,mommass(1),mommass(2),
     +                   momid(i,1), momid(i,2), rho0,
     +                   j0,j01,betlrfx,betlrfy,betlrfz,betacm,
     +                   u3, u4, pxc, pyc,pzc,negp,nhw2,momid(i,4),nhw1,
     +                   hw1,hw2,hw3,  nhw0,iseed,hw4,rhap,rhap2)
c             write(*,*) consen,npot,mommass(1),mommass(2),
c     +                   momid(i,1), momid(i,2), rho0,
c     +                   j0,j01,betlrfx,betlrfy,betlrfz,betacm,
c     +                   u3, u4, pxc, pyc,pzc,negp,nhw2,momid(i,4),nhw1,
c     +          hw1,hw2,hw3,  nhw0,iseed,hw4,rhap,rhap2
c             call f77flush()

              if(.not. negp) then
                momp(1)  = pxc
                momp(2)  = pyc
                momp(3)  = pzc
                momp(0)  = sqrt(pxc**2 + pyc**2 +pzc**2)
                channel  = .true.
                if(momid(i,4).eq.1) then
                  momscapot(momid(i,4)) = u3
                  momscapot(momid(i,3)) = u4
                else
                  momscapot(momid(i,4)) = u4
                  momscapot(momid(i,3)) = u3
                end if


                pinabs = sqrt(pxc**2+pyc**2+pzc**2)
                massba = mommass(momid(i,3))
                massme = mommass(momid(i,4))

                consenout=sqrt((massba+ momscapot(momid(i,3)))**2+
     +                    pinabs**2)
     +             +sqrt((massme+ momscapot(momid(i,4)))**2 +
     +                pinabs**2)


                if(abs(consen-consenout).gt.5.0e-05) then
                  write(79,*)'schon prob in pionem oben.',i
                  write(*,*)' prob in pionem oben 1 ', i
                  write(*,*)e(jj), upot(jj), consen, consenout
                  stop
                end if

              end if

              if(nptest.gt.100) then
                write(*,*)'nptest .gt. 100 '
                write(*,*)consen, e(jj), upot(jj), i
              end if

            end do

          end if

            if(nptest1.gt.10) then
              write(*,*)'nptest1 .gt. 10 '
              write(*,*)consen, e(jj), upot(jj), i
            end if

        end do   !end of while loop, nptest1

c        write(*,*)'negp, channel ', negp,channel
c        call f77flush()
c        write(*,*)'gams ', gamst
            if(.not.channel) then
              gamst(imomch)=0.0
              goto 19
            end if
*--------------------------------------------------------------------*
*             store final properties                                 *


c              if (id1.eq.4) then
c                call dec_1520(e(jj),bmass,wmass,bapot,mespot,
c     &             pxc,pyc,pzc,pabsh)
c                write(*,*) "after dec_1520 - rho is created",ire
c                consenout = sqrt((bmass+ bapot )**2+pabsh**2)
c     +             +sqrt((wmass+ mespot)**2+pabsh**2)
c                goto 77702
c              end if
              bmass  = mommass(momid(imomch,3))
              wmass  = mommass(momid(imomch,4))
              bapot  = momscapot(momid(imomch,3))
              mespot = momscapot(momid(imomch,4))
              pxc    = momp(1)
              pyc    = momp(2)
              pzc    = momp(3)
              pabsh  = momp(0)

          consenout = sqrt((bmass+ bapot )**2+pabsh**2)
     +               +sqrt((wmass+ mespot)**2+pabsh**2)

          if(abs(consen-consenout).gt.1.0e-05) then
            write(79,*)'schon prob in pionem mitte.',ichannel
            write(*,*)'schon prob in pionem mitte.',ichannel,
     +       consen-consenout
            stop
          end if


          vecpotb = -sqrt(pabsh**2+bmass**2)+sqrt(pabsh**2+bmass**2+
     +               2.0*bmass*bapot+bapot**2)
          vecpotm = -sqrt(pabsh**2+wmass**2)+sqrt(pabsh**2+wmass**2+
     +               2.0*wmass*mespot+mespot**2)

          srtit=e(jj)+ upot(jj)- vecpotb- vecpotm
            if( srtit .lt. bmass+wmass)then
              write(*,*)'hiba '
              write(*,*)'pionem ', srtit, e(jj),upot(jj),
     1                   vecpotb, vecpotm,bmass,wmass, i
             stop
c             goto 800
            end if

c***********************************************************************
c77702       continue
            ib =inp
  10        ib= ib+ 1
            if(ipi(1,ib).ne.0 .and. ib .lt.(inp+maxp)  ) go to 10
c            write(*,*) "a meson - ib ",ib
c            call f77flush()
            if(ib .gt. (inp+maxp) ) then
              write(isum,'(10x,''***  too many pions  ***'')')
              write(*,'(''warning in pionem *** too many pions  ***'')')
            else



***   store stuff
              pdx      = p(1,jj)
              pdy      = p(2,jj)
              pdz      = p(3,jj)

c      write(*,*)'pioem 12'
              pinabs = sqrt(pxc*pxc+pyc*pyc+pzc*pzc)

              consenout=sqrt((bmass+ bapot)**2+ pinabs**2)
              consenout=consenout+sqrt((wmass+ mespot)**2+pinabs**2)

             if(abs(consen-consenout).gt.1.0e-05) then
                 write(*,*)'test vor boo :',e(jj),upot(jj),consen,
     &            consenout,ema,bmass,wmass
c                 stop
             end if

             bm2      = (bmass + bapot)**2
             eba      = sqrt(bm2 + pxc**2 + pyc**2 + pzc**2 )
             call lorentz(-betacm(1), -betacm(2), -betacm(3),
     +                      pxc, pyc, pzc, eba)
             p(1,jj)  = pxc
             p(2,jj)  = pyc
             p(3,jj)  = pzc
*
              ntag = 0
              if(ipauli.eq.1.and.
     &          (ichannel.ne.6.and.ichannel.ne.7.and.ichannel.ne.8.and.
     &            ichannel.ne.9)) then
                call pauli(jj,ntag,iseed,phase,r(1,jj),r(2,jj),r(3,jj),
     &                                    p(1,jj),p(2,jj),p(3,jj))
                ibin=int(phase*10.0)
                iphapt(ibin)=iphapt(ibin)+1
              endif
              if(flag) ntag = 0
              if (ntag .eq. -1) then
                p(1,jj)= pdx
                p(2,jj)= pdy
                p(3,jj)= pdz
              else
c      write(*,*)'pioem 14'
************************************************************************
c           meson will now be emitted
                ixx = nint(r(1,jj))
                iyy = nint(r(2,jj))
                izz = nint(r(3,jj))
                if(abs(ixx).lt.maxx .and. abs(iyy).lt.maxx
     &                 .and.abs(izz).lt.maxz) then
                  gridte = .true.
                else
                  gridte = .false.
                endif
************************************************************************
c           meson collision number

                id7 = id(7,jj)
                if(id7 .lt.0 .or. id7.gt.6) then
                  write(*,*) 'hiba in pionem, id7 is out of range', id7
                  id7 = 1
                endif
                if(id7.eq.ire) then
                  ipi(4,ib)=id(4,jj)+1
                else
                  idj = min(id(4,jj) , 50)
                  idj = max(idj,1)
                  id7 = max(id7,1)
                  mlife(id7,idj) = mlife(id7,idj) + 1
                  ipi(4,ib) = 1
                  if(id7.ge.2 .and. ire.eq.1) then
                    mestopi  = mestopi + 1
                  else if(id7.eq.1 .and. ire.gt.1) then
                    pitomes  = pitomes + 1
                  else
                    mestomes  = mestomes + 1
                  endif
                endif
************************************************************************
                dti = - log(rnxx) * hbc / gam * gamma
*                write(*,*)'vor resabs '
*        if(nthw .eq. 15)   write(*,*)jj, id(1,jj),e(jj), ire
                call resabs(jj,dti,pdx,pdy,pdz,ire+1)
************************************************************************
c****        store the position, density, mass, and time at creation *
                rpie(1,ib)= r(1,jj)
                rpie(2,ib)= r(2,jj)
                rpie(3,ib)= r(3,jj)
                if(gridte) then
                  rpie(4,ib)= rhb(ixx,iyy,izz)
                else
                  rpie(4,ib) = 0.0
                end if
c               rpie(5,ib)= e(jj)
                rpie(6,ib)= time
                rpie(7,ib) = density

c***********************************************************************

                rpi(1,ib)= r(1,jj)
                rpi(2,ib)= r(2,jj)
                rpi(3,ib)= r(3,jj)
                ppi(1,ib)= pdx - p(1,jj)
                ppi(2,ib)= pdy - p(2,jj)
                ppi(3,ib)= pdz - p(3,jj)
                epi(ib)  = wmass
                rpie(5,ib)= wmass
                mpot(ib) = mespot
                ipi(1,ib)= ire
                ipi(2,ib)= 0
                ipi(3,ib)= jj
                ipi(5,ib)= id(1,jj)
                ipi(6,ib)= id(6,jj)
                ipi(7,ib)= id(4,jj) ! ???? nemertem id(3,jj) volt ????
                ipi(8,ib)= 0
c                if (ire.eq.5) write(*,*)
c     1            'omega mass: ',jj,ib,wmass, rpi(1,ib)
c                if (ire.eq.3) write(*,*)
c     1               'rho mass: ',jj,ib,wmass, rpi(1,ib)

c                if(id(8,jj).ne.0) then  ! ???? nemertem
c                  if(ipi(3,id(8,jj)).eq.jj) then
c                     ipi(6,ib) = id(8,jj)  ! ???? nemertem
c                     ipi(6,id(8,jj)) = ib  ! ???? nemertem
c                  end if
c                end if
c           if (ire.eq.3 .or. ire .eq. 5)
c           if (ire .eq. 5)
c    1      write(mspfpri,*) 'now a rho-omega has a number: ', id1,ib,
c    2      epi(ib), ire, jj, ppi(1,ib), ppi(2,ib), ppi(3,ib)
            if (ire.eq.3 .or. ire .eq. 5) then
c           if (ire.eq.2 ) then
c             write(49,730)  ib, ire, nthw, epi(ib), gameta, gam
c 730         format( ' meson created ', i8,i3,i5, 3f10.4)
              iii = (ire-1)/2
c             iii = 1
c              m_birth(iii) = m_birth(iii) + epi(ib)
c              m2_birth(iii) = m2_birth(iii) + epi(ib)**2
c              n_birth(iii) = n_birth(iii) + 1
c              t_birth(iii) = t_birth(iii) + dt * nthw
            endif
                if(ichannel.lt.3)then
                  cres(resin,ichannel)= cres(resin,ichannel) +1
                else
                  cres(resin,ichannel-1)= cres(resin,ichannel-1) +1
                end if

***          save the lrf properties

                 betlrfboo(jj,1) = stlrfprop(1)
                 betlrfboo(jj,2) = stlrfprop(2)
                 betlrfboo(jj,3) = stlrfprop(3)
                 betlrfboo(jj,4) = stlrfprop(4)

                 betlrfbom(ib,1) = stlrfprop(1)
                 betlrfbom(ib,2) = stlrfprop(2)
                 betlrfbom(ib,3) = stlrfprop(3)
                 betlrfbom(ib,4) = stlrfprop(4)



*--------------------------------------------------------------------*
*              densitiy of meson creation                            *
                if(gridte) then
                  dendel = rhb(ixx,iyy,izz)/rho0
                else
                  dendel = 0.0
                end if
                iendel = nint(5.0*dendel)
                if(iendel .gt. 50) iendel = 50
                mdens(ire,iendel) = mdens(ire,iendel) + 1

c               write(*,*) ' pionem - lmesc ',id(1,jj)-1,nres+imomch
                lmesc(id(1,jj)-1)  = lmesc(id(1,jj)-1) + 1
                lmesc(nres+imomch) = lmesc(nres+imomch) + 1

***********************************************************************
*               charge assignment


                ihelp = 1
                if(ichannel.eq.6) ihelp = 2
                if(ichannel.eq.7) ihelp = 3

               qtest = id(2,jj)

*     charge assignment

                  ipi(2,ib) = 0
                  xx = rn(iseed)

c---:   D -> N pi, N rho, N* pi
                  if((resprop2(id(1,jj)-1,1) .eq. 3).and.
     &                 ((ichannel.eq.1).or.(ichannel.eq.5).or.
     &                 (ichannel.eq.7))) then
                    if(id(2,jj) .eq. 2) then
                      id(2,jj)  = 1
                      ipi(2,ib) = 1

                      if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                        write(*,*)'1  ', ihelp,id(1,jj), id(2,jj)
                      end if

                    else if(id(2,jj) .eq.-1) then
                      id(2,jj)  = 0
                      ipi(2,ib) =-1

                      if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                        write(*,*)'2  ', ihelp,id(1,jj), id(2,jj)
                      end if

                    else
                      if(xx .gt. 0.66666667) then
                        ipi(2,ib) = 2*id(2,jj)-1
                        id(2,jj)  = 1-id(2,jj)
                        if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                          write(*,*)'3  ', ihelp,id(1,jj), id(2,jj)
                        end if
                      end if
                    end if

c                    if(id(1,jj).eq.2) then
c                      npion(ipi(2,ib))=npion(ipi(2,ib))+1
c                    end if

c---:   D -> D pi
                  else if((resprop2(id(1,jj)-1,1).eq.3).and.
     &                   (ichannel.eq.6)) then
                    if((id(2,jj).eq.2).or.(id(2,jj).eq.-1)) then
                      if(xx.gt.0.4) then
                        ipi(2,ib)=(2*id(2,jj)-1)/3
                        id(2,jj)=(id(2,jj)+1)/3
                        if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                          write(*,*)'4  ', ihelp,id(1,jj), id(2,jj)
                        end if
                      end if
                    else
                      if((xx.gt.1./15.).and.(xx.lt.0.6)) then
                        ipi(2,ib)=2*id(2,jj)-1
                        id(2,jj)=1-id(2,jj)
                      if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                         write(*,*)'5  ', ihelp,id(1,jj), id(2,jj)
                      end if


                      else if(xx.gt.0.6) then
                        ipi(2,ib)=-2*id(2,jj)+1
                        id(2,jj)=3*id(2,jj)-1
                        if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                          write(*,*)'6  ', ihelp,id(1,jj), id(2,jj)
                        end if
                      end if
                    end if

c---:   N -> N pi, N rho, N* pi
                  else if((resprop2(id(1,jj)-1,1).eq.1) .and.
     &                   ((ichannel.eq.1).or.(ichannel.eq.5).or.
     &                   (ichannel.eq.7))) then
                    if(xx .gt. 0.33333) then
                      ipi(2,ib) = 2 * id(2,jj) - 1
                      id(2,jj)  = 1 - id(2,jj)
                      if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                         write(*,*)'7  ', ihelp,id(1,jj), id(2,jj)
                      end if
                    end if

c---:   N -> D pi
                  else if((resprop2(id(1,jj)-1,1).eq.1).and.
     &                   ichannel.eq.6) then
                    if(xx .gt. 0.33333 .and. xx.lt. 0.5) then
                      ipi(2,ib) = 2 * id(2,jj) - 1
                      id(2,jj)  = 1 - id(2,jj)
                      if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                         write(*,*)'8  ', ihelp,id(1,jj), id(2,jj)
                      end if


                    else if(xx .gt. 0.5) then
                      ipi(2,ib) = 1 - 2 * id(2,jj)
                      id(2,jj)  = 3 * id(2,jj) - 1
                      if(ihelp.eq.1 .and. id(2,jj).eq.-1) then
                         write(*,*)'9  ', ihelp,id(1,jj), id(2,jj)
                      end if
                    endif

c---:   N -> L K
                  else if(resprop2(id(1,jj)-1,1).eq.1 .and.
     &                        ichannel.eq.9) then
                    ipi(2,ib)=id(2,jj)
                    id(2,jj)=0

c---:   D -> Si K
                  else if(resprop2(id(1,jj)-1,1).eq.3 .and.
     &                     ichannel.eq.8)then
                    if(id(2,jj) .eq. 2) then
                      id(2,jj)  = 1
                      ipi(2,ib) = 1
                    else if(id(2,jj) .eq.-1) then
                      id(2,jj)  = -1
                      ipi(2,ib) = 0
                    else
                      if(xx .gt. 0.66666667) then
                        ipi(2,ib)= 1-id(2,jj)
                        id(2,jj) = 2*id(2,jj)-1
                      else
                        ipi(2,ib)= id(2,jj)
                        id(2,jj) = 0
                      end if
                    end if

c---:   N -> Si K
                  else if(resprop2(id(1,jj)-1,1).eq.1 .and.
     &                         ichannel.eq.8)then
                    if(xx .gt. 0.33333) then
                      ipi(2,ib)= 1-id(2,jj)
                      id(2,jj) = 2*id(2,jj)-1
                    else
                      ipi(2,ib) = id(2,jj)
                      id(2,jj) = 0
                    end if
                  end if

************** store coulomb pot
                crx = rpi(1,ib)
                cry = rpi(2,ib)
                crz = rpi(3,ib)

                cpx = ppi(1,ib)
                cpy = ppi(2,ib)
                cpz = ppi(3,ib)

                cid2 = ipi(2,ib)
                etot = sqrt(pmass**2+cpx**2+cpy**2+cpz**2)
                call emfoca(crx,cry,crz,cid2,
     &                       cfox,cfoy,cfoz,ncont,cpot)
c               rpie(7,ib) = cpot

*********************************************



                if(qtest.ne. ipi(2,ib)+id(2,jj)) then
                  write(*,*)'chrge prb in pionem ',ichannel
                  write(*,*)id(1,jj), qtest, ipi(2,ib),id(2,jj)
                  stop
                end if
                if(abs(ipi(2,ib)).gt.1 .and. ipi(1,ib).eq.1) then
           write(*,*)'chrge prb in pionem with pi charge',ichannel
                  write(*,*)id(1,jj), qtest, ipi(2,ib),id(2,jj)
                  stop
                end if
               if((ipi(1,ib).eq.2.or.ipi(1,ib).eq.4).and.ipi(2,ib).ne.0)
     &            write(*,*) 'hiba:pionem',ichannel,ipi(1,ib),ipi(1,ib)
     &            ,id(1,jj), idn, gam1pi, gameta, gam


                id(1,jj) = 1
                if(ichannel.eq.6) id(1,jj) = 2
                if(ichannel.eq.7) id(1,jj) = 3
                if(ichannel.eq.9) id(1,jj) = nres+2
                if(ichannel.eq.8) id(1,jj) = nres+3
                id(3,jj) = 0
                if(ichannel.ne.6 .and. ichannel.ne.7) id(4,jj)=0
                if(ichannel.ne.6 .and. ichannel.ne.7) id(5,jj)=0
                if(ichannel.ne.6 .and. ichannel.ne.7) id(7,jj)=0
                id(8,jj) = ib
c                write(*,*)'end of pionem particle id ', id(1,jj)
                e(jj)    = bmass
                upot(jj) = bapot
c                write(*,*) 'pionem upot',bmass,bapot
                pinabs = sqrt(p(1,jj)**2+p(2,jj)**2+p(3,jj)**2)

                stot   = (sqrt((e(jj)+ upot(jj))**2+pinabs**2)
     +                   +sqrt((epi(ib)+ mpot(ib))**2+ppi(1,ib)**2
     +                   +ppi(2,ib)**2+ppi(3,ib)**2))**2-
     +                   (p(1,jj)+ppi(1,ib))**2-(p(2,jj)+ppi(2,ib))**2
     +                   - (p(3,jj)+ppi(3,ib))**2

                if(abs(ema2-stot)/ema2  .gt.1.e-05 ) then
                write(79,*) 'hiba:pionem energy conservation'
                 write(79,*)'ire = ', ire, stot, ema2
                 write(79,*)'ichanel = ', ichannel
                end if
              endif
            endif
          endif
  800   continue
 1000 continue

      write(*,*)'end of pionem '
c     do i = 1,nres
c       write(*,*)i,lmesc(i)
c      write(*,*)(cres(i,j),j= 1,9)
c     end do


      write(*,*)'total decays : ', totdec, ichannel
c      call f77flush()
c      tout = secnds(t0)
c     tin = tout - tin
c     write(*,*)'time ellapsed in pionem = ',tin,'  sec.'
c      write(*,*)'ellapsed in pionem = mass ', mass
c      call f77flush()
      return
      end

