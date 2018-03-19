***********************************************************************
      subroutine crosw(pcm,srt,em1,em2,id1,id2,iblock,sig,
     &                 cutoff,prcm,i1,i2,iz1,iz2,
     &                 u1, u2,u3, u4, j01, j02,
     &                 betacm, betlrfx,betlrfy,betlrfz,testflag,nt,
     &                 sstate, ppi3,pnucl1,pnucl2,upi,izpi,rpi3,
     &                 rhap1,rhap2 )
*
*
*           n-n cross section ; cugnon parametrization
*           pcm(3)   - momentum coordinates of one particle in cm frame
*           srt      - sqrt of s
*           i1,i2    - identificator of particle 1 and 2
*           em1,em2  - rest masses of particle 1 and 2
*           id1,id2  - identificator of particle 1 and 2(n=1,d=2,nr=3)
*           iz1,iz2  - charges of particle 1 and 2
*           potndcm  - selfenergy correction in CMS
*           prcm     - absolute value of CMS-momentum
*           iseed    - seed for random number generator
*           iblock   - the information back
*           sig      - max cross section at cutoff energy
*           cutoff   - cutoff energy = em1 + em2 + 0.02 gev for nucl.
*                                      em1 + em2            for delt.
***********************************************************************
      implicit none

      real*8 gamr, bet2, qqr2
      real*8 m2,m3
c      parameter (m2=18.,m3=8.)!!!! richtig

******************************!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter(gamr = 0.11)
      parameter(bet2 = 0.090)
      parameter(qqr2 = 0.051936)
      include"common"
      include"com_pert"
      include"cominput"

      real*8    srt, em1, em2,  rn
      real*8    cutoff, prcm, sig
      integer id1, id2, iblock,  i1, i2, iz1, iz2
      real*8    pcm(3)
      real*8    u1, u2, u3, u4, j01, j02
      real*8    betacm(3), betlrfx(2), betlrfy(2), betlrfz(2)
      logical testflag, flagdd, ddcoll, momflag, sumtest,sstate
c      logical sflag
      integer nt
      real*8    rm2, pm2, c2, pr
      real*8    srtfree, cutfree, asrtfree, csrt
      integer idi, ida, indexs, indexslor
      real*8    pinitial, integral
      integer is, i, ilauf, iend, im, istart, j, indexm
      integer idres, icount, ifin
      real*8    resmasse, isofac, spinfac
      real*8    matritilde, fac
      real*8    pfinal,  xt
      integer totcharge, reslad, idhelp
      real*8    em3, em4
      real*8    rhap1(1:3), rhap2(1:3)

      real*8    sigel, signrnn, signnnr(1:nres),sigspion(2)
      real*8    signndd, sigddnn, signrnrp(1:nres), sumnnnr
      real*8    sumsstate, sumnrnrp, testsum, test, sigtot
      real*8    sigtt, sigh1, sigh2

      real*8    ppi3(1:3), pnucl1(1:3),pnucl2(1:3), upi, rpi3(1:3)
      real*8    mat
      real*8   srtin, help, matrix
      real*8    pfin(3), x1
      integer id3, id4, izpi, testmom, iwinkel
      real*8    finmul, matin, matout, help7,pnperpp0,pnperpp

      real*8 st

      integer iii
      integer id1inp, id2inp

      parameter (pnperpp0=1.3408)  ! new        after 09Nov2006
c      include"resdata1"

c      parameter (m2=18.,m3=8.)!!!! richtig
******************************!!!!!!!!!!!!!!!!!!!!!!!!!

c      write(*,*) 'crosw: srtfree = ',srtfree

      iblock = 0
      if(id1.eq.-1 .or. id2.eq.-1) return
      id2inp = max(id1,id2)

      m2 =  18.
      m3 = 8.
      help7 = 0.0

      if(inotwopi.eq.1) then
        m2 = 9.
      end if

      srtin = srt
      iblock = 0
      rm2 = rmass**2
      pm2 = pmass**2
      pr  = prcm
      c2  = pcm(3) / pr
      testflag =.false.
      sstate = .false.
      testmom = 0

*------- set sigma-variables to zero
      sigel    = 0.0
      signrnn  = 0.0
      do icount = 1, nres
        signnnr(icount) = 0.0
        signrnrp(icount)= 0.0
      end do
      do icount = 1,2
        sigspion(icount) = 0.0
      end do
      signndd  = 0.0
      sigddnn  = 0.0

*---------------------------------------------------------------------
*      define free quantities for calculation of free X-sections
*

      srtfree = sqrt(em1**2 + pr**2) + sqrt(em2**2 + pr**2)


      cutfree = srtfree - em1 - em2
      asrtfree= srtfree - em1 - em2
      indexs = nint( (srtfree-sigs0)/delsigs)
      indexslor = nint( (srtfree-sigs0)/delsiglor)
      if(indexs.gt.nsmax) then
c        write(*,*)'crosw vigyazz maxsrt index ',indexs,indexslor,srtfree
        indexs = nsmax - 1
      end if
      if(indexslor.gt.nsmaxlor) then
c        write(*,*)'crosw vigyazz max srt index 2 ', indexslor, srtfree
        indexslor = nsmaxlor - 1
      end if

      totcharge = iz1 + iz2
      ddcoll = .false.
      if(id1.eq.2 .and. id2.eq.2) ddcoll = .true.

      if((id1+id2).gt.2) then
        if(id1.gt.1) then
          resmasse = em1
          idres    = id1
          reslad   = iz1
        else
          resmasse = em2
          idres    = id2
          reslad   = iz2
        end if
        indexm = nint((resmasse-dimimass0)/delmassdi)

        if(indexm.lt.0) then
          write(*,*)'crosw falsche masse : ', indexm
          stop
        end if
        if(indexm.gt.nmassmax) then
          write(*,*)'warning!! crosw mass zu hoch : ', indexm
          indexm = nmassmax
        end if

      end if

*---------------------------------------------------------------------*
*      evaluate elastic X-sections or mass change X-sections          *

      idi=min(id1,id2)
      ida=max(id1,id2)

      csrt = srtfree - 2.0*rmass
      if(((id1.eq.2).and.(id2.eq.2)) .or.
     +              ((idi.eq.1).and.(ida.ge.2.and.ida.lt.nres+2))) then
***    NR-NR, ND-ND, DD-DD
        sigel = 35.0 / (1. + csrt * 100.0)  +  20.0


        pinitial = (srtfree**2+em1**2-em2**2)/(2.0*srtfree)
        pinitial = sqrt(pinitial**2 - em1**2)

        if(id1.eq.2 .and. id2.eq.2) then
          integral = ddmassmax(1,indexs)
        else
          integral = nrmassmax(1,ida-1,indexslor)
        end if

        if(pinitial.lt.1.0e-04) pinitial = 1.0e-04

        sigel = sigel * integral / pinitial


        else if((id1.eq.1.or.id1.ge.nres+2) .and.
     &          (id2.eq.1.or.id2.ge.nres+2) ) then
***       NN - NN, NL-NL, ...
        sigel = 35.0 / (1. + csrt * 100.0)  +  20.0
      end if
*
*---------------------------------------------------------------------*
*     is there enough energy available for inelastic scattering ?     *
*     if the available energy is less than the pion-mass, nothing     *
*     can happen any more ==> return (2.0146= 2*rmass +  pi-mass)     *
*     if there were 2 delta and they didn't scatt. elast., return     *



      if(((iz1+iz2) .eq. -2) .or. ((iz1+iz2) .eq. 4))       goto 77
      if(srtfree .le. 2.0146) then
        if(id1+id2.ne.2) then
          write(*,*)'prob in crosw 1 : ', id1, id2,srtfree
        end if
        goto 77
      end if
      if(ida.eq.nres+2 .or. ida.eq.nres+3)                       goto 77

*********** chek if NR - NN possible
      flagdd = .false.
      if(((iz1+iz2) .eq. -1) .or. ((iz1+iz2) .eq. 3)) then
        flagdd = .true.
      end if
**************************************************

*--------------------------------------------------------------------*
*                                                                    *
*          N + R -> N + N    LOOP                                    *

      if((id1+id2).ne.2 .and. idi .eq.1) then
************* N R - N N
         if( (id1+id2).eq.3 .and. .not.flagdd)then
****         N D - N N

*-----------Isospinfaktoren
*-----------(everything normalized to pp-nD++)
            if (iz1.eq.iz2) then
               isofac  =1./3.
               spinfac = 1.0
               fac     = isofac*spinfac !!! spinfac = identische teilchen
            else if((min0(iz1,iz2).eq.-1).or.max0(iz1,iz2).eq.2) then
               fac=1.
            else if(totcharge.eq.3 .or. totcharge.eq.-1) then
               fac = 0.
            else
               isofac = 2./3.
               spinfac= 1.0
               fac = isofac*spinfac
            end if
            signrnn = fac*dmdimi(2,indexs,indexm)
         else if((id1+id2).gt.3 .and.
     &           (id1.eq.1 .or. id2.eq.1).and. .not.flagdd)then
***   N R - N N (all res. except for D(1232))
            if(resprop2(idres-1,1).eq.3) then

               if (iz1.eq.iz2) then
                  isofac  =1./4.
                  spinfac = 0.5
                  fac     = isofac*spinfac
               else if((min0(iz1,iz2).eq.-1).or.max0(iz1,iz2).eq.2) then
                  fac=3./8.
               else
                  isofac = 1./4.
                  spinfac= 1.0
                  fac = isofac*spinfac
               end if


            else  if(resprop2(idres-1,1).eq.1) then
*&&&&&&&&&&&&& das 1535 muss verstaerkt kommen !!! im np stoss
               if (iz1.eq.iz2) then
                  isofac  = 1.0
                  spinfac = 0.5
                  fac     = isofac*spinfac
               else
                  isofac =  1.0
                  spinfac=  0.5
                  pnperpp=pnperpp0
*---------- special treatement for N*(1535)------------------------
                  if(idres .eq. 5) then
                     pnperpp=5.0
                     if((srtfree.le.2.51724).and.
     &                    (srtfree.gt.(2.0*rmass+emass)))
     &                    pnperpp=10.0
                  end if
                  fac = isofac*spinfac*pnperpp
               end if
            end if
            pinitial = (srtfree**2+em1**2-em2**2)/(2.0*srtfree)
            pinitial = sqrt(pinitial**2 - em1**2)
            pfinal   = (srtfree**2)/(2.0*srtfree)
            pfinal   = sqrt(pfinal**2 - rmass**2)
            finmul = 4.0
            matritilde = m2d16pi(idres-1,1)
            if(inotwopi.eq.1) then
               matritilde = m2d16pi(idres-1,2)
            end if

            if(pinitial.lt.1.0e-04) pinitial = 1.0e-04
            signrnn  = matritilde*finmul*pfinal/pinitial
     &           /srtfree**2*fac
         end if
      end if
*                                                                    *
*         END OF N + R - > N + N   LOOP                              *
*--------------------------------------------------------------------*
*         BEGIN OF  NN, NR', DD - > N + R  LOOP                      *
*                                                                    *
*        check for security
      if(iblock .ne. 0) then
        write(*,*)'program terminated in crosw: IBLOCK = ',iblock
        stop
      end if


      if(irescol.eq.0 .and. ida.gt.1)                   goto 77
        if(ida.eq.1 ) then
****      N N - N R (all res, incl. D(1232))

              if (iz1.eq.iz2) then
                 fac = 1.0 + 1./3.
              else
                 fac = 2./3.
              end if
             signnnr(1) = dimixsec(indexs)*fac
          icount = 0
          do icount = 2, nres

            if(resprop2(icount,1).eq.3) then

              if (iz1.eq.iz2) then
                 fac = 1.0
              else
                 fac = 0.5
              end if

            else  if(resprop2(icount,1).eq.1) then

              if (iz1.eq.iz2) then
                fac     = 1.0
              else
                fac     = pnperpp0
              end if

            end if
            if(icount.eq.4.and.totcharge.eq.1)then
             fac = fac*5.0
             if(srtfree.le.2.51724.and.srtfree.gt.(2.0*rmass+emass))
     &                 fac = fac*2.0
            end if

           pinitial = (srtfree**2+em1**2-em2**2)/(2.0*srtfree)
           pinitial = sqrt(pinitial**2 - em1**2)

           finmul = 2.0*(resprop2(icount,3)+1.0)

           mat = m2d16pi(icount,1) * finmul
           if(inotwopi.eq.1) then
             mat = m2d16pi(icount,2) * finmul
           end if


            if(pinitial.lt.1.0e-04) pinitial = 1.0e-04
            signnnr(icount) = nrmassmax(1,icount,indexslor)*fac
     &                        /pinitial/srtfree**2*mat
          end do

c          write(*,*) 'NN->NR cross sections'
c          do icount = 1,nres
c             write(*,*) icount, signnnr(icount)
c          end do

************* s- state pions
          call spion(srtfree,sigspion,iz1,iz2)
******    NN - DD
          if(inotwopi.eq.0) then
          call sddaich(srtfree,iz1,iz2, signndd,1,em1,em2,sigh1,sigh2)
          end if
        else if(idi.eq.1 .and. ida.gt.1) then
****     N R - N R'
****     with the famous isospin correction factors.
          do icount=1,nres
            if(resprop2(idres-1,1).eq.1) then
              if(resprop2(icount,1).eq.1) then
                  isofac = 1.
              else
                if (totcharge.eq.1) then
                  isofac = 1./2.
                else
                  isofac = 1.
                end if
              end if
            else  if(resprop2(idres-1,1).eq.3) then
              if (resprop2(icount,1).eq.1) then
                if(totcharge.gt.2) then
                  isofac = 0.0
                else if(totcharge.le.-1) then
                  isofac = 0.0
                else if(totcharge.eq.1) then
                  isofac  = 1./2.
                else if((reslad.eq.2).or.(reslad.eq.-1)) then
                  isofac = 3./4.
                else
                  isofac = 1./4.
                end if
              else
                isofac =1.
              end if
            end if
c             if ((idres-1.eq.1).or.(icount.eq.1))isofac =8.*isofac
             pinitial = (srtfree**2+em1**2-em2**2)/(2.0*srtfree)
             pinitial = sqrt(pinitial**2 - em1**2)
*********
********** achtung im integeral muss faktor pi/2 verarbeitet sein.

            integral = nrmassmax(1,icount,indexslor)

            finmul = 2.0*(resprop2(icount,3)+1.0)


********** berechne das matrix-element fuer die res. im initial channel
            if(idres.eq.5) then

              m3 = m2d16pi(4,1)
              if(inotwopi.eq.1) then
                m3 = m2d16pi(4,2)
              end if

              if(totcharge.eq.2.or.totcharge.eq.0) then
                matin = m3
              else if(totcharge.eq.1) then
                matin = 5.0*m3
                if(srtfree.le.2.51724.and.
     &             srtfree.gt.(2.0*rmass+emass)) then
                  matin = 2.0*matin
                end if
              else
                write(*,*)'crosw matrix warning !!!!!!!!!!'
                matin = m3
              end if

            else
              matin = m2d16pi(idres-1,1)
              if(inotwopi.eq.1) then
                matin = m2d16pi(idres-1,2)
              end if
              if(totcharge.eq.1) matin = matin*pnperpp0
            end if



********** berechne das matrix-element fuer die res. im final channel
            if(icount.eq.4) then

              m3 = m2d16pi(4,1)
              if(inotwopi.eq.1) then
                m3 = m2d16pi(4,2)
              end if

              if(totcharge.eq.2.or.totcharge.eq.0) then
                matout = m3
              else if(totcharge.eq.1) then
                matout = 5.0*m3
                if(srtfree.le.2.51724.and.
     &             srtfree.gt.(2.0*rmass+emass)) then
                  matout = 2.0*matout
                end if
              else
                matout = m3 !!!!!! nur damit was hier steht
              end if

            else
              matout = m2d16pi(icount,1)
              if(inotwopi.eq.1) then
                matout = m2d16pi(icount,2)
              end if
              if(totcharge.eq.1) matout = matout*pnperpp0
            end if

            matrix = 0.5*(matin+matout)


            signrnrp(icount)=
     &        matrix/pinitial/srtfree**2*isofac*integral*finmul

          end do
*         nicht nr - nr ! schon bei elast abgefangen
          signrnrp(idres-1) = 0.0

c zm: only delta and n1520 exists!
c          signrnrp(2) = 0.
c          do iii = 4,nres
c            signrnrp(iii) = 0.
c          end do

         else if(ddcoll .and. inotwopi.eq.0) then
****     DD - NN
          call sddaich(srtfree,iz1,iz2, sigddnn,2,em1,em2,sigh1,sigh2)
        end if


*-------------------------------------------------------------------*
*       sumation of cross sections needed for MC-decisions          *

 77      continue


      sumnrnrp = 0.0
      sumnnnr  = 0.0
      sumsstate= 0.0

      do icount =  1, nres
        sumnrnrp = sumnrnrp + signrnrp(icount)
        sumnnnr  = sumnnnr  + signnnr(icount)
      end do
      do icount = 1,2
        sumsstate = sumsstate + sigspion(icount)
      end do

      if(isnndd.eq.1) signndd = 0.0
      if(isddnn.eq.1) sigddnn = 0.0
      if(isnrnrp.eq.1) sumnrnrp = 0.0
      if(isstatep.eq.1) sumsstate = 0.0
      if(isdeltaonly.eq.1) then
        sumnnnr = 0.0
        sumnnnr =  signnnr(1)
      end if
      if(isdeltaall.eq.1) then
        signnnr(1)=  sumnnnr
        sumnnnr   =  signnnr(1)
      end if


      sigtot = sigel+signrnn+sigddnn+signndd+sumnrnrp+sumnnnr+
     &         sumsstate



c      if(sigtot.gt.sig)
*     &     write(*,*) 'vigyazz sigtot',sigtot,sigdn+sigrn+sigqn
*     &     ,signd+signr+signq,sigdd,sigddnd
      sigtt = max(sig, sigtot)
      x1    = rn(iseed)




      if(x1.le.sigtot/sigtt) testmom = 1

c      write(*,*)'sig = ', sig
c      write(*,*)'sigtot = ', sigtot
c      write(*,*)'sigtt = ', sigtt
c      write(*,*)'sigel = ', sigel
c      write(*,*)'signrnn = ', signrnn
c      write(*,*)'sigddnn = ', sigddnn
c      write(*,*)'signndd = ', signndd
c      write(*,*)'sumnrnrp ', sumnrnrp
c      write(*,*)'sumnnnr ', sumnnnr
c      write(*,*)'sumsstate =', sumsstate , x1, testflag



      if(x1 .le. sigel/sigtt) then
*--------------------- elastic scattering

        id3 = id1
        id4 = id2

        em3 = em1
        em4 = em2
        call momkin(pcm,pr,srt,em1,em2,id3,id4,j01,j02,
     &              betacm,betlrfx,betlrfy,betlrfz,momflag,nt,
     &              pfin,em3,em4,u3,u4,1,srtfree,help7,rhap1,rhap2)

        if((id3+id4).eq.2) then
          iblock = 1
        else
          iblock = idres
        end if

c        write(*,*) 'id1,id2,id3,id4,iz1,iz2',id1,id2,id3,id4,iz1,iz2

        if(momflag) then

          if((id3+id4).eq.2) then   ! N + N
            iblock = 1
          else if(id3.eq.nres+2 .or. id4.eq.nres+2) then ! L + ...
            iblock = 3*nres+4+nres**2
          else if(id3.eq.nres+3 .or. id4.eq.nres+3) then ! S + ...
            iblock = 3*nres+5+nres**2
          end if

          if((id3+id4).ne.2) then   ! other than N + N
            if(id1.eq.id4) then
              idhelp = iz2
              iz2    = iz1
              iz1    = idhelp
            end if
          end if
          testflag =.true.
        end if

      else if(x1.le.( (sigel+signrnn)/sigtt)) then
*---------------------N R - NN
        iblock = nres+1 + (idres-1) ! max iblock up to now 2nres+1
        id3 = 1
        id4 = 1
        if(idres.eq.2) then
          iwinkel = 2
        else
          iwinkel = 0
        end if
        call momkin(pcm,pr,srt,em1,em2,id3,id4,j01,j02,
     &              betacm,betlrfx,betlrfy,betlrfz,momflag,nt,
     &   pfin,em3,em4,u3,u4,iwinkel,srtfree,resmasse,rhap1,rhap2)

        if(momflag) then
          testflag = .true.
          iblock = nres+1 + (idres-1)   ! max iblock up to now 2nres+1

*             assign charges
          if(totcharge.eq.2) then
            iz1 = 1
            iz2 = 1
          else if(totcharge.eq.0) then
            iz1 = 0
            iz2 = 0
          else if(totcharge.eq.1) then
            xt = rn(iseed)
            if(xt.le.0.5) then
              iz1 = 0
              iz2 = 1
            else
              iz1 = 1
              iz2 = 0
            end if
          else
            write(*,*)'warning crosw A : wrong qtot: '
            write(*,*)id1, id2, iz1, iz2
            write(*,*)id3, id4, totcharge
            stop
          end if
        else
          signrnn = 0.0
          goto 77
        end if
      else if(x1.le.( (sigel+signrnn+sumnnnr)/sigtt)) then
*---------------------N N - N R
*        chose the specific channel
        testsum = 0.0
        sumtest = .true.
        test    = rn(iseed)
        icount = 0
        do while(sumtest)
          icount = icount + 1
          testsum = testsum + signnnr(icount)
          if(test.le.(testsum/sumnnnr)) then
            sumtest = .false.
            ifin    = icount
          end if
          if(icount .gt. nres) then
            write(*,*)'hiba, error from crosw B '
            write(*,*)'icount = ', icount
            stop
          end if
        end do

        iblock = 2*nres+1+ ifin        !max iblock 3*nres+1
        id3 = 1
        id4 = ifin + 1
        if(ifin.eq.1) then
          iwinkel = 2
        else
          iwinkel = 0
        end if
        call momkin(pcm,pr,srt,em1,em2,id3,id4,j01,j02,
     &              betacm,betlrfx,betlrfy,betlrfz,momflag,nt,
     &             pfin,em3,em4,u3,u4,iwinkel,srtfree,help7,rhap1,rhap2)
        if(momflag) then
          testflag = .true.
          iblock = 2*nres+1+ ifin        !max iblock 3*nres+1
**        charge assignment
          test = rn(iseed)
          if(resprop2(ifin,1).eq.3) then

            if(totcharge.eq.2) then
              if(test.le. 3./4.) then
                if(id3.eq.1) then
                  iz1 =  0
                  iz2 =  2
                else if(id4.eq.1) then
                  iz1 =  2
                  iz2 =  0
                end if
              else
                iz1 =  1
                iz2 =  1
              end if
            else if(totcharge.eq.0) then
              if(test.le. 3./4.) then
                if(id3.eq.1) then
                  iz1 =  1
                  iz2 = -1
                else if(id4.eq.1) then
                  iz1 = -1
                  iz2 =  1
                end if
              else
                  iz1 =  0
                  iz2 =  0
              end if
            else if(totcharge.eq.1) then
              if(test.le. 0.5) then
                if(id3.eq.1) then
                  iz1 =  1
                  iz2 =  0
                else if(id4.eq.1) then
                  iz1 =  0
                  iz2 =  1
                end if
              else
                if(id3.eq.1) then
                  iz1 =  0
                  iz2 =  1
                else if(id4.eq.1) then
                  iz1 =  1
                  iz2 =  0
                end if
              end if
            end if
          else if(resprop2(ifin,1).eq.1) then
            if(totcharge.eq.2) then
              iz1 = 1
              iz2 = 1
            else if(totcharge.eq.0) then
              iz1 = 0
              iz2 = 0
            else if(totcharge.eq.1) then
              if(test.le.0.5) then
                iz1 = 0
                iz2 = 1
              else
                iz1 = 1
                iz2 = 0
              end if
            end if
          end if
        else
          signnnr(ifin)= 0.0
          goto 77
*** end of momflag
        end if
      else if(x1.le.((sigel+signrnn+sumnnnr+sumnrnrp)/sigtt))then
*---------------------N R - N R'

*        chose the specific channel
        testsum = 0.0
        sumtest = .true.
        test    = rn(iseed)
        icount = 0
c        write(*,*) 'crosw 769: ',signrnrp(1),signrnrp(3),sumnrnrp,test
        do while(sumtest)
          icount = icount + 1
          testsum = testsum + signrnrp(icount)
          if(test.le.(testsum/sumnrnrp)) then
            sumtest = .false.
            ifin    = icount
          end if
          if(icount .gt. nres) then
            write(*,*)'warning from crosw c '
            write(*,*)'icount = ', icount
            stop
          end if
        end do
c        if (ifin.lt.3 .and. id2inp.eq.4) write(*,*)
c     &     "N1520 lost in collision"

        id3 = 1
        id4 = ifin + 1
        iblock = 3*nres+1 +(idres-2)*nres+ifin !max ibl=3nres+1+nres**2
        call momkin(pcm,pr,srt,em1,em2,id3,id4,j01,j02,
     &              betacm,betlrfx,betlrfy,betlrfz,momflag,nt,
     &              pfin,em3,em4,u3,u4,0,srtfree,help7,rhap1,rhap2)
        if(momflag) then
          testflag =.true.
          iblock = 3*nres+1 +(idres-2)*nres+ifin !max ibl=3nres+1+nres**2
          if( (id3+id4).eq.(id1+id2)) then
            write(*,*)'crosw nr-nrp'
            write(*,*)'id1, id2 ', id1, id2
            write(*,*)'id3, id4 ', id3, id4
            stop
          end if

          if(resprop2(idres-1,1).eq.1) then
            if(resprop2(ifin,1).eq.1) then
              if(totcharge.eq.2) then
                iz1 = 1
                iz2 = 1
              else if(totcharge.eq.0) then
                iz1 = 0
                iz2 = 0
              else if(totcharge.eq.1) then
                if(test.le.0.5) then
                  iz1 = 0
                  iz2 = 1
                else
                  iz1 = 1
                  iz2 = 0
                end if
              end if
            else if(resprop2(ifin,1).eq.3) then
              if(totcharge.eq.2) then
                if(test.le. 3./4.) then
                  if(id3.eq.1) then
                    iz1 =  0
                    iz2 =  2
                  else if(id4.eq.1) then
                    iz1 =  2
                    iz2 =  0
                  end if
                else
                  iz1 =  1
                  iz2 =  1
                end if
              else if(totcharge.eq.0) then
                if(test.le. 3./4.) then
                  if(id3.eq.1) then
                    iz1 =  1
                    iz2 = -1
                  else if(id4.eq.1) then
                    iz1 = -1
                    iz2 =  1
                  end if
                else
                   iz1 =  0
                   iz2 =  0
                end if
              else if(totcharge.eq.1) then
                if(test.le. 0.5) then
                  if(id3.eq.1) then
                    iz1 =  1
                    iz2 =  0
                  else if(id4.eq.1) then
                    iz1 =  0
                    iz2 =  1
                  end if
                else
                  if(id3.eq.1) then
                    iz1 =  0
                    iz2 =  1
                  else if(id4.eq.1) then
                    iz1 =  1
                    iz2 =  0
                  end if
                end if
              end if
            end if
          else if(resprop2(idres-1,1).eq.3) then
            if(resprop2(ifin,1).eq.1) then
              if(totcharge.eq.2) then
                iz1 =  1
                iz2 =  1
              else if(totcharge.eq.0) then
                iz1 =  0
                iz2 =  0
              else if(totcharge.eq.1) then
                if(test.le. 0.5) then
                  if(id3.eq.1) then
                    iz1 =  1
                    iz2 =  0
                  else if(id4.eq.1) then
                    iz1 =  0
                    iz2 =  1
                  end if
                else
                  if(id3.eq.1) then
                    iz1 =  0
                    iz2 =  1
                  else if(id4.eq.1) then
                    iz1 =  1
                    iz2 =  0
                  end if
                end if
              end if
            else if(resprop2(ifin,1).eq.3) then
              if(totcharge.eq.3) then
                if(id3.eq.1) then
                  iz1 =  1
                  iz2 =  2
                else if(id4.eq.1) then
                  iz1 =  2
                  iz2 =  1
                end if
              else if(totcharge.eq.2) then
                if(reslad.eq.2) then
                  if(test.le.5./8.) then
                    if(id3.eq.1) then
                      iz1 =  0
                      iz2 =  2
                    else if(id4.eq.1) then
                      iz1 =  2
                      iz2 =  0
                    end if
                  else
                    iz1 = 1
                    iz2 = 1
                  end if
                else if(reslad.eq.1) then
                  if(test.le.5./8.) then
                      iz1 =  1
                      iz2 =  1
                  else
                    if(id3.eq.1) then
                      iz1 = 0
                      iz2 = 2
                    else if(id4.eq.1) then
                      iz1 =  2
                      iz2 =  0
                    end if
                  end if
                end if
              else if(totcharge.eq.0) then
                if(reslad.eq.0) then
                  if(test.le.5./8.) then
                    iz1 = 0
                    iz2 = 0
                  else
                    if(id3.eq.1)then
                      iz1 =  1
                      iz2 = -1
                    else if(id4.eq.1)then
                      iz1 = -1
                      iz2 =  1
                    end if
                  end if
                else if(reslad.eq.-1)then
                  if(test.le.5./8.) then
                    if(id3.eq.1)then
                      iz1 =  1
                      iz2 = -1
                    else if(id4.eq.1)then
                      iz1 = -1
                      iz2 =  1
                    end if
                  else
                    iz1 = 0
                    iz2 = 0
                  end if
                else
                  write(*,*)'waning crosw D ', iz1, iz2
                end if
              else if(totcharge.eq.1) then
                if(test.le.0.5) then
                  iz1 = 0
                  iz2 = 1
                else
                  iz1 = 1
                  iz2 = 0
                end if
              else if(totcharge.eq.-1) then
                if(id3.eq.1)then
                  iz1 =  0
                  iz2 = -1
                else if(id4.eq.1)then
                  iz1 = -1
                  iz2 =  0
                end if
              end if
            end if
          end if
        else
          signrnrp(ifin)= 0.0
          goto 77
******* end of momflag if-statement
        end if
      else if(x1.le.
     &      ((sigel+signrnn+sumnnnr+sumnrnrp+signndd)/sigtt))then
**------------- N N - D D

        id3 = 2
        id4 = 2
          iblock = 3*nres +1 +nres**2  + 1

      call momkin(pcm,pr,srt,em1,em2,id3,id4,j01,j02,
     &              betacm,betlrfx,betlrfy,betlrfz,momflag,nt,
     &              pfin,em3,em4,u3,u4,0,srtfree,help7,rhap1,rhap2)

        if(momflag) then
          testflag = .true.
          iblock = 3*nres +1 +nres**2  + 1
          call sddaich(srtfree,iz1,iz2, sigddnn,1,em1,em2,sigh1,sigh2)
          test = rn(iseed)

*            charge assignment

*           nucleon- nucleon-collision
          if(totcharge.eq.2)then
*              pp-collision
            if(test.le. sigh1/sigddnn) then
              xt   = rn(iseed)
              if(xt.le.0.5) then
                iz1 =  0
                iz2 = +2
              else
                iz2 =  0
                iz1 = +2
              end if
            else
                iz1 = 1
                iz2 = 1
            end if
          else if(totcharge.eq.0) then
*              nn-collision
            if(test.le. sigh1/sigddnn) then
              xt   = rn(iseed)
              if(xt.le.0.5) then
                iz1 =  1
                iz2 = -1
              else
                iz2 = -1
                iz1 =  1
              end if
            else
                iz1 = 0
                iz2 = 0
            end if
          else if(totcharge.eq.1) then
*              pn-collision
            if(test.le. sigh1/sigddnn) then
              xt   = rn(iseed)
              if(xt.le.0.5) then
                iz1 =  1
                iz2 =  0
              else
                iz2 =  0
                iz1 =  1
              end if
            else if(test .le. (sigh1+sigh2)/sigddnn) then
              xt   = rn(iseed)
              if(xt.le.0.5) then
                iz1 = +2
                iz2 = -1
              else
                iz2 = -1
                iz1 = +2
              end if
            else
              write(*,*)'something wrong with the NN-DD channel'
              write(*,*)id1, id2, iz1,iz2
              stop
            end if
          end if
        else
         signndd = 0.0
         goto 77
*        end of momflag
        end if
      else if(x1.le.
     & ((sigel+signrnn+sumnnnr+sumnrnrp+signndd+sigddnn))/sigtt)then
**------------- D D - N N

        id3 = 1
        id4 = 1
        iblock = 3*nres +1 +nres**2  + 2
        call momkin(pcm,pr,srt,em1,em2,id3,id4,j01,j02,
     &              betacm,betlrfx,betlrfy,betlrfz,momflag,nt,
     &              pfin,em3,em4,u3,u4,0,srtfree,help7,rhap1,rhap2)

        if(momflag) then
          testflag = .true.
            iblock = 3*nres +1 +nres**2  + 2

*           nucleon- nucleon-collision
          if(totcharge.eq.2)then
            iz1 = 1
            iz2 = 1
          else if(totcharge.eq.0) then
            iz1 = 0
            iz2 = 0
          else if(totcharge.eq.1) then
            xt   = rn(iseed)
            if(xt.le.0.5) then
              iz1 =  1
              iz2 =  0
            else
              iz2 =  0
              iz1 =  1
            end if
          else
              write(*,*)'something wrong with the NN-DD channel'
              write(*,*)id1, id2, iz1,iz2
          end if
        else
          sigddnn = 0.0
          goto 77
*     end of momflag
        end if

      else if(x1.le. ((sigel+signrnn+sumnnnr+sumnrnrp+
     &                 signndd+sigddnn+sumsstate))/sigtt)then
**-------------N N - N N pi (s-state scattering)

        id3 = 1
        id4 = 1
            iblock = 3*nres +1 +nres**2  + 3
c            sflag =.false.
c          if(nt.eq.1000 .and. i2.ge.0) then
c            write(*,*)'cw 60'
c            sflag = .true.
c          end if
        momflag =.false.
        call sstatekin(srt,u3, u4, upi, ppi3,pnucl1,pnucl2,momflag,
     +               i1,i2,rpi3,betacm,srtfree,rhap1,rhap2)

        if(momflag) then
          testflag = .true.
          em3 = rmass
          em4 = rmass
            iblock = 3*nres +1 +nres**2  + 3
          sstate = .true.

*         charge assignment
          test = rn(iseed)
          if(totcharge.eq.1)then
            if(test.le. sigspion(1)/sumsstate) then
*            create charge pion
              xt = rn(iseed)
              if(xt.le.0.5) then
                iz1  =  1
                iz2  =  1
                izpi = -1
              else
                iz1  =  0
                iz2  =  0
                izpi =  1
              end if
            else
*             neutral pion
               xt = rn(iseed)
               izpi = 0
               if(xt.le.0.5) then
                 iz1 = 1
                 iz2 = 0
               else
                 iz1 = 0
                 iz2 = 1
               end if
            end if
          else if(totcharge.eq.2) then
            if(test.le. sigspion(1)/sumsstate) then
*            create charged pion
              izpi = 1
              xt = rn(iseed)
               if(xt.le.0.5) then
                 iz1 = 1
                 iz2 = 0
               else
                 iz1 = 0
                 iz2 = 1
               end if
            else
              izpi = 0
              iz1  = 1
              iz2  = 1
            end if
          else if(totcharge.eq.0) then
            if(test.le. sigspion(1)/sumsstate) then
*            create charged pion
              izpi = -1
              xt = rn(iseed)
               if(xt.le.0.5) then
                 iz1 = 1
                 iz2 = 0
               else
                 iz1 = 0
                 iz2 = 1
               end if
            else
              izpi = 0
              iz1  = 0
              iz2  = 0
            end if
          end if

        else
          sumsstate = 0.0
          goto 77
        end if
*         end of sigma loop
      end if

      if( .not.testflag) then
        iblock = 0 !!!!!!!!!!!!!!!!! change upper stuff later!!!!!
      end if


      if(.not.sstate .and. testflag) then
        do i = 1,3
          ppi3(i) = 0.0
        end do
        xt = rn(iseed)
        if(xt.le.0.5) then
          do i = 1, 3
            pnucl1(i) =  pfin(i)
            pnucl2(i) = -pfin(i)
            pcm(i)    =  pfin(i)
          end do
        else
          do i = 1, 3
            pnucl2(i) =  pfin(i)
            pnucl1(i) = -pfin(i)
            pcm(i)    =  pfin(i)
          end do
        end if
      end if
      if(testflag)then
        em1 = em3
        em2 = em4
        id1 = id3
        id2 = id4

        if(.not.sstate) then

          if(em1.le.0.1 .or. em2.le.0.1) then
            write(*,*)'incrosw falsche masse '
            write(*,*)'em1, em2 ', em1, em2
            write(*,*)'iblock = ', iblock
          end if


        help= (em1+u3)**2+pnucl1(1)**2+pnucl1(2)**2+pnucl1(3)**2
        help = sqrt(help)+sqrt((em2+u4)**2+pnucl2(1)**2
     &                     +pnucl2(2)**2+pnucl2(3)**2)
        if(abs(help-srt).gt.1.0e-05) then
        write(*,*)'am ende von crosw '
        write(*,*)'em1 = ', em1
        write(*,*)'em2 = ', em2
        write(*,*)'u3  = ', u3
        write(*,*)'u4  = ', u4
        write(*,*)'pnucl1 = ', pnucl1
        write(*,*)'pnucl2 = ', pnucl2
        write(*,*)'srt = ', srt, srtin, iblock , sstate,help
        stop
        end if

        else if(sstate) then
c          write(*,*)'sstate reaction'

          help= (em1+u3)**2+pnucl1(1)**2+pnucl1(2)**2+pnucl1(3)**2
          help = sqrt(help)+sqrt((em2+u4)**2+pnucl2(1)**2
     &                     +pnucl2(2)**2+pnucl2(3)**2)

          help = help +
     &         sqrt((pmass+upi)**2 +ppi3(1)**2+ppi3(2)**2+ppi3(3)**2)

        help = help**2 - (pnucl1(1)+pnucl2(1)+ppi3(1))**2
     &                 - (pnucl1(2)+pnucl2(2)+ppi3(2))**2
     &                 - (pnucl1(3)+pnucl2(3)+ppi3(3))**2
        help = sqrt(help)

        if(abs(srt-help).gt.1.0e-05) then
        write(*,*)'am ende von crosw '
        write(*,*)'em1 = ', em1
        write(*,*)'em2 = ', em2
        write(*,*)'u3  = ', u3
        write(*,*)'u4  = ', u4
        write(*,*)'upi  = ', upi
        write(*,*)'pnucl1 = ', pnucl1
        write(*,*)'pnucl2 = ', pnucl2
        write(*,*)'ppi3   = ', ppi3
        write(*,*)'srt = ', srt, srtin, iblock , sstate,help
             stop
        end if
        end if

      end if
c      write(*,*)   ' crosw  - end '
      return


************************************************************************
      entry croswread
*
*      reads the cross sections needed for resonance population
*      from file

 10   format(5(1x,e16.8))
******** read Dimitriev dsigma/dm und tot sigma for ND-NN
      open(mdimiread,file='buuinput/DIMIDM.dat',status='old')

      do is = 0, nsmax
        do im = 0, (nmassmax)/5
          istart = im*5
          iend   = istart + 4
          read(mdimiread,10) (dmdimi(1,is,ilauf),ilauf = istart,iend)
        end do
      end do

      do is = 0, nsmax
        do im = 0, (nmassmax)/5
          istart = im*5
          iend   = istart + 4
          read(mdimiread,10) (dmdimi(2,is,ilauf),ilauf = istart,iend)
        end do
      end do

      close(mdimiread)


***** read the total cross sections NN - ND as a function
***** of sqrt(s)
c
      open(mxresread,file='buuinput/DIMIXSEC.dat',status='old')

        do is = 0, (nsmax)/5
          istart = is*5
          iend   = istart + 4
          read(mxresread,10) (dimixsec(ilauf),ilauf = istart,iend)
        end do
      close(mxresread)

****************** integral over mass dist. and corresp max
c   OLD VERSION:
c      open(mmaxresread,file='buuinput/NRINT.dat',status='old')
c      do j = 1,2
c        do i = 1 , nres
c          do is = 0, (nsmax)/5
c            istart = is*5
c            iend   = istart + 4
c         read(mmaxresread,10)(nrmassmax(j,i,ilauf),ilauf = istart,iend)
c          end do
c        end do
c      end do
c      close(mmaxresread)

****************** integral over mass dist. and corresp max
c   NEW VERSION:
      write(*,*) 'before reading lorinteg.dat'
      open(mmaxresread,file='buuinput/lorinteg.dat',status='old')
      read(mmaxresread,'(///)')
      do is = 1, nsmaxlor
         read(mmaxresread,'(f6.3,30e10.3)')
     $        st,(nrmassmax(1,i,is),i=1,nres)
      end do
      close(mmaxresread)
      write(*,*) 'after reading lorinteg',(nrmassmax(1,i,1),i=1,nres)
      open(mmaxresread,file='buuinput/lormax.dat',status='old')
      read(mmaxresread,'(///)')
      do is = 1, nsmaxlor
         read(mmaxresread,*) st,(nrmassmax(2,i,is),i=1,nres)
      end do
      close(mmaxresread)
      write(*,*) 'after reading lormax.dat'

****************** integral over mass dist. and corresp max for D(1232)
      open(mmaxddread,file='buuinput/DDMASSV.dat',status='old')
      do i = 1,2
        do is = 0, (nsmax)/5
          istart = is*5
          iend   = istart + 4
          read(mmaxddread,10) (ddmassmax(i,ilauf),ilauf = istart,iend)
        end do
      end do
      close(mmaxddread)


      open(msstread,file='buuinput/SSTATE.dat',status='old')
      do i = 1, 3
        do is = 0, (nsmax)/5
          istart = is*5
          iend   = istart + 4
          read(msstread,10) (matsstate(i,ilauf),ilauf = istart,iend)
        end do
      end do
      close(msstread)
c
c      open(mxddread,file='buuinput/DDXSEC.dat',status='old')
c      do i = 1, 3
c        do is = 0, (nsmax)/5
c          istart = is*5
c          iend   = istart + 4
c          read(mxddread,10) (secdd(i,ilauf),ilauf = istart,iend)
c        end do
c      end do
c      close(mxddread)

c------------------------------------------------------------------
c----------------------by Peter Kovacs (2017)---------------------

c      write(*,*) 'before reading p_pbar_to_p_pbar_JPsi.dat'
c      open(mmaxresread,
c     $     file='buuinput/p_pbar_to_p_pbar_JPsi.dat',status='old')
c      read(mmaxresread,'(///)')
c      read(mmaxresread,*) NoL_1
c      do is = 1, NoL_1
c         read(mmaxresread,*)       
c     $        sig_ppbar_ppbar_JPsi(is,1), sig_ppbar_ppbar_JPsi(is,2) 
c      end do
c      close(mmaxresread)
c      write(*,*) 'after reading p_pbar_to_p_pbar_JPsi', 
c     $     sig_ppbar_ppbar_JPsi(1,1), sig_ppbar_ppbar_JPsi(1,2)
c      write(*,*) sig_ppbar_ppbar_JPsi(NoL_1,1),
c     $     sig_ppbar_ppbar_JPsi(NoL_1,2)
      
c------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_JPsi_pi0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__JPsi_pi0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_2
      do is = 1, NoL_2
         read(mmaxresread,*)       
     $  sig_ppbar_JPsi_pi0(is,1), sig_ppbar_JPsi_pi0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_JPsi_pi0.dat', 
     $  sig_ppbar_JPsi_pi0(1,1), sig_ppbar_JPsi_pi0(1,2)
      write(*,*) sig_ppbar_JPsi_pi0(NoL_2,1),sig_ppbar_JPsi_pi0(NoL_2,2)

c------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_JPsi_rho0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__JPsi_rho0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_3
      do is = 1, NoL_3
         read(mmaxresread,*)       
     $  sig_ppbar_JPsi_rho0(is,1), sig_ppbar_JPsi_rho0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_JPsi_rho0.dat', 
     $  sig_ppbar_JPsi_rho0(1,1), sig_ppbar_JPsi_rho0(1,2)
      write(*,*)sig_ppbar_JPsi_rho0(NoL_3,1),
     &   sig_ppbar_JPsi_rho0(NoL_3,2)

c------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_Psi1_pi0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__Psi1_pi0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_4
      do is = 1, NoL_4
         read(mmaxresread,*)       
     $  sig_ppbar_Psi1_pi0(is,1), sig_ppbar_Psi1_pi0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_Psi1_pi0.dat', 
     $  sig_ppbar_Psi1_pi0(1,1), sig_ppbar_Psi1_pi0(1,2)
      write(*,*) sig_ppbar_Psi1_pi0(NoL_4,1),sig_ppbar_Psi1_pi0(NoL_4,2)

c------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_Psi1_rho0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__Psi1_rho0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_5
      do is = 1, NoL_5
         read(mmaxresread,*)       
     $  sig_ppbar_Psi1_rho0(is,1), sig_ppbar_Psi1_rho0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_Psi1_pi0.dat', 
     $  sig_ppbar_Psi1_pi0(1,1), sig_ppbar_Psi1_pi0(1,2)
      write(*,*) sig_ppbar_Psi1_pi0(NoL_5,1),sig_ppbar_Psi1_pi0(NoL_5,2)

c------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_Psi2_pi0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__Psi2_pi0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_6
      do is = 1, NoL_6
         read(mmaxresread,*)       
     $  sig_ppbar_Psi2_pi0(is,1), sig_ppbar_Psi2_pi0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_Psi2_pi0.dat', 
     $     sig_ppbar_Psi2_pi0(1,1), sig_ppbar_Psi2_pi0(1,2)
      write(*,*) sig_ppbar_Psi2_pi0(NoL_6,1),sig_ppbar_Psi2_pi0(NoL_6,2)

c-----------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_Psi2_rho0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__Psi2_rho0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_7
      do is = 1, NoL_7
        read(mmaxresread,*)       
     $       sig_ppbar_Psi2_rho0(is,1), sig_ppbar_Psi2_rho0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_Psi2_pi0.dat', 
     $  sig_ppbar_Psi2_pi0(1,1), sig_ppbar_Psi2_pi0(1,2)
      write(*,*) sig_ppbar_Psi2_pi0(NoL_7,1),sig_ppbar_Psi2_pi0(NoL_7,2)

c-----------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_Dp_Dm.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__Dp_Dm.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_DD
      do is = 1, NoL_DD
        read(mmaxresread,*)       
     $       sig_ppbar_DpDm(is,1), sig_ppbar_DpDm(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_DpDm.dat', 
     $  sig_ppbar_DpDm(1,1), sig_ppbar_DpDm(1,2)
      write(*,*) sig_ppbar_DpDm(NoL_DD,1),sig_ppbar_DpDm(NoL_DD,2)
c-----------------------------------------------------------------------
      write(*,*) 'before reading p_pbar_to_D0_D0.dat'
      open(mmaxresread,
     $     file='buuinput/p_pbar__D0_D0.dat',status='old')
      read(mmaxresread,'(///)')
      read(mmaxresread,*) NoL_DD
      do is = 1, NoL_DD
        read(mmaxresread,*)       
     $       sig_ppbar_D0D0(is,1), sig_ppbar_D0D0(is,2) 
      end do
      close(mmaxresread)
      write(*,*) 'after reading p_pbar_to_D0D0.dat', 
     $  sig_ppbar_D0D0(1,1), sig_ppbar_D0D0(1,2)
      write(*,*) sig_ppbar_D0D0(NoL_DD,1),sig_ppbar_D0D0(NoL_DD,2)

c-----------------------------------------------------------------------

c----------------------by Peter Kovacs and Wolf (2017)------------------
c-----------------------------------------------------------------------
      return

      entry croswwrite
      write(isum,*)'in crosw ************'
      write(isum,*)' ohne die querschnitte '

      write(*,*)'in crosw ************'
      write(*,*)' ohne die querschnitte    '


      end


******      subroutine
      subroutine spion(srts,sigspion,iz1,iz2)

      implicit none
      include "cominput"
      include "common"

      real*8 x(5,3) , srts, xfit, s0,sigspion(2)
      integer i
      integer  iz1,iz2


      data (x(i,1), i=1,5) / 61.32865, 6.1777919, 1.520738,3.481812,
     &     2.503206/
      data (x(i,2), i=1,5) / 24.94817,1.928099,3.300971,2.4472713e-03,
     &     0.8496530/
      data (x(i,3), i=1,5) / 7.25069,2.306925, 0.883237,
     &                     3.641924, 6.6390312e-05 /

      s0 = 2.015
      xfit = (srts-s0)*5
      if(xfit.gt.0.0) then
        if(iz1+iz2.eq.1) then
*    index 1 : charged pions, index 2: neutral pions
          sigspion(1) =2.*
     &         x(1,2)* xfit**x(2,2)*exp(-(x(3,2)*xfit**x(4,2)
     &         +x(5,2)*xfit))
          sigspion(2) =
     &         x(1,3)* xfit**x(2,3)*exp(-(x(3,3)*xfit**x(4,3)
     &         +x(5,3)*xfit))
        else
          sigspion(1) =
     &         x(1,1)* xfit**x(2,1)*exp(-(x(3,1)*xfit**x(4,1)
     &         +x(5,1)*xfit))
          sigspion(2) =2.*sigspion(1)
        end if
      else
        sigspion(1) = 0.0
        sigspion(2) = 0.0
      end if
      return
      end

*******************************************************************
***************** N N - D D routines ******************************

      subroutine sddaich(srts,iz1,iz2,sig,imode,em1,em2,sig1,sig2)

      implicit none

      real*8    srts, sig, sdppd0, sdpd0, sdppdm, sdpdp
      real*8    sig_15, sig_16, sig_17, sigtemp, em1, em2
      integer imode, totc, imax,   iz1, iz2
      real*8    multipl, pfinal, pinitial, rmass, sig1, sig2
      parameter(rmass = 0.9383)

*************** NN - DD cross sections from Aichelin/Huber


      if(imode.eq.1) then
********** NN - DD
        if(iz1.eq.iz2) then
          sdppd0 = sig_15(srts)      !pp-d0dpp, nn-dpd-
          sig1   = sdppd0
          sdpdp  = sig_15(srts)/1.5   !pp-dpdp, nn-d0d0
          sig2   = sdpdp
          sig    = sdppd0 + sdpdp
        else
          sdpd0  = sig_16(srts) !pn-d0dp
          sig1   = sdpd0
          sdppdm = sig_17(srts) !pn-d-dpp
          sig2   = sdppdm
          sig    = sdpd0 + sdppdm
        end if
      else if(imode.eq.2) then
********** DD - NN
        totc = iz1+iz2
        imax = max0(iz1,iz2)
        sigtemp = 0.0

        if(totc .eq. 1) then
          if(imax.eq.2) then
***          D++ D- - pn
            sigtemp = sig_17(srts)
          else
***          D+ D0 - pn
            sigtemp = sig_16(srts)
          end if
        else if(totc.eq.0) then
          if(imax.eq.0) then
***          D0 D0 - nn
            sigtemp = sig_15(srts)/1.5*2.0
          else
***          D+ D- - nn
            sigtemp = sig_15(srts)
          end if
          sigtemp = sigtemp*0.5 !identical particles
        else if(totc.eq.2) then
          if(imax.eq.2) then
***          D++ D0 - pp
            sigtemp = sig_15(srts)
          else
***          D+ Dp - pp
            sigtemp = sig_15(srts)/1.5*2.0  !! 2.0 because id part in
                                            !! initial and final
          end if
          sigtemp = sigtemp*0.5 !identical particles
        end if


          pinitial = (srts**2+em1**2-em2**2)/(2.0*srts)
          pinitial = sqrt(pinitial**2 - em1**2)

          pfinal   = (srts**2)/(2.0*srts)
          pfinal   = sqrt(pfinal**2 - rmass**2)
          multipl  = 0.25
          sig = sigtemp*multipl*pfinal**2/pinitial**2

      end if

      return
      end

      real*8 function sig_15(x)
      implicit none
      real*8  x, output

      if(x.le. 2.152)  then
        output = 0.0
      else if(x.ge.2.152 .and. x.le. 2.44) then
        output = 150.4*(x - 2.152)**5.299
      else if(x.gt.2.44  .and. x.le. 2.63) then
        output = 2.794*x - 6.626
      else if(x.gt.2.63  .and. x.le. 2.98) then
        if(x.le.2.898) then
        output = 0.9205 - 2.763*(2.898 - x)**2
        else
        output = 0.9205 - 2.763*(x - 2.898)**2
        end if
      else if(x.gt.2.98 .and. x.le.3.24) then
        output = 2.192 - 0.4320*x
      else if(x.gt.3.24 .and. x.le.4.50) then
        output = 6.898*x**(-1.417) - 0.5140
      else
c        write(*,*)'sig_15 vigyazz sqrt(s) is big',x
        output = 6.898*4.5**(-1.417) - 0.5140
      end if

      sig_15  = output
      return
      end

      real*8 function sig_16(x)

      implicit none

      real*8  x, output

      if(x.le. 2.152)  then
        output = 0.0
      else if(x.ge.2.152 .and. x.le. 2.412) then
        output = 144.8*(x - 2.152)**5.083
      else if(x.gt.2.412 .and. x.le. 2.58) then
        output = 2.841*x - 6.707
      else if(x.gt.2.58  .and. x.le. 2.98) then
        if(x.le.2.898) then
        output = 0.7819 - 1.777*(2.896 - x)**2
        else
        output = 0.7819 - 1.777*(x - 2.896)**2
        end if
      else if(x.gt.2.98 .and. x.le.3.26) then
        output = 1.894 - 0.3798*x
      else if(x.gt.3.24 .and. x.le.4.50) then
        output = 8.371*x**(-1.926) - 0.2035
      else
c        write(*,*)'sig_16 vigyazz sqrt(s) is big',x
        output = 8.371*4.5**(-1.926) - 0.2035
      end if

      sig_16  = output

      return
      end


      real*8 function sig_17(x)

      implicit none

      real*8  x, output

      if(x.le. 2.152)  then
        output = 0.0
      else if(x.ge.2.152 .and. x.le. 2.412) then
        output = 239.6*(x - 2.152)**5.148
      else if(x.gt.2.412 .and. x.le. 2.60) then
        output = 4.738*x - 11.21
      else if(x.gt.2.60  .and. x.le. 2.98) then
        if(x.le.2.884) then
        output = 1.395 - 3.615*(2.884 - x)**2
        else
        output = 1.395 - 3.615*(x - 2.884)**2
        end if
      else if(x.gt.2.98 .and. x.le.3.69) then
        output = 3.517 - 0.7197*x
      else if(x.gt.3.69 .and. x.le.4.50) then
        output = 9.358*x**(-1.242) - 0.9878
      else
c        write(*,*)'sig_17 vigyazz sqrt(s) is big',x
        output = 9.358*4.5**(-1.242) - 0.9878
      end if

      sig_17  = output

      return
      end


******************** END OF NN- DD ******************************
*


      subroutine crosw1(pra,asrt,c2,srt,pcm,iblock,pstore,rnstore,
     &                  mdel,nt)
*-----------------------------------------------------------------------
*           angular distribution
*     com:  parametrisation is taken from the cugnon-paper (elastic)
*                              and   from Wolf90
*----------------------------------------------------------------------
      implicit none
c      include"common"
      include"cominput"

      integer ntheta
      parameter(ntheta = 20)
      real*8    pcm(3)
      real*8    srt, asrt, pra, c2
      real*8    as, ta, x, t1, t2, c1
      real*8    a, s1, s2, ct1, st1, ct2, st2, ss
      real*8    pstore(1:3), rnstore(1:5)
      real*8    sinth, srts, mdel,mn, mpi, fp, fps, lamda
      real*8    u,t,gp, c, pf2, pi2, dsdm, cost,qr2t, q2t
      real*8    zt, qr2u, q2u, zu, m1, m2, m12, aa1,aa2,aa3,aa4
      real*8    ft, fu, pt, pu,pi,rntheta
      integer iblock,nt,i
      real*8    forback, dsdodm(0:ntheta-1)
      real*8    sum,dmass,rn,gdelt
      logical flag


******************************************************************
      parameter(fp = 1.008)
      parameter(fps = 2.202)
      parameter(lamda = 0.63)
      parameter(mpi = 0.14)
      parameter(mn = 0.9383)
      parameter(dmass=1.232,gdelt=0.12,pi=3.141592654)



******************************************************************

c      write(*,*)'in crosw1', iblock


      if(iblock.eq.1) then
        as  = ( 3.65 * asrt )**6
        a   = 6.0 * as / (1.0 + as)
        ta  = -2.0 * pra**2
        x   = rn(iseed)
        if(a.gt.0.0) then
          t1  = log( (1.0-x) * exp(2.*a*ta) + x )  /  a
        else
          t1  = 0.0
        end if
        c1  = 1.0 - t1/ta
        if(abs(c1).gt.1.0) c1=2.0*x-1.0
c       write(*,*) ' crosw1 - c1 ',c1,asrt,t1,ta,x
      else if(iblock .eq. 2) then

       gp=fp*2.*mn/mpi
       c = 0.3
       srts = srt
c      write(*,*)'crosw 1  1a',mdel,srts,mn
      pi2=((srts**2-mdel**2-mn**2)**2-4.*mdel**2*mn**2)/
     &     (4.*srts**2)
      if(pi2.le.0.0) then
        write(*,*)'pi2 lt 00',srts, mdel
        pi2 = 0.0
      end if
      pf2=srts**2/4.-mn**2
      dsdm =0.
c       write(*,*)'crosw 1  1 ',pf2,pi2,mdel,srts,mn
      do i=0,ntheta-1

          cost=(float(i)+0.5)/float(ntheta)
          t=mdel**2+mn**2-2*sqrt(mdel**2+pi2)*sqrt(mn**2+pf2)
     &         +2*cost*sqrt(pi2*pf2)
          u=2.*mn**2-2*sqrt(mn**2+pf2)*sqrt(mn**2+pi2)
     &         -2*cost*sqrt(pi2*pf2)

c        write(*,*)'crosw 1  2 '
          qr2t=(dmass**2+mn**2-t)**2/(4.*dmass**2)-mn**2
          q2t=(mdel**2+mn**2-t)**2/(4.*mdel**2)-mn**2
          zt=((qr2t+c**2)/(q2t+c**2))**2
          qr2u=(dmass**2+mn**2-u)**2/(4.*dmass**2)-mn**2
          q2u=(mdel**2+mn**2-u)**2/(4.*mdel**2)-mn**2
          zu=((qr2u+c**2)/(q2u+c**2))**2
          if(zu.lt.0.) write(*,*)'zu',zu
          if(zt.lt.0.) then
            write(*,*)'zt',zt,t,mdel,srts
            stop
          end if
          m1=t*(t-(mdel-mn)**2)*((mdel+mn)**2-t)**2/(3.*mdel**2)
          m2=u*(u-(mdel-mn)**2)*((mdel+mn)**2-u)**2/(3.*mdel**2)
        aa1   =  (t*u+(mdel**2-mn**2)*(t+u) -mdel**4+mn**4)
        aa2   =  (t*u+mn*(mn+mdel)*(mdel**2-mn**2))
        aa3   =  (t*u-(mdel+mn)**2*(t+u)+(mn+mdel)**4)
        aa4   =  (t*u-mn*(mdel-mn)* (mdel**2-mn**2))
        m12 = (aa1*aa2 - aa3*aa4/3.0) / (2.*mdel**2)
c       m12=1./(2.*mdel**2)*((t*u+(mdel**2-mn**2)*(t+u)
c    &      -mdel**4+mn**4)*
c    &       (t*u+mn*(mn+mdel)*(mdel**2-mn**2))-1./3.*(t*u-
c    &       (mdel+mn)**2*(t+u)+(mn+mdel)**4)*(t*u-mn*(mdel-mn)*
c    &       (mdel**2-mn**2)))
        if(abs(m12).gt.(m1+m2)) then
c               write(*,*)' crosw1 3',m1,m2,m12,srts
c               write(*,*)' mdel ', mdel,mn, t, u
        endif
        ft=(lamda**2-mpi**2)/(lamda**2-t)
        fu=(lamda**2-mpi**2)/(lamda**2-u)
        pt=1./(t-mpi**2)
        pu=1./(u-mpi**2)
c       write(*,*)'crosw 1  4'
************** pp - n D++, do = dcos(theta)
c        dsdodm(i)=1./(32.*pi*srts**2)*sqrt(pf2/pi2)*(fps*gp/mpi)**2*
c     &       (ft**4*zt*pt**2*m1+fu**4*zu*pu**2*m2+
c     &       ft**2*fu**2*sqrt(zu*zt)*pt*pu*m12)*
c     &       2.*mdel*fz*0.389*2.0

        dsdodm(i)=  (ft**4*zt*pt**2*m1+fu**4*zu*pu**2*m2+
     &       ft**2*fu**2*sqrt(zu*zt)*pt*pu*m12)
c       write(*,*) ' dsdodm in crosw ',i,dsdodm(i), ft, pt, m1, zt

c        sigma=sigma+dsdodm/float(ntheta)*dmdelp/float(nm)
        dsdm=dsdm+dsdodm(i)/float(ntheta)
c        dsdo(i)=dsdo(i)+dsdodm*dmdelp/float(nm)
c        write(*,*)'crosw 1  5',i
        end do
c        write(*,*)'crosw 1  6'

        rntheta = rn(iseed)
        forback = 2.0*(nint(rn(iseed))-0.5)

        c1  = .3
        sum = 0.0
        flag = .true.
        i = -1
        do while(flag)
          i = i + 1
          if(i.gt.ntheta-1) then
            write(*,*)'problem in crosw ', i
            write(*,*)rntheta, dsdodm(i-1), sum , dsdm
            stop
          end if
          sum = sum + dsdodm(i)/float(ntheta)
c        write(*,*)'crosw 1  7',i, ntheta, c1, sum, rntheta*dsdm
          if(sum.ge.rntheta*dsdm) then
c        write(*,*)'crosw 1  7a',i, ntheta, c1, forback
            flag = .false.
            c1 = (float(i)+0.5) /float(ntheta)
            c1 = forback * c1
          end if
        end do

c        write(*,*)'crosw 1  8'

      else
        write(*,*)'crosw1 iblock ', iblock
        stop
      end if
      t1  = 2.0 * pi * rn(iseed)
*
*     com: set the new momentum coordinates
*
      if(pcm(1) .eq. 0.0 .and. pcm(2) .eq. 0.0) then
        t2 = 0.0
      else
        t2=atan2(pcm(2),pcm(1))
      end if
      if( (1.0-c1**2).lt.0.0) write(*,*)'crosw1 ,1 ', c1
      s1        = sqrt( 1.0 - c1**2 )
      if( (1.0-c2**2).lt.0.0) write(*,*)'crosw1 ,2 ', c2
      s2        =  sqrt( 1.0 - c2**2 )

      if(ireapl .eq.0) then
        ct1       = cos(t1)
        st1       = sin(t1)
      else if(ireapl .eq. 1) then
        if(pra.gt.0.0) then
          sinth = sqrt( 1-(pcm(3)/pra)**2)
          ct1   = pcm(1)/(pra*sinth)
          st1   = pcm(2)/(pra*sinth)
        else
          ct1   = 0.0
          st1   = 1.0
        end if
      end if
      ct2       = cos(t2)
      st2       = sin(t2)
      ss        = c2 * s1 * ct1  +  s2 * c1
      pstore(1) = ( ss*ct2 - s1*st1*st2 )
      pstore(2) = ( ss*st2 + s1*st1*ct2 )
      pstore(3) = ( c1*c2 - s1*s2*ct1 )
*-----------------------------------------------------------------------
c     write(*,*) ' crosw1 end ', pstore
      return

      end

