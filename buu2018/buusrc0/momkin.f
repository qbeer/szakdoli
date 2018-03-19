***********************************************************************
      subroutine momkin(pcm,pr,srt,em1, em2,id3,id4,j01, j02,
     &                  betacm, betlrfx,betlrfy,betlrfz,testflag,nt,
     &                  pcmfin,em3,em4,u3,u4,iwinkel,srtinfree,deltain,
     &                  rhap1,rhap2)

      implicit none
      include"common"
      include"cominput"

      real*8    pcm(1:3), srt, srtinfree, j01, j02
      real*8    betacm(3), betlrfx(2), betlrfy(2), betlrfz(2)
      real*8    pcmfin(3), em3, em4, u3, u4, em1, em2
      integer id3, id4, nt, nhw0, nhw1
      logical testflag
      integer idhelp1, idhelp2

      real*8    testmass, rm2, pm2, pr, c2
      real*8    cutpi, x1, rn
      real*8    pin, potanal, vecpot, maxnucpot2, maxnucpot1
      real*8    remass3, remass4
      integer iwinkel
      real*8    phx, phy, phz, eh, pabsh, vecmax1, vecmax2
      real*8    srtit, newmass3, newmass4, pxc, pyc, pzc
      integer nptest, npot
      logical massflag, nope, negp
      real*8    srtstart , deltain
      real*8    plrf(1:3), rhap1(1:3), rhap2(1:3)


c      include"resdata"
c      include"resdata1"

      testmass = 5.0*rmass
      rm2 = rmass**2
      pm2 = pmass**2
      if(pr.lt.1.0e-08) then
        write(*,*)'momkin pr<<', pr
      end if
      c2  = pcm(3) / pr
      if(pr.lt.1.0e-08) then
        write(*,*)'after pr', pr
      end if

c      write(*,*)'momk 1'
      cutpi     =2.0*rmass+pmass


      testflag = .false. ! this flag decides if procedure succsessful

*     randomize  ids
      x1 = rn(iseed)
      if(x1.gt. 0.5) then
        idhelp1 = id3
        idhelp2 = id4
      else
        idhelp1 = id4
        idhelp2 = id3
      end if
      id3 = idhelp1
      id4 = idhelp2

      if(id3.lt.1 .or. id3.gt.nres+3) write(*,*) 'hiba momkin1',id3
      remass3 = rmass
      if(id3.ge.2 .and. id3.le.nres+1) then
        remass3   = resprop1(id3-1,1)
      else if(id3.eq.nres+2) then
        remass3 = xlmas
      else if(id3.eq.nres+3) then
        remass3 = xsmas
      end if

      if(id4.lt.1 .or. id4.gt.nres+3) write(*,*) 'hiba momkin2',id4
      remass4 = rmass
      if(id4.ge.2 .and. id4.le.nres+1) then
        remass4   = resprop1(id4-1,1)
      else if(id4.eq.nres+2) then
        remass4 = xlmas
      else if(id4.eq.nres+3) then
        remass4 = xsmas
      end if


c      write(*,*)'momk 2'

*--------------------------------------------------------------------*
*     evaluate the kinematics for the final channel                  *
*                                                                    *

      if(id3+id4.gt.2 .and. id3.le.nres+1 .and. id4.le.nres+1) then
*       needed only for he final channel NR, DD , NOT for NN
*       and not for anything+(S or L) !!! zm




c        write(*,*)'momk 3'
*
*     determine the scalar potentials for p = 0 in the lrf in order  *
*     to be able to determine the thresholds                         *

        pin         = 0.0
        plrf(1) = 0.0
        plrf(2) = 0.0
        plrf(3) = 0.0

        vecpot      = potanal(rho0, j01, plrf, id3,rhap1)
        maxnucpot1  = -remass3 + sqrt(remass3**2 +
     +                2.0*remass3*vecpot + vecpot**2)

        pin         = 0.0
        plrf(1) = 0.0
        plrf(2) = 0.0
        plrf(3) = 0.0

        vecpot      = potanal(rho0, j02, plrf, id4, rhap2)
        maxnucpot2  = -remass4 + sqrt(remass4**2 +
     +                2.0*remass4*vecpot + vecpot**2)

*     for the upper bound in the mass determination one has to express
*     the max. scalar pot. in terms of the vector pot.. therefor one
*     needs the mom. in this frame that corresponds to p=0 in the LRF

        phx     = 0.0
        phy     = 0.0
        phz     = 0.0
        eh      = testmass + maxnucpot1
        call lorentz(-betlrfx(1),-betlrfy(1),-betlrfz(1),phx,phy,phz,eh)
        call lorentz(betacm(1), betacm(2),betacm(3),phx,phy,phz,eh)
        pabsh   = sqrt(phx**2+phy**2+phz**2)
        vecmax1 = -sqrt(pabsh**2+testmass**2)+sqrt(pabsh**2+
     +            testmass**2+2.0*testmass*maxnucpot1+maxnucpot1**2)

        phx     = 0.0
        phy     = 0.0
        phz     = 0.0
        eh      = testmass + maxnucpot2
        call lorentz(-betlrfx(2),-betlrfy(2),-betlrfz(2),phx,phy,phz,eh)
        call lorentz(betacm(1), betacm(2),betacm(3),phx,phy,phz,eh)
        pabsh   = sqrt(phx**2+phy**2+phz**2)
        vecmax2 = -sqrt(pabsh**2+testmass**2)+sqrt(pabsh**2+
     +            testmass**2+2.0*testmass*maxnucpot2+maxnucpot2**2)


*     maximumal possible free sqrt(s) for the final channel
        srtit = srt - vecmax1 - vecmax2

        if(srtit.le.0.0) then
          write(*,*)'warning 1 in crosw: srtit = ', srtit
          stop
        end if

c      write(*,*)'momk 4'
        nope     =.false.
        negp     =.true.
        nptest   = 0
        massflag =.false.
        srtit = srtit + 0.01
        do while(negp .and. .not.nope. and. nptest.lt.20)

          srtit  = srtit - 0.005
          nptest = nptest + 1
          if(id3.eq.2 .and. id4.eq.2) then
            if(srtit .le. 2.0*rmass + 2.0*pmass) then
              nope = .true.
            end if
          else
            if(srtit.le. 2.0*rmass+pmass) then
              write(*,*)'momkin prob. 111'
              nope =.true.
            end if
          end if
c      write(*,*)'momk 5'
c--------------------------------------------------
c                              statement by gyuri for  pp collisions
        if (massta + masspr .eq. 2)  srtit = srt
c-------------------------------------------------------------
          if(.not.nope) then
c      write(*,*)'momk 6'
            newmass3 = 0.0
            newmass4 = 0.0
            call mommass(srtit,id3, id4, newmass3, newmass4, massflag,
     &                    srtinfree)
c      write(*,*)'momk 7', newmass3, newmass4
            if(massflag) then
              em3     = newmass3
              em4     = newmass4
            else
              nope =.false.
            end if
          end if
c      write(*,*)'momk 8', massflag
*          do the momenta                                            *
          if(massflag) then
c      write(*,*)'momk 9'
            npot = 1
            pxc = pcm(1)
            pyc = pcm(2)
            pzc = pcm(3)
            srtstart = sqrt(em1**2+pr**2)+sqrt(em2**2+pr**2)
c      write(*,*)'momk 10'
            if(id3.gt.1) deltain = em3
            if(id4.gt.1) deltain = em4
            nhw1 = 1
            nhw0 = 0
c           write(*,*) ' momkin calls momiter 1'
            call momiter(srt,npot,em3,em4,
     +                   id3, id4,rho0,j01,j02,
     +                   betlrfx,betlrfy,betlrfz,betacm,
     +              u3, u4, pxc, pyc,pzc,negp,nhw1,nhw0,iwinkel,
     +                   pr, c2,srtstart, nt,iseed,deltain,
     &                   rhap1,rhap2)
c           write(*,*)'momk 11'

            if(.not. negp) then
              pcmfin(1)  = pxc
              pcmfin(2)  = pyc
              pcmfin(3)  = pzc
              testflag   =.true.
              nope = .true.
            else
              massflag = .false.
            end if
          end if
c      write(*,*)'momk 12'
        end do   !end of while loop

c      write(*,*)'momk 13'
        if(nptest.gt.20) then
          write(*,*)'in momkin.f laeuft noch was schief '
          stop
        end if
c             write(*,*)'momk 14'



      else if(id3+id4 .eq. 2) then
c       write(*,*)'momk 15'
*          do the momenta  for NN final state                                *
            negp =. true.
            npot = 1
            pxc = pcm(1)
            pyc = pcm(2)
            pzc = pcm(3)
            em3 = rmass
            em4 = rmass
            srtstart = sqrt(em1**2+pr**2)+sqrt(em2**2+pr**2)
c      write(*,*)'momk 16' , iwinkel
c           write(*,*) ' momkin calls momiter 2'
            call momiter(srt,npot,em3,em4,
     +                   id3, id4,rho0,j01,j02,
     +                   betlrfx,betlrfy,betlrfz,betacm,
     +                   u3, u4, pxc, pyc,pzc,negp,1,0,iwinkel,
     +                   pr, c2,srtstart,nt,iseed,deltain,rhap1,rhap2)
c            write(*,*)srt,npot,em3,em4,
c     +                   id3, id4,rho0,j01,j02,
c     +                   betlrfx,betlrfy,betlrfz,betacm,
c     +                   u3, u4, pxc, pyc,pzc,negp,1,0,iwinkel,
c     +                   pr, c2,srtstart,nt

c      write(*,*)'momk 17'
            if(.not. negp) then
              pcmfin(1)  = pxc
              pcmfin(2)  = pyc
              pcmfin(3)  = pzc
              testflag   =.true.
              nope = .true.
            else
              massflag = .false.
            end if
          end if

c            write(*,*)'momk 18'



        return
        end
