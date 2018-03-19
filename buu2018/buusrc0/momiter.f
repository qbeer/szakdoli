      subroutine momiter( srts,npot, mass3, mass4, id3, id4,
     +                   rho0, rho3, rho4,
     +                   betlrfx, betlrfy, betlrfz, betacm,
     +                   u3, u4, pxc, pyc, pzc, negp, switch,imes,
     +                   iwinkel,pra,c2,srt0in, nt,ranval,deltain,
     &                   rhap1, rhap2)

c call crosw1(pra, asrt,c2,srt,pcm,pr,iblock,pstore))
*----------------------------------------------------------------------*
*     srts    : sqrt(s) including the potentials                       *
*     npot    : include potentials or not                              *
*     mass3   : restmass of particle 3                                 *
*     mass4   : restmass of particle 4                                 *
*     rho0    : nuclear density at saturation                          *
*     rho3    : density at creation point of particle 3 in LRF         *
*     rho4    : density at creation point of particle 4 in LRF         *
*     betlrfi : betas for the boost from the calcframe into LRF        *
*     beta()  : betas for the boost from the CMS into the calc.frame   *
*     u3      : scalar pot of particle 3                       (output)*
*     u4      : scalar pot of particle 4                       (output)*
*     pic     : momenta in CMS                                 (output)*
*     switch  : = 1  BB - BB                                           *
*               = 2  B  - b'                                           *
*               = 3  m m                                               *
*      imes   : = 1  particle 3 is assumed to be the mes               *
*             : = 2  particle 4 is assumed to be the mes               *
*     deltain incomming delta mass                                     *
*----------------------------------------------------------------------*

      implicit none

      common /nthwhw/  nthw
      integer nthw

      integer npot, id3, id4, switch, imes,  i
      real*8    srts, mass3, mass4, rho0, rho3, rho4
      real*8    betlrfx(1:2), betlrfy(1:2) , betlrfz(1:2)
      real*8    betacm(1:3), rnstore(1:5)
      real*8    u3, u4, pxc, pyc, pzc, pcm(1:3), rnx
      real*8    betacmx, betacmy, betacmz
      logical negp, numsrt
      integer numsrtc, ranval, nit,  nt,iwinkel

      real*8    srtstart, c2,  pra , srtiter, pbetr, deltain
      real*8    sinth, costh, sinphi, cosphi, srttest, phi, pi
      real*8    asrt, pstore(1:3), rn, srt0in, mdel
      integer numsrtcm
      real*8    rhap1(1:3), rhap2(1:3)

      parameter(pi = 3.141592654)
      parameter(numsrtcm = 5)

      betacmx = betacm(1)
      betacmy = betacm(2)
      betacmz = betacm(3)

      u3 = 0.0
      u4 = 0.0
      negp = .true.
      srtstart = srt0in
c     write(*,*) ' start momiter', nthw,id3,id4,mass3,mass4
*----------------------------------------------------------------------*
*        the following loop has to be done for the ang. dist.          *


      if(switch .eq.2 .or. (switch.eq.1 .and. iwinkel.eq.0).or.
     &   switch.eq.3 ) then
*----------------------------------------------------------------------*
*           determine the angles in the CMS randomly                   *
*           assume isotropic distribution                              *
*        ranval = 179

        rnx =  rn(ranval)
        costh = -1.0 + 2.0*rnx
        rnx =  rn(ranval)
        phi = 2.0*pi*rnx
        sinth   = sqrt(1.0 - costh**2)
        sinphi  = sin(phi)
        cosphi  = cos(phi)

c     write(*,*)  ' momiter 1  , vor call iterp  ', switch, iwinkel
c     write(*,*) srts, npot, mass3, mass4, id3, id4,
c    +           rho0, rho3, rho4,
c    +           betlrfx, betlrfy, betlrfz,
c    +           betacmx, betacmy, betacmz,
c    +           u3, u4, pxc, pyc, pzc,
c    +           negp, switch, imes,
c    +           costh, sinth, sinphi, cosphi, nit,rhap1, rhap2
      call iterp( srts, npot, mass3, mass4, id3, id4,
     +           rho0, rho3, rho4,
     +           betlrfx, betlrfy, betlrfz,
     +           betacmx, betacmy, betacmz,
     +           u3, u4, pxc, pyc, pzc,
     +           negp, switch, imes,
     +           costh, sinth, sinphi, cosphi, nit,rhap1, rhap2)
c     write(*,*) '  result of iterp ',
c    +            srts, npot, mass3, mass4, id3, id4,
c    +           rho0, rho3, rho4,
c    +           betlrfx, betlrfy, betlrfz,
c    +           betacmx, betacmy, betacmz,
c    +           u3, u4, pxc, pyc, pzc,
c    +           negp, switch, imes,
c    +           costh, sinth, sinphi, cosphi, nit,rhap1, rhap2,
c    +           ' end result iterp '
*----------------------------------------------------------------------*
      else if (switch .eq.1 .and. iwinkel .ne. 0) then
c       write(*,*)  ' momiter  1 0  ', nthw
        pcm(1)  = pxc
        pcm(2)  = pyc
        pcm(3)  = pzc

        do i = 1,5
         rnstore(i) =  rn(ranval)
        end do


        numsrtc = 0
        numsrt  = .true.
        srtiter = srt0in
        do while(numsrt)

          numsrtc =  numsrtc + 1

          if(numsrtc.eq.numsrtcm) then
            if(srt0in.ge.srts) then
              srtiter = srt0in
            else
              srtiter = srts
            end if
          end if
            asrt    = srtiter - mass3 - mass4
            mdel = deltain

c      write(*,*)  ' momiter  1 0 = vor crosw1  ',pra,asrt,c2,srtiter,
c    1             pcm,iwinkel,pstore,rnstore,mdel, nt
          call crosw1(pra,asrt,c2,srtiter,pcm,iwinkel,pstore,rnstore,
     &                 mdel, nt)
c       write(*,*)  ' momiter  1 0 = nach crosw1  '
          costh   = pstore(3)
          if((1.0-costh**2).lt.0.0) write(*,*)'momiter 1',costh
          sinth  = sqrt(1.0 - costh**2)
          if(abs(sinth).gt.1.0e-08) then
            cosphi = pstore(1)/sinth
            sinphi = pstore(2)/sinth
          else
            cosphi = 1.0
            sinphi = 0.0
          end if

c       write(*,*)  ' momiter  1 0 = vor iterp 2  '
          call iterp( srts, npot, mass3, mass4, id3, id4,
     +               rho0, rho3, rho4,
     +               betlrfx, betlrfy, betlrfz,
     +               betacmx, betacmy, betacmz,
     +               u3, u4, pxc, pyc, pzc,
     +               negp, switch, imes,
     +               costh, sinth, sinphi, cosphi, nit,rhap1,rhap2)
           pbetr    = pxc**2 + pyc**2 +pzc**2
           srttest  = sqrt(mass3**2 +pbetr) + sqrt(mass4**2+pbetr)
c       write(*,*)  ' momiter  1 0 nach iterp  ', mass3, mass4

           srtstart = srttest

*           if(abs(srtiter-srtstart).lt.1.0e-04.and.(.not.negp)) then
           if(.not.negp) then
             numsrt = .false.
           end if
           if( numsrtc.gt.numsrtcm) then
c              write(*,*)'die 5 winkel sind abgearbeitet '
              numsrt = .false.
              negp   = .true.
           end if
        end do
      end if


*     now everything is given in the cms of the particles

           pbetr   = sqrt(pxc**2 + pyc**2 +pzc**2)
           srttest = sqrt((mass3+u3)**2 +pbetr**2) +
     +               sqrt((mass4+u4)**2 +pbetr**2)

          if(.not.negp.and.abs(srttest-srts).gt.1.0e-05) then
            write(123,*)'prob in momiter ', srts, srttest, nit,iwinkel
            write(*,*)'prob in momiter ', srts, srttest, nit,iwinkel
            write(*,*)mass3, mass4, u3, u4
            write(*,*)rho3, rho4, id3, id4
            stop
          end if

      return
      end
***********************************************************************
*---------------------------------------------------------------------*
*         header routine for the iteration                            *

      subroutine iterp( srts, npot, mass3, mass4, id3, id4,
     +                 rho0, rho3, rho4,
     +                 betlrfx, betlrfy, betlrfz,
     +                 betacmx, betacmy, betacmz,
     +                 u3, u4, px, py, pz,
     +                 negp, switch, imes,
     +                 costh, sinth, sinphi, cosphi, nit,rhap1, rhap2)

      implicit none
      real*8 pi
      parameter(pi = 3.1415)
      common /nthwhw/  nthw
      integer nthw

      real*8  srts, pstart, mass3, mass4, rho3, rho4,rho0
      real*8  epsilon, plast, pnew, p2, pinp2 , delta, psq
      integer icount, icmax, npot
      real*8    betlrfx(1:2), betlrfy(1:2), betlrfz(1:2)
      real*8    betacmx, betacmy, betacmz, srttest
      real*8    costh, sinphi, cosphi, sinth
      integer id3, id4, switch , imes, nit
      logical flag, negp
      real*8    u3, u4,px,py,pz
      real*8    rhap1(1:3), rhap2(1:3)

      epsilon = 1.0e-03
      icmax   = 100

c     p2help = (s + (mass3+u3)**2 - (mass4+u4)**2)**2
c     p2help = p2help/(4.0*s)
      pstart = 0.001
      plast  = pstart
      icount = 0
      u3 = 0.0
      u4 = 0.0
      flag   = .true.

      do while(flag)
        icount = icount + 1
        nit  = icount
        pinp2 = plast

c       write(*,*) ' iterp:  while loop', icount, pinp2,
c    1               ' new ', p2, pnew, plast
        call psquar( srts, pinp2, mass3, mass4, id3, id4,
     +              rho3, rho4,rho0,npot,switch,imes,
     +              betlrfx, betlrfy, betlrfz, betacmx, betacmy,
     +              betacmz,costh, sinth,cosphi,sinphi,u3,u4,psq,
     +              px, py, pz, rhap1, rhap2)

        p2    = psq
        if(p2 .ge. 0.0) then
          pnew = sqrt(p2)
          px = pnew*cosphi*sinth
          py = pnew*sinphi*sinth
          pz = pnew*costh
          negp =.false.
          delta = abs(pnew - plast)
          if(delta.lt.epsilon) then
             srttest = sqrt((mass3+u3)**2 +pnew**2) +
     +       sqrt((mass4+u4)**2 +pnew**2)
             if(abs(srttest-srts).gt.1.0e-05) then
               negp = .true.
             end if
          end if
        else
          negp = .true.
          pnew =  pstart  + float(icount)*0.001
          delta = 100000.0*epsilon
        end if

        if(delta .lt. epsilon) then
          flag = .false.
        else
          plast = pnew
        end if

        if(icount. gt. icmax) then
c   write(*,*)'mimiter fatal error ',switch
c       write(*,*)  '  momiter called iterp  ',
c     +           srts, npot, mass3, mass4, id3, id4,
c     +                 rho0, rho3, rho4,
c     +                 betlrfx, betlrfy, betlrfz,
c     +                 betacmx, betacmy, betacmz,
c     +                 u3, u4, px, py, pz,
c     +                 negp, switch, imes,
c     +                 costh, sinth, sinphi, cosphi, nit,rhap1, rhap2
          flag = .false.
          negp = .true.
        end if
      end do             !  while-flag
      return
      end

*----------------------------------------------------------------------*
*                     do the iteration                                 *
      subroutine psquar( srts, pinp2, mass3, mass4, id3, id4,
     +                  rho3, rho4,rho0,npot,switch,imes,
     +                  betlrfx, betlrfy, betlrfz, betacmx, betacmy,
     +                  betacmz,costh, sinth, cosphi, sinphi,u3,u4,psq,
     +                  px, py, pz, rhap1, rhap2)

      implicit none


      real*8  srts, pinp2, mass3, mass4, rho0, rho3, rho4
      real*8  s, u3, u4 , pin,  potanal, v3, v4, potmes
      real*8  paux(1:3), px, py, pz,  en3, en4 ,psq
      integer npot, id3, id4, switch, imes
      real*8  betlrfx(1:2), betlrfy(1:2), betlrfz(1:2)
      real*8  betacmx, betacmy, betacmz
      real*8  costh, sinth
      real*8  p2help , pabss, cosphi, sinphi
      real*8  veccm4, veccm3, plrf
      logical test
      real*8  delta, invv, invn
      integer nit
      real*8  rhap1(1:3), rhap2(1:3), plrfp(1:3)
      delta = 1.0e-04

*--------------------------------------------------------------------*
*         build explicit vector - needed for boost back to the LRF   *


      paux(1) = pinp2*cosphi*sinth
      paux(2) = pinp2*sinphi*sinth
      paux(3) = pinp2*costh
c     write(*,*)  ' start  psquare ', npot, switch, paux

*--------------------------------------------------------------------*

      if(npot .eq. 1) then
*--------------------------------------------------------------------*
*      evaluate the abs. value of momentum of particle 3 in the LRF  *
        nit = 0
        test = .true.
        do while(test)
          nit = nit + 1
          px  = paux(1)
          py  = paux(2)
          pz  = paux(3)
          en3  = sqrt((mass3+u3)**2 + px**2 + py**2 + pz**2)
          invv = en3**2 - px**2 - py**2 -pz**2
          call lorentz(-betacmx ,-betacmy ,-betacmz ,px,py,pz,en3)
          call lorentz(betlrfx(1),betlrfy(1),betlrfz(1),px,py,pz,en3)
          pin   = sqrt(px**2 + py**2 + pz**2)
          plrf  = pin
          pabss = pin
          plrfp(1) = px
          plrfp(2) = py
          plrfp(3) = pz

          v3=0.0
          if(switch .eq. 1) then
            v3    = potanal(rho0,rho3,plrfp,id3,rhap1)
          else if((switch .eq.2) .and. (imes .eq. 1)) then
            v3    = potmes( rho3,pin,id3)
          else if((switch .eq.2) .and. (imes .eq. 2)) then
            v3    = potanal(rho0,rho3,plrfp,id3,rhap1)
          else if(switch .eq.3) then
            v3    = potmes( rho3,pin,id3)
          end if
          u3   = -mass3 + sqrt(mass3**2 +
     +           2.0*sqrt(pabss**2+mass3**2)*v3 + v3**2)
          invn = (mass3+u3)**2
          if(abs(invn-invv).lt.delta)then
            test = .false.
          end if
c         if (nit .lt. 3)
c    1   write(*,*) 'in momiter: square ', nit, invv, invn,plrfp
          if(nit.gt.200) then
            write(*,*)'problems in momiter 1-3 ',nit, mass3,invv,invn
            write(*,*)u3,en3,paux
            psq = -10.000
            goto 99
*            stop
          end if
        end do
*        the corresponding vectorpot in the cms frame
         pin    = sqrt(paux(1)**2 + paux(2)**2 + paux(3)**2)
         veccm3 = -sqrt(pin**2 + mass3**2)+ sqrt(pin**2 + mass3**2 +
     +                2.0*mass3*u3 + u3**2)


*      evaluate the abs. value of momentum of particle 4 in the LRF  *
        nit = 0
        test = .true.
        do while(test)
          nit = nit + 1
          px  = -paux(1)
          py  = -paux(2)
          pz  = -paux(3)
          en4  = sqrt((mass4+u4)**2 + px**2 + py**2 + pz**2)
          invv = en4**2 - px**2 - py**2 -pz**2
          call lorentz(-betacmx ,-betacmy ,-betacmz ,px, py, pz, en4)
          call lorentz(betlrfx(2),betlrfy(2),betlrfz(2),px, py, pz, en4)
          pin = sqrt(px**2 + py**2 + pz**2)
          pabss = pin
          plrf = pin
          plrfp(1) = px
          plrfp(2) = py
          plrfp(3) = pz

          v4=0.0
          if(switch .eq. 1) then
           v4    = potanal(rho0,rho4,plrfp,id4,rhap2)
          else if((switch .eq.2) .and. (imes .eq. 2)) then
           v4    = potmes(rho4,pin,id4)
          else if((switch .eq.2) .and. (imes .eq. 1)) then
           v4    = potanal(rho0,rho4,plrfp,id4,rhap2)
          else if(switch.eq.3) then
           v4    = potmes(rho4,pin,id4)
          end if
          u4   = -mass4 + sqrt(mass4**2 +
     +           2.0*sqrt(pabss**2+mass4**2)*v4 + v4**2)
          invn = (mass4+u4)**2
          if(abs(invn-invv).lt.delta)then
            test = .false.
          end if
          if(nit.gt.200) then
            write(*,*)'problems in momiter 1-4 ',nit, mass4,invv,invn
            write(*,*)u4,en4,paux
            psq = -10.000
            goto 99
          end if
        end do

        pin    = sqrt(paux(1)**2 + paux(2)**2 + paux(3)**2)
        veccm4 = -sqrt(pin**2 + mass4**2)+ sqrt(pin**2 + mass4**2 +
     +                2.0*mass4*u4 + u4**2)

      else if(npot .eq. 0) then
        u3     = 0.0
        u4     = 0.0
        veccm3 = 0.0
        veccm4 = 0.0

      end if


      s = srts**2
      p2help = (s + (mass3+u3)**2 - (mass4+u4)**2)**2
      p2help = p2help/(4.0*s)
      psq    = p2help - (mass3+u3)**2
 99   continue

      return
      end
