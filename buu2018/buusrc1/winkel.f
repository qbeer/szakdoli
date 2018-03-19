      subroutine sstatekin(srt,un1,un2,upi,ppi3,pn1,pn2,testflag,
     +                i1,i2,rpi3,betacm,srtfreei,rhap1,rhap2)
      implicit none
      include"common"
      include"cominput"



      real*8    srt, upi, un1, un2
      real*8    ppi3(3),pn1(3),pn2(3),rpi3(3)
      integer i1, i2

      real*8    m1, m2, m3, m12

      integer i, icount
      real*8  m12max, m12min, rn, m12del
      real*8  absp3, absp1st, phi3, cos3, phi1st, cos1st
      real*8    betacm(3),deriv(0:4)
      real*8  betarf(1:3), x2
      real*8  x1,   etot, check, test
      real*8  m12ppmin, m12ppmax,   m12ppdel
      real*8  j0, j1, j2, j3, betlrfx, betlrfy, betlrfz
      logical flag, testflag
      real*8  pb1, pb2, pb3, pb0, pabs, pin, potmes
      real*8  ps0, ps1, ps2, ps3, x, y, z,p3(0:3)
      real*8  blrf1(1:4) , blrf2(1:4), prf(1:3), psinv, pbinv
      real*8  vecpo, scapo, srtfree, ppiabs
      real*8  srtfreei, rhap1(1:3), rhap2(1:3)
      integer ntest
*
      testflag = .false.
*
*----------------------------------------------------------------------
*      determine space coordinate of pion
*
      x = 0.5*(r(1,i1)+r(1,i2))
      y = 0.5*(r(2,i1)+r(2,i2))
      z = 0.5*(r(3,i1)+r(3,i2))

c      write(*,*)'flagg = ', flagg
*     det boost prop. for the pion from calc. frame into the LRF

      call linint1(x,y,z,deriv)
      j0    = deriv(0)
      j1    = deriv(1)
      j2    = deriv(2)
      j3    = deriv(3)

      if(j0 .gt. 1.0e-6) then
        betlrfx = j1/j0
        betlrfy = j2/j0
        betlrfz = j3/j0
      else
        betlrfx = 0.0
        betlrfy = 0.0
        betlrfz = 0.0
      end if

cc      if(flagg) write(*,*)' winkel 1'
              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba winkel1 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
      call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)
cc      if(flagg) write(*,*)' winkel 2 '

***** j0 contains the density of the LRF where the pion is going to sit

*
*-----------------------------------------------------------------------
*     do the Monte-Carlo decision for the kinematics
*

      m3 = pmass
      m1 = rmass
      m2 = rmass

      if(srt.lt.srtfreei .and. srt.gt.(2.0*rmass+pmass) ) then
        srtfree = srt
      else
        srtfree = srtfreei
      end if

       m12max = srtfree - m3
       m12min = m1 + m2
      m12del = m12max - m12min
*     boundaries


*     calculate p3 in the restframe of the dec. particle
        absp3 =  (srtfree**2 -(m12min+m3)**2)*
     &           (srtfree**2 -(m12min-m3)**2)
        absp3 = sqrt(absp3)/2./srtfree
*     calculate p1 in the resframe of particle 1 and 2
        absp1st = (m12min**2 -(m1+m2)**2)*(m12min**2 -(m1-m2)**2)
        absp1st = sqrt(absp1st)/2./m12min
        m12ppmin = absp3*absp1st*m12min
*     calculate p3 in the restframe of the dec. particle
        absp3 =  (srtfree**2 -(m12max+m3)**2)*
     &           (srtfree**2 -(m12max-m3)**2)
        absp3 = sqrt(absp3)/2./srtfree
*     calculate p1 in the resframe of particle 1 and 2
        absp1st = (m12max**2 -(m1+m2)**2)*(m12max**2 -(m1-m2)**2)
        absp1st = sqrt(absp1st)/2./m12max
        m12ppmax = absp3*absp1st*m12max
        m12ppdel = m12ppmax - m12ppmin


      ntest = 0
 10   continue
      ntest = ntest + 1
      testflag = .false.

*      determine  m12
        flag = .true.
        icount = 0
        do while(flag)
          x1 = rn(iseed)
          icount = icount + 1
          m12= m12min + x1*m12del
*
*     calculate p3 in the restframe of the dec. particle
          absp3 =  (srtfree**2 -
     &          (m12+m3)**2)*(srtfree**2 -(m12-m3)**2)
          absp3 = sqrt(absp3)/2./srtfree
*     calculate p1 in the resframe of particle 1 and 2
          absp1st = (m12**2 -(m1+m2)**2)*(m12**2 -(m1-m2)**2)
          absp1st = sqrt(absp1st)/2./m12
          test = m12*absp3*absp1st
          x2 = rn(iseed)
          check = m12ppmin + x2 * m12ppmax
          if(test.gt.check) then
            flag     =.false.
            testflag = .true.
          end if
          if(icount.gt.1000) then
            flag =.false.
            write(*,*)'1000   erreicht in sstate kin '
            testflag =.false.
          end if
        end do
*     we have found a m12 (lucky)
*     now do the angles
        if(testflag) then

          x1 = rn(iseed)
          phi3 = 2*pi*x1
          x1 =rn(iseed)
          cos3 = -1.0 + 2.0*x1
          x1 = rn(iseed)
          phi1st = 2*pi*x1
          x1 =rn(iseed)
          cos1st  = -1.0 + 2.0*x1
          p3(1) = absp3 * cos(phi3) * sqrt(1.0-cos3**2)
          p3(2) = absp3 * sin(phi3) * sqrt(1.0-cos3**2)
          p3(3) = absp3 * cos3
          p3(0) = sqrt(m3**2 + p3(1)**2 +p3(2)**2 +p3(3)**2)

****  determine the pion optical pot.
cc      if(flagg) write(*,*)'  winkel 11'
          ps0 = p3(0)
          ps1 = p3(1)
          ps2 = p3(2)
          ps3 = p3(3)

          psinv  = ps0**2 - ps1**2 - ps2**2 - ps3**2
          upi    = 0.0
          icount = 0
          flag   = .true.
          do while(flag)
            icount = icount + 1
            pb1 = ps1
            pb2 = ps2
            pb3 = ps3
            pb0 = sqrt((pmass+upi)**2+ps1**2+ps2**2+ps3**2)
            psinv  = pb0**2 - ps1**2 - ps2**2 - ps3**2
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba winkel2 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
c            stop
          end if
         call lorentz(-betacm(1),-betacm(2),-betacm(3),pb1,pb2,pb3,pb0)
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba winkel3 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
c            stop
          end if
          call lorentz(betlrfx,betlrfy,betlrfz,pb1,pb2,pb3,pb0)
          pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
          pin    = pabs
          vecpo  = potmes(j0,pin ,1)

          scapo  = -m3 + sqrt(m3**2 +
     +              2.0*sqrt(pabs**2+m3**2)*vecpo + vecpo**2)

          upi    = scapo
          pb0 = sqrt((pmass+upi)**2+ps1**2+ps2**2+ps3**2)
          pbinv  = pb0**2 - ps1**2 - ps2**2 - ps3**2
          if(abs(pbinv-psinv).lt.1.0e-04) then
            flag = .false.
          end if
          if(icount.gt.100) then
            write(*,*)'hier laeuft was schief '
            stop
          end if
        end do
        ppi3(1) = p3(1)
        ppi3(2) = p3(2)
        ppi3(3) = p3(3)
        ppiabs = sqrt(p3(1)**2+p3(2)**2+p3(3)**2)
********* det. boost params. for boost into rest-frame of 1,2

         etot = srt - sqrt( (m3+upi)**2 + ppiabs**2 )
         betarf(1) = -p3(1)/etot
         betarf(2) = -p3(2)/etot
         betarf(3) = -p3(3)/etot
         do i = 1,4
           blrf1(i) =  betlrfboo(i1,i)
           blrf2(i) =  betlrfboo(i2,i)
         end do
         testflag =.false.
         m12 = (srt - sqrt((pmass+upi)**2+ppiabs**2) )**2
         m12 = sqrt(m12 - ppiabs**2)
         call momkinsst(m12,betacm,betarf,blrf1,blrf2,prf,
     &                un1,un2,testflag,iseed,rho0, rmass,rhap1,rhap2)



      if(.not.testflag.and.ntest.lt.5) then
        goto 10
      end if

      if(.not.testflag) then
        write(*,*)'problem in winkel '
        write(*,*)srt, srtfreei, srtfree , ntest
      end if
****   boost momenta of nucleons back to the CMS of incomming
****   nucleons


         ps1 = prf(1)
         ps2 = prf(2)
         ps3 = prf(3)
         ps0 = sqrt((rmass+un1)**2+ps1**2+ps2**2+ps3**2)
          if(ps1**2+ps2**2+ps3**2.gt.ps0**2) then
            write(*,*) "hiba winkel4 lorentz, negative mass",
     &                ps1,ps2,ps3,ps0
c            stop
          end if
         call lorentz(-betarf(1),-betarf(2),-betarf(3),
     &               ps1,ps2,ps3,ps0)
         pn1(1) = ps1
         pn1(2) = ps2
         pn1(3) = ps3

         ps1 = -prf(1)
         ps2 = -prf(2)
         ps3 = -prf(3)
         ps0 = sqrt((rmass+un2)**2+ps1**2+ps2**2+ps3**2)
          if(ps1**2+ps2**2+ps3**2.gt.ps0**2) then
            write(*,*) "hiba winkel5 lorentz, negative mass",
     &                ps1,ps2,ps3,ps0
c            stop
          end if
         call lorentz(-betarf(1),-betarf(2),-betarf(3),
     &               ps1,ps2,ps3,ps0)
         pn2(1) = ps1
         pn2(2) = ps2
         pn2(3) = ps3

***     end of testflag-if
       end if
         rpi3(1) = x
         rpi3(2) = y
         rpi3(3) = z


        return
        end


*************************************************************************
*************************************************************************

      subroutine momkinsst(srts,betacm,betarf,blrf1,blrf2,prf,
     &             un1,un2,testflag,ranval,rho0,rmass,rhap1,rhap2)


      implicit none

      real*8     pi
      parameter(pi = 3.141592654)
      real*8     rmass
      real*8      betacm(3),betarf(1:3),blrf1(1:4),blrf2(1:4)
      real*8    prf(1:3), un1, un2, srts
      logical testflag
      real*8    rnx, costh, phi, sinth, cosphi, sinphi
      real*8    pbetr, srttest, rho0, rn
      integer ranval, nit
      real*8    rhap1(1:3), rhap2(1:3)
      un1 = 0.0
      un2 = 0.0

*--------------------------------------------------------------------*
*           determine the angles in the rest frame of the two        *
*           nucleons randomly                                        *
*           assume isotropic distribution                            *

        rnx =  rn(ranval)
        costh = -1.0 + 2.0*rnx
        rnx =  rn(ranval)
        phi = 2.0*pi*rnx
        sinth   = sqrt(1.0 - costh**2)
        sinphi  = sin(phi)
        cosphi  = cos(phi)


*--------------------------------------------------------------------*
      call iterps(srts,betacm,betarf,blrf1,blrf2,prf,un1,un2,testflag,
     &           costh,sinth,sinphi,cosphi,rho0,rmass,nit,rhap1,rhap2)
*--------------------------------------------------------------------*

*     now everything is given in the restframe of 1 and 2

           pbetr   = sqrt(prf(1)**2 + prf(2)**2 +prf(3)**2)
           srttest = sqrt((rmass+un1)**2 +pbetr**2) +
     +               sqrt((rmass+un2)**2 +pbetr**2)

          if(testflag .and. abs(srttest-srts).gt.1.0e-05) then
            write(123,*)'prob in sstate ', srts, srttest, nit
            write(*,*)'stopped in sstate '
            stop
          end if


      return
      end
***********************************************************************
*---------------------------------------------------------------------*
*         header routine for the iteration                            *

      subroutine iterps(srts,betacm,betarf,blrf1,blrf2,prf,un1,un2,
     &            testflag,costh,sinth,sinphi,cosphi,rho0,rmass,nit,
     &            rhap1, rhap2)

      implicit none
      real*8 pi
      parameter(pi = 3.1415)
      real*8    srts, pstart, rho0
      real*8    epsilon, plast, pnew, p2, pinp2 , delta, psq
      integer icount, icmax,   nit
      real*8    costh, sinphi, cosphi, sinth
      logical flag, negp, testflag
      real*8    betacm(1:3), betarf(1:3), blrf1(1:4), blrf2(1:4)
      real*8    un1, un2,rmass
      real*8    u3, u4,px,py,pz, prf(1:3), rhap1(1:3), rhap2(1:3)

      epsilon = 1.0e-03
      icmax   = 100
      pstart = 0.0
      plast  = pstart
      icount = 0
      u3 = 0.0
      u4 = 0.0
      flag   = .true.

      do while(flag)
        icount = icount + 1
        nit  = icount
        pinp2 = plast

        call psquars(srts, pinp2, costh, sinth,cosphi,sinphi,
     &              betacm, betarf, blrf1,blrf2,rho0,rmass,
     &              u3,u4,psq,rhap1, rhap2)


        p2    = psq
        if(p2 .ge. 0.0) then
          pnew = sqrt(p2)
          px = pnew*cosphi*sinth
          py = pnew*sinphi*sinth
          pz = pnew*costh
          negp =.false.
          delta = abs(pnew - plast)
        else
          negp = .true.
          pnew = pstart + float(icount)*0.0001
          delta = 10000.0*epsilon
          px = 0.0
          py = 0.0
          pz = 0.0
        end if

        if(delta .lt. epsilon) then
          flag = .false.
        else
          plast = pnew
        end if

        if(icount. gt. icmax) then
c           write(*,*)'s state iteration error '
c           write(*,*) srts, pnew
          flag = .false.
          negp = .true.
        end if
      end do

      if(negp) then
        testflag = .false.
        un1 = 0.0
        un2 = 0.0
      else if(.not. negp) then
        testflag = .true.
        un1 = u3
        un2 = u4
        prf(1) = px
        prf(2) = py
        prf(3) = pz
      end if

      return
      end

*----------------------------------------------------------------------*
*                     do the iteration                                 *
      subroutine psquars(srts, pinp2, costh, sinth,cosphi,sinphi,
     &                  betacm, betarf, blrf1,blrf2,rho0,rmass,
     &                  u3,u4,psq, rhap1, rhap2)



      implicit none

      integer npot
      parameter(npot = 1)

      real*8    srts, pinp2, mass3, mass4, rho0, rho3, rho4
      real*8    s, u3, u4 , pin,  potanal, v3, v4,rmass
      real*8    paux(1:3), px, py, pz,  en3, en4 ,psq
      real*8    betlrfx(1:2), betlrfy(1:2), betlrfz(1:2)
      real*8    betacmx, betacmy, betacmz, betarfx, betarfy, betarfz
      real*8    costh, sinth
      real*8    p2help , pabss, cosphi, sinphi
      real*8    veccm4, veccm3
      real*8    blrf1(1:4), blrf2(1:4), betacm(1:3), betarf(1:3)
      integer nit
      real*8    invv, invn, delta
      logical test
      real*8    rhap1(1:3), rhap2(1:3), plrf(1:3)
      delta = 1.0e-03
*--------------------------------------------------------------------*
*     store the variables in order to be compatible to the old       *
*     version used in momiter                                        *

      betacmx  = betacm(1)
      betacmy  = betacm(2)
      betacmz  = betacm(3)

      betarfx  = betarf(1)
      betarfy  = betarf(2)
      betarfz  = betarf(3)

      betlrfx(1) =  blrf1(1)
      betlrfy(1) =  blrf1(2)
      betlrfz(1) =  blrf1(3)


      betlrfx(2) = blrf2(1)
      betlrfy(2) = blrf2(2)
      betlrfz(2) = blrf2(3)

      rho3     = blrf1(4)
      rho4     = blrf2(4)

      mass3    = rmass
      mass4    = rmass


*--------------------------------------------------------------------*
*         build explicit vector - needed for boost back to the LRF   *


      paux(1) = pinp2*cosphi*sinth
      paux(2) = pinp2*sinphi*sinth
      paux(3) = pinp2*costh

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
          if(px**2+py**2+pz**2.gt.en3**2) then
            write(*,*) "hiba winkel5a lorentz, negative mass",
     &                px,py,pz,en3
c            stop
          end if
          call lorentz(-betarfx ,-betarfy ,-betarfz ,px,py,pz,en3)
          if(px**2+py**2+pz**2.gt.en3**2) then
            write(*,*) "hiba winkel6 lorentz, negative mass",
     &                px,py,pz,en3
c            stop
          end if
          call lorentz(-betacmx ,-betacmy ,-betacmz ,px,py,pz,en3)
          if(px**2+py**2+pz**2.gt.en3**2) then
            write(*,*) "hiba winkel7 lorentz, negative mass",
     &                px,py,pz,en3
c            stop
          end if
          call lorentz(betlrfx(1),betlrfy(1),betlrfz(1),px,py,pz,en3)
          pin   = sqrt(px**2 + py**2 + pz**2)
          pabss = pin
          plrf(1) = px
          plrf(2) = py
          plrf(3) = pz

          v3    = potanal(rho0,rho3,plrf,1,rhap1)
          u3    = -mass3 + sqrt(mass3**2 +
     +           2.0*sqrt(pabss**2+mass3**2)*v3 + v3**2)
          invn = (mass3+u3)**2
          if(abs(invn-invv).lt.delta)then
            test = .false.
          end if
          if(nit.gt.200) then
            write(*,*)'problems in momiter 1 ',nit, mass3,invv,invn
            write(*,*)u3,en3,paux
            stop
          end if
        end do


*        the corresponding vectorpot in the rest frame
        pin    = sqrt(paux(1)**2 + paux(2)**2 + paux(3)**2)
        veccm3 = -sqrt(pin**2 + mass3**2)+ sqrt(pin**2 + mass3**2 +
     +                2.0*mass3*u3 + u3**2)

        nit = 0
        test = .true.
        do while(test)
          nit = nit + 1

*      evaluate the abs. value of momentum of particle 4 in the LRF  *
          px  = -paux(1)
          py  = -paux(2)
          pz  = -paux(3)
          en4  = sqrt((mass4+u4)**2 + px**2 + py**2 + pz**2)
          invv = en4**2 - px**2 - py**2 -pz**2
          if(px**2+py**2+pz**2.gt.en4**2) then
            write(*,*) "hiba winkel8 lorentz, negative mass",
     &                px,py,pz,en4
c            stop
          end if
          call lorentz(-betarfx ,-betarfy ,-betarfz ,px, py, pz, en4)
          if(px**2+py**2+pz**2.gt.en4**2) then
            write(*,*) "hiba winkel9 lorentz, negative mass",
     &                px,py,pz,en4
            stop
          end if
          call lorentz(-betacmx ,-betacmy ,-betacmz ,px, py, pz, en4)
          if(px**2+py**2+pz**2.gt.en4**2) then
            write(*,*) "hiba winkel10 lorentz, negative mass",
     &                px,py,pz,en4
            stop
          end if
          call lorentz(betlrfx(2),betlrfy(2),betlrfz(2),px, py, pz,en4)
          pin = sqrt(px**2 + py**2 + pz**2)
          pabss = pin
          plrf(1) = px
          plrf(2) = py
          plrf(3) = pz

           v4    = potanal(rho0,rho4,plrf,1,rhap2)
           u4    = -mass4 + sqrt(mass4**2 +
     +           2.0*sqrt(pabss**2+mass4**2)*v4 + v4**2)
          invn = (mass4+u4)**2
          if(abs(invn-invv).lt.delta)then
            test = .false.
          end if
          if(nit.gt.200) then
            write(*,*)'problems in momiter 1 ',nit, mass4,invv,invn
            write(*,*)u4,en4,paux
            stop

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

      return
      end





