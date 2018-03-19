************************************************************************
      function bwdist(s,ires,idec,imode,igain,igaout,igam,ishgam,
     &     ratio)
************************************************************************
*       This function calculates the Breit-Wigner-distribution for     *
*       a resonance.                                                   *
*                                                                      *
*        input  :                                                      *
*          s    :  square of inv. energy (dynamical mass squared)      *
*          ires :  what kind of resonance                     (integer)*
*                    1 = Delta(1232)                                   *
*                    2 = N(1440)                                       *
*                    3 = N(1520)                                       *
*                    4 = N(1535)                                       *
*                    5=  N(1650)                                       *
*                    6=  N(1675)                                       *
*                    7=  N(1680)                                       *
*                    8=  N(1720)                                       *
*                    9=  Delta(1620)                                   *
*                   10=  Delta(1700)                                   *
*          idec :   what kind of decay is allowed for the resonace     *
*                     0 = without 2-pion decay (nicht sinnvoll!!!!)    *
*                     1 = with 2-pion decay                            *
*                     2 = with 2-pion decay (constant width)           *
*                                                                      *
*          imode:   specifies the functional form of the BW-dist.      *
*                    1 = nonrelativistic form                          *
*                    2 = relativistic form using m                     *
*                    3 = relativistic form using s                     *
*                                                                      *
*          igain:   gamma in ::  values see igaout                     *
*                                                                      *
*          igaout:   0 = total gamma                                   *
*                    1 = 1-pi-decay width                              *
*                    2 = eta-decay-width                               *
*                    3 = 2-pi-decay-width: total                       *
*                    4 = 2-pi-decay-width: Delta pi                    *
*                    5 = 2-pi-decay-width: N(1440) pi                  *
*                    6 = 2-pi-decay-width: N rho                       *
*                    7 = 2-pi-decay-width: N sigma                     *
*                    8 = Lambda K-decay-width                          *
*                    9 = Sigma K-decay-width                           *
*                    10= omega-decay-width                             *
*                                                                      *
*          igam:     specifies the kind of gamma to be used            *
*                    0 = gamma = const                                 *
*                    1 = Monitz (for Delta) + standard gamma (else)    *
*                    2 = Kitazoe (for Delta) + standard gamma (else)   *
*                    3 = standard gamma (mosel phys. rep ) rossz fit   *
*                    4 = Manley gamma                                  *
*                                                                      *
*          ishgam    0 = display bw-dist gamma                         *
*                    1 = display bw-dist nominator squared             *
*                    2 = display bw-dist                               *
*                    3 = display bw-dist normalized the max. to 1      *
*                                                                      *
*          ratio     branching ratios                                  *
*                                                                      *
*          ibreite   0 = Vakuumbreiten                                 *
*                    1 = Mediumbreiten                                 *
************************************************************************

      implicit none

      integer nres
      parameter(nres = 24)

      real*8  resprop1(1:nres,1:11)
      integer resprop2(1:nres,3)
      real*8  m2d16pi(1:nres,1:2)
      common /resproperties / resprop1, resprop2, m2d16pi

      real*8    dimimass0
      parameter(dimimass0= 1.076)
c      include "resdata1"


      real*8 bwdist,s
c      real*8 delpi, n14pi, nrho, nsig
c      external delpi, n14pi, nrho, nsig
*     parameters needed for kinematics
      real*8 msig, gamsig0, mrho, gamrho0, rmass, pmass, emass, omass
      real*8 kmass, lmass, simas
      parameter(rmass = 0.9383)
      parameter(pmass = 0.138)
      parameter(emass = 0.548, omass = 0.782)
      parameter(kmass = 0.494, lmass = 1.116,  simas = 1.189)
      parameter(msig = 0.800 , gamsig0 = 0.800)
      parameter(mrho = 0.770 , gamrho0 = 0.118 )

*     parameters for cutoff factor
      real*8 beta0, beta1
      common/betacutoff/ beta0, beta1
      real*8 beta2, beta4
c      parameter(beta2 = 0.09)
c      parameter(beta4 = 0.16)

      integer ires, imode, igam, ishgam, igaout, igain
      integer idec, ang
      real*8 resm, resm0, gam0
      real*8 rm2 , pm2, fact1, fact2, fact3, br1pi, br2pi
      real*8 gam1pi, gam2pi, gamtot, bwdis
      real*8 stest, vs, vsr, cutoff, fact4
      real*8 breta, gameta, q2, kpi2, kpi, kr2, kr, gamin, gamout
      real*8 brdpi, brnro, brnsi, brnrp, brlak, brsik, brnom

      real*8 intnsig, intnrho, intnsi0, intnrh0, intn14p0, intn14pi

      real*8 intdelp0, intdelpi, gamnsi, gamnrp, gamnro, gamdpi
      real*8 hatmin, hatmax, stot, check, gamlak, gamsik, gamnom
      integer i

      real*8 gamdel0, betacut2
      real*8 gamn140, br1pi14, br2pi14,smass,dmass, test1, test0

      real*8 ratio(9)


      integer index
c      real*8 iintdelp0, iintdelpi
c      real*8 iintnsig, iintnrho, iintnsi0, iintnrh0, iintn14p0, iintn14pi

      real*8 delint
      integer nsrts
      parameter(nsrts = 4000)
      parameter(delint = 0.001)
      real*8 intresm0(1:nres,1:4)
      real*8 intresm(0:nsrts,1:4)
      common /bwfelder/ intresm0, intresm


c      common /param/ stot, betacut2
c      common /deldata/ dmass, gamdel0
c      common /n1440data/ smass, gamn140, br1pi14, br2pi14
c      include"resdata1"

      beta0=0.1325
      beta1=0.2025
      beta2=beta0
      beta4=beta1


      dmass   = resprop1(1,1)
      smass   = resprop1(2,1)
      gamdel0 = resprop1(1,2)
      gamn140 = resprop1(2,2)
      br1pi14 = resprop1(2,3)
      br2pi14 = resprop1(2,5)+resprop1(2,6)+resprop1(2,7)+resprop1(2,8)
c      betacut2 = beta2
c
*       det. mass of resonance
      resm = sqrt(s)

*----------------------------------------------------------------------*
*             store particle properties                                *

      resm0  = resprop1(ires,1)
      gam0   = resprop1(ires,2)
      ang    = resprop2(ires,2)
      br1pi  = resprop1(ires,3)
      breta  = resprop1(ires,4)
      brnsi  = resprop1(ires,5)
      brnro  = resprop1(ires,6)
      brdpi  = resprop1(ires,7)
      brnrp  = resprop1(ires,8)
      brsik  = resprop1(ires,9)
      brlak  = resprop1(ires,10)
      brnom  = resprop1(ires,11)
      br2pi  = brdpi+brnrp+brnro+brnsi
c      write(*,*)'bwstart',ires,resm,resm0,gam0,ang,br1pi,breta,br2pi
      if(idec.eq.0) then
        br1pi = 1.0 - breta
        brdpi = 0.0
        brnrp = 0.0
        brnro = 0.0
        brnsi = 0.0
        br2pi = 0.0
      end if
*                                                                     *
      if(ishgam.eq.0) then
        if(igain.ne.0 .and. igaout.ne.0) then
          if(igain.ne.1 .and. igaout.ne.1) br1pi  = 0.0
          if(igain.ne.2 .and. igaout.ne.2) breta  = 0.0
          if(igain.ne.3 .and. igaout.ne.3) brnsi  = 0.0
          if(igain.ne.4 .and. igaout.ne.4) brnro  = 0.0
          if(igain.ne.5 .and. igaout.ne.5) brdpi  = 0.0
          if(igain.ne.6 .and. igaout.ne.6) brnrp  = 0.0
          if(igain.ne.7 .and. igaout.ne.7) brsik  = 0.0
          if(igain.ne.8 .and. igaout.ne.8) brlak  = 0.0
          if(igain.ne.9 .and. igaout.ne.9) brnom  = 0.0
        endif
      endif
      br2pi  = brdpi + brnrp + brnro + brnsi
*                                                                     *
*---------------------------------------------------------------------*
*       determine gamma according to the igam choice                  *

      gam1pi = 0.0
      gameta = 0.0
      gam2pi = 0.0
      gamnsi = 0.0
      gamnro = 0.0
      gamdpi = 0.0
      gamnrp = 0.0
      gamsik = 0.0
      gamlak = 0.0
      gamnom = 0.0
      if(igam .eq. 0) then
*         use constant gamma
          gam1pi = gam0*br1pi
          gam2pi = gam0*br2pi
          gameta = gam0*breta
          gamdpi = gam0*brdpi
          gamnro = gam0*brnro
          gamnsi = gam0*brnsi
          gamnrp = gam0*brnrp
          gamlak = gam0*brlak
          gamsik = gam0*brsik
          gamnom = gam0*brnom

*       momentum dep. parametrization for gamma
      else if(igam .gt. 0) then
        if(br1pi .gt. 1.d-5) then
          stest = resm - rmass - pmass
          if( stest .lt. 0.0) then
c            write(*,*)"warning from bwdist: 1-pi decay not possible",
c     &           resm
            gam1pi = 0.0
          else
            rm2   = rmass**2
            pm2   = pmass**2
            kpi2 = (s-(rmass+pmass)**2)*(s-(rmass-pmass)**2)/(4.0*s)
            kr2  = 0.25*(resm0**2-rm2+pm2)**2/resm0**2-pm2
            kpi  = sqrt(kpi2)
            kr   = sqrt(kr2)
            if((igam .eq. 1) .or. (igam/2.eq.1 .and. ires.eq.1)) then
*             the Moniz-parameterization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta2 + kr2)/(beta2 + kpi2))**(ang+1.0)
              gam1pi= fact1 * fact2 * fact3 * gam0*br1pi
            else if((igam .eq. 2) .and. (ires .ne. 1) ) then
*             the Manley-parametrization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta4 + kr2)/(beta4 + kpi2))**ang
              gam1pi= fact1 * fact2 * fact3 * gam0*br1pi
            else if((igam .eq. 5) .or. (igam.eq.3.and.ires.ne.1)) then
ccc           igam=3 nem ad jo fittet talan a fact2 hianyzik
*             use standard-momentum dep. parametrization for gamma
              q2   = (resm0 - rmass - pmass)**2 + gam0**2/4.0
              if(ires.eq.4) q2 = 0.25
              gam1pi = gam0*br1pi*(kpi/kr)**(2*ang+1)*
     +                  ((kr2+q2)/(kpi2+q2))**(ang+1.0)
*              if(ires.eq.4) gam1pi=gam0*br1pi*(kpi/kr)**(2*ang+1)
            else if(igam .eq. 4) then
*             the Manley-parametrization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta4 + kr2)/(beta4 + kpi2))**ang
              gam1pi= fact1 * fact2 * fact3 * gam0*br1pi
            end if
          end if
        end if
c        write(*,*)'bw1',gam1pi,igam,beta2,kpi,kr,resm0,resm,gam0,br1pi

        if(breta.gt.1.d-5) then
*       do eta-width
          stest = resm - rmass - emass
          if( stest .le. 0.0) then
*           write(*,*)"warning from bwdist: eta decay not possible"
            gameta = 0.0
          else
            kpi2=(s-(rmass+emass)**2)*(s-(rmass-emass)**2)/(4.0*s)
            kr2 =(resm0**2-(rmass+emass)**2)*(resm0**2-(rmass-emass)**2)
     +            /(4.0*resm0**2)
c            if(kpi2.le.0.0 .or. kr2.le.0.0)write(*,*)'hibaeta',kpi2,kr2
            kpi = sqrt(kpi2)
            kr  = sqrt(abs(kr2))

            if(igam .le. 3) then
*             use standard-momentum dep. parametrization for gamma
              q2  =(resm0 - rmass - emass)**2 + gam0**2/4.0
*Achtung Aenderung
              if(ires.eq.4) q2 = 0.25
*sonst macht cutoff keinen Sinn
              gameta = gam0*breta*(kpi/kr)**(2.0*ang+1.0)*
     +                    ((kr2+q2)/(kpi2+q2))**(ang+1.0)
*              if(ires.eq.4) gameta=gam0*breta*(kpi/kr)**(2.0*ang+1.0)
            else if(igam .eq. 4) then
*             the Manley-parametrization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta4 + kr2)/(beta4 + kpi2))**ang
              gameta= fact1 * fact2 * fact3 * gam0*breta
            end if
          end if
        end if

        if(brlak.gt.1.d-5) then
*       do Lambda K-width
          stest = resm - lmass - kmass
          if( stest .lt. 0.0) then
*           write(*,*)"warning from bwdist: Lambda K decay not possible"
            gamlak = 0.0
          else
            kpi2=(s-(lmass+kmass)**2)*(s-(lmass-kmass)**2)/(4.0*s)
            kr2 =(resm0**2-(rmass+pmass)**2)*(resm0**2-(rmass-pmass)**2)
     +            /(4.0*resm0**2)
            kpi = sqrt(kpi2)
            kr  = sqrt(kr2)

            if(igam .le. 3) then
*             use standard-momentum dep. parametrization for gamma
              q2  =(resm0 - lmass - kmass)**2 + gam0**2/4.0
              gamlak = gam0*brlak*(kpi/kr)**(2.0*ang+1.0)*
     +                    ((kr2+q2)/(kpi2+q2))**(ang+1.0)
            else if(igam .eq. 4) then
*             the Manley-parametrization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta4 + kr2)/(beta4 + kpi2))**ang
              gamlak= fact1 * fact2 * fact3 * gam0*brlak
            end if
          end if
        end if

        if(brsik.gt.1.d-5) then
*       do Sigma K-width
          stest = resm - simas - kmass
          if( stest .lt. 0.0) then
*           write(*,*)"warning from bwdist: Sigma K decay not possible"
            gamlak = 0.0
          else
            kpi2=(s-(simas+kmass)**2)*(s-(simas-kmass)**2)/(4.0*s)
            kr2 =(resm0**2-(rmass+pmass)**2)*(resm0**2-(rmass-pmass)**2)
     +            /(4.0*resm0**2)
            kpi = sqrt(kpi2)
            kr  = sqrt(kr2)

            if(igam .le. 3) then
*             use standard-momentum dep. parametrization for gamma
              q2  =(resm0 - simas - kmass)**2 + gam0**2/4.0
              gamsik = gam0*brsik*(kpi/kr)**(2.0*ang+1.0)*
     +                    ((kr2+q2)/(kpi2+q2))**(ang+1.0)
            else if(igam .eq. 4) then
*             the Manley-parametrization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta4 + kr2)/(beta4 + kpi2))**ang
              gamsik= fact1 * fact2 * fact3 * gam0*brsik
            end if
          end if
        end if

        if(brnom.gt.1.d-5) then
*       do omega-width
          stest = resm - rmass - omass
          if( stest .le. 0.0) then
*           write(*,*)"warning from bwdist: omega decay not possible"
            gamnom = 0.0
          else
            kpi2=(s-(rmass+omass)**2)*(s-(rmass-omass)**2)/(4.0*s)
            kr2 =(resm0**2-(rmass+pmass)**2)*(resm0**2-(rmass-pmass)**2)
     +            /(4.0*resm0**2)
            kpi = sqrt(kpi2)
            kr  = sqrt(kr2)

            if(igam .le. 3) then
*             use standard-momentum dep. parametrization for gamma
              q2  =(resm0 - rmass - omass)**2 + gam0**2/4.0
              gamnom = gam0*brnom*(kpi/kr)**(2.0*ang+1.0)*
     +                    ((kr2+q2)/(kpi2+q2))**(ang+1.0)
            else if(igam .eq. 4) then
*             the Manley-parametrization
              fact1 = (kpi/kr)**(2.0*ang+1.0)
              fact2 = resm0/resm
              fact3 = ((beta4 + kr2)/(beta4 + kpi2))**ang
              gamnom= fact1 * fact2 * fact3 * gam0*brnom
            end if
          end if
        end if

        if(br2pi .gt.1.d-5) then
*       calculate the 2-pion-width
          stest = resm - rmass - 2.0 * pmass
          if(stest.gt.1.e-3) then
            if(idec.eq.2) then
              gamdpi = gam0*brdpi
              gamnro = gam0*brnro
              gamnsi = gam0*brnsi
              gamnrp = gam0*brnrp
              gam2pi = gam0*br2pi
            else if(idec.eq.1) then
              vs =1.0/(beta4+(resm  -(rmass+2.*pmass))**2)
              vsr=1.0/(beta4+(resm0 -(rmass+2.*pmass))**2)
*              cutoff= (vs/vsr)**(ang+1.0)
              cutoff=(vs/vsr)**2.0
              hatmin = (rmass+pmass)**2
              if(brdpi.gt.1.e-5) then
c                stot = s
c                hatmax = (resm-pmass)**2
c                CALL QRomb(delpi,hatmin,HATmax,test1)
c                stot = resm0**2
c                hatmax = (resm0-pmass)**2
c                CALL QRomb(delpi,hatmin,HATmax,test0)

                index  = min(nsrts,nint((resm-dimimass0)/delint))
                intdelpi =  sngl(intresm(index,1))
                intdelp0 =  sngl(intresm0(ires,1))
                gamdpi = resm0/resm*intdelpi/intdelp0*gam0*brdpi*cutoff

              endif
              if(brnrp.gt.1.d-5) then
c                stot = s
c                hatmax = (resm-pmass)**2
c                CALL QRomb(n14pi,hatmin,HATmax,test1)
c                stot = resm0**2
c                hatmax = (resm0-pmass)**2
c                CALL QRomb(n14pi,hatmin,HATmax,test0)

                index  = min(nsrts,nint((resm-dimimass0)/delint))
                intn14pi =  sngl(intresm(index,2))
                intn14p0 =  sngl(intresm0(ires,2))
                gamnrp = resm0/resm*intn14pi/intn14p0*gam0*brnrp*cutoff

              endif
              hatmin = (2.0*pmass)**2
              if(brnro.gt.1.d-5) then
c                stot = s
c                hatmax = (resm-rmass)**2
c                CALL QRomb(nrho,hatmin,HATmax,test1)
c                stot = resm0**2
c                hatmax = (resm0-rmass)**2
c                CALL QRomb(nrho,hatmin,HATmax,test0)

                index  = min(nsrts,nint((resm-dimimass0)/delint))
                intnrho =  sngl(intresm(index,3))
                intnrh0 =  sngl(intresm0(ires,3))

                gamnro = resm0/resm*intnrho/intnrh0*gam0*brnro*cutoff
              endif
              if(brnsi.gt.1.d-5) then
c                stot=s
c                hatmax=(resm-rmass)**2
c                call qromb(nsig,hatmin,hatmax,test1)
c                stot=resm0**2
c                hatmax=(resm0-rmass)**2
c                call qromb(nsig,hatmin,hatmax,test0)

                index  = min(nsrts,nint((resm-dimimass0)/delint))

                intnsig =  sngl(intresm(index,4))
                intnsi0 =  sngl(intresm0(ires,4))

                gamnsi=resm0/resm*intnsig/intnsi0*gam0*brnsi*cutoff

*                gamnsi=brnsi*gam0*(sqrt(s)-rmass-2.*pmass)/
*     &               (resm0-rmass-2.*pmass)*cutoff
              endif
              gam2pi = gamdpi+gamnro+gamnsi+gamnrp
            end if
          end if
        end if
      end if

*----------------------------------------------------------------------*
*           total gamma                                                *

        gamtot = gam1pi + gam2pi + gameta + gamlak + gamsik + gamnom
c      write(*,*)'bwdistgam',gamtot,gam1pi,gam2pi,gameta,gamlak,gamsik,
c     &         gamnom,gam0,resm,resm0,ires

        if(gamtot.ne.0) then
          ratio(1)=gam1pi/gamtot
          ratio(2)=gameta/gamtot
          ratio(3)=gamnsi/gamtot
          ratio(4)=gamnro/gamtot
          ratio(5)=gamdpi/gamtot
          ratio(6)=gamnrp/gamtot
          ratio(7)=gamsik/gamtot
          ratio(8)=gamlak/gamtot
          ratio(9)=gamnom/gamtot

********************** check
          check = 0.0
          do i = 1, 9
            check = check + ratio(i)
          end do
          check = abs(check-1.0)
          if(check .gt. 1.0e-03) then
            write(*,*)'check warning in bw :  '
            write(*,*)'check = ', check
            write(*,*)'ires = ', ires
            write(*,*)'idec = ', idec
            write(*,*)'resm = ', resm
            write(*,*)'ratio ', ratio
          end if

        end if
*                                                                      *
*----------------------------------------------------------------------*
*           gamma to be used in nominator of bw-dist .                 *

      if(igain .eq. 0) then
         gamin  = gamtot
      else if(igain .eq. 1) then
         gamin  = gam1pi
      else if(igain .eq. 2) then
         gamin  = gameta
      else if(igain .eq. 3) then
         gamin  = gamnsi
      else if(igain .eq. 4) then
         gamin  = gamnro
      else if(igain .eq. 5) then
         gamin  = gamdpi
      else if(igain .eq. 6) then
         gamin  = gamnrp
      else if(igain .eq. 7) then
         gamin  = gamsik
      else if(igain .eq. 8) then
         gamin  = gamlak
      else if(igain .eq. 9) then
         gamin  = gamnom
      else
        write(*,*)"wrong choice of igain !!!!!"
      end if


      if(igaout .eq. 0) then
         gamout  = gamtot
      else if(igaout .eq. 1) then
         gamout  = gam1pi
      else if(igaout .eq. 2) then
         gamout  = gameta
      else if(igaout .eq. 3) then
         gamout  = gamnsi
      else if(igaout .eq. 4) then
         gamout  = gamnro
      else if(igaout .eq. 5) then
         gamout  = gamdpi
      else if(igaout .eq. 6) then
         gamout  = gamnrp
      else if(igaout .eq. 7) then
         gamout = gamsik
      else if(igaout .eq. 8) then
         gamout = gamlak
      else if(igaout .eq. 9) then
         gamout = gamnom
      else
        write(*,*)"wrong choice of igaout !!!!!"
      end if

*                                                                      *
      if(ishgam.eq. 0) then
        if(igain .eq. igaout) then
          bwdist = gamin
          return
        else
          write(*,*) "In order to display a width the parameters: "
          write(*,*) "igain and igaout must have identical values !!"
          stop
        end if
      end if

*----------------------------------------------------------------------*
*            determine Breit-Wigner distribution                       *
      if(imode .eq. 1) then
*      non-relativistic version
        fact1 = 0.25 * gamin*gamout
        fact2 = sqrt(fact1)
        fact3 = (resm - resm0)**2
        fact4 = 0.25 *gamtot*gamtot
        bwdist= fact1/(fact3+fact4)
        bwdis = fact2/(fact3+fact4)
      else if(imode .eq. 2) then
*      relativistic version using the resonance mass
        fact1 = resm0**2 * gamin*gamout
        fact2 = sqrt(fact1)
        fact3 = (resm**2 - resm0**2)**2
        fact4 = resm0**2*gamtot*gamtot
        bwdist= fact1/(fact3 + fact4)
        bwdis = fact2/(fact3+fact4)
      else if(imode .eq. 3) then
*      relativistic version using s
        fact1 = resm**2* gamin*gamout
        fact2 = sqrt(fact1)
        fact3 = (resm**2 - resm0**2)**2
        fact4 = resm*resm * gamtot*gamtot
        bwdist= fact1/(fact3+fact4)
        bwdis = fact2/(fact3+fact4)
      end if

      if(ishgam .eq. 2) bwdist = bwdis
      if(ishgam .eq. 3) bwdist = bwdis*gam0*resm
*----------------------------------------------------------------------*
*
c      write(*,*) bwdist, gamin,gamout,resm,resm0,ang,imode
      return
      end

c**********************************************************************c
      function delpi(s)
      implicit none
      real*8 delpi, s
      real*8 stot, p32, beta2
      real*8 rmass, pmass, emass, dmass, gamdel0
      real*8 gamdels, srt, kpi2, kr2, dm2, fact3
      parameter(rmass = 0.9383)
      parameter(pmass = 0.138 )
      parameter(emass = 0.548 )
      common /deldata/ dmass, gamdel0
      common /param/ stot, beta2
      delpi = 0.0
      dm2 = dmass**2
      srt = dsqrt(s)
      p32= 0.25*(stot-s+pmass**2)**2/stot-pmass**2
      if(p32.lt.1.d-16) return
      kpi2 =max(0.d0,(s-(rmass+pmass)**2)*(s-(rmass-pmass)**2)/(4.d0*s))
      if(kpi2.lt.1.d-16) return
      kr2  = (dm2-(rmass+pmass)**2)*(dm2-(rmass-pmass)**2)/(4.*dm2)
      if(kr2.lt.1.d-16) write(*,*) 'hiba in delpi:kr2',kr2
      fact3 = ((beta2 + kr2)/(beta2 + kpi2))**2
      gamdels= gamdel0*dmass * (dsqrt(kpi2/kr2))**3 * fact3
      delpi=dsqrt(p32)*gamdels/((s-dm2)**2+gamdels**2)
      RETURN
      END

      function n14pi(s)
      implicit none
      real*8 n14pi, s
      real*8 stot,   p32, beta2
      real*8 rmass, pmass, emass, smass, gamn140, br1pi14, br2pi14
      real*8 gamn14, srt, kpi2, kr2, sm2, fact3
      parameter(rmass = 0.9383)
      parameter(pmass = 0.138 )
      parameter(emass = 0.548 )
      common /n1440data/ smass, gamn140, br1pi14, br2pi14
      common /param/ stot, beta2
      n14pi = 0.0
      sm2 = smass**2
      srt = dsqrt(s)
      p32= 0.25*(stot-s+pmass**2)**2/stot-pmass**2
      if(p32.lt.1.d-16) return
      kpi2 =max(0.d0,(s-(rmass+pmass)**2)*(s-(rmass-pmass)**2)/(4.d0*s))
      kr2  = (sm2-(rmass+pmass)**2)*(sm2-(rmass-pmass)**2)/(4.*sm2)
      if(kpi2.lt.1.d-16) return
      fact3 = ((beta2 + kr2)/(beta2 + kpi2))**2
      gamn14 = (dsqrt(kpi2/kr2))**3 * fact3 * gamn140 * br1pi14
      if(srt .gt. rmass+2.*pmass) gamn14=gamn14+br2pi14*gamn140
      n14pi=dsqrt(p32)*gamn14*srt/((s-sm2)**2+s*gamn14**2)
      RETURN
      END

      function nrho(s)
      implicit none
      real*8 nrho, s
      real*8 stot, p32, beta2
      real*8 rmass, pmass, emass, mrho, gamrho0
      real*8 rm2, pm2, gamrhos, kro2, krro2, srt, fact3
      parameter(rmass = 0.9383)
      parameter(pmass = 0.138 )
      parameter(emass = 0.548 )
      parameter(mrho = 0.770 , gamrho0 = 0.15 )
      common /param/ stot, beta2
      nrho = 0.0
      pm2= pmass**2
      rm2= rmass**2
      p32= 0.25*(stot-s+rm2)**2/stot-rm2
      if(p32.lt.1.d-16) return
      srt = dsqrt(s)
      kro2 = max(0.d0,0.25*s - pm2)
      krro2 = 0.25*mrho - pm2
      if(kro2.lt.1.d-16) return
      fact3 = ((beta2 + krro2)/(beta2 + kro2))
      gamrhos = gamrho0 * (kro2/krro2)**1.5 * mrho * fact3
      nrho=dsqrt(p32)*gamrhos/((s-mrho**2)**2+gamrhos**2)
      RETURN
      END

      function nsig(s)
      implicit none
      real*8 nsig, s
      real*8 stot,  p32, beta2
      real*8 rmass, pmass, emass, msig, gamsig0
      real*8 rm2, pm2, gamsigs, kro2, krro2, srt, fact3
      parameter(rmass = 0.9383)
      parameter(pmass = 0.138 )
      parameter(emass = 0.548 )
      parameter(msig = 0.800 , gamsig0 = 0.800)
      common /param/ stot, beta2
      nsig = 0.0
      pm2= pmass**2
      rm2= rmass**2
      p32= 0.25*(stot-s+rm2)**2/stot-rm2
      if(p32.lt.1.d-16) return
      srt = dsqrt(s)
      kro2 = max(0.d0,0.25*s - pm2)
      krro2 = 0.25*msig - pm2
      if(kro2.lt.1.d-16) return
      fact3 = ((beta2 + krro2)/(beta2 + kro2))
      gamsigs = gamsig0 * (kro2/krro2)**0.5 * msig * fact3
      nsig=dsqrt(p32)*gamsigs/((s-msig**2)**2+gamsigs**2)
      RETURN
      END


************************************************************************
      subroutine   bwini
*-----------------------------------------------------------------------
      implicit none

      integer ilauf,nres
      integer nsrts
      parameter(nsrts = 4000)
      parameter(nres = 24)
      real*8 intresm0(1:nres,1:4)
      real*8 intresm(0:nsrts,1:4)
      common /bwfelder/ intresm0, intresm

      real*8 resprop1(1:nres,1:11)
      integer resprop2(1:nres,3)
      real*8 m2d16pi(1:nres,1:2)
      common /resproperties / resprop1, resprop2, m2d16pi
c      include"resdata1"
      real*8 dimimass0
      parameter(dimimass0= 1.076)

      real*8 rmass, pmass, emass, msig,mrho
      parameter(rmass = 0.9383)
      parameter(pmass = 0.138)
      parameter(emass = 0.548)
      parameter(msig = 0.800)
      parameter(mrho = 0.770 )

      real*8 delpi,intdelp0,n14pi,intn14p0,nsig,intnsi0,nrho,intnrh0
      external delpi, n14pi, nrho, nsig

*     parameters for cutoff factor
      real*8 beta0, beta1, beta2, beta4
      common/betacutoff/ beta0, beta1
c      parameter(beta2 = 0.09)
c      parameter(beta4 = 0.16)

      real*8 stot, hatmin, hatmax,betacut2
      real*8 intdelpi, intn14pi, intnrho, intnsig
      integer isrt

      real*8  delint
      parameter(delint = 0.001)

      real*8 resm, resm0, s
      real*8 gamn140, br1pi14, br2pi14,smass,dmass,gamdel0

      common /param/ stot, betacut2
      common /deldata/ dmass, gamdel0
      common /n1440data/ smass, gamn140, br1pi14, br2pi14

      beta0=0.1325
      beta1=0.2025
      beta2=beta0
      beta4=beta1

      dmass   = resprop1(1,1)
      smass   = resprop1(2,1)
      gamdel0 = resprop1(1,2)
      gamn140 = resprop1(2,2)
      br1pi14 = resprop1(2,3)
      br2pi14 = resprop1(2,5)+resprop1(2,6)+resprop1(2,7)+resprop1(2,8)
      betacut2 = beta4

c      write(6,*)'in bwini',dmass,smass,gamdel0,gamn140,br1pi14,br2pi14
*       det. mass of resonance
c      resm = sqrt(s)

*----------------------------------------------------------------------*
*             store particle properties                                *

      do ilauf = 1, nres
        do isrt = 1, 4
          intresm0(ilauf,isrt) = 0.0
        end do
      end do

      do ilauf = 0, nsrts
        do isrt = 1, 4
          intresm(ilauf,isrt) = 0.0
        end do
      end do

c      write(*,*)'vor loop'
      do ilauf = 2, nres

        resm0  = resprop1(ilauf,1)
        stot   = dble(resm0**2)
c        write(*,*)'resm0 =', resm0

        hatmin = dble((rmass+pmass)**2)
        hatmax = dble((resm0-pmass)**2)
c        write(*,*)'delpi0 elott'
        CALL dQRomb(delpi,hatmin,HATmax,intdelp0)
        intresm0(ilauf,1) = intdelp0
c        write(*,*)'delpi0 utan',intdelp0

c        write(*,*)'n14pi0 elott'
        CALL dQRomb(n14pi,hatmin,HATmax,intn14p0)
        intresm0(ilauf,2) = intn14p0
c        write(*,*)'n14pi0 utan',intn14p0

        hatmin = dble((2.0*pmass)**2)
        hatmax = dble((resm0-rmass)**2)
c        write(*,*)'nrho0 elott'
        CALL dQRomb(nrho,hatmin,HATmax,intnrh0)
        intresm0(ilauf,3) = intnrh0
c        write(*,*)'nrho0 utan',intnrh0

c        write(*,*)'nsig0 elott'
        call dqromb(nsig,hatmin,hatmax,intnsi0)
        intresm0(ilauf,4) = intnsi0
c        write(*,*)'nsig0 utan',intnsi0

      end do
      open(88,file='bwinres0.dat')
      do isrt=1,nres
        write(88,'(4e12.4)') (sngl(intresm0(isrt,ilauf)),ilauf=1,4)
      end do
      close(88)
      do ilauf = 0, nsrts
        resm = dimimass0 + dble(ilauf)*delint
        s    = resm**2
        stot = dble(s)

        if(resm.gt.rmass+2.0*pmass) then

          hatmin = dble((rmass+pmass)**2)
          hatmax = dble((resm-pmass)**2)
          CALL dQRomb(n14pi,hatmin,HATmax,intn14pi)
          intresm(ilauf,2) = intn14pi

          CALL dQRomb(delpi,hatmin,HATmax,intdelpi)
          intresm(ilauf,1) = intdelpi

          hatmin = dble((2.0*pmass)**2)
          hatmax = dble((resm-rmass)**2)
          CALL dQRomb(nrho,hatmin,HATmax,intnrho)
          intresm(ilauf,3) = intnrho

          call dqromb(nsig,hatmin,hatmax,intnsig)
          intresm(ilauf,4) = intnsig

        end if
      end do

      open(89,file='bwinresm.dat')
      do isrt=0,nsrts
        write(89,'(4e12.4)') (sngl(intresm(isrt,ilauf)),ilauf=1,4)
      end do
      close(89)

      return
*----------------------------------------------------------------------*

      end
c***********************************************************************
      subroutine dqromb(func,a,b,ss)
*
*      return ss as the integral of the function func from a to b.
*      integration is performed by romberg's method of order 2k,
*      where k=2 is simpson's rule
***********************************************************************
      implicit none
      real*8 func,a,b,ss,s,h,eps,dss
      external func
      integer jmax,jmaxp,k,km,j,l
      parameter(eps=1.d-6,jmax=90,jmaxp=jmax+1,k=10,km=k-1)
*      here eps is the fractional accuracy desired, as determined by
*      the extrapolation error estimate;
*      jmax limits the total number of steps;
*      k is the number of points used in the extrapolation
      dimension s(jmaxp),h(jmaxp)
      h(1)=1.
      do 11 j=1,jmax
        call trapzd(func,a,b,s(j),j)
        if (j.ge.k) then
	  l=j-km
          call polint(h(l),s(l),k,0.d0,ss,dss)
          if (dabs(dss).lt.eps*dabs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      write(*,*)'qrom:too many steps hiba:',abs(dss),eps*abs(ss)
      return
      end
************************************************************************
      subroutine trapzd(func,a,b,s,n)
      implicit none
      real*8 func,a,b,s,del,tnm,x,sum
      integer N,it,j
      save it
C     external func line may be wrong, but without it the compiler cries
c      external func
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
        it=1
      else
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
        it=2*it
      endif
      return
      end
***********************************************************************
      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      real*8 xa,ya,x,y,dy,c,d,dif,dift,ho,hp,w,den
      integer N,nmax,ns,I,m
      parameter (nmax=20)
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0) write(*,*) 'pause in qromb.f'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end

