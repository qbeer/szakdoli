*-----------------------------------------------------------------------
      function bwmes(s,ires,idec, imode,igain,igaout,igam,ishgam)
************************************************************************
*       This function calculates the Breit-Wigner-distribution for     *
*       a resonance.                                                   *
*                        3pi decay should be recalculated via rho pi   *
*        input  :                                                      *
*          s    :  square of inv. energy (dynamical mass squared)      *
*          ires :  what kind of resonance                     (integer)*
*                    1 = rho                                           *
*                    2 = sigma                                         *
*                    3 = omega                                         *
*                    4 = phi                                         *
*                                                                      *
*          idec : dummy here                                           *
*          imode:   specifies the functional form of the BW-dist.      *
*                    1 = nonrelativistic form                          *
*                    2 = relativistic form using m                     *
*                    3 = relativistic form using s                     *
*                                                                      *
*          igain:   gamma in ::  values see igaout                     *
*                                                                      *
*          igaout:   0 = total gamma                                   *
*                                                                      *
*          igam:     specifies the kind of gamma to be used            *
*                    0 = gamma = const                                 *
*                    1 = Monitz (for Delta) + standard gamma (else)    *
*                    2 = without cut off                               *
*                                                                      *
*          ishgam    0 = display bw-dist gamma                         *
*                    1 = display bw-dist nominator squared             *
*                    2 = display bw-dist                               *
*                    3 = display bw-dist normalized the max. to 1      *
*                                                                      *
************************************************************************

      implicit none

      real*8 bwmes
*     parameters needed for kinematics
      real*8 pmass,kmass
      parameter(pmass = 0.138, kmass = 0.496)
      real*8 mrho, gamrho0
      real*8 msig, gamsig0
      real*8 mome, gamome0
      real*8 mphi, gamphi0
      parameter(msig  = 0.800, gamsig0 = 0.800 )
      parameter(mrho  = 0.770, gamrho0 = 0.118 )
      parameter(mome  = 0.782, gamome0 = 0.0085 )
      parameter(mphi  = 1.019, gamphi0 = 0.0043 )

*     parameters for cutoff factor
      real*8   beta2, beta4
      parameter(beta2 = 0.09)
      parameter(beta4 = 0.16)

      integer ires, imode, igam, ishgam, igaout, igain
      integer ang , idec
      real*8    s, resm, resm0,       gam0
      real*8    pm2, fact1, fact2, fact3,   br2pi
      real*8    gam2pi, gamtot, bwdis
      real*8    stest,    fact4
      real*8    kpi2, kpi, kr2, kr, gamin, gamout
      real*8    br3pi,br2k,gam3pi,gam2k

      resm = sqrt(s)

*----------------------------------------------------------------------*
      bwmes  =  .0
      if (ires .gt. 4)  return
*             store particle properties                                *
      if( ires .eq. 1) then
*       rho(770)
        resm0  = mrho
        gam0   = gamrho0
        ang    = 1
        br2pi  = 1.0
        br3pi  = 0.0
        br2K   = 0.0
      else if( ires .eq. 2) then
*       sigma(800)
        resm0  = msig
        gam0   = gamsig0
        ang    = 0
        br2pi  = 1.0
        br3pi  = 0.0
        br2K   = 0.0
      else if( ires .eq. 3) then
*       omega(782)
        resm0  = mome
        gam0   = gamome0
        ang    = 0
        br2pi  = 1.0
        br3pi  = 0.0
        br2K   = 0.0
      else if( ires .eq. 4) then
*       phi(1019)
        resm0  = mphi
        gam0   = gamphi0
        ang    = 1
        br2pi  = 0.0
        br3pi  = 0.15
        br2K   = 0.85
      endif
*                                                                     *
      gam2pi = 0.0
      gam3pi = 0.0
      gam2k  = 0.0

c      write(*,*)'bwmes 0',ires,resm,gam2pi,gam3pi,gam2k
      if(igam.eq.0.and.(br2pi.gt.0.000001.or.br3pi.gt.0.000001)) then
        gam2pi = gam0 * br2pi
        gam3pi = gam0 * br3pi
      else
        stest = resm - 2.* pmass
        if( stest .lt. 0.0) then
          write(*,*)"warning from bwmes: 2-pi decay not possible"
          gam2pi = 0.0
        else
          pm2   = pmass**2
          kpi2 = 0.25*s-pm2
          kr2  = 0.25*resm0**2-pm2
          kpi  = sqrt(kpi2)
          kr   = sqrt(kr2)
          fact1 = (kpi/kr)**(2.0*ang+1.0)
          fact2 = resm0/resm
          fact3 = ((beta2 + kr2)/(beta2 + kpi2))**(ang+1.0)
          if(igam.eq.1) gam2pi= fact1 * fact2 * fact3 * gam0
          if(igam.eq.2) gam2pi= fact1 * fact2 * gam0
        end if
        if(resm-3.*pmass.gt.0.0) gam3pi = br3pi*gam2pi
        gam2pi = br2pi * gam2pi
      end if
c      write(*,*)'bwmes1',ires,resm,gam2pi,gam3pi,gam2k

      if(ires.eq.4.and.br2k.gt.0.000001) then
        if(igam.eq.0) then
          gam2k = gam0 * br2k
        else
          stest = resm - 2.* kmass
          if( stest .lt. 0.0) then
            write(*,*)"warning bwmes:2K decay not possible",resm,kmass
            gam2k = 0.0
          else
            pm2   = kmass**2
            kpi2 = 0.25*s-pm2
            kr2  = 0.25*resm0**2-pm2
            kpi  = sqrt(kpi2)
            kr   = sqrt(kr2)
            fact1 = (kpi/kr)**(2.0*ang+1.0)
            fact2 = resm0/resm
            fact3 = ((beta2 + kr2)/(beta2 + kpi2))**(ang+1.0)
            if(igam.eq.1) gam2k= fact1 * fact2 * fact3 * gam0
            if(igam.eq.2) gam2k= fact1 * fact2 * gam0
          end if
          gam2k = br2k * gam2k
        end if
      end if

*-----------------------------------------------------------------------*
*           total gammma                                                *
c      write(*,*)'bwmes2',ires,resm,gam2pi,gam3pi,gam2k
      gamtot = gam2pi+gam3pi+gam2k
*                                                                       *
*-----------------------------------------------------------------------*
*           gammma to be used in nominator of bw-dist .                 *

      if(igain .eq. 0) then
         gamin  = gamtot
      else if(igain .eq. 1) then
         gamin  = gam2pi
      else if(igain .eq. 2) then
         gamin  = gam3pi
      else if(igain .eq. 3) then
         gamin  = gam2k
      end if

      if(igaout .eq. 0) then
         gamout  = gamtot
      else if(igaout .eq. 1) then
         gamout  = gam2pi
      else if(igaout .eq. 2) then
         gamout  = gam3pi
      else if(igaout .eq. 3) then
         gamout  = gam2k
      else
        write(*,*)"wrong choice of igaout !!!!!"
      end if

*                                                                      *
      if(ishgam .eq. 0) then
        if(igain .eq. igaout) then
          bwmes = gamin
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
        bwmes= fact1/(fact3+fact4)
        bwdis = fact2/(fact3+fact4)
      else if(imode .eq. 2) then
*      relativistic version using the resonance mass
        fact1 = resm0**2 * gamin*gamout
        fact2 = sqrt(fact1)
        fact3 = (s - resm0**2)**2
        fact4 = resm0**2*gamtot*gamtot
        bwmes= fact1/(fact3 + fact4)
        bwdis = fact2/(fact3+fact4)
      else if(imode .eq. 3) then
*      relativistic version using s
        fact1 = s * gamin*gamout
        fact2 = sqrt(fact1)
        fact3 = (s - resm0**2)**2
        fact4 = s * gamtot*gamtot
        bwmes= fact1/(fact3+fact4)
        bwdis = fact2/(fact3+fact4)
      end if
      if(ishgam .eq. 2) bwmes = bwdis
      if(ishgam .eq. 3) bwmes = bwdis*gam0*resm

*----------------------------------------------------------------------*
*
      return
      end
