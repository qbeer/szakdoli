************************************************************************
*                                                                      *
      subroutine kinema(massta,masspr,mstapr,msprpr,mpion,elab,b,ipou,
     &          icoll,insys,rdist,radta,radpr,rmax,rxta,rzta,rxpr,rzpr,
     &       yref,yta,ypr,pxta,pzta,pxpr,pzpr,gammta,gammpr,betata,
     &          betapr, surfta, surfpr)
*                                                                      *
*       purpose:     providing initial kinematics for the heavy-ion col*
*                                                                      *
*       variables:   (all input)                                       *
*         massta  - mass of the target                       (integer) *
*         masspr  - mass of the projectile                   (integer) *
*         mstapr  - charge of the target                     (integer) *
*         msprpr  - charge of the projectile                 (integer) *
*         elab    - bombarding energy in gev                    (real) *
*         b       - impact parameter  (fm)                      (real) *
*         ipou    - 0-> no coulomb; 1-> with coulomb;        (integer) *
*         icoll   - 0-> buu; 1-> vlasov; -1-> cascade;       (integer) *
*         insys   - 0-> lab. system; 1-> cm system;          (integer) *
*         rdist   - initial distance of projectile and target   (real) *
*        output                                                        *
*         radta(pr) radius of target(projectile);               (real) *
*         rmax    - initial distance of projectile and target   (real) *
*         rxpr-rzta-initial displacement of nuclei (fm)         (real) *
*         yref    - initial rapidity of the target              (real) *
*         pxpr-pzta-initial average mom. of particles(gev/c)    (real) *
*         gammta  - gamma factor of target                      (real) *
*         gammpr  - gamma factor of projectile                  (real) *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      real*8 radta,radpr,rmax,elab,b,rdist,rxta,rzta,aas
      real*8 rxpr,rzpr,yref,yta,ypr,pxta,pzta,pxpr,pzpr,gammta,gammpr
      real*8 betata,betapr, surfta, surfpr,cammpr,cammta,cetapr,cetata
      real*8 cost,cpr,cta,cypr,czpr,cyta,czta,betacm,denspr,densta,eccm
      real*8 aas1,aat1,aat2,bbs,epr,eta,gamacm,pccf,sint,theta,zeroz
      integer massta,masspr,mstapr,msprpr,mpion,ipou,icoll,insys
*
*-----------------------------------------------------------------------
      radta  = 1.124 * float(massta)**(1./3.)
      radpr  = 1.124 * float(abs(masspr))**(1./3.)

      write(*,*)'vor nuclfit',radta,radpr

      call  nuclfit(massta,mstapr,radta,surfta,densta)
      if(abs(masspr).gt.0) then
        call  nuclfit(abs(masspr),abs(msprpr),radpr,surfpr,denspr)
      else
        radpr = 0.0
        surfpr = 0.0
        denspr = 0.0
      end if


!       write(20,*)'radius of target     = ', radta
!       write(20,*)'surf   of target     = ', surfta
!       write(20,*)'dens0  of target     = ', densta
!
!
!
!       write(20,*)'radius of projectile     = ', radpr
!       write(20,*)'surf   of projectile    = ', surfpr
!       write(20,*)'dens0  of projectile     = ', denspr



      rmax   = rdist + radta + radpr

      if(ipou.eq.1)then
      rmax   = sqrt(b**2 + rmax**2)
      end if
c      rmax   = sqrt(rmax**2 - b**2)
      rzta   = rmax*float(abs(masspr))/float(massta+abs(masspr))
      rzpr   =-rmax*float(massta)/float(massta+masspr)
      rxta   =-b/2.0
      rxpr   = b/2.0
*
*   if the calculation in lab.,then the system is shifted to left the
*   most to be in the region where the density is calculated (20 fm)
*     such a way the system stay the longest in the region where the
*      density and pauli is calculated
*
      zeroz  = 0.0
      if(insys.eq.0) zeroz = rmax*float(masspr)/float(massta+masspr)
     &                       + 20.0-radpr-rmax
      if(insys.eq.0 .and. masspr.le.1) then
        rzta=0.0
        rxta=0.0
        rzpr=-rmax
        rxpr=b
        zeroz=0.0
      end if
      if(masspr.eq.0 .and. mpion.eq.0) then
      rmax   = 0.0
      rzta   = 0.0
      rzpr   = 0.0
      rxta   = 0.0
      rxpr   = 0.0
      zeroz  = 0.0
      b      = 0.0
      end if
      
*
*  relativistic kinematics
*
*  1) labsystem
*
        eta    = float(massta) * rmass
        pzta   = 0.0
        betata = 0.0
        gammta = 1.0
        yta    = 0.0
*
        if(masspr.ne.0) then
          epr    = float(abs(masspr)) * (rmass + elab)
          pzpr   = sqrt( epr**2 - (rmass * float(abs(masspr)))**2 )
          betapr = pzpr / epr
          gammpr = 1.0 / sqrt( 1.0 - betapr**2 )
          ypr    = 0.5 * log( (epr+pzpr)/(epr-pzpr) )
        else if(masspr.eq.0 .and. mpion.ne.0 .and.mpion.le.2) then
          epr    = pmass + elab
          pzpr   = sqrt( epr**2 - pmass**2 )
          betapr = pzpr / epr
          gammpr = 1.0 / sqrt( 1.0 - betapr**2 )
          ypr    = 0.5 * log( (epr+pzpr)/(epr-pzpr) )
        else
          epr=0.0
          pzpr=0.0
          betapr=0.0
          gammpr=1.0
          return
        end if
*-----------------------------------------------------------------------
*
        betacm=pzpr/(eta+epr)
        gamacm=1./sqrt(1.-betacm**2)
*
*-----------------------------------------------------------------------
*
*  2) c.m. system
*
        czta   = pzta+betacm*gamacm*(gamacm/(1.+gamacm)*pzta*betacm-eta)
        cta    = gamacm*(eta-betacm*pzta)
        cetata = czta / cta
        cammta = 1.0 / sqrt( 1.0 - cetata**2 )
        cyta   = 0.5 * log( (cta+czta)/(cta-czta) )
*
        if(masspr.ne.0 .or. mpion.ne.0) then
        czpr   = pzpr+betacm*gamacm*(gamacm/(1.+gamacm)*pzpr*betacm-epr)
        cpr    = gamacm*(epr-betacm*pzpr)
        cetapr = czpr/ cpr
        cypr   = 0.5 * log( (cpr+czpr)/(cpr-czpr) )
        else
        czpr   = 0.0
        cpr    = 0.0
        cetapr = 0.0
        end if
        cammpr = 1.0 / sqrt( 1.0 - cetapr**2 )
*-----------------------------------------------------------------------
*
*  3) correction of coulomb trajectory (non-relativistic)
*
        write(*,*)'vor coulomb'
        pxta   = 0.0
        pxpr   = 0.0
        if((ipou.eq.1).and.(mstapr.ne.0.and.msprpr.ne.0).
     &      and.(icoll.ne.-1) .and. (mpion.eq.0)) then
        eccm   = cta+cpr-float(massta+abs(masspr))*rmass
        write(*,*)'vor kiritical '
        pccf   = sqrt(1.-float(mstapr*abs(msprpr))*.00144/eccm/rmax
     &                -(b/rmax)**2)
        write(*,*)'nach kritical'
        if(b.ne.0.0) then
          if(b.gt.rmax-2.0) rmax=b+3.0
          aas    = 2.*eccm*b/float(mstapr*abs(msprpr))/.00144
          bbs    = 1./sqrt(1.+aas**2)
          aas1   = (1.+aas*b/rmax)*bbs
          aat1   = aas1/sqrt(1.-aas1**2)
          aat2   = bbs/sqrt(1.-bbs**2)
          theta  = atan(aat1)-atan(aat2)
          cost   = cos(theta)
          sint   = sin(theta)
        else
          cost   = 1.0
          sint   = 0.0
        end if
*
        czpr   = czpr*( cost*pccf+sint*b/rmax)
        pxpr   = czpr*(-sint*pccf+cost*b/rmax)
        cpr    = sqrt(czpr**2+pxpr**2+(float(abs(masspr))*rmass)**2)
        cetapr = czpr/cpr
        czta   =-czpr
        pxta   =-pxpr
        cta    = sqrt(czta**2+pxta**2+(float(massta)*rmass)**2)
        cetata = czta/cta
        cammta = 1.0/sqrt(1.0-cetata**2)
*
        pzta   = czta+betacm*gamacm*(gamacm/(1.+gamacm)*czta*betacm+cta)
        eta    = gamacm*(cta+betacm*czta)
        pzpr   = czpr+betacm*gamacm*(gamacm/(1.+gamacm)*czpr*betacm+cpr)
        epr    = gamacm*(cpr+betacm*czpr)
        betata = pzta / eta
        gammta = 1.0 / sqrt( 1.0 - betata**2 )
        betapr = pzpr / epr
        gammpr = 1.0 / sqrt( 1.0 - betapr**2 )
*
        rzta   = rmax*cost*float(abs(masspr))/float(massta+abs(masspr))
        rzpr   =-rmax*cost*float(massta)/float(massta+abs(masspr))
        rxpr   = rmax/2.0*sint
        rxta   =-rxpr
        if(insys.eq.0 .and. masspr.le.1) then
          rzta=0.0
          rxta=0.0
          rzpr=-rmax*cost
          rxpr= rmax*sint
        end if
      end if
*


        rzta   = rzta - zeroz
        rzpr   = rzpr - zeroz


*
*-----------------------------------------------------------------------
*   write kinematic parameters
        write(isum,'(/''c:'',
     &          43x,''======= kinematical parameters ========''/)')
        write(isum,'(''c:'',
     &               43x,''1) lab-frame:      target   projectile'')')
        write(isum,'(''c:'',
     &               43x,''   etotal "gev" '',2f11.4)') eta, epr
        write(isum,'(''c:'',
     &               43x,''   pz "gev/c"   '',2f11.4)') pzta, pzpr
        write(isum,'(''c:'',
     &               43x,''   beta         '',2f11.4)') betata, betapr
        write(isum,'(''c:'',
     &               43x,''   gamma        '',2f11.4)') gammta, gammpr
        write(isum,'(''c:'',
     &               43x,''   rapidity     '',2f11.4)') yta, ypr
*
        write(isum,'(''c:'',
     &               43x,''2) c.m.-frame:  '')')
        write(isum,'(''c:'',
     &               43x,''   etotal "gev" '',2f11.4)') cta, cpr
        write(isum,'(''c:'',
     &               43x,''   pz "gev/c"   '',2f11.4)') czta, czpr
        write(isum,'(''c:'',
     &               43x,''   beta         '',2f11.4)') cetata, cetapr
        write(isum,'(''c:'',
     &               43x,''   gamma        '',2f11.4)') cammta, cammpr
        write(isum,'(''c:'',
     &               43x,''   rapidity     '',2f11.4)') cyta, cypr
        write(isum,'(''c:'',
     &               43x,''3) correction of coulomb'')')
        write(isum,'(''c:'',
     &               43x,''   px "gev/c"   '',2f11.4)')  pxta,pxpr
*
      if (insys .eq. 0) then
        write(isum,'(/''c:'',
     &               43x,''==== calculation done in lab-frame ===='')')
      else
        eta=cta
        pzta=czta
        betata=cetata
        gammta=cammta
        epr=cpr
        pzpr=czpr
        betapr=cetapr
        gammpr=cammpr
        ypr   =cypr
        yta   =cyta
        write(isum,'(/''c:'',
     &               43x,''==== calculation done in cm-frame ====='')')
      end if
      yref = - yta
        write(isum,'(''c:'',48x,''impact parameter b = '',
     &               f7.3,'' fm'')') b
        if(ipou.eq.0) then
        write(isum,'(''c:'',48x,''initial distance   = '',
     &               f7.3,'' fm'')') sqrt(rmax**2+b**2)
        else
        write(isum,'(''c:'',48x,''initial distance   = '',
     &               f7.3,'' fm'')') rmax
        end if
*
*-----------------------------------------------------------------------
*   momentum per particle
      pzta = pzta / float(massta)
      pxta = pxta / float(massta)
      if(masspr.ne.0) then
      pzpr = pzpr / float(abs(masspr))
      pxpr = pxpr / float(abs(masspr))
      else if(masspr.eq.0 .and. mpion.eq.0) then
      pzpr=0.0
      end if
*
      write(*,*)'end of kinema'

        write(*,*)'******************************+'
        write(*,*)'test test test test '
        write(*,*)'rzta = ', rzta
        write(*,*)'rzpr = ', rzpr
        write(*,*)'pzta = ', pzta
        write(*,*)'pzpr = ', pzpr
        write(*,*)'pxpr = ', pxpr
        write(*,*)'pxta = ', pxta
        write(*,*)'rxta = ', rxta
        write(*,*)'rxpr = ', rxpr

        write(*,*)'end of test *   ********++'

      return
      end

**********************************************************************
************* fit nuclear properties *********************************


      subroutine nuclfit(ma,mz,mradi,msurf,mdens)
      implicit none
      real*8    msurf, mradi, mdens,rd,ad,rhod,rms,rnp,anp,xmsq,a1
      real*8 z,pi,fpi,an,anz,dnz,a3,x,rmsp,pcor,rc
      integer ma, mz,i
      COMMON/RADPN/rd(3),ad(3),rhod(3),RMS(3)
      dimension rnp(3,2),anp(2,2)
      data rnp/1.2490,-0.5401,-0.9582,
     &         1.2131,-0.4415, 0.8931/
      data anp/0.4899,-0.1236, 0.4686,0.0741/
      data xmsq/0.64d0/

**    index 1 = protons
**    index 2 = neutron
**    index 3 = together

      write(*,*)'in nuclfit ', ma,mz,mradi,msurf,mdens

      a1 = dfloat(ma)
      z  = dfloat(mz)
      pi=atan(1.d0)*4.d0
      fpi=4.d0*pi
      an=a1-z
      rd(3)=0.d0
      ad(3)=0.d0
      anz=z
      dnz=(an-z)/a1
      a3=a1**(1./3.)
      do 10 i=1,2
        if(anz.eq.0)    goto 10
        rd(i) =a3*rnp(1,i)+rnp(2,i)+rnp(3,i)*dnz
        ad(i) =anp(1,i)+anp(2,i)*dnz
        rhod(i)=0.75*anz/(rd(I)**3*pi
     &        *(1.+(pi*ad(i)/rd(i))**2))
c <r**2>
        x=(pi*ad(i)/rd(i))**2
        rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*
     &       fpi*rhod(i)/anz
c correct to point-particle distribution (xmsq is the rms of the nucleon)
        rmsp=rms(i)-xmsq
        pcor=sqrt(rmsp/rms(i))
        rd(i)=pcor*rd(i)
        rms(i)=rmsp
        rhod(i)=0.75*anz/(rd(i)**3*pi
     &        *(1.+(pi*ad(i)/rd(i))**2))
c
        rd(3)=rd(3)+anz*rd(i)
        ad(3)=ad(3)+anz*ad(i)
        anz=an
  10  continue
      RC=sqrt(rms(1)/0.6d0)
      rd(3)=rd(3)/a1
      ad(3)=ad(3)/a1
      rhod(3)=0.75*a1/(rd(3)**3*(1.+(pi*ad(3)/rd(3))**2)*pi)
      i=3
      x=(pi*ad(i)/rd(i))**2
      rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*
     &       fpi*rhod(i)/a1
      do 50 i=1,3
   50 rms(i)=sqrt(rms(i))

      mradi = rd(3)
      msurf = ad(3)
      mdens = rhod(3)

      return
      end

