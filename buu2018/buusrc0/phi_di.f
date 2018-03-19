***********************************************************************
*                                                                      *
      subroutine phi_din
*                                                                      *
*       purpose:    calculating phi_ (and omega)production from:       *
*                                n+n  collision      (entry phi_dbb)   *
*                                n+pi collision      (entry phi_dpi)   *
*--------------------------------------------------------------------  *
*                   initialization by subroutine phi_in                *
*--------------------------------------------------------------------  *
*                   final cross section by entry kaonout               *
*--------------------------------------------------------------------  *
*                                                                      *
*         yref    - -1 * target rapidity in the frame                  *
*         nyk     - number of rapidities               (integer,input) *
*         nqtk    - number of transverse momenta       (integer,input) *
*         nfk     - number of angles                   (integer,input) *
*                                                                      *
*         ika(1,  ) = 1 (kaon+),      =2 (phi)                         *
*                                                                      *
*         ireac: 1-> nn      ->  nn phi                                *
*                2-> nD      ->  nn phi                                *
*                3-> DD      ->  nn phi                                *
*                4-> RR      ->  nn phi                                *
*                5-> n pi    ->  n  phi                                *
*                6-> D pi    ->  n  phi                                *
*                7-> n1440 pi->  n  phi                                *
*                8-> n1520 pi->  n  phi                                *
*                9-> R pi    ->  n  phi                                *
*               10-> N rho   ->  n  phi                                *
*               11-> D rho   ->  n  phi                                *
*               12-> R rho   ->  n  phi                                *
*               13-> pi pi   ->  pi phi                                *
*               14-> pi rho  ->     phi                                *
*               15-> K+K-    ->     phi                                *
*               18-> N ome   ->  n  phi                                *
*               19-> D ome   ->  n  phi                                *
*               20-> R ome   ->  n  phi                                *
*               21-> N eta   ->  n  phi                                *
*               22-> D eta   ->  n  phi                                *
*               23-> R eta   ->  n  phi                                *
*               38-> n eta   ->  n  phi  etacoll (non perturbative eta *
*               39-> D eta   ->  n  phi  etacoll (non perturbative eta *
*               28-> N ome   ->  n  phi  mesphi  (non perturbative ome *
*               29-> D ome   ->  n  phi  mesphi  (non perturbative ome *
*               30-> R ome   ->  n  phi  mesphi  (non perturbative ome *
*               38-> n eta   ->  n  phi  mesphi  (non perturbative eta *
*               39-> D eta   ->  n  phi  mesphi  (non perturbative eta *
*               40-> R eta   ->  n  phi  mesphi  (non perturbative eta *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
      include 'com_pert'
*----------------------------------------------------------------------*
      real*8   phi_b, del, xz, zig,zzig
      integer  inkrun, inkmin,inkrunk
      integer  bin_dens,bin_time,binsqrts
      integer  p_dbb(0:999,0:999),p_dpi(0:999,0:999),p_rho(0:999,0:999)
      integer  p_pipi(0:999,0:999),p_pirho(0:999,0:999)
      integer  c_pn(0:999)
      integer  c_pp(0:999),c_DN(0:999),c_DD(0:999)
      integer  iy,i1,i2,ireac,ix,iz,maxde,id2
      integer  i61,id6,iz1,iz2,id1,id12,id22,id21,j2
      integer  irun,i_pert,ntag,kl,numprodd,numprod
      real*8   gamma,srt, meff1,meff2, xxx,yyy,zzz,distr
      real*8   facreac,em12,em22,tmass2,rmass2,emm2,s,pmax2,pmax,denst
      real*8   szig0, szig, szig_pn, szig1, traf, spr, sig0
      real*8   xx,yy,zz,rr,transf,qqx,qqy,qqz, xxy, xn
      real*8   e3,pbeta,ede,pdx,pdy,pdz, srt0,phase
      real*8   gamm,betax,betay,betaz,q0,qq2,qqabs
      real*8   rn,facphi_di2
      real*8   srtmin,etotal,valkabb,valkapi,xx2,yy2,zz2
      real*8   ranp,pka,eka,enupr,pnupr,factn1,factn2, pphi
      real*8   vx,vy,vz, vxx,vyy,vzz, pxx, tmass, mmass
      real*8   phi_rhoDtab,ds12,As1,As2,As3,As4,As5,As6,As7,Bs1,Bs2,fak
      real*8   Ah1,Bh1,Ch1,Dh1,Eh1,gradm,bwmes
      real*8   em1,em2,px1,py1,pz1,e1,px2,py2,pz2,e2,dx,dy,dz,rsqare
      real*8   p12,p1dr,p2dr,a12,b12,c12,brel,b21,t1,t2,ddlt,qk2,sigkk
c      real*8   bi(0:999)
c      integer  iop,iwp
*----------------------------------------------------------------------*
      parameter (maxde=20)
c      data phi_mass /1.020/
*----------------------------------------------------------------------*
      real*8  p3(3),beta(3)
*----------------------------------------------------------------------*
      save pphi,numprodd,numprod!
!      save p_dbb,p_dpi,p_rho,p_pipi,p_pirho
*----------------------------------------------------------------------*
      pphi = 1./real(pphi_num)
      return

************************************************************************
      entry phi_omegaeta(ede,pdx,pdy,pdz,srt,meff1,meff2,xxx,yyy,zzz,
     &         szig0,i2,id2,id6,iz1,iz2,id1,irun)
* purpose: phi production in omega+N -> N+phi           zm             *
*                     and in omega+D -> N+phi                          *
*       variables:    1 = eta,omega       2 = baryon                       *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*----------------------------------------------------------------------*
*        the  following  call should be in subroutine   pionab.f
*----------------------------------------------------------------------*

      if(iz1+iz2 .gt. 1)                                       return
      if(iz1+iz2 .lt. 0)                                       return
      if(id1.ne.2 .and. id1.ne.5 .and. id1.ne.25)              return !eta omega
c      if(id1.ne.25)                               return  !omega
c      if((id2.ne.1).and.(id2.ne.2))                            return

      pxx =  .0
      if (iphi_pot .ge. 1) then
c        write(*,*) 'phi_omegaeta 2', xxx,yyy,zzz, pxx,pxx,pxx,mmass,
c     1                 vx,vy,vz,vxx,vyy,vz
         call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                      vx,vy,vz,vxx,vyy,vzz)
      else
         mmass = xphimas
      endif
      tmass  = mmass
      srt0  = rmass + tmass
      if(srt.le.srt0)                                             return

c      write(*,*) 'srt is above threshold'
      tmass2= tmass**2
      rmass2= rmass**2
      s     = srt**2
c       q0  : energy of the phi_
      q0      = 0.5 * (s + tmass2 - rmass2) / srt
      if(q0 .le. tmass)                                    return
c      write(*,*) 'phi energy is enough'
      qq2     = q0**2 - tmass2
      qqabs   = sqrt(qq2)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gamm  = ede / srt
      traf  = gamm / (gamm+1.0)
*----------------------------------------------------------------------*
      facreac = 1.0
      xx = srt - srt0
      if(id1.eq.5 .or. id1.eq.25) then
        if(id2.eq.1) then         ! omega+n -> phi+n      henry
          ireac = 18
          szig = sqrt(xx) * 0.05275/(1.-0.8732*meff1+1.467*xx) *
     &        exp(1. - 1.507*meff1 - 1.523*xx +
     &        1.051*meff1**2 + 0.4646*xx**2 - 0.8308*meff1*xx)
          if (iz1 .ne. 0) facreac = 2.0
        else if(id2.eq.2) then    ! omega+D -> phi+n      zm
          ireac = 19
          szig = phi_rhoDtab(meff2, meff1, xx)
c          call  sig_del_rho(meff2, meff1, xx, szig)
          if (iz1 .eq. 0)  then
            facreac = 2.0
          else
            xz = abs (iz2 -0.5) - 1.0
            if (xz .gt. .0)  then
              facreac = 3.0
            else
              facreac = 1.0
            endif
          endif
        else ! omega+R -> phi+n  taken from delta
          ireac = 20
          szig = phi_rhoDtab(meff2, meff1, xx)
c          call  sig_del_rho(meff2, meff1, xx, szig)
          if (iz1 .eq. 0)  then
            facreac = 2.0
          else
            xz = abs (iz2 -0.5) - 1.0
            if (xz .gt. .0)  then
              facreac = 3.0
            else
              facreac = 1.0
            endif
          endif
cc         szig = 1.116*(1-1.142*xx)/(1+47.47*xx**0.7894)
        endif
      else if(id1.eq.2) then ! eta cross sections are stollen from omega
        if(id2.eq.1) then         ! eta+n -> phi+n
          ireac = 21
          szig = sqrt(xx) * 0.05275/(1.-0.8732*meff1+1.467*xx) *
     &        exp(1. - 1.507*meff1 - 1.523*xx +
     &        1.051*meff1**2 + 0.4646*xx**2 - 0.8308*meff1*xx)
          if (iz1 .ne. 0) facreac = 2.0
        else if(id2.eq.2) then    ! eta+D -> phi+n      zm
          ireac = 22
          szig = phi_rhoDtab(meff2, meff1, xx)
c          call  sig_del_rho(meff2, meff1, xx, szig)
          if (iz1 .eq. 0)  then
            facreac = 2.0
          else
            xz = abs (iz2 -0.5) - 1.0
            if (xz .gt. .0)  then
              facreac = 3.0
            else
              facreac = 1.0
            endif
          endif
        else ! eta+R -> phi+n
          ireac = 23
          szig = phi_rhoDtab(meff2, meff1, xx)
c          call  sig_del_rho(meff2, meff1, xx, szig)
          if (iz1 .eq. 0)  then
            facreac = 2.0
          else
            xz = abs (iz2 -0.5) - 1.0
            if (xz .gt. .0)  then
              facreac = 3.0
            else
              facreac = 1.0
            endif
          endif
cc         szig = 1.116*(1-1.142*xx)/(1+47.47*xx**0.7894)
        endif
      end if
      szig   = szig/szig0 * facreac

c***********************************************************************
      inkrun = max_pert / num
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
* find a free place for the new phi, or the phi with the smallest prob.
      do 210 kl=1,inkrun                    !        iphi_bb       hw
        i_pert = i_pert + 1
        if (nx_pert(id_phi,0,i_pert) .eq. 0) goto 212
        if (valkapi .gt. p_pert(id_phi,4,i_pert)) then
          inkmin = i_pert
          valkapi  = p_pert(id_phi,4,i_pert)
        endif
  210 continue
      i_pert = inkmin
c      write(*,*) 'vigyazz phi all positions filled'
      if (szig .lt. valkapi) goto 297
  212 continue
c      write(*,*) 'prob. of phi is big enough'

  270 continue
      xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 270

      qqx= xx*qqabs/rr
      qqy= yy*qqabs/rr
      qqz= zz*qqabs/rr
*   phi_ momentum in observable system
      pbeta  = betax*qqx + betay*qqy + betaz*qqz
      transf = gamm * (traf * pbeta + q0)

c      if (nx_omega(0,i_pert).ne.1) iop = iop +1
c      if (nx_omega(0,i_pert).eq.1) iwp = iwp +1
c      if (nx_omega(0,i_pert).ne.1) p_omega(4,i_pert) = 1.
c      write(20,*) p_omega(4,i_pert)

      p_pert(id_phi,1,i_pert)= qqx + betax * transf
      p_pert(id_phi,2,i_pert)= qqy + betay * transf
      p_pert(id_phi,3,i_pert)= qqz + betaz * transf
      p_pert(id_phi,4,i_pert)= szig * p_omega(4,i_pert)
      r_pert(id_phi,1,i_pert)= xxx
      r_pert(id_phi,2,i_pert)= yyy
      r_pert(id_phi,3,i_pert)= zzz
      nx_pert(id_phi,0,i_pert) = 1
      nx_pert(id_phi,2,i_pert) = ireac
      nx_pert(id_phi,3,i_pert) = i2
      nx_pert(id_phi,4,i_pert) = 0
      nx_pert(id_phi,5,i_pert) = 1
      nx_pert(id_phi,6,i_pert) = id1
      nx_pert(id_phi,7,i_pert) = id2
      nx_pert(id_phi,8,i_pert) = iz1
      nx_pert(id_phi,9,i_pert) = iz2

c          write(20,*) iop,iwp

      goto  297
c-----------------------------------------------------------------------
 297  continue
      return


*----------------------------------------------------------------------*
      entry omega_dbb(id1,id2,beta,gamma,srt,xxx,yyy,zzz,sig0,i61
     &  ,irun,i1,i2,etotal,id21,id22)
c  perturbative omega prod.
c  B+B -> omega + x,  B < 3
*----------------------------------------------------------------------*
      if(id1.gt.2 .or. id2.gt.2)                      return
cc
      pxx =  .0

c      if (iphi_pot .ge. 1) then
c         call gradupi(xxx,yyy,zzz, pxx,pxx,pxx,tmass,5,0,i_pert,
c     1                vx,vy,vz,vxx,vyy,vzz,gradm,1.0,dt,1)
c       else
      tmass = omass
c      endif

      srtmin = 2.*rmass + tmass
      s     = srt**2
      pmax2=.25*(s-(rmass+rmass+tmass)**2)*(s-(rmass+rmass-tmass)**2)/s
      if(pmax2 .le. 0.0)                                          return
c     pmax  = the maximal omega_ momentum
      pmax  = sqrt(pmax2)
      traf  = gamma / (gamma+1.0)
***                                                                  ***
      ireac = 0
      if(id1+id2 .eq. 2) then
        ireac = 21
      else
!      return      !ANKE  nur NN kollisionen, die ein omega erzeugen
      endif
      if(id1+id2 .eq. 3) ireac = 22
      if(id1.eq.2 .and. id2.eq.2) ireac = 23

      facreac = 1.0
      xx = log(srt - srtmin + .000005)
      ds12 = exp(xx)

      if (ireac .eq. 22) facreac = 6.0
      if (ireac .eq. 23) facreac = 2.0

!   kaptari & kaempfer with FSI, 2004

*************** constants for sig_pp **************************
      As1 = -28.16974e-3
      As2 = 328.1566e-3
      As3 = 0.473162
      As4 = 0.189525

      szig = As2 + (As1-As2)/(1.+exp((ds12-As3)/As4))
*************** constants for sig_pn **************************
      Bs1 = 1307.62e-3
      Bs2 = 1.35106

      if ( ((id22.eq.0 .and. id21.ne.0) .or.
     &        (id21.eq.0 .and. id22.ne.0)) .and. ireac.eq.21) then

        szig_pn = Bs1*exp(Bs2*xx)
        szig   = szig_pn/sig0 * facreac         !sig(pn) = sig(np)
      else
        szig   = szig/sig0 * facreac            !sig(pp) = sig(nn)
      endif

      szig = 1.*szig
c        write(20,*) "omega[mb]: ",szig*sig0

c***********************************************************************
      inkrun = (max_omega/num)
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkabb  = 1000.0
      do 110 kl=1,inkrun
        i_pert = i_pert + 1
        if (nx_omega(0,i_pert) .eq. 0) goto 112
c        if (ipi(1,i_pert).eq.5 .or. ipi(1,i_pert).eq.25) goto 199
        if (valkabb .gt. p_omega(4,i_pert)) then
           inkmin = i_pert
           valkabb  = p_omega(4,i_pert)
        endif
  110 continue

      i_pert = inkmin
      if (szig .lt. valkabb) goto 199

  112 continue

      phi_b  = 1.+2.*rmass/tmass
  115 continue
      ranp     = rn(iseed)
      xx = ranp**2 * sqrt(1.00001-ranp**2) / 0.385
      if(xx .lt. rn(iseed)) go to 115

      pka       = pmax*ranp
  120 xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 120
      eka      = sqrt(pka**2 + tmass**2)
      spr      = s - 2.0*srt*eka + tmass**2
      enupr    = 0.5* sqrt(spr)
      pnupr    = sqrt(enupr**2-rmass**2)
      qqx      = pka * xx / rr
      qqy      = pka * yy / rr
      qqz      = pka * zz / rr                   !   omega -meson

 130  xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 130
      factn1   = spr + sqrt(spr)*(srt-eka)
      factn2=pnupr*(qqx*xx+qqy*yy+qqz*zz)/factn1/rr-enupr/sqrt(spr)
*   p3:                    nucleon momentum in i1-i2-c.m. system
      p3(1)    = pnupr*xx/rr + factn2*qqx
      p3(2)    = pnupr*yy/rr + factn2*qqy
      p3(3)    = pnupr*zz/rr + factn2*qqz
      e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p3:                    nucleon momentum in observable system
      pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
      transf = gamma * (traf * pbeta + e3)
      p3(1)  = p3(1) + beta(1) * transf
      p3(2)  = p3(2) + beta(2) * transf
      p3(3)  = p3(3) + beta(3) * transf
      if(ipauli.eq.1)
     &  call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
      facphi_di2 = szig*(1.0-phase)
*
*   p4:                    nucleon momentum in i1-i2-c.m. system
      p3(1)    = -pnupr*xx/rr + factn2*qqx
      p3(2)    = -pnupr*yy/rr + factn2*qqy
      p3(3)    = -pnupr*zz/rr + factn2*qqz
      e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p4:                    nucleon momentum in observable system
      pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
      transf = gamma * (traf * pbeta + e3)
      p3(1)  = p3(1) + beta(1) * transf
      p3(2)  = p3(2) + beta(2) * transf
      p3(3)  = p3(3) + beta(3) * transf
      if(ipauli.eq.1)
     &  call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))

      facphi_di2 = facphi_di2 * (1.0-phase)
*
*   omega_ momentum in observable system
      pbeta  = beta(1)*qqx + beta(2)*qqy + beta(3)*qqz
      transf = gamma * (traf * pbeta + eka)

c      ppi(1,i_pert)= qqx + beta(1) * transf
c      ppi(2,i_pert)= qqy + beta(2) * transf
c      ppi(3,i_pert)= qqz + beta(3) * transf

      p_omega(1,i_pert)= qqx + beta(1) * transf
      p_omega(2,i_pert)= qqy + beta(2) * transf
      p_omega(3,i_pert)= qqz + beta(3) * transf
c      ipi(1,i_pert) = 25
c      rpi(1,i_pert)= xxx
c      rpi(2,i_pert)= yyy
c      rpi(3,i_pert)= zzz
      r_omega(1,i_pert)= xxx
      r_omega(2,i_pert)= yyy
      r_omega(3,i_pert)= zzz
      p_omega(4,i_pert)= facphi_di2
      zig = p_omega(4,i_pert)

      nx_omega(0,i_pert) = 1
      nx_omega(2,i_pert) = ireac
      nx_omega(3,i_pert) = i1
      nx_omega(4,i_pert) = i2
      nx_omega(5,i_pert) = 0              !el. charge
      nx_omega(6,i_pert) = id1
      nx_omega(7,i_pert) = id2
      nx_omega(8,i_pert) = id21
      nx_omega(9,i_pert) = id22
      numprod = numprod + 1
c      write(20,*) "#omega prod: ",i_pert,i1,i2,irun,numprod

 199  continue

      return
************************************************************************

*----------------------------------------------------------------------*
      entry eta_dbb(id1,id2,beta,gamma,srt,xxx,yyy,zzz,sig0,i61
     &  ,irun,i1,i2,etotal,id21,id22)
c  perturbative eta prod. not used
c  B+B -> eta + x,  B < 3
*----------------------------------------------------------------------*

      if(id1.gt.2 .or. id2.gt.2)                      return
cc
      pxx =  .0
      tmass = emass

      srtmin = 2.*rmass + tmass
      s     = srt**2
      pmax2=.25*(s-(rmass+rmass+tmass)**2)*(s-(rmass+rmass-tmass)**2)/s
      if(pmax2 .le. 0.0)                                          return
c      pmax  = the maximal eta_momentum
      pmax  = sqrt(pmax2)
      traf  = gamma / (gamma+1.0)
***                                                                  ***
      ireac = 0
      if(id1+id2 .eq. 2) then
        ireac = 31
      else
!        return	         !ANKE  nur NN kollisionen, die ein eta erzeugen
      endif
      if(id1+id2 .eq. 3) ireac = 32
      if(id1.eq.2 .and. id2.eq.2) ireac = 33

      facreac = 1.0
      xx = log(srt - srtmin + .000005)
      ds12 = exp(xx)

      if (ireac .eq. 32) facreac = 6.0
      if (ireac .eq. 33) facreac = 2.0

!       kaptari & kaempfer with FSI, 2004

*************** constants for sig_pp **************************
      As1 = 1.1064
      As2 = 1.7167
      szig = As1*exp(As2*xx)
*************** constants for sig_pn **************************
      Bs1 = 2.624
      Bs2 = 1.1892

      if ( ((id22.eq.0 .and. id21.ne.0) .or.
     &        (id21.eq.0 .and. id22.ne.0)) .and. ireac.eq.31) then
        szig_pn = Bs1*exp(Bs2*xx)
        szig   = szig_pn/sig0 * facreac         !sig(pn) = sig(np)
      else
        szig   = szig/sig0 * facreac            !sig(pp) = sig(nn)
      endif

c***********************************************************************
      inkrun = (maxeta/num)
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkabb  = 1000.0
      do 510 kl=1,inkrun
        i_pert = i_pert + 1
        if (nx_eta(0,i_pert) .eq. 0) goto 512
        if (valkabb .gt. p_eta(4,i_pert)) then
           inkmin = i_pert
           valkabb  = p_eta(4,i_pert)
        endif
  510 continue

      i_pert = inkmin
      if (szig .lt. valkabb) goto 599

  512 continue
      phi_b  = 1.+2.*rmass/tmass
  515 continue
      ranp     = rn(iseed)
      xx = ranp**2 * sqrt(1.00001-ranp**2) / 0.385
      if(xx .lt. rn(iseed)) go to 515

      pka       = pmax*ranp
 520  xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 520
      eka      = sqrt(pka**2 + tmass**2)
      spr      = s - 2.0*srt*eka + tmass**2
      enupr    = 0.5* sqrt(spr)
      pnupr    = sqrt(enupr**2-rmass**2)
      qqx      = pka * xx / rr
      qqy      = pka * yy / rr
      qqz      = pka * zz / rr                   !   eta -meson

 530  xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 530
      factn1   = spr + sqrt(spr)*(srt-eka)
      factn2=pnupr*(qqx*xx+qqy*yy+qqz*zz)/factn1/rr-enupr/sqrt(spr)
*   p3:                    nucleon momentum in i1-i2-c.m. system
      p3(1)    = pnupr*xx/rr + factn2*qqx
      p3(2)    = pnupr*yy/rr + factn2*qqy
      p3(3)    = pnupr*zz/rr + factn2*qqz
      e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p3:                    nucleon momentum in observable system
      pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
      transf = gamma * (traf * pbeta + e3)
      p3(1)  = p3(1) + beta(1) * transf
      p3(2)  = p3(2) + beta(2) * transf
      p3(3)  = p3(3) + beta(3) * transf
      if(ipauli.eq.1)
     &  call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
      facphi_di2 = szig*(1.0-phase)
*
*   p4:                    nucleon momentum in i1-i2-c.m. system
      p3(1)    = -pnupr*xx/rr + factn2*qqx
      p3(2)    = -pnupr*yy/rr + factn2*qqy
      p3(3)    = -pnupr*zz/rr + factn2*qqz
      e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p4:                    nucleon momentum in observable system
      pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
      transf = gamma * (traf * pbeta + e3)
      p3(1)  = p3(1) + beta(1) * transf
      p3(2)  = p3(2) + beta(2) * transf
      p3(3)  = p3(3) + beta(3) * transf
      if(ipauli.eq.1)
     &  call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))

      facphi_di2 = facphi_di2 * (1.0-phase)
*
*   eta_ momentum in observable system
      pbeta  = beta(1)*qqx + beta(2)*qqy + beta(3)*qqz
      transf = gamma * (traf * pbeta + eka)
      p_eta(1,i_pert)= qqx + beta(1) * transf
      p_eta(2,i_pert)= qqy + beta(2) * transf
      p_eta(3,i_pert)= qqz + beta(3) * transf

      r_eta(1,i_pert)= xxx
      r_eta(2,i_pert)= yyy
      r_eta(3,i_pert)= zzz
      p_eta(4,i_pert)= facphi_di2
      zzig = p_eta(4,i_pert)

      nx_eta(0,i_pert) = 1
      nx_eta(2,i_pert) = ireac
      nx_eta(3,i_pert) = i1
      nx_eta(4,i_pert) = i2
      nx_eta(5,i_pert) = 0              !el. charge
      nx_eta(6,i_pert) = id1
      nx_eta(7,i_pert) = id2
      nx_eta(8,i_pert) = id21
      nx_eta(9,i_pert) = id22

 599  continue
      return
************************************************************************

*----------------------------------------------------------------------*
      entry phi_dbb(id1,id2,beta,gamma,srt,xxx,yyy,zzz,sig0,i61,irun
     &  ,i1,i2,etotal,id21,id22)
*       variables:                                                     *
*         ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*     cross-sections are taken from sibirtsev   or  Chung PL B401(97)1 *
*----------------------------------------------------------------------*

c      if(id1.gt.2 .or. id2.gt.2)                      return
c      if(id1.lt.0 .or. id2.lt.0)                      return
      pphi = 1./real(pphi_num)
      if(rn(iseed) .ge. pphi)                    return
cc
      ix = nint(xxx)            ! nint: closest integer
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &     denst = rhb(ix,iy,iz)

      distr = sqrt(xxx**2 + yyy**2 + zzz**2)

      pxx =  .0
cc      write(*,*)  '  phi_ mass phi_dbb1'
      if (iphi_pot .ge. 1) then
c        write(*,*) 'phi_dbb 2', xxx,yyy,zzz, pxx,pxx,pxx,mmass,
c     1                 vx,vy,vz,vxx,vyy,vz
        call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                   vx,vy,vz,vxx,vyy,vzz)
      else
        mmass = xphimas
      endif
      tmass  = mmass
      srtmin = 2.*rmass + tmass
cc      write(*,*)  '  phi_ mass phi_dbb ',tmass, srtmin, srt
      s     = srt**2
      pmax2=.25*(s-(rmass+rmass+tmass)**2)*(s-(rmass+rmass-tmass)**2)/s
      if(pmax2 .le. 0.0)                                          return
c     pmax  = the maximal phi_ momentum
      pmax  = sqrt(pmax2)
      traf  = gamma / (gamma+1.0)
***                                                                  ***
      ireac = 0
      if(id1+id2 .eq. 2) ireac = 1
      if(id1+id2 .eq. 3) ireac = 2
      if(id1.eq.2 .and. id2.eq.2) ireac = 3
      if(id1.gt.2 .or. id2.gt.2) ireac = 4
      if(id1.lt.0 .or. id2.lt.0) ireac = 4
      facreac = 1.0
      xx = log(srt - srtmin + .000005)
      szig1 = 0.002
      if (iphi_cr .eq. 2) then                 !    sibirtsev
        facreac = 1.0
        if (xx  .gt. 0.) then
          szig = szig1
        else
          szig = szig1 * exp(-3.7 * (abs(xx)/2.66)**1.64)
        endif
      elseif (iphi_cr .eq. 3) then        !  titov, kaempfer, PR C59
        facreac = 1.0
        if (xx .lt. -2.3)  then
          szig   = .015 * exp ( 2. * xx)
        elseif (xx .lt. 0.26) then
          szig   = .015 * exp(-1.770 + xx * (.7112 - .2258*xx))
        else
          szig = .003
        endif
      elseif (iphi_cr .eq. 4) then   !   kaptari & kaempfer with FSI, 2004
*************** constants for sig_pn/sig_pp **************************
        As1 = 6.040e-1
        As2 = 1.583e-1
        As3 = -1.217e-4
        As4 = 7.320
        As5 = -3.137
        As6 = 1.305
        As7 = 1.564
*************** constants for sig_pp **************************
        Bs1 = 5.360e-3
        Bs2 = 1.298
****************constants for sig_pp from Maeda/Hartmann******
        Ah1 = -6.815e-6
        Bh1 = 2.224e-4
        Ch1 = 0.084938
        Dh1 = 0.030879
        Eh1 = 1./0.491
**************************************************************

        facreac = 1.0
        szig = Bs1*exp(Bs2*xx)
        ds12 = exp(xx)
!        szig = Eh1*(Bh1+(Ah1-Bh1)/(1.+exp((ds12-Ch1)/Dh1)))
!        write(54,*) ds12*1000.
c        binsqrts = nint(ds12/0.01)
        szig_pn = 0.0
        fak = 0.0
        if ( ((id22.eq.0 .and. id21.ne.0) .or.
     &        (id21.eq.0 .and. id22.ne.0)) .and. ireac.eq.1) then
!        if (id22.eq.0 .and. ireac.eq.1) then

          fak =   As1*log(As2*ds12-As3)+As4+As5*ds12+As6*ds12**As7
!          fak =   1.                            !ANKE test
          szig_pn = szig * fak
        endif
        if (ireac .eq. 2) then
          facreac = 6.0
c          c_DN(binsqrts) = c_DN(binsqrts) + 1
        endif
        if (ireac .eq. 3) then
          facreac = 2.0
c          c_DD(binsqrts) = c_DD(binsqrts) + 1
        endif
        if (ireac .eq. 4) then
          facreac = 1.0 ! why not		
c          c_DD(binsqrts) = c_DD(binsqrts) + 1
        endif

      else   !   chung + ... PL B401
        szig = szig1 * exp(-3.9 + 1.9 * (xx +2.3) )
        szig_pn = szig

        facreac = 1.0
        if (ireac .eq. 2) facreac = 6.0
        if (ireac .eq. 3) facreac = 2.0
      endif

      if ( ((id22.eq.0 .and. id21.ne.0) .or.
     &        (id21.eq.0 .and. id22.ne.0)) .and. ireac.eq.1) then
!          if (id22.eq.0 .and. ireac.eq.1) then
        szig   = szig_pn/sig0 * facreac / pphi         !sig(pn) = sig(np)
c        c_pn(binsqrts) = c_pn(binsqrts) + 1
      else
        szig   = szig/sig0 * facreac / pphi            !sig(pp) = sig(nn)
c        if (ireac .eq. 1) c_pp(binsqrts) = c_pp(binsqrts) + 1
      endif

c                write(20,*) "phi[mb]: ",szig*sig0
c***********************************************************************
      inkrun = (max_pert/num)
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkabb  = 1000.0
      do 10 kl=1,inkrun
        i_pert = i_pert + 1
        if (nx_pert(id_phi,0,i_pert) .eq. 0) goto 12
        if (valkabb .gt. p_pert(id_phi,4,i_pert)) then
           inkmin = i_pert
           valkabb  = p_pert(id_phi,4,i_pert)
        endif
   10 continue
      i_pert = inkmin
      if (szig .lt. valkabb) goto 99
   12 continue
cc      write(*,*)   '  phi_bb   phi found ',i_pert
cc    if (szig .lt. valkabb) goto 99                    !  error
      phi_b  = 1.+2.*rmass/xphimas
   15 continue
      ranp     = rn(iseed)
c      xx       = ranp * (ranp-sqrt(phi_b * (1.-ranp**2)))**2 /0.53
      xx = ranp**2 * sqrt(1.00001-ranp**2) / 0.385
      if(xx .lt. rn(iseed)) go to 15
c
!          write(54,*) rhb(nint(xxx),nint(yyy),nint(zzz))
      pka       = pmax*ranp

  20  xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 20
      eka      = sqrt(pka**2 + tmass**2)
      spr      = s - 2.0*srt*eka + tmass**2
      enupr    = 0.5* sqrt(spr)
      pnupr    = sqrt(enupr**2-rmass**2)
      qqx      = pka * xx / rr
      qqy      = pka * yy / rr
      qqz      = pka * zz / rr                   !    Phi -meson

c          write(20,*) srt,beta(3),pka

  30  xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 30
      factn1   = spr + sqrt(spr)*(srt-eka)
      factn2=pnupr*(qqx*xx+qqy*yy+qqz*zz)/factn1/rr-enupr/sqrt(spr)
*   p3:                    nucleon momentum in i1-i2-c.m. system
      p3(1)    = pnupr*xx/rr + factn2*qqx
      p3(2)    = pnupr*yy/rr + factn2*qqy
      p3(3)    = pnupr*zz/rr + factn2*qqz
      e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p3:                    nucleon momentum in observable system
      pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
      transf = gamma * (traf * pbeta + e3)
      p3(1)  = p3(1) + beta(1) * transf
      p3(2)  = p3(2) + beta(2) * transf
      p3(3)  = p3(3) + beta(3) * transf
      if(ipauli.eq.1)
     &  call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
      facphi_di2 = szig*(1.0-phase)
*
*   p4:                    nucleon momentum in i1-i2-c.m. system
      p3(1)    = -pnupr*xx/rr + factn2*qqx
      p3(2)    = -pnupr*yy/rr + factn2*qqy
      p3(3)    = -pnupr*zz/rr + factn2*qqz
      e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*   p4:                    nucleon momentum in observable system
      pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
      transf = gamma * (traf * pbeta + e3)
      p3(1)  = p3(1) + beta(1) * transf
      p3(2)  = p3(2) + beta(2) * transf
      p3(3)  = p3(3) + beta(3) * transf
      if(ipauli.eq.1)
     &  call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
      facphi_di2 = facphi_di2 * (1.0-phase)
*
*   phi_ momentum in observable system
      pbeta  = beta(1)*qqx + beta(2)*qqy + beta(3)*qqz
      transf = gamma * (traf * pbeta + eka)

      p_pert(id_phi,0,i_pert)= tmass
      p_pert(id_phi,1,i_pert)= qqx + beta(1) * transf
      p_pert(id_phi,2,i_pert)= qqy + beta(2) * transf
      p_pert(id_phi,3,i_pert)= qqz + beta(3) * transf
      p_pert(id_phi,4,i_pert)= facphi_di2
      r_pert(id_phi,1,i_pert)= xxx
      r_pert(id_phi,2,i_pert)= yyy
      r_pert(id_phi,3,i_pert)= zzz
      r_pert(id_phi,4,i_pert)= denst/rho0
      r_pert(id_phi,5,i_pert)= time
      nx_pert(id_phi,0,i_pert) = 1
      nx_pert(id_phi,1,i_pert) = 1
      nx_pert(id_phi,2,i_pert) = ireac
      nx_pert(id_phi,3,i_pert) = i1
      nx_pert(id_phi,4,i_pert) = i2
      nx_pert(id_phi,5,i_pert) = 1
      nx_pert(id_phi,6,i_pert) = id1
      nx_pert(id_phi,7,i_pert) = id2
      nx_pert(id_phi,8,i_pert) = id21
      nx_pert(id_phi,9,i_pert) = id22
      numprodd = numprodd + 1

          if(p_pert(id_phi,0,i_pert).lt.0.05) then
            write(*,*) 'phi_dbb phimass',p_pert(id_phi,0,i_pert)
            stop
          end if
c          write(20,*) "#phi prod: ",i_pert,i1,i2,irun,numprodd
c-----------------------------------------------------------------------
c           density dependence
c!      if (zzz.le.0.0) then                  !only front-hemisphere z<0
c      zzz = zzz + 18.876 - rmax	             !=radta+radpr+rdist
c
c      if (zzz.lt.-10.)        bi(0) = bi(0) + 1. !NOT weightened with b (MB)
c      if (zzz.ge.-10. .and. zzz.lt.-9.) bi(1) = bi(1) + 1.		
c      if (zzz.ge.-9. .and. zzz.lt.-8.)  bi(2) = bi(2) + 1.
c      if (zzz.ge.-8. .and. zzz.lt.-7.)  bi(3) = bi(3) + 1.
c      if (zzz.ge.-7. .and. zzz.lt.-6.)  bi(4) = bi(4) + 1.
c      if (zzz.ge.-6. .and. zzz.lt.-5.)  bi(5) = bi(5) + 1.
c      if (zzz.ge.-5. .and. zzz.lt.-4.)  bi(6) = bi(6) + 1.
c      if (zzz.ge.-5. .and. zzz.lt.-4.)  bi(6) = bi(6) + 1.
c      if (zzz.ge.-5. .and. zzz.lt.-4.)  bi(6) = bi(6) + 1.
c      if (zzz.ge.-4. .and. zzz.lt.-3.)  bi(7) = bi(7) + 1.
c      if (zzz.ge.-3. .and. zzz.lt.-2.)  bi(8) = bi(8) + 1.
c      if (zzz.ge.-2. .and. zzz.lt.-1.)  bi(9) = bi(9) + 1.
c      if (zzz.ge.-1. .and. zzz.lt.0.)   bi(10) = bi(10) + 1.
c      if (zzz.ge.0. .and. zzz.lt.1.)    bi(11) = bi(11) + 1.
c      if (zzz.ge.1. .and. zzz.lt.2.)    bi(12) = bi(12) + 1.
c      if (zzz.ge.2. .and. zzz.lt.3.)    bi(13) = bi(13) + 1.
c      if (zzz.ge.3. .and. zzz.lt.4.)    bi(14) = bi(14) + 1.
c      if (zzz.ge.4. .and. zzz.lt.5.)    bi(15) = bi(15) + 1.
c      if (zzz.ge.5. .and. zzz.lt.6.)    bi(16) = bi(16) + 1.
c      if (zzz.ge.6. .and. zzz.lt.7.)    bi(17) = bi(17) + 1.
c      if (zzz.ge.7. .and. zzz.lt.8.)    bi(18) = bi(18) + 1.
c      if (zzz.ge.8. .and. zzz.lt.9.)    bi(19) = bi(19) + 1.
c      if (zzz.ge.9. .and. zzz.lt.10.)   bi(20) = bi(20) + 1.
c      if (zzz.ge.10.)                   bi(21) = bi(21) + 1.
c!      endif
      pkao(8,i_pert)= denst
      pkao(13,i_pert)= time
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      p_dbb(bin_dens,bin_time)=p_dbb(bin_dens,bin_time)+1

!           write(54,*) time, denst,ireac, xxx,yyy,zzz
!       write(20,*) "#phis: ", p_pert(id_phi,4,i_pert)
!      co = 1234
!       write(52,233) time, denst,ireac
cccc      write(*,*) 'end of phi_dbb'
   99 continue
      return
************************************************************************
      entry phi_dpi(ede,pdx,pdy,pdz,srt,xxx,yyy,zzz,szig0,i2,id2,id6,
     &        iz1,iz2,id1,irun)
*       variables:    1 = pion       2 = baryon                        *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*     cross-sections are taken from Chung et al., Phys. lett. B401(97)1*
*----------------------------------------------------------------------*
c      write(*,*) 'in phi_dpi',srt

      if(iz1+iz2 .gt. 1)                                          return
      if(iz1+iz2 .lt. 0)                                          return
      if(id1.ne.1)                                                return
      if(id2.gt.nres+1)                                           return
      pxx =  .0

      ix = nint(xxx)
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz)
      distr = sqrt(xxx**2 + yyy**2 + zzz**2)

      if (iphi_pot .ge. 1) then
c        write(*,*) 'phi_dpi 2', xxx,yyy,zzz, pxx,pxx,pxx,mmass,
c     1                 vx,vy,vz,vxx,vyy,vz
        call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                      vx,vy,vz,vxx,vyy,vzz)
      else
        mmass = xphimas
      endif
      tmass  = mmass
      srt0  = rmass + tmass
      if(srt.le.srt0)                                             return
      tmass2= tmass**2
      rmass2= rmass**2
      s     = srt**2
c       q0  : energy of the phi_
      q0      = 0.5 * (s + tmass2 - rmass2) / srt
      if(q0 .le. tmass)                                    return
      qq2     = q0**2 - tmass2
      qqabs   = sqrt(qq2)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gamm  = ede / srt
      traf  = gamm / (gamm+1.0)
*----------------------------------------------------------------------*
c...  cross sections:
      facreac = 0.0

      if(srt-srt0 .le. 0.00001)                                   return
      if(id2.eq.1) then         ! N+pi
        ireac = 5
        xx = log (srt - srt0 + .00001)
        szig = 1.*0.015 * exp (0.4292 + 0.3497*xx - 0.00664*xx*xx) !faktor 1
            write(20,*) szig,xx,srt
      else if(id2.eq.2) then    ! D+pi
        ireac = 6
        xx = log (srt - srt0 + .00001)
        szig =1.*(0.0018 - 0.00195*xx) * exp(0.3408*xx - 0.03285*xx*xx)
      else if(id2 .eq. 3) then  ! N(1440)
        ireac   = 7                       ! N(1440)+pi
        xx = srt - srt0 + .00001
c        szig = .0
        szig = .028 * xx / (xx**1.69  + .1170) !!war bisher auf 0.0 gesetzt
      else if(id2 .eq. 4) then                 ! N(1520)
        ireac   = 8                           ! N(1520)+pi
        xx = srt - srt0 + .00001
c        xxs = sqrt(xx)
c        szig = .0030 * xxs / (xxs*xx + .0223)
        xxy = xx / .05                     !  .05  gives maximum
        xn  = 1.75
        szig = 2.*xn*.036 * sqrt(xxy) / (xxy**xn+2.*xn-1.)

      else if(id2 .gt. 4) then                 ! R > N(1520)
        ireac   = 9                            ! R+pi
        xx = srt - srt0 + .00001
c        xxs = sqrt(xx)
c        szig = .0030 * xxs / (xxs*xx + .0223)
        xxy = xx / .05                     !  .05  gives maximum
        xn  = 1.75
        szig = 2.*xn*.036 * sqrt(xxy) / (xxy**xn+2.*xn-1.)
      endif
c...  isospin factors:
      if(id2.eq.1 .or. id2.eq.3 .or. id2.ge.4) then
c      if(id2.eq.1 .or. id2.eq.4) then
        facreac = 1.0
        if(iz1.ne.0) facreac=2.0
      else if(id2.eq.2) then
        facreac = 1.0
        if(iz1.eq.0) facreac = 2.0
        if(iz2.eq.-1 .or. iz2.eq.2) facreac = 3.0
      end if
      szig   = szig/szig0 * facreac
c***********************************************************************
      inkrun = max_pert / num
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
      do 710 kl=1,inkrun                    !        iphi_bb       hw
        i_pert = i_pert + 1
        if (nx_pert(id_phi,0,i_pert) .eq. 0) goto 712
        if (valkapi .gt. p_pert(id_phi,4,i_pert)) then
           inkmin = i_pert
           valkapi  = p_pert(id_phi,4,i_pert)
        endif
  710 continue
      i_pert = inkmin
      if (szig .lt. valkapi) goto 98
  712 continue
   60 continue
      xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 60
cc       write(*,*) 'phi_dpi',szig,traf,q0, i_pert, inkm, irun, kl
      qqx= xx*qqabs/rr
      qqy= yy*qqabs/rr
      qqz= zz*qqabs/rr
*   phi_ momentum in observable system
      pbeta  = betax*qqx + betay*qqy + betaz*qqz
      transf = gamm * (traf * pbeta + q0)
      p_pert(id_phi,0,i_pert)= tmass
      p_pert(id_phi,1,i_pert)= qqx + betax * transf
      p_pert(id_phi,2,i_pert)= qqy + betay * transf
      p_pert(id_phi,3,i_pert)= qqz + betaz * transf
      p_pert(id_phi,4,i_pert)= szig
      r_pert(id_phi,1,i_pert)= xxx
      r_pert(id_phi,2,i_pert)= yyy
      r_pert(id_phi,3,i_pert)= zzz
      r_pert(id_phi,4,i_pert)= denst/rho0
      r_pert(id_phi,5,i_pert)= time
      nx_pert(id_phi,0,i_pert) = 1
      nx_pert(id_phi,1,i_pert) = 1
      nx_pert(id_phi,2,i_pert) = ireac
      nx_pert(id_phi,3,i_pert) = i2
      nx_pert(id_phi,4,i_pert) = 0
      nx_pert(id_phi,5,i_pert) = 1
      nx_pert(id_phi,6,i_pert) = id1
      nx_pert(id_phi,7,i_pert) = id2
      nx_pert(id_phi,8,i_pert) = iz1
      nx_pert(id_phi,9,i_pert) = iz2
          if(p_pert(id_phi,0,i_pert).lt.0.05) then
            write(*,*) 'phi_dpi phimass',p_pert(id_phi,0,i_pert)
            stop
          end if
!      write(54,*) 'id1,id2',id1,id2
!      write(54,*) 'iz1,iz2',iz1,iz2
c      write(*,*) 'new phi written to array'
c-----------------------------------------------------------------------
c           density dependence
!      if (zzz.le.0.0) then                  !only front-hemisphere z<0
      numprodd = numprodd + 1
c      if (zzz.lt.-10.) bi(0) = bi(0) + 1.       !NOT weightened with b (MB)
c      if (zzz.ge.-10. .and. zzz.lt.-9.)     bi(1) = bi(1) + 1.		
c      if (zzz.ge.-9. .and. zzz.lt.-8.)      bi(2) = bi(2) + 1.
c      if (zzz.ge.-8. .and. zzz.lt.-7.)      bi(3) = bi(3) + 1.
c      if (zzz.ge.-7. .and. zzz.lt.-6.)      bi(4) = bi(4) + 1.
c      if (zzz.ge.-6. .and. zzz.lt.-5.)      bi(5) = bi(5) + 1.
c      if (zzz.ge.-5. .and. zzz.lt.-4.)      bi(6) = bi(6) + 1.
c      if (zzz.ge.-4. .and. zzz.lt.-3.)      bi(7) = bi(7) + 1.
c      if (zzz.ge.-3. .and. zzz.lt.-2.)      bi(8) = bi(8) + 1.
c      if (zzz.ge.-2. .and. zzz.lt.-1.)      bi(9) = bi(9) + 1.
c      if (zzz.ge.-1. .and. zzz.lt.0.)       bi(10) = bi(10) + 1.
c      if (zzz.ge.0. .and. zzz.lt.1.)        bi(11) = bi(11) + 1.
c      if (zzz.ge.1. .and. zzz.lt.2.)        bi(12) = bi(12) + 1.
c      if (zzz.ge.2. .and. zzz.lt.3.)        bi(13) = bi(13) + 1.
c      if (zzz.ge.3. .and. zzz.lt.4.)        bi(14) = bi(14) + 1.
c      if (zzz.ge.4. .and. zzz.lt.5.)        bi(15) = bi(15) + 1.
c      if (zzz.ge.5. .and. zzz.lt.6.)        bi(16) = bi(16) + 1.
c      if (zzz.ge.6. .and. zzz.lt.7.)        bi(17) = bi(17) + 1.
c      if (zzz.ge.7. .and. zzz.lt.8.)        bi(18) = bi(18) + 1.
c      if (zzz.ge.8. .and. zzz.lt.9.)        bi(19) = bi(19) + 1.
c      if (zzz.ge.9. .and. zzz.lt.10.)       bi(20) = bi(20) + 1.
c      if (zzz.ge.10.)                       bi(21) = bi(21) + 1.
!      endif
c      pkao(8,i_pert)= denst
c      pkao(13,i_pert)= time
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      p_dpi(bin_dens,bin_time)=p_dpi(bin_dens,bin_time)+1

!           write(54,*) time, denst,ireac, xxx,yyy,zzz
!        co = co + 1
!           write(52,233) time, denst,ireac
      goto  98
c-----------------------------------------------------------------------
 98   continue
c***********************************************************************
c      write(*,*) 'end of phi_dpi'
      return
************************************************************************
      entry phi_rho(ede,pdx,pdy,pdz, srt,meff1,meff2, xxx,yyy,zzz,szig0,
     &              i2,id2,id6,iz1,iz2,id1,irun)
* purpose: phi production in rho+N -> N+phi           zm               *
*                     and in rho+D -> N+phi                            *
*       variables:    1 = pion       2 = baryon                        *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*----------------------------------------------------------------------*
*        the  following  call should be in subroutine   pionab.f
*
*      call phi_rho(ede,pdx,pdy,pdz,srt, meff1,meff2, rdx,rdy,rdz,
*    &     sig0,i2,id2,id6,ipi(2,i1),id(2,i2),ipi(1,i1),irun)
*
*----------------------------------------------------------------------*
cc      write(*,*) 'in phi_rho'
cc      write(*,*)  ede,pdx,pdy,pdz,srt,xxx,yyy,zzz,szig0,i1,i2,id2,id6,
cc     &        iz1,iz2,id1,irun, epi(i1), e(i2)

      if(iz1+iz2 .gt. 1)                                          return
      if(iz1+iz2 .lt. 0)                                          return
      if(id1.ne.3)                                                return
      if((id2.ne.1).and.(id2.ne.2))                               return
c      write(*,*) 'rho+n has been found'

      ix = nint(xxx)
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz)
      distr = sqrt(xxx**2 + yyy**2 + zzz**2)

      pxx =  .0
      if (iphi_pot .ge. 1) then
c        write(*,*) 'phi_rho 2', xxx,yyy,zzz, pxx,pxx,pxx,mmass,
c     1                 vx,vy,vz,vxx,vyy,vz
         call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                      vx,vy,vz,vxx,vyy,vzz)
       else
         mmass = xphimas
      endif
      tmass  = mmass
      srt0  = rmass + tmass

      if(srt.le.srt0)                                             return
c      write(*,*) 'srt is above threshold'
      tmass2= tmass**2
      rmass2= rmass**2
      s     = srt**2
c       q0  : energy of the phi_
      q0      = 0.5 * (s + tmass2 - rmass2) / srt
      if(q0 .le. tmass)                                    return
c      write(*,*) 'phi energy is enough'
      qq2     = q0**2 - tmass2
      qqabs   = sqrt(qq2)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gamm  = ede / srt
      traf  = gamm / (gamm+1.0)
*----------------------------------------------------------------------*
      facreac = 1.0
      xx = srt - srt0
      if(id2.eq.1) then         ! rho+n -> phi+n      zm
        ireac = 10
        szig = sqrt(xx) * 0.05275/(1.-0.8732*meff1+1.467*xx) *
     &        exp(1. - 1.507*meff1 - 1.523*xx +
     &        1.051*meff1**2 + 0.4646*xx**2 - 0.8308*meff1*xx)

c        szig = szig * 5.0                  !!!! TEST: Faktor 5

c        szig = sqrt(xx) * 0.09490/(1.+4.490*xx) *
c     &        exp(1. - 2.195*xx + 0.4853*xx*xx)
        if (iz1 .ne. 0) facreac = 2.0
      else if(id2.eq.2) then    ! rho+D -> phi+n      zm
        ireac = 11
        szig = phi_rhoDtab(meff2, meff1, xx)
c        szig = szig * 5.0                  !!!! TEST: Faktor 5
c         call  sig_del_rho(meff2, meff1, xx, szig)
        if (iz1 .eq. 0)  then
          facreac = 2.0
        else
          xz = abs (iz2 -0.5) - 1.0
          if (xz .gt. .0)  then
            facreac = 3.0
          else
            facreac = 1.0
          endif
        endif
      else if(id2.ge.3) then    ! rho+R -> phi+n      zm
        ireac = 12
        szig = phi_rhoDtab(meff2, meff1, xx)
c        szig = szig * 5.0                  !!!! TEST: Faktor 5
c         call  sig_del_rho(meff2, meff1, xx, szig)
        if (iz1 .eq. 0)  then
          facreac = 2.0
        else
          xz = abs (iz2 -0.5) - 1.0
          if (xz .gt. .0)  then
            facreac = 3.0
          else
            facreac = 1.0
          endif
        endif
cc       szig = 1.116*(1-1.142*xx)/(1+47.47*xx**0.7894)
      else
        szig= .0
      end if
c      write(*,*) 'phi from rhos (mdelta, mrho, srt-mn-mphi, sig): ',
c     &      id2, ireac, meff2, meff1, xx, szig, facreac
      szig   = szig/szig0 * facreac

c***********************************************************************
      inkrun = max_pert / num
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
* find a free place for the new phi, or the phi with the smallest prob.
      do 810 kl=1,inkrun                    !        iphi_bb       hw
        i_pert = i_pert + 1
        if (nx_pert(id_phi,0,i_pert) .eq. 0) goto 812
        if (valkapi .gt. p_pert(id_phi,4,i_pert)) then
          inkmin = i_pert
          valkapi  = p_pert(id_phi,4,i_pert)
        endif
  810 continue
      i_pert = inkmin
      if (szig .lt. valkapi) goto 97
  812 continue
c      write(*,*) 'prob. of phi is big enough'
   70 continue
      xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 70
cc     write(*,*) 'phi_rho',szig,traf,q0, i_pert, inkm, irun, kl
      qqx= xx*qqabs/rr
      qqy= yy*qqabs/rr
      qqz= zz*qqabs/rr
*   phi_ momentum in observable system
      pbeta  = betax*qqx + betay*qqy + betaz*qqz
      transf = gamm * (traf * pbeta + q0)
      p_pert(id_phi,0,i_pert)= tmass
      p_pert(id_phi,1,i_pert)= qqx + betax * transf
      p_pert(id_phi,2,i_pert)= qqy + betay * transf
      p_pert(id_phi,3,i_pert)= qqz + betaz * transf
      p_pert(id_phi,4,i_pert)= szig
      r_pert(id_phi,1,i_pert)= xxx
      r_pert(id_phi,2,i_pert)= yyy
      r_pert(id_phi,3,i_pert)= zzz
      r_pert(id_phi,4,i_pert)= denst/rho0
      r_pert(id_phi,5,i_pert)= time
      nx_pert(id_phi,0,i_pert) = 1
      nx_pert(id_phi,1,i_pert) = 1
      nx_pert(id_phi,2,i_pert) = ireac
      nx_pert(id_phi,3,i_pert) = i2
      nx_pert(id_phi,4,i_pert) = 0
      nx_pert(id_phi,5,i_pert) = 1
      nx_pert(id_phi,6,i_pert) = id1
      nx_pert(id_phi,7,i_pert) = id2
      nx_pert(id_phi,8,i_pert) = iz1
      nx_pert(id_phi,9,i_pert) = iz2
          if(p_pert(id_phi,0,i_pert).lt.0.05) then
            write(*,*) 'phi_rho phimass',p_pert(id_phi,0,i_pert)
            stop
          end if
c       write(*,*) 'new phi from rho+n!'
c-----------------------------------------------------------------------
c       density dependence
!      if (zzz.le.0.0) then   !only front-hemisphere z<0
      numprodd = numprodd + 1
c      if (zzz.lt.-10.) bi(0) = bi(0) + 1. !NOT weightened with b (MB)
c      if (zzz.ge.-10. .and. zzz.lt.-9.) bi(1) = bi(1) + 1.		
c      if (zzz.ge.-9. .and. zzz.lt.-8.)  bi(2) = bi(2) + 1.
c      if (zzz.ge.-8. .and. zzz.lt.-7.)  bi(3) = bi(3) + 1.
c      if (zzz.ge.-7. .and. zzz.lt.-6.)  bi(4) = bi(4) + 1.
c      if (zzz.ge.-6. .and. zzz.lt.-5.)  bi(5) = bi(5) + 1.
c      if (zzz.ge.-5. .and. zzz.lt.-4.)  bi(6) = bi(6) + 1.
c      if (zzz.ge.-4. .and. zzz.lt.-3.)  bi(7) = bi(7) + 1.
c      if (zzz.ge.-3. .and. zzz.lt.-2.)  bi(8) = bi(8) + 1.
c      if (zzz.ge.-2. .and. zzz.lt.-1.)  bi(9) = bi(9) + 1.
c      if (zzz.ge.-1. .and. zzz.lt.0.)   bi(10) = bi(10) + 1.
c      if (zzz.ge.0. .and. zzz.lt.1.)    bi(11) = bi(11) + 1.
c      if (zzz.ge.1. .and. zzz.lt.2.)    bi(12) = bi(12) + 1.
c      if (zzz.ge.2. .and. zzz.lt.3.)    bi(13) = bi(13) + 1.
c      if (zzz.ge.3. .and. zzz.lt.4.)    bi(14) = bi(14) + 1.
c      if (zzz.ge.4. .and. zzz.lt.5.)    bi(15) = bi(15) + 1.
c      if (zzz.ge.5. .and. zzz.lt.6.)    bi(16) = bi(16) + 1.
c      if (zzz.ge.6. .and. zzz.lt.7.)    bi(17) = bi(17) + 1.
c      if (zzz.ge.7. .and. zzz.lt.8.)    bi(18) = bi(18) + 1.
c      if (zzz.ge.8. .and. zzz.lt.9.)    bi(19) = bi(19) + 1.
c      if (zzz.ge.9. .and. zzz.lt.10.)   bi(20) = bi(20) + 1.
c      if (zzz.ge.10.)                   bi(21) = bi(21) + 1.
!      endif
      pkao(8,i_pert)= denst
      pkao(13,i_pert)= time
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      p_rho(bin_dens,bin_time)=p_rho(bin_dens,bin_time)+1
!      write(54,*) time, denst,ireac, xxx,yyy,zzz
!      co = co + 1
!      write(52,233) time, denst,ireac
      goto  97
c-----------------------------------------------------------------------
 97   continue
c      write(*,*) 'end of phi_rho', i_pert, nx_pert(id_phi,0,i_pert),nx_pert(id_phi,5,i_pert),szig
      return
c***********************************************************************
*                                                                      *
      entry phi_pi_pipi(ede,pdx,pdy,pdz, srt, xxx,yyy,zzz,szig0,
     &              i1,i2, iz1,iz2,irun)
* purpose: phi production in pi+pi  -> pi+phi                          *
*       szig0   -  referenz cross section in mb units                  *
*       variables:    1 = pion       2 = pion                          *
*         ede,pdx,- energy, momentum of the pi+pi     sys.(real,input) *
*         srt     - cms energy (mass of the pi+pi     sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*----------------------------------------------------------------------*

c      write(*,*) 'in phi_pipi', srt

c      write(*,*) 'srt is above threshold'
      if(iz1+iz2 .gt. 1)                                          return
      if(iz1*iz2 .eq. 0)                                          return
      if(iz1+iz2 .lt. -1)                                         return

      ix = nint(xxx)
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz)
      distr = sqrt(xxx**2 + yyy**2 + zzz**2)

      pxx =  .0
      if (iphi_pot .ge. 1) then
c        write(*,*) 'phi_pi_pipi 2', xxx,yyy,zzz, pxx,pxx,pxx,mmass,
c     1                 vx,vy,vz,vxx,vyy,vz
        call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                 vx,vy,vz,vxx,vyy,vzz)
      else
        mmass = xphimas
      endif
      tmass  = mmass
      tmass2 = tmass*tmass
      srt0  =  xphimas + pmass
      if(srt.le.srt0)    return
      s     = srt**2
c       q0  : energy of the phi_
      q0      = 0.5 * (s + xphimas**2 - pmass**2) / srt
      if(q0 .le. tmass)                                    return
      qq2     = q0**2 - tmass2
      qqabs   = sqrt(qq2)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
      gamm  = ede / srt
      traf  = gamm / (gamm+1.0)
*----------------------------------------------------------------------*
      szig = 5.0e-3 * (10.*(srt-srt0))**1.5     !  mb
      szig   = szig/szig0
c---------------------------------------------
      inkrun = max_pert / num
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
* find a free place for the new phi, or the phi with the smallest prob.
      do 910 kl=1,inkrun                    !        iphi_bb       hw
        i_pert = i_pert + 1
        if (nx_pert(id_phi,0,i_pert) .eq. 0) goto 912
        if (valkapi .gt. p_pert(id_phi,4,i_pert)) then
          inkmin = i_pert
          valkapi  = p_pert(id_phi,4,i_pert)
        endif
  910 continue
      i_pert = inkmin
      if (szig .lt. valkapi) goto 997
  912 continue
c      write(*,*) 'prob. of phi is big enough'
   71 continue
      xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 71
cc      write(*,*) 'phi_rho',szig,traf,q0, i_pert, inkm, irun, kl
      qqx= xx*qqabs/rr
      qqy= yy*qqabs/rr
      qqz= zz*qqabs/rr
*   phi_ momentum in observable system
      pbeta  = betax*qqx + betay*qqy + betaz*qqz
      transf = gamm * (traf * pbeta + q0)
      ireac = 13
      p_pert(id_phi,0,i_pert)= tmass
      p_pert(id_phi,1,i_pert)= qqx + betax * transf
      p_pert(id_phi,2,i_pert)= qqy + betay * transf
      p_pert(id_phi,3,i_pert)= qqz + betaz * transf
      p_pert(id_phi,4,i_pert)= szig
      r_pert(id_phi,1,i_pert)= xxx
      r_pert(id_phi,2,i_pert)= yyy
      r_pert(id_phi,3,i_pert)= zzz
      r_pert(id_phi,4,i_pert)= denst/rho0
      r_pert(id_phi,5,i_pert)= time
      nx_pert(id_phi,0,i_pert) = 1
      nx_pert(id_phi,1,i_pert) = 1
      nx_pert(id_phi,2,i_pert) = ireac
      nx_pert(id_phi,3,i_pert) = i1
      nx_pert(id_phi,4,i_pert) = i2
      nx_pert(id_phi,5,i_pert) = 1
c      write(*,*) 'new phi from pi + pi!',tmass
c-----------------------------------------------------------------------
c         density dependence
!      if (zzz.le.0.0) then          !only front-hemisphere z<0
      numprodd = numprodd + 1
c      if (zzz.lt.-10.) 	bi(0) = bi(0) + 1.   !NOT weightened with b (MB)
c      if (zzz.ge.-10. .and. zzz.lt.-9.) bi(1) = bi(1) + 1.		
c      if (zzz.ge.-9. .and. zzz.lt.-8.)  bi(2) = bi(2) + 1.
c      if (zzz.ge.-8. .and. zzz.lt.-7.)  bi(3) = bi(3) + 1.
c      if (zzz.ge.-7. .and. zzz.lt.-6.)  bi(4) = bi(4) + 1.
c      if (zzz.ge.-6. .and. zzz.lt.-5.)  bi(5) = bi(5) + 1.
c      if (zzz.ge.-5. .and. zzz.lt.-4.)  bi(6) = bi(6) + 1.
c      if (zzz.ge.-4. .and. zzz.lt.-3.)  bi(7) = bi(7) + 1.
c      if (zzz.ge.-3. .and. zzz.lt.-2.)  bi(8) = bi(8) + 1.
c      if (zzz.ge.-2. .and. zzz.lt.-1.)  bi(9) = bi(9) + 1.
c      if (zzz.ge.-1. .and. zzz.lt.0.)   bi(10) = bi(10) + 1.
c      if (zzz.ge.0. .and. zzz.lt.1.)    bi(11) = bi(11) + 1.
c      if (zzz.ge.1. .and. zzz.lt.2.)    bi(12) = bi(12) + 1.
c      if (zzz.ge.2. .and. zzz.lt.3.)    bi(13) = bi(13) + 1.
c      if (zzz.ge.3. .and. zzz.lt.4.)    bi(14) = bi(14) + 1.
c      if (zzz.ge.4. .and. zzz.lt.5.)    bi(15) = bi(15) + 1.
c      if (zzz.ge.5. .and. zzz.lt.6.)    bi(16) = bi(16) + 1.
c      if (zzz.ge.6. .and. zzz.lt.7.)    bi(17) = bi(17) + 1.
c      if (zzz.ge.7. .and. zzz.lt.8.)    bi(18) = bi(18) + 1.
c      if (zzz.ge.8. .and. zzz.lt.9.)    bi(19) = bi(19) + 1.
c      if (zzz.ge.9. .and. zzz.lt.10.)   bi(20) = bi(20) + 1.
c      if (zzz.ge.10.)                   bi(21) = bi(21) + 1.
!      endif
      pkao(8,i_pert)= denst
      pkao(13,i_pert)= time
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      p_pipi(bin_dens,bin_time)=p_pipi(bin_dens,bin_time)+1
!           write(54,*) time, denst,ireac, xxx,yyy,zzz

!      co = co + 1
!           write(52,233) time, denst,ireac
      goto  997
c-----------------------------------------------------------------------
 997  continue
c      write(*,*) 'end of phi_pi_pipi'
      return
c***********************************************************************
*                                                                      *
      entry phi_KK(i1,irun)
* purpose: phi production in pi+rho  -> phi                            *
*       szig0   -  referenz cross section in mb units                  *
*       variables:    1 = pion       2 = pion                          *
*         srt     - cms energy (mass of the pi+pi     sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*----------------------------------------------------------------------*

c      write(*,*) 'in phi_KK' , i1,irun,ipi(1,i1), max_kminu,num,
c     &     rpi(1,i1),rpi(2,i1),rpi(3,i1),ppi(1,i1),ppi(2,i1),ppi(3,i1)
      xx  = rpi(1,i1)
      yy  = rpi(2,i1)
      zz  = rpi(3,i1)
      px1 = ppi(1,i1)
      py1 = ppi(2,i1)
      pz1 = ppi(3,i1)
      em1 = epi(i1)
      e1  = sqrt( em1**2 + px1**2 + py1**2 + pz1**2 )
      pxx = 0.0
      if (ikaonpot .gt. 0) then
        call gradukaon(-1, xx,yy,zz, pxx,pxx,pxx,em1,1,
     1         vx,vy,vz,vxx, vyy, vzz)
      else
        em1 = xkmas
      endif
      e1  = sqrt( em1**2 + px1**2 + py1**2 + pz1**2 )

c      write(*,*) 'phi_kk e1', max_kminu,num
      call f77flush()
      inkrunk = (max_kminu/num)
      do 410 j2=(irun-1) * inkrunk+1,irun*inkrunk !        iphi_KK
        if (nx_kminu(0,j2) .eq. 0) goto 410
c        write(*,*) 'phi_kk e1',j2,nx_kminu(0,j2),nx_kminu(1,j2),
c     &   p_kminu(1,j2),p_kminu(2,j2),p_kminu(3,j2),p_kminu(4,j2)
        dx  = xx - r_kminu(1,j2)
        if (abs(dx) .gt. delpi)                                goto 410
        dy  = yy - r_kminu(2,j2)
        if (abs(dy) .gt. delpi)                                goto 410
        dz  = zz - r_kminu(3,j2)
        if (abs(dz) .gt. delpi)                                goto 410
        rsqare = dx**2 + dy**2 + dz**2
        if (rsqare .gt. delpi**2)                              goto 410
*         now particles are close enough to each other !
        px2 = p_kminu(1,j2)
        py2 = p_kminu(2,j2)
        pz2 = p_kminu(3,j2)
        xx2 = r_kminu(1,j2)
        yy2 = r_kminu(2,j2)
        zz2 = r_kminu(3,j2)
        if(px2.gt.1.d20 .or. py2.gt.1.d20 .or. pz2.gt.1.d20) 
     &    write(*,*)'phi_KK K- imp nagy',px2,py2,pz2,nx_kminu(1,j2)        
c        write(*,*)'phi_KK K- imp nagy',px2,py2,pz2,nx_kminu(1,j2)        
        em2 = xkmas
        pxx = 0.0
        if (i_kminu_pot  .gt. 0) then
          call gradukao2(-1, xx2,yy2,zz2, pxx,pxx,pxx,em2,-1,
     1         vx,vy,vz,vxx, vyy, vzz)
        else
          em2 = xkmas
        endif
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
c        write(*,*) 'phi_kk brel', brel,dispi
*   is their impact parameter small enough?
*            write(*,*)'small enough 2'
        if (brel .gt. dispi)                               goto 410
        b21    = - p2dr / em2 + p1dr * em2 / p12
        t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
        t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
        ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
*            write(*,*)'closest approach '
        if ( abs(t1+t2) .gt. dt )                          goto 410
*   now  the kaons may annihilate
        xx2     = 0.5*(xx  + xx2)
        yy2     = 0.5*(yy  + yy2)
        zz2     = 0.5*(zz  + zz2)
        s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2

        ix = nint(xx2)            ! nint: closest integer
        iy = nint(yy2)
        iz = nint(zz2)
        denst = 0.0
        if(iabs(ix).le.maxx .and.iabs(iy).le.maxx .and.iabs(iz).le.maxz)
     &     denst = rhb(ix,iy,iz)

        pxx =  .0
        mmass = xphimas
        if (iphi_pot .ge. 1) then
c        write(*,*) 'phi_KK 2', xx2,yy2,zz2, pxx,pxx,pxx,mmass,
c     1                 vx,vy,vz,vxx,vyy,vz
          call graduphi(-1, xx2,yy2,zz2, pxx,pxx,pxx,mmass,
     1                 vx,vy,vz,vxx,vyy,vzz)
        endif
        del = .0044*((mmass-em1-em2)/(xphimas-2.0*xkmas))                          !    e-width / scale
        if(abs(sqrt(s) - mmass) .gt. del)  goto 410
        qk2 = 0.25 * s - xkmas**2
        sigkk=120.*pi/qk2*hbc**2*
     &               bwmes(s,4,idec2pi,iresmode,3,0,iwidth,1)
        szig =  sigkk / (10.*pi*dispi**2)
c        write(*,*) 'phi_kk szig', szig,sigkk,qk2,dispi,pi
c---------------------------------------------
        inkrun = max_pert / num
        i_pert = (irun-1) * inkrun
        inkmin = 0
        valkapi  = 1000.0
* find a free place for the new phi, or the phi with the smallest prob.
        do 415 kl=1,inkrun                    !        iphi_bb       hw
        i_pert = i_pert + 1
        if (nx_pert(id_phi,0,i_pert) .eq. 0) goto 417
c        write(*,*) 'phi_KK i_pert', i_pert,nx_pert(id_phi,1,i_pert),valkapi,p_pert(id_phi,4,i_pert)
        if (valkapi .gt. p_pert(id_phi,4,i_pert)) then
           inkmin = i_pert
           valkapi  = p_pert(id_phi,4,i_pert)
        endif
 415    continue
        i_pert = inkmin
        if (szig .lt. valkapi) goto 410
 417    continue
c        write(*,*) 'prob. of phi is big enough',i_pert,irun,inkrun,
c     &      max_pert,valkapi
        ireac = 15
        p_pert(id_phi,0,i_pert)= mmass
        p_pert(id_phi,1,i_pert)= px1+px2
        p_pert(id_phi,2,i_pert)= py1+py2
        p_pert(id_phi,3,i_pert)= pz1+pz2
        p_pert(id_phi,4,i_pert)= szig
        r_pert(id_phi,1,i_pert)= xx2
        r_pert(id_phi,2,i_pert)= yy2
        r_pert(id_phi,3,i_pert)= zz2
        r_pert(id_phi,4,i_pert)= denst/rho0
        r_pert(id_phi,5,i_pert)= time
        nx_pert(id_phi,0,i_pert) = 1
        nx_pert(id_phi,1,i_pert) = 1
        nx_pert(id_phi,2,i_pert) = ireac
        nx_pert(id_phi,3,i_pert) = i1
        nx_pert(id_phi,4,i_pert) = j2
        nx_pert(id_phi,5,i_pert) = 1
        nx_pert(id_phi,6,i_pert) = 6   !id1
        nx_pert(id_phi,7,i_pert) = 7   !id2
        nx_pert(id_phi,8,i_pert) = 1   !ipi(2,i1)
        nx_pert(id_phi,9,i_pert) = -1  !iz2
          if(p_pert(id_phi,0,i_pert).lt.0.05) then
            write(*,*) 'phi_KK phimass',p_pert(id_phi,0,i_pert)
            stop
          end if
c        write(*,*) 'new phi from K+K-!'
c-----------------------------------------------------------------------
c           density dependence
        ix = nint(xx2)
        iy = nint(yy2)
        iz = nint(zz2)
c-----------------------------------------------------------------------
  410 continue


*----------------------------------------------------------------------*
c      write(*,*) 'end of phi_KK'
      return

c***********************************************************************
*                                                                      *
      entry phi_pirho(ede,pdx,pdy,pdz, srt, xxx,yyy,zzz,szig0,
     &                i1,i2, iz1,iz2,irun)
* purpose: phi production in pi+rho  -> phi                            *
*       szig0   -  referenz cross section in mb units                  *
*       variables:    1 = pion       2 = pion                          *
*         ede,pdx,- energy, momentum of the pi+pi     sys.(real,input) *
*         srt     - cms energy (mass of the pi+pi     sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*----------------------------------------------------------------------*

c      write(*,*) 'in phi_pirho' , srt

      if(iz1+iz2 .ne. 0)   return
      if(iz1.eq.0)         return
      del = 4.*3.14 * .0044                          !    e-width / scale
      if(abs(srt - xphimas) .gt. del)  return

      ix = nint(xxx)
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz)
      distr = sqrt(xxx**2 + yyy**2 + zzz**2)

*----------------------------------------------------------------------*
      szig = .125 * 0.13 * 73. / szig0   !  scale * branching * sigtot (mb)
c---------------------------------------------
      inkrun = max_pert / num
      i_pert = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
* find a free place for the new phi, or the phi with the smallest prob.
      do 920 kl=1,inkrun                    !        iphi_bb       hw
        i_pert = i_pert + 1
        if (nx_pert(id_phi,1,i_pert) .eq. 0) goto 922
        if (valkapi .gt. p_pert(id_phi,4,i_pert)) then
           inkmin = i_pert
           valkapi  = p_pert(id_phi,4,i_pert)
        endif
  920 continue
      i_pert = inkmin
      if (szig .lt. valkapi) goto 927
  922 continue
c      write(*,*) 'prob. of phi is big enough'
      ireac = 14
      p_pert(id_phi,0,i_pert)= srt
      p_pert(id_phi,1,i_pert)= pdx
      p_pert(id_phi,2,i_pert)= pdy
      p_pert(id_phi,3,i_pert)= pdz
      p_pert(id_phi,4,i_pert)= szig
      r_pert(id_phi,1,i_pert)= xxx
      r_pert(id_phi,2,i_pert)= yyy
      r_pert(id_phi,3,i_pert)= zzz
      r_pert(id_phi,4,i_pert)= denst/rho0
      r_pert(id_phi,5,i_pert)= time
      nx_pert(id_phi,0,i_pert) = 1
      nx_pert(id_phi,1,i_pert) = 2
      nx_pert(id_phi,2,i_pert) = ireac
      nx_pert(id_phi,3,i_pert) = i1
      nx_pert(id_phi,4,i_pert) = i2
      nx_pert(id_phi,5,i_pert) = 1
      nx_pert(id_phi,6,i_pert) = 1
      nx_pert(id_phi,7,i_pert) = 3
      nx_pert(id_phi,8,i_pert) = iz1
      nx_pert(id_phi,9,i_pert) = iz2
          if(p_pert(id_phi,0,i_pert).lt.0.05) then
            write(*,*) 'phi_pirho phimass',p_pert(id_phi,0,i_pert)
            stop
          end if
c      write(*,*) 'new phi from pi + rho!'
c-----------------------------------------------------------------------
c           density dependence
!      if (zzz.le.0.0) then         !only front-hemisphere z<0
      numprodd = numprodd + 1
!      endif
      pkao(8,i_pert)= denst
      pkao(13,i_pert)= time
      bin_dens = int(denst*10.)
      bin_time = int(time*1.)
      p_pirho(bin_dens,bin_time)=p_pirho(bin_dens,bin_time)+1

!       write(54,*) time, denst,ireac, xxx,yyy,zzz
!       co = co + 1
!       write(52,233) time, denst,ireac
      goto  997
c-----------------------------------------------------------------------
 927  continue
      write(*,*) 'end of phi_pirho'
      return
c*-*-*_*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      end

************************************************************************
      real*8 function phi_rhoDtab(mdelta,mrho,dss0)
*
*     calculate rho+Delta -> phi+N cross secton
*     by interpolation from table
*-----------------------------------------------------------------------
      implicit none
      real*8 mdelta,mrho,dss0
      real*8 srt
      real*8 wdelta,wrho,ddelta,drho,dleps,xidelta,xirho,xileps,
     &   leps,lepsmin,lepsmax,adel,arho,aeps,sig
      integer nrho,ndelta,nleps,idelta,irho,ileps,idelta2,irho2,ileps2
      parameter(ndelta=40, nrho=40, nleps=80)
      parameter(wdelta=0.4, wrho=0.4)
      real*8 sig_rd(0:nleps, 0:nrho, 0:ndelta)
c                    eps      rho       delta
      integer  opening
      save  opening, sig_rd
      data opening /0/

      if (opening .eq. 0) then
        opening = 1
        open(59, file='buuinput/sig_rho_delzm.dat', status='old')
        read(59,*) sig_rd
        close(59)
      end if
c      open(12, file='sig_rho_deltest.dat', status='unknown')
c      write(12,700) sig_rd
c 700  format(81e12.4)

      lepsmin = log(.004)
      lepsmax = log(1.)
      leps = log(dss0)
      ddelta = 2.*wdelta/ndelta
      drho   = 2.*wrho/nrho
      dleps  = (lepsmax-lepsmin)/nleps
      xidelta = ndelta/2 + (mdelta-1.232)/ddelta
      xirho   = nrho/2 + (mrho-0.770)/drho
      xileps  = (leps - lepsmin)/dleps
      idelta = nint(xidelta-0.5)
      irho   = nint(xirho-0.5)
      ileps  = nint(xileps-0.5)
      idelta2 = idelta+1
      irho2   = irho+1
      ileps2  = ileps+1
      adel = xidelta-idelta
      arho = xirho-irho
      aeps = xileps-ileps
      if(idelta.lt.0) then
        idelta = 0
        idelta2 = 0
        adel = 0
        write(*,*) 'vigyazz rho+D -> phi+N: outside the table !!'
      end if
      if(idelta.ge.ndelta) then
        idelta = ndelta
        idelta2 = ndelta
        adel = 0
        write(*,*) 'vigyazz rho+D -> phi+N: outside the table !!'
      end if
      if(irho.lt.0) then
        irho = 0
        irho2 = 0
        arho = 0
        write(*,*) 'vigyazz rho+D -> phi+N: outside the table !!'
      end if
      if(irho.ge.nrho) then
        irho = nrho
        irho2 = nrho
        arho = 0
        write(*,*) 'vigyazz rho+D -> phi+N: outside the table !!'
      end if
      if(ileps.lt.0) then
        ileps = 0
        ileps2 = 0
        aeps = 0
        write(*,*) 'vigyazz rho+D -> phi+N: outside the table !!'
      end if
      if(ileps.ge.nleps) then
        ileps = nleps
        ileps2 = nleps
        aeps = 0
        write(*,*) 'vigyazz rho+D -> phi+N: outside the table !!'
      end if

      sig =
     &   + (1.-aeps)*(1.-arho)*(1.-adel)*sig_rd(ileps, irho, idelta)
     &   + (1.-aeps)*(1.-arho)*(adel)*sig_rd(ileps, irho, idelta2)
     &   + (1.-aeps)*(arho)*(1.-adel)*sig_rd(ileps, irho2, idelta)
     &   + (1.-aeps)*(arho)*(adel)*sig_rd(ileps, irho2, idelta2)
     &   + (aeps)*(1.-arho)*(1.-adel)*sig_rd(ileps2, irho, idelta)
     &   + (aeps)*(1.-arho)*(adel)*sig_rd(ileps2, irho, idelta2)
     &   + (aeps)*(arho)*(1.-adel)*sig_rd(ileps2, irho2, idelta)
     &   + (aeps)*(arho)*(adel)*sig_rd(ileps2, irho2, idelta2)

      srt = dss0+0.938+1.020
      if(srt.lt.mdelta+mrho) sig=0
      if (sig .gt. 0.5d0)  sig = .5d0

      phi_rhoDtab = sig
      return
      end

      subroutine sig_del_rho(xdelta, xrho, epsilon, sig)
c     NOT USED !!!!(zm)
      implicit none
      integer ndelta,neps
      real*8 xdelta, xrho, epsilon, sig
      parameter(ndelta =10, neps = 20)
      real*8 sig_rd(0:neps, 0:ndelta, 0:ndelta)
c                    eps      rho       delta
      integer  opening
      save  opening, sig_rd
      data opening /0/
      integer ndel,idel,irho,ieps,iborder,id1,ir1,ie1
      real*8 wdelta,wrho,ddelta,drho,xleps1,xleps2,dleps,xleps,adel
      real*8 arho,aeps,sig0
      if (opening .eq. 0) then
        opening = 1
        open(59, file='buuinput/sig_rho_del.dat', status='old')
        read(59,*) sig_rd
        close(59)
      endif
      ndel   = ndelta / 2
      wdelta = 0.4
      wrho   = 0.4
      ddelta = wdelta / ndel
      drho   = wrho   / ndel
      xleps1  = log(.004d0)  !   e_min
      xleps2  = log(1.d0)    !   e_max
      dleps   = (xleps2-xleps1) / neps
c
      idel = nint( (xdelta - 1.232)/ddelta ) + ndel
      irho = nint( (xrho   - 0.770)/drho   ) + ndel
      xleps = log(epsilon)
      ieps = nint( (xleps-xleps1)/dleps )
      iborder = 1
      if (idel .lt. 0) idel = 0
      if (idel .gt. ndelta) idel = ndelta
      if (irho .lt. 0) irho = 0
      if (irho .gt. ndelta) irho = ndelta
      if (ieps .lt. 0) ieps = 0
      if (ieps .gt. neps) ieps = neps
cc      write(*,*) ' idel, irho, ieps ',  idel, irho, ieps
      adel = (xdelta - 1.232) / ddelta - idel+ndel
      arho = (xrho   - 0.77 ) / ddelta - irho+ndel
      aeps = (xleps - xleps1) / dleps  - ieps
      if (adel .gt .0) then
          id1  = idel + 1
      else
          id1 =  idel - 1
      endif
      id1 = max(0,id1)
      id1 = min(ndelta,id1)
c
      if (arho .gt .0) then
          ir1  = irho + 1
      else
          ir1 =  irho - 1
      endif
      ir1 = max(0,ir1)
      ir1 = min(ndelta,ir1)
c
      if (aeps .gt .0) then
          ie1  = ieps + 1
      else
          ie1 =  ieps - 1
      endif
      ie1 = max(0, ie1)
      ie1 = min(neps, ie1)
c
      sig0 = sig_rd(ieps, irho, idel)
      sig  = sig0 + abs(adel) *( sig_rd(ieps, irho, id1 ) -sig0)
     2            + abs(arho) *( sig_rd(ieps, ir1 , idel) -sig0)
     3            + abs(aeps) *( sig_rd(ie1 , irho, idel) -sig0)
      write(*,*)  '  sig_rho_delta ', xdelta, xrho, epsilon, sig
      return
      end
