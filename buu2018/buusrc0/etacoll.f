************************************************************************
*                                                                      *
      subroutine etacoll(isu)
*                                                                      *
*     purpose: eta baryon elastic/inelastic scattering  (henry)      *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"

      real*8 em3,e3,piet,croset,sigetaphi
      real*8 crosphi,s,xx,yy,zz,rr,betax,betay,betaz,gamma
      real*8 x1, y1, z1, px1, py1, pz1,em1,em12,e1,e1cm,transf
      real*8 dx, dy, dz, rsqare, px2, py2, pz2, em2, em22, e2, x2, y2,z2
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1,t2,ddlt,b21,dxm,rn
      real*8 e2cm,dxp,denst,distr
      real*8 alp,eps_abs, decide
      integer isu,i,ix,iy,iz
      integer maxk,irun,ink,ini,i1,i2,ii,inew
      integer inp
      real*8 szig0,ede,pdx,pdy,pdz,xxx,yyy,zzz,vx,vy,vz,pot1,pot2,srt
      real*8 srt_0,pxx,mmass,vxx,vyy,vzz,tmass,srt0,rmass2,tmass2,q0
      real*8 qq2,qqabs,gamm,traf,facreac,szig,meff1,xz,valkapi
      real*8 qqx,qqy,qqz,phi_mass,pbeta

      integer id1,id2,id6,iz1,iz2,id62,ibar,ireac,inkrun,inkmin,kl

      parameter(piet = 1.262)        !sig0 = 50mb
      parameter(sigetaphi = 0.07)
      data phi_mass /1.020/

c      return

      call dens_4
      call potcalc

      maxk = maxeta/num

*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        ink = (irun-1) * maxk
        ini = (irun-1) * maxb

c        do 800 ii  = 1,maxk
c          i1  = ii + ink
        do 800 ii  = 1,maxp
          i1  = ii + inp
          if(ipi(1,i1).ne. 2)  goto 800               !eta
c          if(nx_eta(0,i1).eq. 0)  goto 800

c          x1  = r_eta(1,i1)
c          y1  = r_eta(2,i1)
c          z1  = r_eta(3,i1)
c          px1 = p_eta(1,i1)
c          py1 = p_eta(2,i1)
c          pz1 = p_eta(3,i1)
c          em1 = emass
c          em12= em1**2
c          e1  = sqrt( em12 + px1**2 + py1**2 + pz1**2 )

          x1  = rpi(1,i1)
          y1  = rpi(2,i1)
          z1  = rpi(3,i1)
          px1 = ppi(1,i1)
          py1 = ppi(2,i1)
          pz1 = ppi(3,i1)
          em1        = epi(i1)
          meff1      = em1 +  mpot(i1)
          e1         = sqrt(meff1**2+px1**2+py1**2+pz1**2)

*     look for a scattering pseudonucleon in the same run
          i2  = ini
  600     i2  = i2 + 1
          if(i2 .gt. ini+maxb)                                 goto 800
          if(id(1,i2) .eq. 0) goto 600

          x2 = r(1,i2)
          dx     = x1 - x2
            if(nbound.eq.1) then
              dxp = dx+2.0*boxx
              dxm = dx-2.0*boxx
              if(abs(dx) .gt. abs(dxp)) dx=dxp
              if(abs(dx) .gt. abs(dxm)) dx=dxm
            end if
          if (abs(dx) .gt. delpi)                          goto 600
          y2 = r(2,i2)
          dy     = y1 - y2
            if(nbound.eq.1) then
              dxp = dy+2.0*boxx
              dxm = dy-2.0*boxx
              if(abs(dy) .gt. abs(dxp)) dy=dxp
              if(abs(dy) .gt. abs(dxm)) dy=dxm
            end if
          if (abs(dy) .gt. delpi)                          goto 600
          z2 = r(3,i2)
          dz     = z1 - z2
            if(nbound.eq.1) then
              dxp = dz+2.0*boxz
              dxm = dz-2.0*boxz
              if(abs(dz) .gt. abs(dxp)) dz=dxp
              if(abs(dz) .gt. abs(dxm)) dz=dxm
            end if
          if (abs(dz) .gt. delpi)                          goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. delpi**2)                        goto 600
*         now particles are close enough to each other !

          px2    = p(1,i2)
          py2    = p(2,i2)
          pz2    = p(3,i2)
          em2    = e(i2)
          em22   = e(i2)**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                           - (pz1+pz2)**2
          eps_abs = sqrt(s) - em1 - em2

          croset = 10.0 * pi * piet**2
          srt = sqrt(s)

*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if (brel .gt. piet)                                goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
*   now  the phi will collide or be absorbed in this time step


             szig0 = croset
             ede = e1 + e2
             pdx = px1+px2
             pdy = py1+py2
             pdz = pz1+pz2
             xxx = r(1,i2)
             yyy = r(2,i2)
             zzz = r(3,i2)
c             id1 = nx_eta(0,i1)              !eta from eta_dbb = 1
             id1 = 1
             id2 = id(1,i2)
             id6 = id(6,i2)
c             iz1 = nx_eta(5,i1)              !el. charge of eta0 = 0
             iz1 = 0

             iz2 = id(2,i2)
             id62 = 1             !          not used

               vx = pdx / ede
               vy = pdy / ede
               vz = pdz / ede

      ibar = 0
      call potcalc_i(i1, ibar, vx, vy, vz, pot1)
      ibar = 1
      call potcalc_i(i2, ibar, vx, vy, vz, pot2)
          srt_0 = srt - pot1 - pot2
          srt = srt_0

      if(iz1+iz2 .gt. 1)                                   goto 600
      if(iz1+iz2 .lt. 0)                                   goto 600
      if(id1.ne.1)                                         goto 800  !eta from eta_dbb
      if((id2.ne.1).and.(id2.ne.2))                        goto 600

c      write(20,*) 'eta+n has been found'
      pxx =  .0
      if (iphi_pot .ge. 1) then
        call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                      vx,vy,vz,vxx,vyy,vzz)
      else
        mmass = phi_mass
      endif
      tmass  = mmass
      srt0  = rmass + tmass

      if(srt.le.srt0)                            goto 600
c      write(*,*) 'srt is above threshold'
      tmass2= tmass**2
      rmass2= rmass**2
      s     = srt**2
c       q0  : energy of the phi_
      q0      = 0.5 * (s + tmass2 - rmass2) / srt
      if(q0 .le. tmass)                          goto 600
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
      if(id2.eq.1) then         ! eta+n -> phi+n      zm
         ireac = 38
         szig = sigetaphi* 1.0                        !in  mb
         if (iz1 .ne. 0) facreac = 2.0
      else if(id2.eq.2) then    ! eta+D -> phi+n      zm
         ireac = 39
         szig = sigetaphi*1.0                         !in mb
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
      else
        szig= .0
      end if

      szig   = szig/szig0 * facreac

c***********************************************************************
      ix = nint(xxx)            ! nint: closest integer
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &     denst = rhb(ix,iy,iz)

      inkrun = max_pert / num
      ink = (irun-1) * inkrun
      inkmin = 0
      valkapi  = 1000.0
* find a free place for the new phi, or the phi with the smallest prob.
      do 810 kl=1,inkrun                    !        iphi_bb       hw
        ink = ink + 1
        if (nx_pert(id_phi,0,ink) .eq. 0) goto 812
        if (valkapi .gt. p_pert(id_phi,4,ink)) then
           inkmin = ink
           valkapi  = p_pert(id_phi,4,ink)
        endif
  810  continue
      ink = inkmin
      if (szig .lt. valkapi) goto 97
  812  continue
c      write(*,*) 'prob. of phi is big enough'
   70  continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 70

          qqx= xx*qqabs/rr
          qqy= yy*qqabs/rr
          qqz= zz*qqabs/rr
*   phi_ momentum in observable system

          pbeta  = betax*qqx + betay*qqy + betaz*qqz
          transf = gamm * (traf * pbeta + q0)
          p_pert(id_phi,0,ink)= tmass
          p_pert(id_phi,1,ink)= qqx + betax * transf
          p_pert(id_phi,2,ink)= qqy + betay * transf
          p_pert(id_phi,3,ink)= qqz + betaz * transf
c          p_pert(id_phi,4,ink)= szig * p_eta(4,i1)
          p_pert(id_phi,4,ink)= szig
          r_pert(id_phi,1,ink)= xxx
          r_pert(id_phi,2,ink)= yyy
          r_pert(id_phi,3,ink)= zzz
          r_pert(id_phi,4,ink)= denst
          r_pert(id_phi,5,ink)= time
          nx_pert(id_phi,0,ink) = 1
          nx_pert(id_phi,2,ink) = ireac
          nx_pert(id_phi,3,ink) = i2
          nx_pert(id_phi,4,ink) = 0
          nx_pert(id_phi,5,ink) = 1
          nx_pert(id_phi,6,ink) = id1
          nx_pert(id_phi,7,ink) = id2
          nx_pert(id_phi,8,ink) = iz1
          nx_pert(id_phi,9,ink) = iz2

          if(ink.eq.7105) write(*,*) 'etacoll 7105',ink,tmass
          if(p_pert(id_phi,0,ink).lt.0.05) then
            write(*,*) 'etacoll phimass',p_pert(id_phi,0,ink)
            stop
          end if
          !this eta was absorbed by a nucleon:
          p_eta(4,i1)= 0.
          nx_eta(0,i1) = 0



  800   continue
 1000 continue

c-----------------------------------------------------------------------
 97   continue


c      write(20,*) "phi from dbb: ", tfg

      return
c***********************************************************************
      end


