************************************************************************
*                                                                      *
      subroutine mesphi(irun,i1,i2,id1,id2,iz2,srt0,pdx,pdy,pdz,ede,
     &                     xxx,yyy,zzz,maxcross)
*                                                                      *
*     purpose: omega baryon elastic/inelastic scattering  (henry)      *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
c      include"com_kminu"
      include"com_pert"

      integer i1,i2,maxcross,id2,iz2,maxk,ibar,ireac,inkrun,irun
      integer inkmin,kl,ink,id1,ix,iy,iz
      real*8 pdx,pdy,pdz,ede,xxx,yyy,zzz,sigomegaphi,phi_mass,srt,srt0,s
      real*8 pot1,pot2,vxx,vyy,vzz,pxx,mmass,tmass,srt_th,tmass2,rmass2
      real*8 q0,qq2,qqabs,xz,szig,facreac,valkapi,xx,yy,zz,rr,rn,pbeta
      real*8 qqx,qqy,qqz,transf,gamm,traf,vx,vy,vz,sigetaphi,denst
      parameter(sigomegaphi = 0.07)
      parameter(sigetaphi = 0.07)
      data phi_mass /1.020/

      vx = pdx / ede
      vy = pdy / ede
      vz = pdz / ede

      ibar = 0
      call potcalc_i(i1, ibar, vx, vy, vz, pot1)
      ibar = 1
      call potcalc_i(i2, ibar, vx, vy, vz, pot2)
      srt = srt0 - pot1 - pot2

c      write(20,*) 'omega+n has been found'
      pxx =  .0
      if (iphi_pot .ge. 1) then
         call graduphi(-1, xxx,yyy,zzz, pxx,pxx,pxx,mmass,
     1                      vx,vy,vz,vxx,vyy,vzz)
       else
         mmass = phi_mass
      endif
      tmass  = mmass
      srt_th  = rmass + tmass

c      write(20,*) srt, srt0
      if(srt.le.srt_th)                                  return
c      write(*,*) 'srt is above threshold'
      tmass2= tmass**2
      rmass2= rmass**2
      s     = srt**2
c       q0  : energy of the phi_
      q0      = 0.5 * (s + tmass2 - rmass2) / srt
      if(q0 .le. tmass)                                return
c      write(*,*) 'phi energy is enough'
      qq2     = q0**2 - tmass2
      qqabs   = sqrt(qq2)
      gamm  = ede / srt
      traf  = gamm / (gamm+1.0)
*----------------------------------------------------------------------*
      facreac = 1.0
      xx = srt - srt_th
      if(id1.eq.5) then
        if(id2.eq.1) then         ! omega+n -> phi+n      zm
          ireac = 28
!          szig = sqrt(xx) * 0.05275/(1.-0.8732*meff1+1.467*xx) *
!     &        exp(1. - 1.507*meff1 - 1.523*xx +
!     &        1.051*meff1**2 + 0.4646*xx**2 - 0.8308*meff1*xx)
          szig = sigomegaphi  * 1.0                      !in  mb
          facreac = 2.0
        else if(id2.eq.2) then    ! omega+D -> phi+n      zm
          ireac = 29
!          szig = phi_rhoDtab(meff2, meff1, xx)
          szig = sigomegaphi    *1.0                  !in mb
          xz = abs (iz2 -0.5) - 1.0
          if (xz .gt. .0)  then
            facreac = 3.0
          else
            facreac = 1.0
          endif
cc         szig = 1.116*(1-1.142*xx)/(1+47.47*xx**0.7894)
        else
          ireac = 30
          szig = sigomegaphi    *1.0                  !in mb
          xz = abs (iz2 -0.5) - 1.0
          if (xz .gt. .0)  then
            facreac = 3.0
          else
            facreac = 1.0
          endif
        end if
      else if(id1.eq.2) then
        if(id2.eq.1) then         ! eta+n -> phi+n      zm
          ireac = 38
          szig = sigetaphi* 1.0                        !in  mb
        else if(id2.eq.2) then    ! eta+D -> phi+n      zm
          ireac = 39
          szig = sigetaphi*1.0                         !in mb
          xz = abs (iz2 -0.5) - 1.0
          if (xz .gt. .0)  then
            facreac = 3.0
          else
            facreac = 1.0
          endif
        else
          ireac = 40
          szig = sigetaphi    *1.0                  !in mb
          xz = abs (iz2 -0.5) - 1.0
          if (xz .gt. .0)  then
            facreac = 3.0
          else
            facreac = 1.0
          endif
        end if
      end if

      szig   = szig/maxcross * facreac

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
  810 continue
      ink = inkmin
      if (szig .lt. valkapi) goto 97
  812 continue
c      write(*,*) 'prob. of phi is big enough'
   70 continue
      xx       = 1. - 2. * rn(iseed)
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 70

      qqx= xx*qqabs/rr
      qqy= yy*qqabs/rr
      qqz= zz*qqabs/rr
*   phi_ momentum in observable system

      pbeta  = vx*qqx + vy*qqy + vz*qqz
      transf = gamm * (traf * pbeta + q0)
      p_pert(id_phi,0,ink)= tmass
      p_pert(id_phi,1,ink)= qqx + vx * transf
      p_pert(id_phi,2,ink)= qqy + vy * transf
      p_pert(id_phi,3,ink)= qqz + vz * transf
      p_pert(id_phi,4,ink)= szig
      r_pert(id_phi,1,ink)= xxx
      r_pert(id_phi,2,ink)= yyy
      r_pert(id_phi,3,ink)= zzz
      r_pert(id_phi,3,ink)= denst/rho0
      r_pert(id_phi,5,ink)= time
      nx_pert(id_phi,0,ink) = 1
      nx_pert(id_phi,2,ink) = ireac
      nx_pert(id_phi,3,ink) = i2
      nx_pert(id_phi,4,ink) = 0
      nx_pert(id_phi,5,ink) = 1
      nx_pert(id_phi,6,ink) = id1
      nx_pert(id_phi,7,ink) = id2
      nx_pert(id_phi,8,ink) = 0  ! iz1
      nx_pert(id_phi,9,ink) = iz2

          if(ink.eq.7105) write(*,*) 'mesphi 7105',ink,tmass
          if(p_pert(id_phi,0,ink).lt.0.05) then
            write(*,*) 'mesphi phimass',p_pert(id_phi,0,ink)
            stop
          end if
c-----------------------------------------------------------------------
 97   continue

      return
c***********************************************************************
      end
