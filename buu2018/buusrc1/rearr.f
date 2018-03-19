      subroutine rearr(icomp,rho,plrf,rhap,rearsky, rearmd)

*        this function provides the rearragement terms
c       for the threefluid model, uncomment the lines in the if
c       loop ithree.eq.1


      implicit none

      include"common"
      include"commonthreef"
      include"potparam"

      integer icomp, ithree

      real*8   as, bs, cs, rho, taus, lam

      real*8   rearsky, rearmd, uproprho, urear
      real*8   uges,plrf(1:3), rhap(1:3)
      real*8   mdpart, pin

      if(icomp.gt.maxset) then
        write(*,*) 'problem in potanal '
        write(*,*)'icomp = ',icomp
        stop
      end if


        as   = a(icomp)
        bs   = b(icomp)
        taus = tau(icomp)
        cs   = c(icomp)
        lam = la(icomp)

        rearsky = 0.0
        rearmd  = 0.0


*     the constnts a,b,c convert then the potential in gev



*************** do Skyrme part **************************************
          uproprho = as*rho/rho0 + 2.*bs/(taus+1)*rho**taus/rho0**taus

          urear    = bs*(taus-1)/(taus+1)*rho**taus/rho0**taus
          rearsky  = 0.5*uproprho + urear

**********************************************************************
          if(icomp.ne.3 .and. icomp.ne.4 .and. icomp.ne.6) then
**************** evaluate mom-dep part ******************************
*------ the funcion potanal needs the momentumin GEV!!!!!

            uges = 0.0
c            if(ithree.eq.0) then
              pin = sqrt(plrf(1)**2 + plrf(2)**2 + plrf(3)**2)
              uges = mdpart(icomp,rho0,rho,pin)
c            else if(ithree.eq.1) then
c
c
c*         add the momentum-dependent part according to the 3-fl.-mod *
c*         using the funcion mdpart
c*         mdpart needs inputs in GeV and gives back the md-part in GeV
c
c***       locations where the potential is to be calculated
c
c              ix = nint(rhap(1))
c              iy = nint(rhap(2))
c              iz = nint(rhap(3))
c
c
c***       eval. potential stemming from target-like nucleons
c
c             pxtar = pprota(1,1,ix,iy,iz)
c             pytar = pprota(1,2,ix,iy,iz)
c             pztar = pprota(1,3,ix,iy,iz)

c             dentar= rho*nsphere(1,ix,iy,iz)
c             pintar= (plrf(1)-pxtar)**2 + (plrf(2)-pytar)**2 +
c     &               (plrf(3)-pztar)**2
c             pintar= sqrt(pintar)
c
c             pottar= mdpart(icomp,rho0, dentar,pintar)
c
c             if(dentar.gt.1.0e-12) then
c               write(*,*)'rearr2 ', dentar
c               stop
c             end if
c
c
c***       eval. potential stemming from projectile-like nucleons
c
c             pxpro = pprota(2,1,ix,iy,iz)
c             pypro = pprota(2,2,ix,iy,iz)
c             pzpro = pprota(2,3,ix,iy,iz)
c
c             pinpro= (plrf(1)-pxpro)**2 + (plrf(2)-pypro)**2 +
c     &               (plrf(3)-pzpro)**2
c             pinpro= sqrt(pinpro)
c
c             denpro= rho*nsphere(2,ix,iy,iz)
c             potpro= mdpart(icomp,rho0, denpro,pinpro)
c
c             if(denpro.gt.1.0e-12) then
c               write(*,*)'rearr1 ', denpro
c               stop
c             end if
c
c
c***       eval. potential stemming from lrf-like nucleons
c
c             pinlrf = sqrt(plrf(1)**2 + plrf(2)**2 + plrf(3)**2)
c
c             denlrf= rho*nsphere(3,ix,iy,iz)
c
c             potlrf= mdpart(icomp,rho0, denlrf,pinlrf)
c
c**************** add up contributions
c             uges = potpro + pottar + potlrf
c
c           end if

c          write(*,*)'rearr 3', icomp
              rearmd = 0.5*uges
          end if

      return
      end


