
************************************************************************
*                                                                      *
      subroutine gradupi(rx, ry, rz, px, py, pz, etot, id1, id2, mass,
     &                   gradxr,gradyr,gradzr, gradxp, gradyp,gradzp,
     &                   gradm, ddt,dtprime,ipropag)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))                          *
*       variables:                                                     *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy o meson       (real, input)     *
*         id2                 - charge of meson      (integer, input)  *
*         gradx, grady, gradz - gradient of u        (real,output)     *
*         ddt    - fraction of sub-timestep                            *
*         ipropag         - 1: spectral func propag; else Hamilton eq  *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"
      include"com_pert"

      integer id1, id2, ipropag,idjpsi0,idjpsi1,idjpsi2

      real*8 rx, ry, rz, px, py, pz, etot,ddt
      real*8 emfox, emfoy, emfoz, cpot

      real*8 gradxr, gradyr, gradzr, gradxp, gradyp, gradzp,gradm
      real*8 step,mstp,mass,dgdx,dgdy,dgdz,dgdpx,dgdpy,dgdpz,dgde,dgdt
      real*8 dsdx,dsdy,dsdz,dsdpx,dsdpy,dsdpz,dsde,dsdt
      real*8 mass0

      real*8 dddt,dtprime,p_max,mstp_max,mstpp
      real*8 m1,m2,en,cfac,dfac, fac, cfac0
      complex dcdx,dcdy,dcdz,dcdpx,dcdpy,dcdpz,dcde,dcdt,self
      complex propa_self_mes

c      write(*,*) 'in gradupi, id1: ',id1,rx, ry, rz, px, py, pz, etot, 
c     &      id1, id2, mass,
c     &                   gradxr,gradyr,gradzr, gradxp, gradyp,gradzp,
c     &                   gradm, ddt,dtprime,ipropag
      call f77flush()

c      i1hw = idn

      gradxr = 0.0
      gradyr = 0.0
      gradzr = 0.0
      gradxp = 0.0
      gradyp = 0.0
      gradzp = 0.0
      gradm  = 0.0

      idjpsi0= 100+id_JPsi(1)
      idjpsi1= 100+id_JPsi(2)
      idjpsi2= 100+id_JPsi(3)
      if (id1.eq.3 .or. id1.eq.5 .or. ! rho or omega
     &     id1.eq.idjpsi0 .or. id1.eq.idjpsi1 .or. id1.eq.idjpsi2) then ! JPsi 3 states
        en = sqrt(mass**2+px**2+py**2+pz**2)

        if (ipropag.eq.1) then
          if (id1.eq.3) mass0 = romas
          if (id1.eq.5) mass0 = omass
          if (id1.eq.idjpsi0) mass0 = JPsi_prop(1,1)
          if (id1.eq.idjpsi1) mass0 = JPsi_prop(2,1)
          if (id1.eq.idjpsi2) mass0 = JPsi_prop(3,1)
          step = 1.
          mstp = .025
          dddt = 0.5*dtprime/dt

c      write(*,*)  ' in gradupi self1 ', mass, px,py,pz
          self = propa_self_mes(id1,mass,px,py,pz,rx,ry,rz,ddt)
c      write(*,*)  ' in gradupi self2 ', mass, px,py,pz
          if (abs(imagpart(self)).lt..00001) then
c            write(*,*) 'warning in gradupi: small rho width!!',self
            return
          end if
c          fac = .0
c         else
          fac = (mass**2-mass0**2-realpart(self)) / imagpart(self)
          dcdx = (propa_self_mes(id1,mass,px,py,pz,rx+step,ry,rz,ddt)
     &       - propa_self_mes(id1,mass,px,py,pz,rx-step,ry,rz,ddt)) /
     &       (2.*step)
          dcdy = (propa_self_mes(id1,mass,px,py,pz,rx,ry+step,rz,ddt)
     &       - propa_self_mes(id1,mass,px,py,pz,rx,ry-step,rz,ddt)) /
     &       (2.*step)
          dcdz = (propa_self_mes(id1,mass,px,py,pz,rx,ry,rz+step,ddt)
     &       - propa_self_mes(id1,mass,px,py,pz,rx,ry,rz-step,ddt)) /
     &       (2.*step)

          dsdx = realpart(dcdx)
          dgdx = imagpart(dcdx)
          dsdy = realpart(dcdy)
          dgdy = imagpart(dcdy)
          dsdz = realpart(dcdz)
          dgdz = imagpart(dcdz)

c      E, py, pz fixed, px is constrained by the requirement, that mass>0
          mstp_max = sqrt(mass**2+px**2) - abs(px)
          mstpp = min(mstp,0.2*mstp_max)
          m1 = sqrt(mass**2+px**2-(px+mstpp)**2)
          m2 = sqrt(mass**2+px**2-(px-mstpp)**2)
c------------------------------------------------------------
c         write(*,*) ' in gradupi ', idn, id1, px,py,pz, mass,
c    1                               m1,m2,mstpp
c
c         if (idn .eq. 9201) then
c          cc1 = propa_self_mes(id1,m1,px+mstpp,py,pz,rx,ry,rz,ddt)
c          cc2 = propa_self_mes(id1,m2,px-mstpp,py,pz,rx,ry,rz,ddt)
c          write(45,*)  ' propa_self ', m1, cc1, m2, cc2
c          write(45,*) ' p-xyz ', mstp, mstpp, mass, px,py,pz
c          call f77flush()
c         endif
c-----------------------------------------------------------
          dcdpx = (propa_self_mes(id1,m1,px+mstpp,py,pz,rx,ry,rz,ddt)
     &       - propa_self_mes(id1,m2,px-mstpp,py,pz,rx,ry,rz,ddt)) /
     &       (2.*mstpp)

          mstp_max = sqrt(mass**2+py**2) - abs(py)
          mstpp = min(mstp,0.2*mstp_max)
          m1 = sqrt(mass**2+py**2-(py+mstpp)**2)
          m2 = sqrt(mass**2+py**2-(py-mstpp)**2)
c          write(*,*) 'gradupi 3 y: ',m1,m2, mstpp
c          call f77flush()
          dcdpy = (propa_self_mes(id1,m1,px,py+mstpp,pz,rx,ry,rz,ddt)
     &       - propa_self_mes(id1,m2,px,py-mstpp,pz,rx,ry,rz,ddt)) /
     &       (2.*mstpp)

          mstp_max = sqrt(mass**2+pz**2) - abs(pz)
          mstpp = min(mstp,0.2*mstp_max)
          m1 = sqrt(mass**2+pz**2-(pz+mstpp)**2)
          m2 = sqrt(mass**2+pz**2-(pz-mstpp)**2)
c          write(*,*) 'gradupi 4 z: ',m1,m2, mstpp
c          call f77flush()
          dcdpz = (propa_self_mes(id1,m1,px,py,pz+mstpp,rx,ry,rz,ddt)
     &       - propa_self_mes(id1,m2,px,py,pz-mstpp,rx,ry,rz,ddt)) /
     &       (2.*mstpp)

          dsdpx = realpart(dcdpx)
          dgdpx = imagpart(dcdpx)
          dsdpy = realpart(dcdpy)
          dgdpy = imagpart(dcdpy)
          dsdpz = realpart(dcdpz)
          dgdpz = imagpart(dcdpz)

c         mstp_max = en - sqrt(px**2+py**2+pz**2)
c         mstpp = min(mstp,0.2*mstp_max)
c         m1 = sqrt(mass**2 + 2.*en*mstpp + mstpp**2)
c         m2 = sqrt(mass**2 - 2.*en*mstpp + mstpp**2)

          mstpp = .01 * mass
          m1 = mass + mstpp
          m2 = mass - mstpp
c          write(*,*) 'gradupi 5m: ',m1,m2,mstpp
c          call f77flush()

c          if (id1 .eq.3)   write(*,*) 'gradupi start dcde: '
c
          dcde = (propa_self_mes(id1,m1,px,py,pz,rx,ry,rz,ddt)
     &          - propa_self_mes(id1,m2,px,py,pz,rx,ry,rz,ddt))
     &       / (2.*mstpp)    * en/mass
c           write(*,*) 'gradupi 6end: '
           call f77flush()

          dsde = realpart(dcde)
          dgde = imagpart(dcde)

          dcdt = (propa_self_mes(id1,mass,px,py,pz,rx,ry,rz,ddt+dddt)
     &       - propa_self_mes(id1,mass,px,py,pz,rx,ry,rz,ddt-dddt)) /
     &       dtprime
c          write(*,*) 'gradupi 7end t: '
           call f77flush()

          dsdt = realpart(dcdt)
          dgdt = imagpart(dcdt)

          dfac =  (dsde + fac*dgde)/(2.*en)

ccccccc proba WGy CCCCCCCCCCCCCCC
c          dfac = 0.0
          if (dfac .gt. .5)  then
              cfac = 1.33333 * (1. + dfac) / (2.*en)
          else
              cfac = 1. / ((1. - dfac)*2.*en)
          endif
c---------------------
c          if (dfac .gt. 0.7 )
c     1    write(*,*) 'dfac in gradupi: ',id1,dfac,cfac,en,
c     2                fac, dsde,dgde, mass, px,py,pz
          call f77flush()
c---------------------  mesons below threshold -----------
          if (id1 .eq. 3  .and.  mass .lt. 2.*pmass) then
c               fac = .0
               cfac = 1.0 / (2.*en)
          endif
c-----------------------------------------------------
          gradxp = cfac * (2.*px + dsdpx + fac*dgdpx)
          gradyp = cfac * (2.*py + dsdpy + fac*dgdpy)
          gradzp = cfac * (2.*pz + dsdpz + fac*dgdpz)
          gradxr = cfac * (dsdx + fac * dgdx)
          gradyr = cfac * (dsdy + fac * dgdy)
          gradzr = cfac * (dsdz + fac * dgdz)
c----------
          gradm  = cfac * (dsdt + fac * dgdt + px/en*(dsdx+fac*dgdx) +
     &       py/en*(dsdy+fac*dgdy) + pz/en*(dsdz+fac*dgdz))

c          write(*,111) id1, mass
c     1     gradm, cfac, dfac,fac,dsdt,dsdx,dsdy,dsdz,
c     1     dgdt, dgdx, dgdy,dgdz,
c     1     realpart(self),imagpart(self),px/en, py/en, pz/en,
c     1     gradxp,gradyp,gradzp,gradxr,gradyr,gradzr,gradm
          call f77flush()
 111      format(2i3,2f6.3,' gradm:',4f8.4,' dsdx:',4f8.4,' dgdx ',
     1           4f8.4,' self ',2f8.4,' pp ',3f7.3, ' grad ', 7f10.6)
c-----------
c         gradm  = cfac * (dfac*dsdt + fac * dgdt + px/en*
c    1             (dfac*dsdx+fac*dgdx) + py/en*(dfac*dsdy+fac*dgdy) +
c    2              pz/en*(dfac*dsdz+fac*dgdz))


c         if (idn .eq. 201)
c    1    write(*,*) 'gradupi 201: ',mass, id1,
c    2    px,py,pz,en,cfac,fac, dsdt, dgdt,  gradm
c          write(*,*) 'gradupi 7.: ',dgdx,dgdy,dgdz,dgdt,
c     &       dgdpx,dgdpy,dgdpz,dgde
c          write(*,*) 'gradupi 8.: ',gradxr,gradyr,gradzr,
c     &       gradxp,gradyp,gradzp,gradm
c          call f77flush()
        else
           gradxp = px/en
           gradyp = py/en
           gradzp = pz/en
        end if
      end if

c the following coulomb part is inconsistent with the previous part!!! (zm)
*----------------------------------------------------------------------*
*       call routine emfoca to calculate the electro-magnetic force    *
*       emfox, emfoy, emfoz  ( link to Coulomb-routine )               *

      if(ipou .eq. 1    .and.  id2.ne.0) then
        call emfoca(rx, ry, rz, id2,
     &             emfox, emfoy,emfoz,ncont,cpot)

        gradxr = gradxr - emfox
        gradyr = gradyr - emfoy
        gradzr = gradzr - emfoz
c       if (rx**2+ry**2 .lt. 3.)
c    1      write(*,*)  '   in gradupi - coul ', rz, gradxr, emfox
      end if
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function propa_self_mes(idm,mass,p1,p2,p3,x1,x2,x3,ddt)
      implicit none
      include"common"
      include"com_pert"

      real*8 density, sgam

      integer idm,idjpsi0,idjpsi1,idjpsi2
      real*8 mass,p1,p2,p3,x1,x2,x3,ddt
      complex propa_self_mes

      real*8 en,deriv(0:4),j0,j1,j2,j3,vrel,jj0,jj1,jj2,jj3
      real*8 betacm(1:3),betalrf(1:3)
      real*8 rself,iself
c      write(*,*) 'in propa_self_mes: ',mass,p1,p2,p3,x1,x2,x3,ddt
      call f77flush()
*---------------------------------------------------------------------*
*     beta that boosts from calc.frame into the resframe of resonance *
c      write(*,*) 'propa_self_mes: ',idm,mass,p1,p2,p3,x1,x2,x3,ddt
      en = sqrt(p1**2+p2**2+p3**2+mass**2)
      betacm(1) =   p1 / en
      betacm(2) =   p2 / en
      betacm(3) =   p3 / en

*---------------------------------------------------------------------*
*     determine the denisty in the LRF  and the trafo-variables for   *
*     the boost into the LRF                                          *
      call linint1(x1,x2,x3,deriv)
      j0    = deriv(0)
      j1    = deriv(1)
      j2    = deriv(2)
      j3    = deriv(3)

      call linint_prev(x1,x2,x3,deriv)
      jj0    = deriv(0)
      jj1    = deriv(1)
      jj2    = deriv(2)
      jj3    = deriv(3)

      j0 = jj0*(1.-ddt) + j0*ddt
      j1 = jj1*(1.-ddt) + j1*ddt
      j2 = jj2*(1.-ddt) + j2*ddt
      j3 = jj3*(1.-ddt) + j3*ddt

c     call f77flush()

      density = sqrt(j0**2-j1**2-j2**2-j3**2)

      if (density.lt.1.0e-3) then
        density=0.
        vrel = 0.
      else
        call f77flush()
c        write(*,*) 'lorentz call gradupi',betacm(1),betacm(2),betacm(3),
c     &     j1,j2,j3,j0
              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba gradupi lorentz, mass<0",j0,j1,j2,j3
                 stop
              end if
        call lorentz(betacm(1),betacm(2),betacm(3),j1,j2,j3,j0)
        vrel = sqrt(j1**2+j2**2+j3**2)/j0
      end if

      if (vrel.gt.1.0) then
        write(*,*) 'vrel>0 in propa_self_mes: ',idm,mass,p1,p2,p3,
     &     x1,x2,x3,ddt
        write(*,*) 'four densities: ',j0,j1,j2,j3,jj0,jj1,jj2,jj3
      end if

      idjpsi0= 100+id_JPsi(1)
      idjpsi1= 100+id_JPsi(2)
      idjpsi2= 100+id_JPsi(3)

      if (idm.eq.3) then
        call self_rho(vrel,density,mass,rself,iself,sgam)
c        write(*,*) 'rho selfen. in propa_self_mes: ',rself,iself,
c    1               vrel,density,mass,rself,iself,sgam
      else if (idm.eq.5) then
        call self_omega(vrel,density,mass,rself,iself, sgam)
c        write(*,*) 'omega selfen. in propa_self_mes: ',rself,iself
      else if (idm.eq.100+id_JPsi(1)) then
        call self_JPsi(1,vrel,density,rself,iself, sgam)
c        write(*,*) 'JPsi selfen. in propa_self_mes: ',rself,iself
      else if (idm.eq.100+id_JPsi(2)) then
        call self_JPsi(2,vrel,density,rself,iself, sgam)
c        write(*,*) 'Psi1 selfen. in propa_self_mes: ',rself,iself
      else if (idm.eq.100+id_JPsi(3)) then
        call self_JPsi(3,vrel,density,rself,iself, sgam)
c        write(*,*) 'Psi2 selfen. in propa_self_mes: ',rself,iself
      else
        write(*,*) 'unknown meson in propa_self_mes, idm=',idm
c        stop
      end if

      propa_self_mes =  cmplx(rself, -iself)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
