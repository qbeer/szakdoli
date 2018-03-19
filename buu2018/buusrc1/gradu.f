************************************************************************
*                                                                      *
      subroutine gradu(rx, ry, rz, px, py, pz, etot, id1, id2, idn,
     &                 gradxr, gradyr, gradzr,
     &                 gradxp, gradyp, gradzp, nt, isw )
*                                                                      *
*       purpose:     determine grad(u(x,y,z))                          *
*       variables:                                                     *
*         ipou                - choice of coulomb      (integer,common)*
*         rx, ry, rz          - coordinates of grid        (real,input)*
*         gradx, grady, gradz - gradient of u             (real,output)*
*         id1                 - particle type           (integer,input)*
*         id2                 - particle charge         (integer,input)*
*         px,py,pz            - momentum of the particle   (real,input)*
*         etot                - total energy              (real,common)*
*                                                                      *
************************************************************************
      implicit none
      real*8 udel0
      parameter(udel0=-0.03)
      include"common"
      include"cominput"

      real*8    epsi
      parameter(epsi = 1.0e-03)

      integer id1,id2, idn, k, isw
      integer ix, iy, iz, jx, jy, jz, kx, ky, kz

      real*8    rx, ry, rz, px, py, pz, etot
      real*8    dx, dy, dz
      real*8    tden, tdenxp, tdenxm, tdenyp, tdenym, tdenzm, tdenzp
      real*8    rxderv, ryderv, rzderv
      real*8    dxderv, dyderv, dzderv
      real*8    constn, constd
      real*8    gradxr, gradyr, gradzr
      real*8    gradxp, gradyp, gradzp
      real*8    udelt
      real*8    yxplus, yyplus, yzplus, yxmins, yymins, yzmins
      real*8    emfox, emfoy, emfoz
      real*8    j0, j1, j2, j3, rho, potanal
      real*8    p1x, p1y, p1z, pabs
      real*8    pboo(0:3)
      real*8    potgradr(1:12)
      real*8    dp
      integer ib
      real*8    pst
      real*8    deriv(0:4,1:13), dd,sd, deriv1(0:4),ort(1:3,1:13)
      real*8    potn, potd, potnd, vecpot, masse
      integer iloop
      real*8    coord(1:13,1:3), scapo(1:13)
      real*8    cpot, betax, betay, betaz

      real*8    pxneu, pyneu, pzneu
      real*8    pxstat, pystat, pzstat
      real*8 pxtest, pytest, pztest
      real*8    scold, scneu, ehelp
      real*8    plrf(1:3), rhap(1:3)

      logical flagit
      integer itcount

      integer nt
*     delta p needed to evaluate the derivative wrsp. to the momentum
      parameter(dp = 0.01)


*---------------------------------------------------------------------*
*     determing force using the normal buu-grid (1 fm)                *
*                                                                     *
      ix = nint(rx)
      iy = nint(ry)
      iz = nint(rz)
*---------------------------------------------------------------------*
*     set gradients to zero; if particle is out off grid then these   *
*     zeros are the return values                                     *
*                                                                     *
      gradxr = 0.0
      gradyr = 0.0
      gradzr = 0.0
      gradxp = 0.0
      gradyp = 0.0
      gradzp = 0.0
*---------------------------------------------------------------------*
*     set electromagn. force to zero                                  *
      emfox = 0.0
      emfoy = 0.0
      emfoz = 0.0

*


      dx=2.0
      jx=ix+1
      kx=ix-1
      dy=2.0
      jy=iy+1
      ky=iy-1
      dz=2.0
      jz=iz+1
      kz=iz-1

*
c      if(ipronew .eq. 0) then
c        write(*,*)'better date up the old version of propagation '
c        stop
c        if( abs(ix).lt.maxx .and. abs(iy).lt.maxx .and.
c     &    abs(iz).lt.maxz)then
c
c
c          tden   =   rhb(ix,iy,iz)
c          tdenxp =   rhb(jx,iy,iz)
c          tdenxm =   rhb(kx,iy,iz)
c          tdenyp =   rhb(ix,jy,iz)
c          tdenym =   rhb(ix,ky,iz)
c          tdenzp =   rhb(ix,iy,jz)
c          tdenzm =   rhb(ix,iy,kz)
c
c           tden   = tden  / rhob_4(6,ix,iy,iz)
c            tdenxp = tdenxp/ rhob_4(6,jx,iy,iz)
c            tdenxm = tdenxm/ rhob_4(6,kx,iy,iz)
c            tdenyp = tdenyp/ rhob_4(6,ix,jy,iz)
c            tdenym = tdenym/ rhob_4(6,ix,ky,iz)
c            tdenzp = tdenzp/ rhob_4(6,ix,iy,jz)
c            tdenzm = tdenzm/ rhob_4(6,ix,iy,kz)
*
c
c          rxderv = (tdenxp-tdenxm)/dx
c          ryderv = (tdenyp-tdenym)/dy
c          rzderv = (tdenzp-tdenzm)/dz
c
*-----------------------------------------------------------------------
*
c          dxderv = (tdenxp**(1.+rpot)-tdenxm**(1.+rpot))/dx
c          dyderv = (tdenyp**(1.+rpot)-tdenym**(1.+rpot))/dy
c          dzderv = (tdenzp**(1.+rpot)-tdenzm**(1.+rpot))/dz
c          constn= 3./4.*tt0
c          constd= 3./8.*tt3*(2.+rpot)
c          gradxr = constn*rxderv+constd*dxderv
c          gradyr = constn*ryderv+constd*dyderv
c          gradzr = constn*rzderv+constd*dzderv
c*
c          if(idelpot.eq.1 .and. id1.eq.2) then
c            udelt = udel0/u00
c            gradxr = udelt*gradxr
c            gradyr = udelt*gradyr
c            gradzr = udelt*gradzr
c          endif
*
c
c          if(idelpot.eq.2 .and. id1.eq.2) then
c            gradxr = udel0*rxderv
c            gradyr = udel0*ryderv
c            gradzr = udel0*rzderv
c          endif
c
c          call epot(rx,ry,rz,potn,potd,potnd)
c          write(*,*) 'gradu epot',rx,ry,rz,potn,potd,potnd
c          if(id1.eq.2) then
c            upot(idn) = potd
c          else if(id1.ne.2) then
c            upot(idn) = potn
c          end if
c
c        end if
*-------------------------------------------------------------------

*-------------------------------------------------------------------
*

        if(isplipi .eq. 1 ) then
          call splinint(rx,ry,rz,deriv,dd)
        else if(isplipi .eq. 0 ) then
          if(isw.eq.1) then
            call linint(rx,ry,rz,deriv,dd,ort)
          else if(isw.eq.2) then
            call linint1(rx,ry,rz,deriv1)
            dd = 1.0
            do k = 0,4
              deriv(k,1) = deriv1(k)
             end do
           end if
c            write(*,*) 'gradu linint',isw,isplipi,ipronew,rx,ry,rz,
c     &    deriv(0,1),deriv(1,1),deriv(2,1),deriv(3,1),deriv(4,1),dd,ort
        end if


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
*     evaluate the derivatives wrsp to r                              *



        iloop = 0
        if(iorder .eq. 1) then
          iloop = 7
        else if(iorder .eq. 2) then
          iloop = 13
          write(*,*)'for iorder = 2  something has to be changed in
     &                  linint'
          stop
        end if

        if(isw.eq.1) then
          if(isplipi.eq.1) iloop = 13

          rhap(1) = ort(1,1)
          rhap(2) = ort(2,1)
          rhap(3) = ort(3,1)


          do k = 2, iloop

            j0 = deriv(0,k)
            j1 = deriv(1,k)
            j2 = deriv(2,k)
            j3 = deriv(3,k)
            sd = deriv(4,k)
c            write(*,*) 'gradu deriv',j0,j1,j2,j3,sd

************ location where the pot. has to be eval.

            rho = sqrt(j0**2 - j1**2 - j2**2 - j3**2)

            if(j0.lt.1.0e-06) then
              j0 = 0.0
              j1 = 0.0
              j2 = 0.0
              j3 = 0.0
              betax = 0.0
              betay = 0.0
              betaz = 0.0
            else
              betax = j1/j0
              betay = j2/j0
              betaz = j3/j0
            end if

            pboo(1) = px
            pboo(2) = py
            pboo(3) = pz
            pboo(0) = sqrt((e(idn)+upot(idn))**2+px**2+py**2+pz**2)

            pxstat = pboo(1)
            pystat = pboo(2)
            pzstat = pboo(3)
            scneu = upot(idn)
            masse = e(idn)

            flagit = .false.
            itcount = 0
            do while(.not.flagit)

              itcount = itcount + 1

              pboo(1) = pxstat
              pboo(2) = pystat
              pboo(3) = pzstat
              pboo(0) = sqrt( (masse+scneu)**2 + pboo(1)**2 +
     &                         pboo(2)**2 + pboo(3)**2)
              scold   = scneu
              call lorentz(betax, betay, betaz,
     &                    pboo(1),pboo(2),pboo(3),pboo(0))

c             call lorlrf(j0,j1,j2,j3,pboo,rho,flag )

*      now the vector pboo contains the energy and momenta
*      of the particle in the LRF of the cell k.
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

              p1x = pboo(1)
              p1y = pboo(2)
              p1z = pboo(3)
              pabs= sqrt(p1x**2+p1y**2+p1z**2)
              pst = pabs

              plrf(1) = p1x
              plrf(2) = p1y
              plrf(3) = p1z



              pst       = pabs
*org             vecpot    = potanal(rho0,masse,rho,pabs,id1)
              vecpot    = potanal(rho0,rho,plrf,id1,rhap)
c              write(500,*)vecpot

              ehelp = sqrt(masse**2 + pabs**2)
              vecpot = - masse + sqrt((ehelp+vecpot)**2-pabs**2)

              pboo(0) = sqrt( (masse+vecpot)**2 + pabs**2)


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*             boost back to calc. frame                              *

              if(pboo(1)**2+pboo(2)**2+pboo(3)**2.gt.pboo(0)**2) then
                 write(*,*) "hiba gradu lorentz2, negative mass",
     &                pboo(0),pboo(1), pboo(2), pboo(3)
c                 stop
              end if
              call lorentz(-betax, -betay, -betaz,
     &                    pboo(1),pboo(2),pboo(3),pboo(0))

c              call lorlrf(j0,-j1,-j2,-j3,pboo,dummy,flag )

*                                                                    *
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
              pxneu = pboo(1)
              pyneu = pboo(2)
              pzneu = pboo(3)
              scneu = vecpot

*            check iteration condition

              pxtest = abs(pxneu - pxstat)
              pytest = abs(pyneu - pystat)
              pztest = abs(pzneu - pzstat)


              if(pxtest.le.epsi .and. pytest.le.epsi .and.
     &          pztest.le.epsi) flagit = .true.

              if(itcount.ge.100) then
                write(*,*)'in iteration gradu 1a ', itcount
                write(*,*)'old ', pxstat, pystat, pzstat, masse,scold
                write(*,*)'neu ', pxneu, pyneu, pzneu, masse,scneu
ccc_hw               stop
               goto 777
              end if

            end do
  777 continue

            potgradr(k-1)=pboo(0)


          end do
*--------------------------------------------------------------------*
*    potgradr contains the 6/12 pot.values needed for the evaluation *
*    of gradxr, gradyr, gradzr                                       *

          if(iorder .eq. 1) then

            gradxr    = (potgradr(1) - potgradr(2))/(2.0*dd)
            gradyr    = (potgradr(3) - potgradr(4))/(2.0*dd)
            gradzr    = (potgradr(5) - potgradr(6))/(2.0*dd)

          else if(iorder .eq. 2) then

            gradxr = 1.0/(12.0*dd)*(potgradr(8) - 8.0*potgradr(2) +
     +               8.0*potgradr(1) - potgradr(7))
            gradyr = 1.0/(12.0*dd)*(potgradr(10) - 8.0*potgradr(4) +
     +               8.0*potgradr(3) - potgradr(9))
            gradzr = 1.0/(12.0*dd)*(potgradr(12) - 8.0*potgradr(6) +
     +               8.0*potgradr(5) - potgradr(11))


          end if

*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*      now do the same for the derivative wrsp to p                  *
        else if(isw.eq.2) then
          j0 = deriv(0,1)
          j1 = deriv(1,1)
          j2 = deriv(2,1)
          j3 = deriv(3,1)
          sd = deriv(4,1)


          rho = sqrt(j0**2 - j1**2 - j2**2 - j3**2)

          if(j0.lt.1.0e-06) then
            j0 = 0.0
            j1 = 0.0
            j2 = 0.0
            j3 = 0.0
            betax = 0.0
            betay = 0.0
            betaz = 0.0
          else
            betax = j1/j0
            betay = j2/j0
            betaz = j3/j0


          end if


          coord(1,1) = px
          coord(1,2) = py
          coord(1,3) = pz

          coord(2,1) = px + dp
          coord(2,2) = py
          coord(2,3) = pz

          coord(3,1) = px - dp
          coord(3,2) = py
          coord(3,3) = pz

          coord(4,1) = px
          coord(4,2) = py + dp
          coord(4,3) = pz

          coord(5,1) = px
          coord(5,2) = py - dp
          coord(5,3) = pz


          coord(6,1) = px
          coord(6,2) = py
          coord(6,3) = pz + dp


          coord(7,1) = px
          coord(7,2) = py
          coord(7,3) = pz - dp

          coord(8,1) = px + 2.0*dp
          coord(8,2) = py
          coord(8,3) = pz

          coord(9,1) = px - 2.0*dp
          coord(9,2) = py
          coord(9,3) = pz

          coord(10,1) = px
          coord(10,2) = py + 2.0*dp
          coord(10,3) = pz

          coord(11,1) = px
          coord(11,2) = py - 2.0*dp
          coord(11,3) = pz

          coord(12,1) = px
          coord(12,2) = py
          coord(12,3) = pz + 2.0*dp

          coord(13,1) = px
          coord(13,2) = py
          coord(13,3) = pz - 2.0*dp

          rhap(1)     = rx
          rhap(2)     = ry
          rhap(3)     = rz


          do ib = 1, 13

            pboo(1) = coord(ib,1)
            pboo(2) = coord(ib,2)
            pboo(3) = coord(ib,3)
            pboo(0) = sqrt((e(idn)+upot(idn))**2+pboo(1)**2+pboo(2)**2+
     +                       pboo(3)**2)


            pxstat = pboo(1)
            pystat = pboo(2)
            pzstat = pboo(3)
            scneu = upot(idn)
            masse = e(idn)


            flagit = .false.
            itcount = 0
            do while(.not.flagit)

              itcount = itcount + 1

              pboo(1) = pxstat
              pboo(2) = pystat
              pboo(3) = pzstat
              pboo(0) = sqrt( (masse+scneu)**2 + pboo(1)**2 +
     &                         pboo(2)**2 + pboo(3)**2)
              scold   = scneu

              call lorentz(betax, betay, betaz,
     &                    pboo(1),pboo(2),pboo(3),pboo(0))


c              call lorlrf(j0,j1,j2,j3,pboo,rho,flag )

*      now the vector pboo contains the energies and momenta
*      of the particle in the LRF of the cell ix,iy,iz.
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              p1x = pboo(1)
              p1y = pboo(2)
              p1z = pboo(3)
              pabs= sqrt(p1x**2+p1y**2+p1z**2)

              plrf(1) = p1x
              plrf(2) = p1y
              plrf(3) = p1z


              pst       = pabs
              vecpot    = potanal(rho0,rho,plrf,id1,rhap)
c              write(500,*)vecpot


              ehelp = sqrt(masse**2 + pabs**2)
              vecpot = - masse + sqrt((ehelp+vecpot)**2-pabs**2)


              pboo(0) = sqrt( (masse+vecpot)**2 + pabs**2)


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*             boost back to calc. frame                              *

              if(pboo(1)**2+pboo(2)**2+pboo(3)**2.gt.pboo(0)**2) then
                 write(*,*) "hiba gradu lorentz4, negative mass",
     &                pboo(0),pboo(1), pboo(2), pboo(3)
c                 stop
              end if
              call lorentz(-betax, -betay, -betaz,
     &                    pboo(1),pboo(2),pboo(3),pboo(0))



c              call lorlrf(j0,-j1,-j2,-j3,pboo,dummy,flag )

*                                                                    *
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
              pxneu = pboo(1)
              pyneu = pboo(2)
              pzneu = pboo(3)
              scneu = vecpot

*            check iteration condition

              pxtest = abs(pxneu - pxstat)
              pytest = abs(pyneu - pystat)
              pztest = abs(pzneu - pzstat)


              if(pxtest.le.epsi .and. pytest.le.epsi .and.
     &          pztest.le.epsi) then

                flagit = .true.
              end if

              if(itcount.ge.100) then
                write(*,*)'in iteration gradu 1b ', itcount
                write(*,*)'old ', pxstat, pystat, pzstat, masse,scold
                write(*,*)'neu ', pxneu, pyneu, pzneu, masse,scneu
ccc_hw               stop
                goto 778
              end if

            end do
  778       continue

            scapo(ib) = pboo(0)
            if(ib.eq.1) scapo(ib) = vecpot

          end do

          gradxp = 1.0/(12.0*dp)*(scapo(9) - 8.0*scapo(3) +
     +                8.0*scapo(2) - scapo(8))
          gradyp = 1.0/(12.0*dp)*(scapo(11) - 8.0*scapo(5) +
     +                8.0*scapo(4) - scapo(10))
          gradzp = 1.0/(12.0*dp)*(scapo(13) - 8.0*scapo(7) +
     +                8.0*scapo(6) - scapo(12))

c          write(*,*) 'gradu 2',scapo(1)
          upot(idn) = scapo(1)

c          if(idelpot.eq.1 .and. id1.eq.2) then
c            write(*,*)'idelpot = 1 is not implemented'
c            stop
c          endif
*
c          if(idelpot.eq.2 .and. id1.eq.2) then
c            write(*,*)'idelpot = 2 is not implemented'
c            stop
c          endif

c        if(ipronew.eq.0) then
c          dx=2.0
c          jx=ix+1
c          kx=ix-1
c          dy=2.0
c          jy=iy+1
c          ky=iy-1
c          dz=2.0
c          jz=iz+1
c          kz=iz-1
c
c          yxplus   = yup(jx,iy,iz)
c          yxmins   = yup(kx,iy,iz)
c          yyplus   = yup(ix,jy,iz)
c          yymins   = yup(ix,ky,iz)
c          yzplus   = yup(ix,iy,jz)
c          yzmins   = yup(ix,iy,kz)
c          gradxr=gradxr+(yxplus-yxmins)/dx
c          gradyr=gradyr+(yyplus-yymins)/dy
c          gradzr=gradzr+(yzplus-yzmins)/dz
c        end if

      end if


*---------------------------------------------------------------------*
*        electro-magnetic part                                        *

      if(isw.eq.1) then
        if(id2.ne.0) then

          if(ipou.eq.1) then

*         call routine emfoca to calculate the electro-magnetic force
*         emfox, emfoy, emfoz  (link to coulomb-routine)
*     use emfoca1.f !!!!!!!!!!!!!!
            call emfoca(rx, ry, rz, id2,
     &                  emfox, emfoy, emfoz,ncont,cpot)

            gradxr = gradxr - emfox
            gradyr = gradyr - emfoy
            gradzr = gradzr - emfoz
          end if
        end if
      end if
*                                                                    *
*--------------------------------------------------------------------*

      return
      end

***************************************************************************
      subroutine lorlrf(j0,j1,j2,j3,pboo,j0p,flag)
*
*    this subroutine does the lorentz-boost into the LRF of each grid cell
*
*---------------------------------------------------------------------*
      implicit none
      include'common'

      real*8    j0,j1,j2,j3
      real*8    j0p,j1p,j2p,j3p
      real*8    p0p, p1p, p2p, p3p
      real*8    p0,  p1 , p2 , p3
      real*8    pboo(0:3)
      real*8    betax, betay, betaz, beta, gamma
      real*8    trans, rinv, rinvp
      logical flag

      if(j0 .lt. 1e-08) then
        j0 = 0.0
        j1 = 0.0
        j2 = 0.0
        j3 = 0.0
        j0p= 0.0
        flag = .true.
        return
      end if
      rinv = j0**2-j1**2-j2**2-j3**2

      betax = j1/j0
      betay = j2/j0
      betaz = j3/j0
      beta  = sqrt(betax**2 +betay**2 + betaz**2)
      gamma = 1.0/(sqrt(1.0 - beta**2))

      trans = betax*j1 + betay*j2 + betaz*j3
      trans = trans * gamma
      trans = trans /(gamma+1.0)
      trans = trans - j0

      j1p = j1 + gamma*betax*trans
      j2p = j2 + gamma*betay*trans
      j3p = j3 + gamma*betaz*trans
      j0p = j0 - betax*j1 - betay*j2 - betaz*j3
      j0p = gamma*j0p

      rinvp=j0p**2-j1p**2-j2p**2-j3p**2

      if(abs(rinv-rinvp).gt. 1.0e-07) then
       write(115,*)rinv, rinvp
      end if

      rinvp = j1p**2 + j2p**2 +j3p**2
      if(rinvp.gt.1.0e-07) then
        write(115,*) rinvp
        write(*,*)'problem in lrf-trafo '
        stop
      end if

        p0 = pboo(0)
        p1 = pboo(1)
        p2 = pboo(2)
        p3 = pboo(3)

        rinv = sqrt(p0**2 -p1**2 -p2**2-p3**2)

        trans = betax*p1 + betay*p2 + betaz*p3
        trans = trans * gamma
        trans = trans /(gamma+1.0)
        trans = trans - p0

        p1p = p1 + gamma*betax*trans
        p2p = p2 + gamma*betay*trans
        p3p = p3 + gamma*betaz*trans
        p0p = p0 - betax*p1 - betay*p2 - betaz*p3
        p0p = gamma*p0p
        rinvp=sqrt(p0p**2-p1p**2-p2p**2-p3p**2 )

      if(abs(rinv-rinvp).gt. 1.0e-04) then
        write(115,*)rinv, rinvp
        write(*,*)'problem in lrf-trafo '
        stop
      end if

      pboo(0) = p0p
      pboo(1) = p1p
      pboo(2) = p2p
      pboo(3) = p3p
      return
      end
