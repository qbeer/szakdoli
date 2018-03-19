      subroutine emfoca(rx, ry, rz, id2,
     &                  emfox, emfoy, emfoz,ncont,cpot)
*---------------------------------------------------------------------*
*             this routine calculates the electomagnetic-forces       *
*             acting on a prticle (baryon or meson)                   *
*          viables:                                                   *
*               ncont    - only coulomb (=0) or full calculation (=3) *
*                          (cominput)                                 *
*               rx,ry,rz - space coordinates of particle (input,real) *
*               px,py,pz - momentum of particle          (input,real) *
*               etot     - energy or particle            (input,real) *
*               id2      - charge of particle         (input,integer) *
*               emfo..   - electromagn. force           (output,real) *
*                                                                     *
*---------------------------------------------------------------------*
      implicit none
      include"coucom"

      real*8    rx, ry, rz, px, py, pz, etot, emfox, emfoy, emfoz
      real*8    charge
      integer id2, ixc, iyc, izc, ncont
      real*8    relx, rely, relz, relr, emfoab

      integer ixo, iyo, izo, ixu, iyu, izu
      real*8 tlin, ulin, vlin
      integer izr, ixr, iyr, ipkt
      real*8 xx, yy, zz
      real*8 funlin(1:8), deriv(1:7)
      real*8 gradx, grady, gradz, cpot
      real*8 coord(1:7,1:3)

*---------------------------------------------------------------------*
*        set foreces to zero

      emfox = 0.0
      emfoy = 0.0
      emfoz = 0.0

      charge = float(id2)
*---------------------------------------------------------------------*
*        determine the indices corresponding to the particle posiion  *


      ixc = nint(rx/dgrid)
      iyc = nint(ry/dgrid)
      izc = nint(rz/dgrid)

*----------------------------------------------------------------------*
*        check wethe indices stay within the Coulomb-grid or not       *

      if(abs(ixc).le.cmaxx-1 .and. abs(iyc).le.cmaxx-1 .and.
     +   abs(izc).le.cmaxz-1) then


*----------------------------------------------------------------------*
*        electric force (in all cases)                                 *

*     belegen der 7 ortskoordinaten                                  *
        coord(1,1) = rx
        coord(1,2) = ry
        coord(1,3) = rz

        coord(2,1) = rx + dgrid
        coord(2,2) = ry
        coord(2,3) = rz

        coord(3,1) = rx - dgrid
        coord(3,2) = ry
        coord(3,3) = rz

        coord(4,1) = rx
        coord(4,2) = ry + dgrid
        coord(4,3) = rz

        coord(5,1) = rx
        coord(5,2) = ry - dgrid
        coord(5,3) = rz


        coord(6,1) = rx
        coord(6,2) = ry
        coord(6,3) = rz + dgrid


        coord(7,1) = rx
        coord(7,2) = ry
        coord(7,3) = rz - dgrid

        do ipkt = 1,7

          xx = coord(ipkt,1)
          yy = coord(ipkt,2)
          zz = coord(ipkt,3)

          ixr = nint(xx/dgrid)
          iyr = nint(yy/dgrid)
          izr = nint(zz/dgrid)



          ixu = int(xx/dgrid)
          iyu = int(yy/dgrid)
          izu = int(zz/dgrid)


          ixo = ixu+1
          iyo = iyu+1
          izo = izu+1


          funlin(1) = emrpot(0,ixu,iyu,izu)
          funlin(2) = emrpot(0,ixo,iyu,izu)
          funlin(3) = emrpot(0,ixo,iyo,izu)
          funlin(4) = emrpot(0,ixu,iyo,izu)
          funlin(5) = emrpot(0,ixu,iyu,izo)
          funlin(6) = emrpot(0,ixo,iyu,izo)
          funlin(7) = emrpot(0,ixo,iyo,izo)
          funlin(8) = emrpot(0,ixu,iyo,izo)


          tlin = (xx-float(ixu)*dgrid)/dgrid
          ulin = (yy-float(iyu)*dgrid)/dgrid
          vlin = (zz-float(izu)*dgrid)/dgrid


          deriv(ipkt) =  (1.-tlin)*(1.-ulin)*(1-vlin)*funlin(1) +
     +                       tlin *(1.-ulin)*(1-vlin)*funlin(2) +
     +                       tlin *    ulin *(1-vlin)*funlin(3) +
     +                   (1.-tlin)*    ulin *(1-vlin)*funlin(4) +
     +                   (1.-tlin)*(1.-ulin)*   vlin *funlin(5) +
     +                       tlin *(1.-ulin)*   vlin *funlin(6) +
     +                       tlin *    ulin *   vlin *funlin(7) +
     +                   (1.-tlin)*    ulin *   vlin *funlin(8)






c              write(20,*)'deriv: '
c              write(20,*)'ipkt =', ipkt
c              write(20,*)'deriv = ', deriv(ipkt)
c              write(20,*)'tlin = ', tlin
c              write(20,*)'ulin = ', ulin
c              write(20,*)'vlin = ', vlin
c              write(20,*)'r s ',xx,yy,zz
c              write(20,*)'ixu,iyu,izu',ixu,iyu,izu
c              write(20,*)'ixo,iyo,izo',ixo,iyo,izo

c              do ij = 1, 8
c                write(20,*)ij, funlin(ij)
c              end do

        end do


        gradx = charge*(deriv(2)-deriv(3))/2./dgrid
        grady = charge*(deriv(4)-deriv(5))/2./dgrid
        gradz = charge*(deriv(6)-deriv(7))/2./dgrid


        emfox = -gradx
        emfoy = -grady
        emfoz = -gradz
        cpot = charge*deriv(1)






*----------------------------------------------------------------------*


      else if(abs(ixc).ge.cmaxx .or. abs(iyc).ge.cmaxx .or.
     +        abs(izc).ge.cmaxz) then

        relx  = rx - rchar(1)
        rely  = ry - rchar(2)
        relz  = rz - rchar(3)
        relr  = relx**2 + rely**2 + relz**2

        if(relr .lt. 1.0) then
          write(*,*)'problems in emfoca; relr**2 = ', relr
          write(*,*)'particle coordinates :', rx, ry, rz
          write(*,*)'center of charge :',rchar(1),rchar(2),rchar(3)
          		write(50,*) "stop HS 20"
          stop
        end if



        relr   = sqrt(relx**2 + rely**2 + relz**2)
        cpot   = charge * elmcon * chatot/relr
        relr   = relr**3.0
        emfoab = charge * elmcon * chatot/relr
c       if (ncont .eq. 777)
c    1  write(*,*) '  in emfoca ', charge, elmcon, chatot, relz

        emfox  = emfoab*relx
        emfoy  = emfoab*rely
        emfoz  = emfoab*relz


      end if

      return
      end


