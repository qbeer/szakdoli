************************************************************************
************** spline-routine for pion-potential************************
************************************************************************
************************************************************************
      subroutine linint(rx, ry, rz, deriv, dd,ort)
*---------------------------------------------------------------------*
*      diese routine stellt das interface zwischen dem buu-programm    *
*      und der spline routine dar.                                    *
*                                                                     *
*     rx, ry, rz           : ort des teilchens                (input) *
*     gradx, grady, gradz  : kraft, die aufs teilchen wirkt  (output) *
*                                                                     *
*---------------------------------------------------------------------*

      implicit none
      include"common"
*      include"commsp"


      real*8    rx, ry, rz
      real*8    deriv(0:4,1:13)
      integer ix, iy, iz
      real*8    xx, yy, zz
      integer ipkt
      real*8    coord(1:13,1:3),ort(1:3,1:13)

      integer  ivec
      real*8    dd

        dd = 1.0


*--------------------------------------------------------------------*
*     belegen der 7 ortskoordinaten                                  *

        coord(1,1) = rx
        coord(1,2) = ry
        coord(1,3) = rz

        coord(2,1) = rx + dd
        coord(2,2) = ry
        coord(2,3) = rz

        coord(3,1) = rx - dd
        coord(3,2) = ry
        coord(3,3) = rz

        coord(4,1) = rx
        coord(4,2) = ry + dd
        coord(4,3) = rz

        coord(5,1) = rx
        coord(5,2) = ry - dd
        coord(5,3) = rz


        coord(6,1) = rx
        coord(6,2) = ry
        coord(6,3) = rz + dd


        coord(7,1) = rx
        coord(7,2) = ry
        coord(7,3) = rz - dd

c        coord(8,1) = rx + 2.0*dd
c        coord(8,2) = ry
c        coord(8,3) = rz

c        coord(9,1) = rx - 2.0*dd
c        coord(9,2) = ry
c        coord(9,3) = rz

c        coord(10,1) = rx
c        coord(10,2) = ry + 2.0*dd
c        coord(10,3) = rz

c        coord(11,1) = rx
c        coord(11,2) = ry - 2.0*dd
c        coord(11,3) = rz

c        coord(12,1) = rx
c        coord(12,2) = ry
c        coord(12,3) = rz + 2.0*dd

c        coord(13,1) = rx
c        coord(13,2) = ry
c        coord(13,3) = rz - 2.0*dd

*                                                                    *
*--------------------------------------------------------------------*


       do ipkt = 1,7

          xx = coord(ipkt,1)
          yy = coord(ipkt,2)
          zz = coord(ipkt,3)

          ort(1,ipkt) = xx
          ort(2,ipkt) = yy
          ort(3,ipkt) = zz

          ix = nint(xx)
          iy = nint(yy)
          iz = nint(zz)

          if( abs(ix).le.maxx .and. abs(iy).le.maxx .and.
     &        abs(iz).le.maxz)then
            do ivec = 0, 4
              deriv(ivec,ipkt) = rhob_4(ivec,ix,iy,iz)
c              write(*,*) 'linint',rhob_4(0,ix,iy,iz),rhob_4(1,ix,iy,iz),
c     &   rhob_4(2,ix,iy,iz),rhob_4(3,ix,iy,iz),rhob_4(4,ix,iy,iz)
            end do
          else
            do ivec = 0, 4
              deriv(ivec,ipkt) = 0.0
            end do
          end if

       end do

       return

       end


**********************************************************************

      subroutine linint1(rx, ry, rz, deriv)
*---------------------------------------------------------------------*
*      calculates at r the four current stored in deriv               *
*                                                                     *
*     rx, ry, rz           : position of the particle         (input) *
*                                                                     *
*---------------------------------------------------------------------*

      implicit none
      include"common"


      real*8    rx, ry, rz
      real*8    deriv(0:4)

      integer ix, iy,iz
      integer ivec
*--------------------------------------------------------------------*
      if(abs(rx).lt.maxx .and. abs(ry).lt.maxx .and.
     &            abs(rz).lt.maxz) then
        ix = nint(rx)
        iy = nint(ry)
        iz = nint(rz)
        do ivec = 0, 4
          deriv(ivec) = rhob_4(ivec,ix,iy,iz)
        end do
      else
        do ivec = 0, 4
          deriv(ivec) = 0.0
        end do
      end if

      return
      end


**********************************************************************

      subroutine linint_prev(rx, ry, rz, deriv)
*---------------------------------------------------------------------*
*      calculates at r the previous four current stored in deriv      *
*                                                                     *
*     rx, ry, rz           : position of the particle         (input) *
*                                                                     *
*---------------------------------------------------------------------*

      implicit none
      include"common"
      real*8    rx, ry, rz
      real*8    deriv(0:4)

      integer ix, iy,iz
      integer  ivec, ivec2
*--------------------------------------------------------------------*
      if(abs(rx).lt.maxx .and. abs(ry).lt.maxx .and.
     &            abs(rz).lt.maxz) then
        ix = nint(rx)
        iy = nint(ry)
        iz = nint(rz)
        do ivec = 0, 4
          if (ivec.eq.0) then
            ivec2 = 10
          else
            ivec2 = ivec + 6
          end if

          deriv(ivec) = rhob_4(ivec2,ix,iy,iz)
        end do
      else
        do ivec = 0, 4
          deriv(ivec) = 0.0
        end do

      end if

      return
      end

