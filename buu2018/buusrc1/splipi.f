***********************************************************************
************* spline-routine for pion-potential************************
***********************************************************************
***********************************************************************
***********************************************************************
      subroutine splinint(rx, ry, rz, deriv, dd)
*---------------------------------------------------------------------*
*      diese routine stellt das interface zwischen dem buu-programm    *
*      und der spline routine dar.                                    *
*                                                                     *
*     rx, ry, rz           : ort des teilchens                (input) *
*                                                                     *
*---------------------------------------------------------------------*

      implicit none
      include"common"
      include"commsp"


      real*8    ddg, ddk
      parameter(ddk = dg)
      parameter(ddg = 1.0)
      real*8    rx, ry, rz
      real*8    funlin(1:8)
      real*8    deriv(0:idim,1:13)

      logical klgro  !true : abl auf kleinem gitter
      real*8    testx,   testz
      integer ixo, ixu, iyo, iyu, izo,izu
      real*8    tlin, ulin, vlin
      real*8    xx, yy, zz
      integer   itestx, itesty, itestz , ipkt
      real*8    coord(1:13,1:3)

      integer ixr, iyr, izr
      integer  ivec,ij
      real*8    dd
      integer indx, indy, indz

*     feststellen, ob das kleine oder das grosse Gitter verwendet
*     werden soll.

      testx = proz*float(maxx)-5.0*dg
      testz = 1.0*float(maxz)-5.0*dg

      indx  = nint(abs(rx))+ 2
      indy  = nint(abs(ry))+ 2
      indz  = nint(abs(rz))+ 2



      klgro = .false.
      if(abs(rx) .le. testx .and. abs(ry) .le. testx
     +  .and. abs(rz) .le. testz) klgro = .true.

      if(klgro) then
*     kleines gitter

        dd = ddk

*--------------------------------------------------------------------*
*     belegen der 7 ortskoordinaten                                  *
        coord(1,1) = rx
        coord(1,2) = ry
        coord(1,3) = rz

        coord(2,1) = rx + ddk
        coord(2,2) = ry
        coord(2,3) = rz

        coord(3,1) = rx - ddk
        coord(3,2) = ry
        coord(3,3) = rz

        coord(4,1) = rx
        coord(4,2) = ry + ddk
        coord(4,3) = rz

        coord(5,1) = rx
        coord(5,2) = ry - ddk
        coord(5,3) = rz


        coord(6,1) = rx
        coord(6,2) = ry
        coord(6,3) = rz + ddk


        coord(7,1) = rx
        coord(7,2) = ry
        coord(7,3) = rz - ddk

        coord(8,1) = rx + 2.0*ddk
        coord(8,2) = ry
        coord(8,3) = rz

        coord(9,1) = rx - 2.0*ddk
        coord(9,2) = ry
        coord(9,3) = rz

        coord(10,1) = rx
        coord(10,2) = ry + 2.0*ddk
        coord(10,3) = rz

        coord(11,1) = rx
        coord(11,2) = ry - 2.0*ddk
        coord(11,3) = rz

        coord(12,1) = rx
        coord(12,2) = ry
        coord(12,3) = rz + 2.0*ddk

        coord(13,1) = rx
        coord(13,2) = ry
        coord(13,3) = rz - 2.0*ddk
*                                                                  *
*------------------------------------------------------------------*

        do ivec = 0, idim
        do ipkt = 1, 13

          xx = coord(ipkt,1)
          yy = coord(ipkt,2)
          zz = coord(ipkt,3)

          ixr = nint(xx)
          iyr = nint(yy)
          izr = nint(zz)



          itestx = int((xx + proz*float(maxx))/dg) + 1
          itesty = int((yy + proz*float(maxx))/dg) + 1
          itestz = int((zz + 1.0*float(maxz))/dg) + 1


          ixu = itestx
          iyu = itesty
          izu = itestz


          ixo = ixu+1
          iyo = iyu+1
          izo = izu+1


          funlin(1) = resg(ivec,ixu,iyu,izu)
          funlin(2) = resg(ivec,ixo,iyu,izu)
          funlin(3) = resg(ivec,ixo,iyo,izu)
          funlin(4) = resg(ivec,ixu,iyo,izu)
          funlin(5) = resg(ivec,ixu,iyu,izo)
          funlin(6) = resg(ivec,ixo,iyu,izo)
          funlin(7) = resg(ivec,ixo,iyo,izo)
          funlin(8) = resg(ivec,ixu,iyo,izo)


          tlin = (xx-xval(ixu))/dg
          ulin = (yy-yval(iyu))/dg
          vlin = (zz-zval(izu))/dg


          deriv(ivec,ipkt) =  (1.-tlin)*(1.-ulin)*(1-vlin)*funlin(1) +
     +                       tlin *(1.-ulin)*(1-vlin)*funlin(2) +
     +                       tlin *    ulin *(1-vlin)*funlin(3) +
     +                   (1.-tlin)*    ulin *(1-vlin)*funlin(4) +
     +                   (1.-tlin)*(1.-ulin)*   vlin *funlin(5) +
     +                       tlin *(1.-ulin)*   vlin *funlin(6) +
     +                       tlin *    ulin *   vlin *funlin(7) +
     +                   (1.-tlin)*    ulin *   vlin *funlin(8)





          if(ivec.eq.0) then
            if(deriv(ivec,ipkt) .lt. 0.0) then
              write(*,*)'deriv: '
              write(*,*)'ipkt =', ipkt
              write(*,*)'deriv = ', deriv(ivec,ipkt)
              write(*,*)'tlin = ', tlin
              write(*,*)'ulin = ', ulin
              write(*,*)'vlin = ', vlin
              write(*,*)'r s ',xx,yy,zz,xval(ixu),yval(iyu),zval(izu)
              write(*,*)'ixu,iyu,izu',ixu,iyu,izu
              write(*,*)'ixo,iyo,izo',ixo,iyo,izo

              do ij = 1, 8
                write(*,*)ij, funlin(ij)
              end do

              stop
            end if
          end if


        end do


*        gradx = (deriv(ivec,2)-deriv(ivec,3))/2./ddk
*        grady = (deriv(ivec,4)-deriv(ivec,5))/2./ddk
*        gradz = (deriv(ivec,6)-deriv(ivec,7))/2./ddk

        end do

      else if( .not.klgro  .and. indx.lt.maxx .and. indy.lt.maxx
     +        .and. indz .lt. maxz) then
*--------------------------------------------------------------------*
*     belegen der 7 ortskoordinaten                                  *

        dd = ddg
        coord(1,1)  = rx
        coord(1,2)  = ry
        coord(1,3)  = rz

        coord(2,1)  = rx + ddg
        coord(2,2)  = ry
        coord(2,3)  = rz

        coord(3,1)  = rx - ddg
        coord(3,2)  = ry
        coord(3,3)  = rz

        coord(4,1)  = rx
        coord(4,2)  = ry + ddg
        coord(4,3)  = rz

        coord(5,1)  = rx
        coord(5,2)  = ry - ddg
        coord(5,3)  = rz


        coord(6,1)  = rx
        coord(6,2)  = ry
        coord(6,3)  = rz + ddg


        coord(7,1)  = rx
        coord(7,2)  = ry
        coord(7,3)  = rz - ddg

        coord(8,1)  = rx + 2.0*ddg
        coord(8,2)  = ry
        coord(8,3)  = rz

        coord(9,1)  = rx - 2.0*ddg
        coord(9,2)  = ry
        coord(9,3)  = rz

        coord(10,1) = rx
        coord(10,2) = ry + 2.0*ddg
        coord(10,3) = rz

        coord(11,1) = rx
        coord(11,2) = ry - 2.0*ddg
        coord(11,3) = rz

        coord(12,1) = rx
        coord(12,2) = ry
        coord(12,3) = rz + 2.0*ddg

        coord(13,1) = rx
        coord(13,2) = ry
        coord(13,3) = rz - 2.0*ddg
*                                                                    *
*--------------------------------------------------------------------*

       do ivec = 0, idim
       do ipkt = 1,13

          xx = coord(ipkt,1)
          yy = coord(ipkt,2)
          zz = coord(ipkt,3)


          ixu = int(xx)
          if(xx.lt.0.0) ixu = ixu - 1
          iyu = int(yy)
          if(yy.lt.0.0) iyu = iyu - 1
          izu = int(zz)
          if(zz.lt.0.0) izu = izu - 1

          ixo = ixu + 1
          iyo = iyu + 1
          izo = izu + 1


          funlin(1) = rhob_4(ivec,ixu,iyu,izu)
          funlin(2) = rhob_4(ivec,ixo,iyu,izu)
          funlin(3) = rhob_4(ivec,ixo,iyo,izu)
          funlin(4) = rhob_4(ivec,ixu,iyo,izu)
          funlin(5) = rhob_4(ivec,ixu,iyu,izo)
          funlin(6) = rhob_4(ivec,ixo,iyu,izo)
          funlin(7) = rhob_4(ivec,ixo,iyo,izo)
          funlin(8) = rhob_4(ivec,ixu,iyo,izo)

          tlin = (xx-float(ixu))/dg
          ulin = (yy-float(iyu))/dg
          vlin = (zz-float(izu))/dg


          deriv(ivec,ipkt) =  (1.-tlin)*(1.-ulin)*(1-vlin)*funlin(1) +
     +                       tlin *(1.-ulin)*(1-vlin)*funlin(2) +
     +                       tlin *    ulin *(1-vlin)*funlin(3) +
     +                   (1.-tlin)*    ulin *(1-vlin)*funlin(4) +
     +                   (1.-tlin)*(1.-ulin)*   vlin *funlin(5) +
     +                       tlin *(1.-ulin)*   vlin *funlin(6) +
     +                       tlin *    ulin *   vlin *funlin(7) +
     +                   (1.-tlin)*    ulin *   vlin *funlin(8)



        end do


*        gradx = (deriv(ivec,2)-deriv(ivec,3))/2./ddg
*        grady = (deriv(ivec,4)-deriv(ivec,5))/2./ddg
*        gradz = (deriv(ivec,6)-deriv(ivec,7))/2./ddg

        end do
      else
***    particle out of grid

       do ivec = 0, idim
         do ipkt = 1,13
           deriv(ivec,ipkt) =  0.0
         end do
       end do

      end if


      return
      end

**********************************************************************
      subroutine splinint1(rx, ry, rz, deriv)
*---------------------------------------------------------------------*
*      diese routine stellt das interface zwischen dem buu-programm   *
*      und der spline routine dar.                                    *
*                                                                     *
*     rx, ry, rz           : ort des teilchens                (input) *
*     deriv(0:idim)        : dichte und stroeme am ort r     (output) *
*                                                                     *
*---------------------------------------------------------------------*

      implicit none
      include"common"
      include"commsp"


      real*8    ddg, ddk
      parameter(ddk = dg)
      parameter(ddg = 1.0)
      real*8    rx, ry, rz

      real*8    funlin(1:8)
      real*8    deriv(0:idim)

      logical klgro  !true : abl auf kleinem gitter
      real*8    testx,  testz
      integer ixo, ixu, iyo, iyu, izo,izu
      real*8    tlin, ulin, vlin
      real*8    xx, yy, zz
      integer  itestx, itesty, itestz
      integer ixr, iyr, izr
      integer  ivec,ij
      real*8    dd
      integer indx, indy, indz

*     feststellen, ob das kleine oder das grosse Gitter verwendet
*     werden soll.

      testx = proz*float(maxx)-5.0*dg
      testz = 1.0*float(maxz)-5.0*dg

      indx  = nint(abs(rx))+ 2
      indy  = nint(abs(ry))+ 2
      indz  = nint(abs(rz))+ 2



      klgro = .false.
      if(abs(rx) .le. testx .and. abs(ry) .le. testx
     +  .and. abs(rz) .le. testz) klgro = .true.

      if(klgro) then
*     kleines gitter

        dd = ddk

        do ivec = 0, idim

          xx  = rx
          yy  = ry
          zz  = rz

          ixr = nint(xx)
          iyr = nint(yy)
          izr = nint(zz)



          itestx = int((xx + proz*float(maxx))/dg) + 1
          itesty = int((yy + proz*float(maxx))/dg) + 1
          itestz = int((zz + 1.0*float(maxz))/dg) + 1


          ixu = itestx
          iyu = itesty
          izu = itestz


          ixo = ixu+1
          iyo = iyu+1
          izo = izu+1


          funlin(1) = resg(ivec,ixu,iyu,izu)
          funlin(2) = resg(ivec,ixo,iyu,izu)
          funlin(3) = resg(ivec,ixo,iyo,izu)
          funlin(4) = resg(ivec,ixu,iyo,izu)
          funlin(5) = resg(ivec,ixu,iyu,izo)
          funlin(6) = resg(ivec,ixo,iyu,izo)
          funlin(7) = resg(ivec,ixo,iyo,izo)
          funlin(8) = resg(ivec,ixu,iyo,izo)


          tlin = (xx-xval(ixu))/dg
          ulin = (yy-yval(iyu))/dg
          vlin = (zz-zval(izu))/dg


          deriv(ivec) =  (1.-tlin)*(1.-ulin)*(1-vlin)*funlin(1) +
     +                       tlin *(1.-ulin)*(1-vlin)*funlin(2) +
     +                       tlin *    ulin *(1-vlin)*funlin(3) +
     +                   (1.-tlin)*    ulin *(1-vlin)*funlin(4) +
     +                   (1.-tlin)*(1.-ulin)*   vlin *funlin(5) +
     +                       tlin *(1.-ulin)*   vlin *funlin(6) +
     +                       tlin *    ulin *   vlin *funlin(7) +
     +                   (1.-tlin)*    ulin *   vlin *funlin(8)





          if(ivec.eq.0) then
            if(deriv(ivec) .le. 0.0) then
              write(*,*)'deriv: '
c              write(*,*)'ipkt =', ipkt
              write(*,*)'deriv = ', deriv(ivec)
              write(*,*)'tlin = ', tlin
              write(*,*)'ulin = ', ulin
              write(*,*)'vlin = ', vlin
              write(*,*)'r s ',xx,yy,zz,xval(ixu),yval(iyu),zval(izu)
              write(*,*)'ixu,iyu,izu',ixu,iyu,izu
              write(*,*)'ixo,iyo,izo',ixo,iyo,izo

              do ij = 1, 8
                write(*,*)ij, funlin(ij)
              end do

              stop
            end if
          end if



        end do

      else if( .not.klgro  .and. indx.lt.maxx .and. indy.lt.maxx
     +        .and. indz .lt. maxz) then

        dd = ddg

       do ivec = 0, idim

          xx = rx
          yy = ry
          zz = rz

          ixu = int(xx)
          if(rx .lt. 0.0) ixu = ixu - 1
          iyu = int(yy)
          if(ry .lt. 0.0) iyu = iyu - 1
          izu = int(zz)
          if(rz .lt. 0.0) izu = izu - 1

          ixo = ixu + 1
          iyo = iyu + 1
          izo = izu + 1


          funlin(1) = rhob_4(ivec,ixu,iyu,izu)
          funlin(2) = rhob_4(ivec,ixo,iyu,izu)
          funlin(3) = rhob_4(ivec,ixo,iyo,izu)
          funlin(4) = rhob_4(ivec,ixu,iyo,izu)
          funlin(5) = rhob_4(ivec,ixu,iyu,izo)
          funlin(6) = rhob_4(ivec,ixo,iyu,izo)
          funlin(7) = rhob_4(ivec,ixo,iyo,izo)
          funlin(8) = rhob_4(ivec,ixu,iyo,izo)

          tlin = (xx-float(ixu))/dg
          ulin = (yy-float(iyu))/dg
          vlin = (zz-float(izu))/dg


          deriv(ivec) =  (1.-tlin)*(1.-ulin)*(1-vlin)*funlin(1) +
     +                       tlin *(1.-ulin)*(1-vlin)*funlin(2) +
     +                       tlin *    ulin *(1-vlin)*funlin(3) +
     +                   (1.-tlin)*    ulin *(1-vlin)*funlin(4) +
     +                   (1.-tlin)*(1.-ulin)*   vlin *funlin(5) +
     +                       tlin *(1.-ulin)*   vlin *funlin(6) +
     +                       tlin *    ulin *   vlin *funlin(7) +
     +                   (1.-tlin)*    ulin *   vlin *funlin(8)



        end do
      else
***    particle out of grid

       do ivec = 0, idim
           deriv(ivec) =  0.0
       end do

      end if


      return
      end


***********************************************************************
***********************************************************************
      subroutine spline3di

      implicit none
      include"common"
      include"commsp"

      integer i, j,k ,  ivec

      real*8  ypkt, zpkt
      integer ig, jg, kg


      integer roindx, roindy, roindz
*-----------------------------------------------------------------*
*         baryonendichte einlesen                                 *


      do ivec = 0, idim
      do i = 1, l
        roindx = -maxx + (i-1)
        do j = 1,m
          roindy = -maxx + (j-1)
          do k = 1,n
            roindz  = -maxz + (k-1)
            ya(i,j,k) = rhob_4(ivec,roindx, roindy,roindz)
          end do
        end do
      end do



      call splie3

      do kg = 1, ng

        zpkt = zval(kg)
        call splin3(zpkt)


*       ableitung nach y
        call splie2

       do jg = 1,mg

          ypkt = yval(jg)


      call splin2(ypkt)
         do ig = 1,lg

           resg(ivec,ig,jg,kg) = ressp(ig)

         end do
        end do
      end do

      end do
      return
********************************************************************

      entry spline3di1
*-------------------------------------------------------------------*
*     belegen der koordinaten   (grosses gitter)                    *

      do i = 1, l
        x1a(i) = - float(maxx) + float(i-1)
      end do

      do i = 1, m
        x2a(i) = -float(maxx) + float(i-1)
      end do

      do i = 1, n
        x3a(i) = -float(maxz) + float(i-1)
      end do
*                                                                  *
*------------------------------------------------------------------*
*     (kleines gitter)                                             *

      do kg = 1, lg
        xval(kg) = -proz*float(maxx) + float(kg-1)*dg
      end do

      do kg = 1,mg
        yval(kg) = -proz*float(maxx) + float(kg-1)*dg
      end do

      do kg = 1,ng
        zval(kg) = -1.0*float(maxz) + float(kg-1)*dg
      end do
*                                                                 *
*-----------------------------------------------------------------*
      return
      end
********************************************************************
************ subroutines *******************************************
*********************************************************************
      subroutine splie3

*-------------------------------------------------------------------*
*     diese routine berechnet die 2. ableitungen der zu fittenden   *
*     function nach z                                               *
*                                                                   *
*     l,m,n      : dim des gitters auf dem die fundion ya           *
*                  definiert ist                          (input)   *
*     ya(l,m,n)  : zu fittende funktion                   (input)   *
*     x3a(n)     : stuetzpunkte von ya in z-richtung      (input)   *
*     y2az(l,m,n): enthalt die 2. ableitung               (output)  *
*-------------------------------------------------------------------*

      implicit none

      include"common"
      include"commsp"

      integer i, j, k

      integer nmax
      parameter(nmax =n)
c      parameter(nmax = max(l,m,n))
c      parameter(nmax = 17)
      real*8 y2tmp(nmax),ytmp(nmax)

      do i = 1, l
        do j = 1, m
          do k = 1, n
            ytmp(k) = ya(i,j,k)
          end do

          call spline(x3a,ytmp,n,1.d30,1.d30,y2tmp)

          do k = 1, n
            y2az(i,j,k)=y2tmp(k)
          end do

        end do
      end do

      return
      end

*******************************************************************
******************************************************************

      subroutine splin3(zpkt)

*-------------------------------------------------------------------*
*    diese routine berechnet die splinekoeffizienten in der         *
*    x-y-ebene an dem z, das von aussen eingegeben wird.            *
*                                                                   *
*     l,m,n      : dim des gitters auf dem die fundion ya           *
*                  definiert ist                          (input)   *
*     x3a(n)     : stuetzpkte in z-richtung               (input)   *
*     ya(l,m,n)  : zu fittende funktion                   (input)   *
*     y2az(l,m,n): enthalt die 2. ableitung               (input)   *
*     zpkt       : z-koordinate des punktes               (input)   *
*     yaz(l,m)   : neue funktion, def. in der xy-ebene    (output)  *
*-------------------------------------------------------------------*

      implicit none

      include"common"
      include"commsp"

      integer i, j, k

      integer nmax
      parameter(nmax =n)
c      parameter(nmax = max(l,m,n))
c      parameter(nmax = 17)
      real*8 y2tmp(nmax),ytmp(nmax)
      real*8    zpkt

      do i = 1, l
        do j = 1, m
          do k = 1, n
            ytmp(k)  = ya(i,j,k)
            y2tmp(k) = y2az(i,j,k)
          end do

          call splint(x3a, ytmp, y2tmp, n, zpkt, yaz(i,j))

        end do
      end do

      return
      end

*********************************************************************
*********************************************************************

      subroutine splie2

*-------------------------------------------------------------------*
*     berechnung der 2. ableitung bezueglich y fuer die output-     *
*     funktion von splin2 (yaz)                                     *
*                                                                   *
*     x2a(m)    : stuetzstellen in y-richtung           (input)     *
*     l,m       : dimension in x-, y-richtung           (input)     *
*     yaz(l,m)  : output vin splin3                     (input)     *
*     y2ay(l,m) : 2. ableitung nach y von yaz           (output)    *
*-------------------------------------------------------------------*

      implicit none

      include"common"
      include"commsp"

      integer i, j

      integer nmax
      parameter(nmax =n)
c      parameter(nmax = max(l,m,n))
c      parameter(nmax = 17)
      real*8 y2tmp(nmax),ytmp(nmax)

      do i = 1, l
        do j = 1, m
          ytmp(j) = yaz(i,j)
        end do

        call spline(x2a, ytmp, m, 1.d30, 1.d30, y2tmp)

        do j = 1, m
          y2ay(i,j) = y2tmp(j)
        end do
      end do

      return
      end

***********************************************************************
***********************************************************************

      subroutine splin2(ypkt)

*---------------------------------------------------------------------*
*     diese routine berechnet die splines, die zu den entsprechenden  *
*     stuetzstellen in der x-y-ebene.                                 *
*                                                                     *
*     yaz(l,m)  : funktion def. i. d. x-y-ebene; output der routine   *
*                 splin3.                                    (input)  *
*     l,m       : dim von x-y-ebene                          (input)  *
*     x1a(l)    : stuetzstellen in x - richtung              (input)  *
*     x2a(m)    :                  y - richtung              (input)  *
*     y2ay(l,m) : 2. ableitung nach y der funktion yaz, nur in der    *
*                 x-y-ebene                                  (input)  *
*     xval(5)   : x-koordinate der pkte., an denen die funktion       *
*                 bestimmt werden soll.(max = 5)             (input)  *
*     yval(5)   : y-koordinate der pkte., an denen die funktion       *
*                 bestimmt werden soll.(max = 5)             (input)  *
*     npkt      : anzahl der werte, die berechnet werden sollen       *
*                                                            (input)  *
*     ressp(5)  : splinewerte zu den npkt punkten           (output)  *
*---------------------------------------------------------------------*

      implicit none

      include"common"
      include"commsp"
      real*8     ypkt

      integer i, j

      integer nmax
      parameter(nmax =n)
c      parameter(nmax = max(l,m,n))
c      parameter(nmax = 17)
      real*8     y2tmp(nmax),ytmp(nmax), yytmp(nmax), yy2tmp(nmax)

      integer ipkt

        do i =1 , l
          do j = 1, m
            ytmp(j)  = yaz(i,j)
            y2tmp(j) = y2ay(i,j)
          end do

          call splint(x2a,ytmp,y2tmp,m,ypkt,yytmp(i))

        end do

*--------ableitungen nach x bilden ----------------------------------*

        call spline(x1a,yytmp,l,1.d30,1.d30,yy2tmp)

        do ipkt = 1,lg
          call splint(x1a,yytmp,yy2tmp,l,xval(ipkt),ressp(ipkt))
        end do


      return
      end

**********************************************************************
**********************************************************************
      subroutine splint(xa,ya,y2a,n,x,y)

***********************************************************************
*     diese routine berechnet den wert eines splines an der stelle    *
*     x, der returnwert ist y.                                        *
*     xa(1:n), ya(1:n) tabellierte funktion; xa, ya, n stuetzstellen  *
**********************************************************************
      implicit none

      integer n
      real*8 x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      real*8 a,b,h

      klo=1
      khi=n

1     if (khi-klo.gt.1) then
        k=(khi+klo)/2

        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif

      goto 1
      endif

      h=xa(khi)-xa(klo)

      if (h.eq.0.) write(*,*) 'error: bad xa input in splint'

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.

      return
      end
c  (c) copr. 1986-92 numerical recipes software ,4#1-$@!.
*************************************************************************
************************************************************************
      subroutine spline(x,y,n,yp1,ypn,y2)

*************************************************************************
*     diese routine berechnet die 2. Ableitung einer funktion f(x),     *
*     wobei als unabhaengige variablen die eintraege des feldes         *
*     x(1:n) verwendet werden. die funktionswerte selbst sind an die-   *
*     sen  stellen in y(1:n) tabelliert.                                *
*     yp1 und ypn enthalten die ersten ableitungen an den raendern      *
*     der splines. es sind eingabeparameter; falls sie groesser als     *
*     1e30 sind, so wird der natuerliche spline verwendet, sprich die   *
*     ableitung am rand ist null.                                       *
*     die rueckgabe der 2. ableitungen erfolgt in dem feld y2(1:n)      *
************************************************************************
      implicit none
      integer n,nmax
      real*8 yp1,ypn,x(n),y(n),y2(n)

      parameter (nmax=500)

      integer i,k
      real*8 p,qn,sig,un,u(nmax)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue

      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
c  (c) copr. 1986-92 numerical recipes software ,4#1-$@!.
*****************************************************************************
