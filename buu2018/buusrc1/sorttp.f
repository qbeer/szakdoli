************************************************************************
      subroutine entro_coord(num)
*
**  do the double sum

      implicit none
      include'common'

******************************************************************
      integer ppnx, ppnz, ppnmax
      parameter(ppnx = 2 * maxx + 1)
      parameter(ppnz = 2 * maxz + 1)
      parameter(ppnmax = ppnx*ppnx*ppnz)

      integer iv(maxpar), icopy(maxpar), ipoint(maxpar)
      integer igp(ppnmax), ngp(ppnmax)

      common / sortcom / iv, icopy, igp, ngp, ipoint

      real*8   entro(0:3,-maxx:maxx,-maxx:maxx,-maxz:maxz)
      real*8   storeen(0:3), totentro(0:3)
      integer ix0, iy0, iz0
      integer ipx, ipy, ipz
      integer jj, ipa


      integer i, j, icellx, icelly, icellz
      integer idstore, num
      real*8    rad, zahl, count
      integer ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly, lz
      real*8    px, py, pz, ppx, ppy, ppz, sek, phase
      real*8    pad, pad2, mspi, miso
      integer ic, ie, indx, indz
c      real*8    tin, tout, t0

c      t0  = 0.0
c      tin = secnds(t0)

*

      rad    = 3.0
      zahl   = 1./6.
      count  = zahl*float(num)


      indx =  1
      indz =  1

*-------------- loop over grid --------------------------------*
        do 40 icellz = -indz, indz
        do 40 icelly = -indx, indx
        do 40 icellx = -indx, indx



          ix0 = icellx
          iy0 = icelly
          iz0 = icellz

          write(*,*)ix0, iy0, iz0

*------------------- momentum integration ---------------------*

          do i = 0,3
            storeen(i)  = 0.0
            entro(i,icellx,icelly,icellz) = 0.0
          end do

          do 42 ipx = -5, 5
            px = float(ipx)*0.2
          do 42 ipy = -5, 5
            py = float(ipy)*0.2
          do 42 ipz = -5, 5
            pz = float(ipz)*0.2

c          write(*,*)'ips ', ipx, ipy, ipz

            sek = 0.0
*------------- loop over point +/- 1 --------------------------*
            do 43 iz = icellz-1, icellz+1
            do 43 iy = icelly-1, icelly+1
            do 43 ix = icellx-1, icellx+1
              if(abs(iz).gt.maxz) goto 43
              if(abs(iy).gt.maxx) goto 43
              if(abs(iz).gt.maxx) goto 43



            i = 1 + ( maxx+ix ) + ppnx * ( maxx+iy )
     +            + ppnx*ppnx * ( maxz+iz )
           if( ngp(i) .gt. 0) then

            do 50 jj=igp(i)+1, igp(i)+ngp(i)
              j = ipoint(jj)
              idstore = id(1,j)
              ppx = p(1,j)
              ppy = p(2,j)
              ppz = p(3,j)
              if(idstore.eq.1) ipa = 1
              if(idstore.eq.2) ipa = 2
              if(idstore.gt.2) ipa = 3

              if(idstore.eq.1) then
                mspi = 2.0
                miso = 2.0
              else
                mspi = float(resprop2(idstore-1,2)) + 1.0
                miso = float(resprop2(idstore-1,1)) + 1.0
              end if

              pad    = (3.0*zahl/(4.0*pi*mspi*miso))**(1./3.)
     &                *(2.*pi*.197)/rad
              pad2   = pad**2


             if((px-ppx)**2+(py-ppy)**2+(pz-ppz)**2.gt.pad2)
     &                                           goto 50
*
             jx=nint(r(1,j))
             kx=ix0-jx
             if(abs(kx).gt.ir) goto 50
             jy=nint(r(2,j))
             ky=iy0-jy
             if(abs(ky).gt.ir) goto 50
             jz=nint(r(3,j))
             kz=iz0-jz
             if(abs(kz).gt.ir) goto 50
*
             lx=nint(float(2*ip+1)*(r(1,j)-float(jx)))
             if(abs(lx) .eq. ip+1) lx = lx/abs(lx) * ip
             ly=nint(float(2*ip+1)*(r(2,j)-float(jy)))
             if(abs(ly) .eq. ip+1) ly = ly/abs(ly) * ip
             lz=nint(float(2*ip+1)*(r(3,j)-float(jz)))
             if(abs(lz) .eq. ip+1) lz = lz/abs(lz) * ip
*
             ic=1+(lz+ip)+(ly+ip)*(2*ip+1)+(lx+ip)*(2*ip+1)**2
             ie=1+(kz+ir)+(ky+ir)*(2*ir+1)+(kx+ir)*(2*ir+1)**2
*
             if(dm(ic,ie).le.1.0e-12) dm(ic,ie) = 0.0
             sek=sek+dm(ic,ie)*float(num)

 50          continue

           end if
 43          continue


             phase = amin1(1.0,sek/count)
             if(phase.gt.1.0e-03) then
              storeen(0)     = storeen(0)+phase*log(phase)
              storeen(ipa) = storeen(ipa)+phase*log(phase)
             end if
             if(phase.le.0.9999) then
              storeen(0) = storeen(0) + (1.0-phase)*log(1.0-phase)
              storeen(ipa)=storeen(ipa)+(1.0-phase)*log(1.0-phase)
             end if

  42          continue
             do i = 0,3
               entro(i,icellx,icelly,icellz) = storeen(i)*(0.2**3)
             end do

 40          continue


*--------------- space integration ---------------------------------
         do i = 0,3
           totentro(i) = 0.0

           do  icellz = -indz, indz
             do  icelly = -indx, indx
               do  icellx = -indx, indx
                 totentro(i) = totentro(i)+
     &               entro(i,icellx,icelly,icellz)
                end do
              end do
            end do
            totentro(i) = totentro(i)*(-1.0)
          end do

          write(*,*)'entropy = ', totentro(0),totentro(1),totentro(2),
     &                            totentro(3)

c      tout = secnds(t0)
c     write(*,*)'time ellapsed in entropy : ', tout-tin ,'secnds.'

      return

      end



************************************************************************
*                                                                      *
      subroutine sorttp( ntotal)
*                                                                      *
*     purpose: sort testparticle-vectors according to the spatial      *
*              position of the particles on the grid                   *
*                                                                      *
************************************************************************
      implicit none
      include 'common'
******************************************************************

      integer ppnx, ppnz, ppnmax
      parameter(ppnx = 2 * maxx + 1)
      parameter(ppnz = 2 * maxz + 1)
      parameter(ppnmax = ppnx*ppnx*ppnz)

      integer iv(maxpar), icopy(maxpar), ipoint(maxpar)
      integer igp(ppnmax), ngp(ppnmax)

      common / sortcom / iv, icopy, igp, ngp, ipoint

**************************************************************


      integer i, nlast, ntotal
      integer ix, iy, iz, ixc, iyc, izc, i1, i2, i3
      integer itemp(maxpar)
      real*8    zwi, rix, riy, riz
c      real*8    tin, tout, t0

      logical warning
      warning = .false.
c      t0 = 0.0
c      tin = secnds(t0)

      do i = 1, maxpar
        ipoint(i) = i
        if(id(1,i).eq.0)  ipoint(i) = 0
      end do

         warning =.false.

      nlast = ntotal
      do 10 i = 1, maxpar
        if(id(1,i).eq.0) goto 10
        ix = nint(r(1,i))
        iy = nint(r(2,i))
        iz = nint(r(3,i))
        if( abs(ix) .gt. maxx ) then
           warning = .true.
           goto 10
        end if
        if( abs(iy) .gt. maxx ) then
           warning = .true.
           goto 10
        end if
        if( abs(iz) .gt. maxz ) then
            warning = .true.
            goto 10
        end if
        iv(i) = 1 + ( maxx+ix ) + ppnx * ( maxx+iy )
     +            + ppnx*ppnx * ( maxz+iz )
        ixc = iv(i)-1 - maxx
        iyc = iv(i)- ppnx - maxx
        izc = iv(i)- ppnx**2 - maxz
        i1  = iv(i)- ppnx
        i2  = i1   -ppnx**2
        zwi = float(i2)/ppnx**2
        i3  = int(zwi) - maxz + 1
        izc = i3

        i1  = iv(i)- ppnx
        i2  = i1   -ppnx**2*(maxz+izc)
        zwi = float(i2)/ppnx
        i3  = int(zwi) - maxx + 1
        iyc = i3


        i2  = iv(i) -ppnx**2*(maxz+izc)-ppnx*(maxx+iyc)
        zwi = float(i2)
        i3  = int(zwi) - maxx - 1
        ixc = i3

        riz = abs(float(iz -izc))
        riy = abs(float(iy -iyc))
        rix = abs(float(ix -ixc))


c        if(riz.gt.1e-08 .or. riy.gt.1e-08.or. rix.gt.1e-08) then
c          write(115,*)'problem in time step ', nt
c          write(115,*)ix, ixc
c          write(115,*)iy, iyc
c          write(115,*)iz, izc
c        end if


10    continue

      do 20 i = 1, ppnx*ppnx*ppnz
        igp(i) = 0
20    continue


      if( warning) then

        do 30 i = 1, maxpar
          if(id(1,i).eq.0) goto 30

          ix = nint(r(1,i))
          iy = nint(r(2,i))
          iz = nint(r(3,i))

          if(abs(ix).le.maxx.and.abs(iy).le.maxx.and.
     &       abs(iz).le.maxz) then
            igp(iv(i)) = igp(iv(i)) + 1
          end if
30      continue
      else
        do 40 i = 1, maxpar
          if(id(1,i).ne.0) igp(iv(i)) = igp(iv(i)) + 1
40      continue
      end if


      do 50 i = 1, ppnx*ppnx*ppnz
         ngp(i) = igp(i)
50    continue

      do 60 i = 2, ppnx*ppnx*ppnz
        igp(i) = igp(i) + igp(i-1)
60    continue
c       do i = 1, ppnmax
c        write(15,*)i,igp(i)
c       end do

      if( warning ) then
        do 70 i = 1, maxpar
          if(id(1,i).eq.0) goto 70
          ix = nint(r(1,i))
          iy = nint(r(2,i))
          iz = nint(r(3,i))

          if(abs(ix).le.maxx.and.abs(iy).le.maxx.and.
     &       abs(iz).le.maxz) then
            icopy(i) = igp(iv(i))
            igp(iv(i)) = igp(iv(i)) - 1
          else
            icopy(i) = nlast
            nlast = nlast - 1
        end if
70      continue
      else
        do 80 i = 1, maxpar
          if(id(1,i).eq.0) goto 80
          icopy(i) = igp(iv(i))
          igp(iv(i)) = igp(iv(i)) - 1
80      continue
      end if



      do 170 i = 1, maxpar
        if(id(1,i).ne.0) itemp(icopy(i)) = ipoint(i)
170   continue

      do 180 i = maxpar, 1, -1
        ipoint(i) =  itemp(i)
c        if(ipoint(i).eq.0) then
c          write(*,*)'pointer gleich null ', i
c        end if
180   continue

c      tout = secnds(t0)
c     write(*,*)'time ellapsed in sort : ', tout-tin ,'secnds.'

      return
      end

