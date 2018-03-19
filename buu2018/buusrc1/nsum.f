      subroutine thermodyn

      implicit none

      include"common"
      include"cominput"


      real*8 epsi
      parameter(epsi = 1.0e-03)

      integer ppnx, ppnz, ppnmax
      parameter(ppnx = 2 * maxx + 1)
      parameter(ppnz = 2 * maxz + 1)
      parameter(ppnmax = ppnx*ppnx*ppnz)

      integer iv(maxpar), icopy(maxpar), ipoint(maxpar)
      integer igp(ppnmax), ngp(ppnmax)

      common / sortcom / iv, icopy, igp, ngp, ipoint

      integer l, i
      real*8    px, py, pz, rearmd, rearsky, potanal
      real*8    meff
      real*8    en, entot
      integer  jj
      real*8     rx,ry,rz, pabs , gammar
      real*8    calccoul, coulen

      integer ix, iy,iz

      real*8    sd, enpart
      real*8    betax, betay, betaz
      real*8    j0, j1, j2, j3
      real*8    encell, ehelp, plrf(1:3), rhap(1:3)
      real*8    uges
      real*8    pboo(0:3), pxtest, pytest, pztest, pxstat, pystat,pzstat
      real*8    pxneu, pyneu, pzneu, p1x, p1y, p1z
      real*8    masse, pst, scneu, scold

      integer itcount, jx, jy, jz
      logical flagit, raus

      real*8    pxcell, pycell, pzcell

c      real*8    as, bs, taus
c      data as,bs,taus  /-286.99993, 233.6517, 1.225707/

        write(*,*)'vor sorttp'
        call  sorttp( ntotal)
        write(*,*)'nach sorttp'

c        if(icomp.ne.3) then
c          write(*,*)'testversion nur icomp = 3 zugelassen',icomp
c          stop
c        end if

      entot = 0.0

            do 43 iz = -maxz, maxz
            do 43 iy = -maxx, maxx
            do 43 ix = -maxx, maxx

              pxcell = 0.0
              pycell = 0.0
              pzcell = 0.0
              encell = 0.0


            i = 1 + ( maxx+ix ) + ppnx * ( maxx+iy )
     +            + ppnx*ppnx * ( maxz+iz )


           if( ngp(i) .gt. 0) then


            sd   = rhob_4(4,ix,iy,iz)
            j0   = rhob_4(0,ix,iy,iz)
            j1   = rhob_4(1,ix,iy,iz)
            j2   = rhob_4(2,ix,iy,iz)
            j3   = rhob_4(3,ix,iy,iz)
c      if(ix .eq. -4 .and. iy .eq.-2) write(*,*)
c    1              ' upot 10', j0, j1, j2, j3, iz

            if(j0.gt.1.0e-06) then
              betax = j1/j0
              betay = j2/j0
              betaz = j3/j0
            else
              betax = 0.0
              betay = 0.0
              betaz = 0.0
            end if

            do 50 jj=igp(i)+1, igp(i)+ngp(i)
              l = ipoint(jj)

              rx   = r(1,l)
              ry   = r(2,l)
              rz   = r(3,l)

              rhap(1) = rx
              rhap(2) = ry
              rhap(3) = rz

              jx   = nint(rx)
              jy   = nint(ry)
              jz   = nint(rz)

              if(jx.ne.ix) then
                write(*,*)'something wrong with jx :', ix, jx
              end if

              if(jy.ne.iy) then
                write(*,*)'something wrong with jy :', iy, jy
              end if

              if(jz.ne.iz) then
                write(*,*)'something wrong with jz :', iz, jz
              end if

              meff = e(l)+upot(l)
              px   = p(1,l)
              py   = p(2,l)
              pz   = p(3,l)
              en   = sqrt(meff**2 + px**2+py**2+pz**2)

              pxstat = px
              pystat = py
              pzstat = pz
              scneu  = upot(l)
              masse  = e(l)

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
c      if(ix .eq. -4 .and. iy .eq.-2) write(*,*)
c    1              ' upot 12',itcount, scneu

*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*         boost prticle to its LRF                                     *

          call lorentz(betax,betay,betaz,pboo(1),pboo(2),pboo(3),
     &                 pboo(0))
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

              p1x = pboo(1)
              p1y = pboo(2)
              p1z = pboo(3)
              pabs= sqrt(p1x**2+p1y**2+p1z**2)
              pst = pabs
          sd = sqrt(j0**2 - j1**2 - j2**2 - j3**2)
c          uproprho = as*sd/rho0 + 2.*bs/(taus+1)*sd**taus/rho0**taus
c          uproprho = uproprho/1000.

c          urear    = bs*(taus-1)/(taus+1)*sd**taus/rho0**taus
c          urear    = urear/1000.
c          uges     = as*sd/rho0 + bs*sd**taus/rho0**taus
c          uges     = uges/1000.

             plrf(1) = p1x
             plrf(2) = p1y
             plrf(3) = p1z

          uges  = potanal(rho0,sd,plrf,1,rhap)

          ehelp = sqrt(masse**2 + pabs**2)
          uges  = - masse + sqrt((ehelp+uges)**2-pabs**2)

          meff = e(l)+uges

          en = sqrt(pabs**2 + meff**2)
          pboo(0)  = en

*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*             boost back to calc. frame                              *

          if(pboo(1)**2+pboo(2)**2+pboo(3)**2.gt.pboo(0)**2) then
            write(*,*) "hiba nsum2 lorentz, negative mass",
     &                pboo(1),pboo(2),pboo(3),pboo(0)
c            stop
          end if
          call lorentz(-betax,-betay,-betaz,pboo(1),pboo(2),pboo(3),
     &                 pboo(0))

*                                                                    *
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
             pxneu = pboo(1)
             pyneu = pboo(2)
             pzneu = pboo(3)
             scneu = uges

*            check iteration condition

             pxtest = abs(pxneu - pxstat)
             pytest = abs(pyneu - pystat)
             pztest = abs(pzneu - pzstat)

             if(pxtest.le.epsi .and. pytest.le.epsi .and.
     &          pztest.le.epsi) then
                enpart = enpart + encell

               flagit = .true.

             gammar = 1./sqrt(1. - betax**2 - betay**2 - betaz**2)


             plrf(1) = p1x
             plrf(2) = p1y
             plrf(3) = p1z

c             call rearr(icomp,0,sd,plrf,rhap,rearsky0,rearmd0)
c             call rearr(icomp,1,sd,plrf,rhap,rearsky1,rearmd1)

c             if(abs(rearsky0-rearsky1).ne.0.0) then
c               write(*,*)'rearr prob 1 ', rearsky0, rearsky1
c               stop
c             end if
c            if(abs(rearmd0-rearmd1).ne.0.0) then
c               write(*,*)'rearr prob 2 ', rearmd0, rearmd1
c               stop
c             end if

             call rearr(icomp,sd,plrf,rhap,rearsky,rearmd)

               encell =
     &            encell + en -1./gammar*(rearsky+rearmd)

               pxcell = pxcell + p1x
               pycell = pycell + p1y
               pzcell = pzcell + p1z

             end if

             if(itcount.ge.100) then
               write(*,*)'in iteration gradu 1 nsum ',
     1                     itcount, ix, iy, iz
               write(*,*)'old ', pxstat, pystat, pzstat, masse,scold
               write(*,*)'neu ', pxneu, pyneu, pzneu, masse,scneu
c              stop
               goto 50
             end if

            end do

 50          continue

*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*             boost back to calc. frame                              *

          if(pxcell**2+pycell**2+pzcell**2.ge.encell**2) then
            write(*,*) "nsum3 lorentz, negative mass",
     &                encell,pxcell,pycell,pzcell
            stop
          end if
          call lorentz(-betax,-betay,-betaz,pxcell,pycell,pzcell,
     &                 encell)

*                                                                    *
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

          entot = entot + encell

           end if
 43          continue

       coulen = 0.0
***************** add Coulomb energy ********************************
       if(ipou.eq.1) then
       coulen = calccoul()
       end if

*********************************************************************

           do l = 1, maxppar
             if(ipi(1,l).ne.0) then
              meff   = epi(l)
              px     = ppi(1,l)
              py     = ppi(2,l)
              pz     = ppi(3,l)
              en     = sqrt(meff**2 + px**2+py**2+pz**2)
              entot = entot + en
             end if
           end do


************* add particles which left the grid ********************
           raus = .false.
       do i = 1, maxpar
         if(id(1,i).ne.0 .and. (nint(abs(r(1,i))).gt.maxx .or. 
     &   nint(abs(r(2,i))).gt.maxx .or. nint(abs(r(3,i))).gt.maxz)) then

              raus   = .true.
              meff   = e(i)
              px     = p(1,i)
              py     = p(2,i)
              pz     = p(3,i)
              en     = sqrt(meff**2 + px**2+py**2+pz**2)

              entot  = entot + en
          end if
       end do
       entot = entot/float(num)

********** add coulomb

       entot = entot + coulen

       write(*,*) 'thermo(nsum)', time,
     &      (entot)/float(massta+masspr) - rmass ,raus




      return
      end


*********************************************************
      real*8 function calccoul()

      implicit none

      include"coucom"


      real*8 coulen
      integer cix, ciy, ciz


      coulen = 0.0
      do cix = -cmaxx, cmaxx
        do ciy = -cmaxx, cmaxx
          do  ciz = -cmaxz,cmaxz
            coulen = coulen + 0.5* chrho(0,cix,ciy,ciz)*
     &                     emrpot(0,cix,ciy,ciz)*dgrid**3.0
          end do
        end do
      end do


      calccoul = coulen
      return
      end
