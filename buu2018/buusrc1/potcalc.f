************************************************************************
      subroutine potcalc
*                                                                      *
*     This routine calculates the rel.scalar potential (using potanal) *
*      for each particle in the corresponding local rest frame.                *
*      The  s c a l a r !!! potentials are stored in the vector upot(i)*
*                                                                      *
***********************************************************************
      implicit none
      include"common"
      include"cominput"
      include"commsp"

      integer irun, j, i, id1
      real*8    x, y, z, px, py, pz, en, scapo, masse
      real*8    deriv(0:4), j0, j1, j2, j3
      real*8    betlrfx, betlrfy, betlrfz
      real*8    pabs, pin , vecpo, potanal, potmes
c      real*8    tin, tout, t0
      integer ix, iy, iz
      real*8    pb0, pb1, pb2, pb3, psinv, pbinv
      logical flag
      integer icount
      real*8    rhap(1:3), plrf(1:3)

c      t0 = 0.0
c      tin = secnds(t0)

      if(isplipi.eq.1) then
        call spline3di1
        call spline3di
      end if

      write(*,*)'begin potcalc'
      if(isplipi .eq. 1) then
        write(*,*)'spline option is switched off '
        write(*,*)'if YOU want to switch it on make sure that the '
        write(*,*)'potentials are build in correctly'
        stop
*     loop over all parallel runs
        do irun = 1,num

*     loop over all pseudoparticles in the same run
          do  j = 1,maxb
            i  = j + (irun - 1) * maxb
            if(id(1,i).gt.0) then
            x     = r(1,i)
            y     = r(2,i)
            z     = r(3,i)
            px    = p(1,i)
            py    = p(2,i)
            pz    = p(3,i)
            masse = e(i)
            scapo = upot(i)
            en    = sqrt( (masse+scapo)**2 + px**2 + py**2 + pz**2 )
            id1   = id(1,i)

            rhap(1) = x
            rhap(2) = y
            rhap(3) = z

            call splinint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

            if(j0 .gt. 1.0e-6) then
              betlrfx = j1/j0
              betlrfy = j2/j0
              betlrfz = j3/j0
            else
              write(*,*)'warning from potcalc j0 = ', j0
              betlrfx = 0.0
              betlrfy = 0.0
              betlrfz = 0.0
            end if

*         det. density in LRF
              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba potcalc lorentz, negative mass",
     &              j1,j2,j3,j0
                 stop
              end if
          call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)

          icount = 0
          flag   = .true.
          do while(flag)
            icount = icount + 1

            pb1 = px
            pb2 = py
            pb3 = pz
            pb0 = sqrt((masse+scapo)**2+pb1**2+pb2**2+pb3**2)
            pbinv  = pb0**2 - pb1**2 - pb2**2 - pb3**2
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba potcalc2 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
            stop
          end if
            call lorentz(betlrfx, betlrfy, betlrfz, pb1, pb2, pb3, pb0)

            pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
            pin    = pabs
            plrf(1)= pb1
            plrf(2)= pb2
            plrf(3)= pb3

            vecpo  = potanal(rho0, j0, plrf , id1, rhap)

            scapo  = -masse + sqrt(masse**2 +
     +               2.0*sqrt(pabs**2+masse**2)*vecpo + vecpo**2)
            psinv  = (masse+scapo)**2

            if(abs(pbinv-psinv).lt.1.0e-04) then
              flag = .false.
            end if
            if(icount.gt.100) then
              write(*,*)'hiba 1 in potcalc',icount
              stop
            end if
          end do
            upot(i)= scapo
            betlrfboo(i,1) = betlrfx
            betlrfboo(i,2) = betlrfy
            betlrfboo(i,3) = betlrfz
            betlrfboo(i,4) = j0
          end if
          end do
        end do

      else if(isplipi .eq. 0) then

*     loop over all parallel runs
        do irun = 1,num

*     loop over all pseudoparticles in the same run
          do  j = 1,maxb
            i  = j + (irun - 1) * maxb
            if(id(1,i).gt.0) then
            x     = r(1,i)
            y     = r(2,i)
            z     = r(3,i)
            px    = p(1,i)
            py    = p(2,i)
            pz    = p(3,i)
            masse = e(i)
            scapo = upot(i)
            en    = sqrt( (masse+scapo)**2 + px**2 + py**2 + pz**2 )
            id1   = id(1,i)

            rhap(1) = x
            rhap(2) = y
            rhap(3) = z

c            write(*,*) 'potcalc before linint1 en, mass',en,masse
            call linint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

c            write(*,*) 'potcalc after linint1 j',j0,j1,j2,j3
            if(j0 .gt. 1.0e-6) then
              betlrfx = j1/j0
              betlrfy = j2/j0
              betlrfz = j3/j0
            else
              ix = nint(x)
              iy = nint(y)
              iz = nint(z)

              betlrfx = 0.0
              betlrfy = 0.0
              betlrfz = 0.0
            end if

*         det. density in LRF
c            write(*,*) 'potcalc before lorentz beta,j',
c     &              betlrfx,betlrfy,betlrfz,j0,j1,j2,j3
              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba potcalc3 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
          call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)
c         if (i .eq. 3112)  write(*,*)  ' in potcalc 1 ',
c    1               betlrfx, betlrfy, betlrfz, j1, j2, j3, j0

          icount = 0
          flag   = .true.
          do while(flag)
            icount = icount + 1

            pb1 = px
            pb2 = py
            pb3 = pz
            pb0 = sqrt((masse+scapo)**2+pb1**2+pb2**2+pb3**2)
c            pbinv  = pb0**2 - pb1**2 - pb2**2 - pb3**2
            pbinv  = (masse+scapo)**2
c            write(*,*) 'potcalc before lorentz2 beta,j',
c     &              betlrfx,betlrfy,betlrfy,pb0,pb1,pb2,pb3
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba potcalc4 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
            stop
          end if
            call lorentz(betlrfx, betlrfy, betlrfz, pb1, pb2, pb3, pb0)
c         if (i .eq. 3112)  write(*,*)  ' in potcalc 1a ',
c    1               px, py, pz, pb1, pb2, pb3, pb0

            pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
            pin    = pabs
            plrf(1)= pb1
            plrf(2)= pb2
            plrf(3)= pb3

c            write(*,*) 'potcalc before potanal',rho0,j0,plrf,id1,rhap
            vecpo  = potanal(rho0, j0, plrf , id1,rhap)
c            if (i .eq. 3112)  write(*,*) ' in potcalc 2 ', icount,
c     1           rho0, masse, j0, plrf , id1,rhap

            scapo  = -masse + sqrt(masse**2 +
     +               2.0*sqrt(pabs**2+masse**2)*vecpo + vecpo**2)
            psinv  = (masse+scapo)**2

            if(abs(pbinv-psinv).lt.1.0d-03) then
              flag = .false.
            end if
            if(icount.gt.100) then
              write(*,*)'hier laeuft was schief 2 in potcalc',icount
              write(*,*)i,id(1,i), id(2,i),masse
              write(*,*)betlrfx, betlrfy, betlrfz,j0
              write(*,*)pabs, pbinv, psinv
              write(*,*)scapo, vecpo
              goto 11
            end if
          end do
 11         continue
            upot(i)= scapo
            betlrfboo(i,1) = betlrfx
            betlrfboo(i,2) = betlrfy
            betlrfboo(i,3) = betlrfz
            betlrfboo(i,4) = j0

          end if
          end do
        end do
      else
        write(*,*)'potcalc : something wrong with isplipi = ', isplipi
        stop
      end if

      write(*,*)'in potcalc mesons', iseed
*************************************************************************
*       now do the mesons                                               *

      if(isplipi .eq. 1) then

*     loop over all mesons
        do i = 1, maxppar
          if(ipi(1,i) .ne. 0) then
            x     = rpi(1,i)
            y     = rpi(2,i)
            z     = rpi(3,i)
            px    = ppi(1,i)
            py    = ppi(2,i)
            pz    = ppi(3,i)
            masse = epi(i)
            scapo = mpot(i)
            en    = sqrt( (masse+scapo)**2 + px**2 + py**2 + pz**2 )
            id1   = ipi(1,i)

            call splinint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

            if(j0 .gt. 1.0e-6) then
              betlrfx = j1/j0
              betlrfy = j2/j0
              betlrfz = j3/j0
            else
              betlrfx = 0.0
              betlrfy = 0.0
              betlrfz = 0.0
            end if

              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba potcalc5 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
            call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)
            icount = 0
            flag   = .true.
            do while(flag)
              icount = icount + 1

              pb1 = px
              pb2 = py
              pb3 = pz
              pb0 = sqrt((masse+scapo)**2+pb1**2+pb2**2+pb3**2)
              pbinv  = pb0**2 - pb1**2 - pb2**2 - pb3**2
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba potcalc6 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
c            stop
          end if
           call lorentz(betlrfx, betlrfy, betlrfz, pb1, pb2, pb3, pb0)

              pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
              pin    = pabs
              vecpo  = potmes( j0, pin , id1)
              scapo  = -masse + sqrt(masse**2 +
     +                 2.0*sqrt(pabs**2+masse**2)*vecpo + vecpo**2)

              psinv  = (masse+scapo)**2

              if(abs(pbinv-psinv).lt.1.0e-04) then
                flag = .false.
              end if
              if(icount.gt.100) then
                write(*,*)'hier laeuft was schief 3 in potcalc',icount
                stop
              end if
            end do

            mpot(i)        = scapo
            betlrfbom(i,1) = betlrfx
            betlrfbom(i,2) = betlrfy
            betlrfbom(i,3) = betlrfz
            betlrfbom(i,4) = j0
          end if
        end do
      else if(isplipi .eq. 0) then

*     loop over all mesons
        do i = 1, maxppar
          if(ipi(1,i).ne.0) then
            x     = rpi(1,i)
            y     = rpi(2,i)
            z     = rpi(3,i)
            px    = ppi(1,i)
            py    = ppi(2,i)
            pz    = ppi(3,i)
            masse = epi(i)
            scapo = mpot(i)
            en    = sqrt( (masse+scapo)**2 + px**2 + py**2 + pz**2 )
            id1   = ipi(1,i)

c        write(*,*) ' pseudo particles in potcalc ',i,id1,epi(i),en,
c     1               scapo

            call linint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

            if(j0 .gt. 1.0e-6) then
              betlrfx = j1/j0
              betlrfy = j2/j0
              betlrfz = j3/j0
            else
              betlrfx = 0.0
              betlrfy = 0.0
              betlrfz = 0.0
            end if

              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba potcalc7 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
            call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)
            icount = 0
            flag   = .true.
            do while(flag)
              icount = icount + 1

              pb1 = px
              pb2 = py
              pb3 = pz
              pb0 = sqrt((masse+scapo)**2+pb1**2+pb2**2+pb3**2)
              pbinv  = pb0**2 - pb1**2 - pb2**2 - pb3**2
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba potcalc8 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
c            stop
          end if
            call lorentz(betlrfx, betlrfy, betlrfz, pb1, pb2, pb3, pb0)

              pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
              pin    = pabs
              vecpo  = potmes(  j0, pin , id1)
              scapo  = -masse + sqrt(masse**2 +
     +                 2.0*sqrt(pabs**2+masse**2)*vecpo + vecpo**2)

              psinv  = (masse+scapo)**2

              if(abs(pbinv-psinv).lt.1.0e-04) then
                flag = .false.
              end if
              if(icount.gt.100) then
                write(*,*)'hier laeuft was schief 4 in potcalc',
     1          icount, i, id1, x,y,z, px,py,pz, pbinv,psinv
c                stop
                go to 14
              end if
            end do
 14         continue
            mpot(i)= scapo
            betlrfbom(i,1) = betlrfx
            betlrfbom(i,2) = betlrfy
            betlrfbom(i,3) = betlrfz
            betlrfbom(i,4) = j0

          end if
        end do

      else
        write(*,*)'potcalc : something wrong with isplipi = ', isplipi
        stop
      end if

c      write(*,*)'end of potcalc'

c      tout = secnds(t0)
c      tin = tout - tin
c     write(*,*)'time ellapsed in potcalc = ',tin,'  sec.'
      return
      end

************************************************************************
      subroutine potcalc_i(nparticle, ibar, vx, vy, vz, pot)
*                                                                      *
*      This routine calculates the potentials (using potanal) for      *
*      particle nparticle in the corresponding local rest frame.        *
*      The  s c a l a r !!! potentials are stored in  pot
*      nparticle  - number of particle
*      ibar      - 1/0   for baryon / meson                            *
*      vx, vy, vz-       velocity
*      pot               result
***********************************************************************
      implicit none
      include"common"
      include"cominput"
      include"commsp"

      integer  i, id1, ibar, nparticle
      real*8    x, y, z, px, py, pz, en, scapo, masse, vx, vy, vz
      real*8    deriv(0:4), j0, j1, j2, j3
      real*8    betlrfx, betlrfy, betlrfz, pot, gamma
      real*8    pabs, pin , vecpo, potanal, potmes
c      real*8    tin, tout, t0
      integer ix, iy, iz
      real*8    pb0, pb1, pb2, pb3, psinv, pbinv
      logical flag
      integer icount
      real*8    rhap(1:3), plrf(1:3)

c      t0 = 0.0
c      tin = secnds(t0)

c      write(*,*)'begin potcalc'
      if(isplipi .eq. 1) then

        write(*,*)'spline option is switched off '
        write(*,*)'if YOU want to switch it on make sure that the '
        write(*,*)'potentials are build in correctly'

        stop
      endif
      IF (ibar .eq. 1)   then
            i  = nparticle
            x     = r(1,i)
            y     = r(2,i)
            z     = r(3,i)
            masse = e(i)
            scapo = upot(i)
            gamma = 1. / sqrt(1. - vx*vx - vy*vy - vz*vz)
            px    = gamma * masse * vx
            py    = gamma * masse * vy
            pz    = gamma * masse * vz
            en    = sqrt( (masse+scapo)**2 + px**2 + py**2 + pz**2 )
            id1   = id(1,i)

            rhap(1) = x
            rhap(2) = y
            rhap(3) = z

            call linint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

c            write(*,*) "potcalc_i 0",id1,x,y,z,deriv,px,py,pz,en,masse
            if(j0 .gt. 1.0e-6) then
              betlrfx = j1/j0
              betlrfy = j2/j0
              betlrfz = j3/j0
            else
              ix = nint(x)
              iy = nint(y)
              iz = nint(z)

              betlrfx = 0.0
              betlrfy = 0.0
              betlrfz = 0.0
            end if

*         det. density in LRF
              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba potcalc9 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
          call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)
c          write(*,*)  ' in potcalc_i 1 ',
c     1      nparticle,masse, upot(nparticle), vx,vy,vz,
c     2             betlrfx, betlrfy, betlrfz, j1, j2, j3, j0

          icount = 0
          flag   = .true.
          do while(flag)
            icount = icount + 1

            pb1 = px
            pb2 = py
            pb3 = pz
            pb0 = sqrt((masse+scapo)**2+pb1**2+pb2**2+pb3**2)
            pbinv  = pb0**2 - pb1**2 - pb2**2 - pb3**2
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba potcalc10 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
c            stop
          end if
            call lorentz(betlrfx, betlrfy, betlrfz, pb1, pb2, pb3, pb0)
c            write(*,*)  ' in potcalc_i 1a ',
c     1               px, py, pz,  pb1, pb2, pb3, pb0


            pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
            pin    = pabs
            plrf(1)= pb1
            plrf(2)= pb2
            plrf(3)= pb3

            vecpo  = potanal(rho0, j0, plrf , id1,rhap)
c           if (nparticle .eq. 3112)  write(*,*) ' in potcalc_i 2 ',
c     1          icount, rho0, masse, j0, plrf , id1,rhap

            scapo  = -masse + sqrt(masse**2 +
     +               2.0*sqrt(pabs**2+masse**2)*vecpo + vecpo**2)
            psinv  = (masse+scapo)**2

            if(abs(pbinv-psinv).lt.1.0e-03) then
              flag = .false.
            end if
            if(icount.gt.100) then
              write(*,*)'hier laeuft was schief 2 in potcalc_i',icount
              write(*,*)i,id(1,i), id(2,i)
              write(*,*)betlrfx, betlrfy, betlrfz,j0
              write(*,*)pabs, pbinv, psinv
              write(*,*)scapo, vecpo
              stop
              goto 11
            end if
          end do
 11         continue
            pot= scapo
            betlrfboo(i,1) = betlrfx
            betlrfboo(i,2) = betlrfy
            betlrfboo(i,3) = betlrfz
            betlrfboo(i,4) = j0

      ELSEIF (ibar.eq. 0) then

            i     = nparticle
            if (ipi(1,i) .eq. 0)  then
              write(*, *) ' no meson in potcalc_i '
              stop
            endif
            x     = rpi(1,i)
            y     = rpi(2,i)
            z     = rpi(3,i)
            masse = epi(i)
            scapo = mpot(i)
            gamma = 1. / sqrt(1. - vx*vx - vy*vy - vz*vz)
            px    = gamma * masse * vx
            py    = gamma * masse * vy
            pz    = gamma * masse * vz
            en    = sqrt( (masse+scapo)**2 + px**2 + py**2 + pz**2 )
            id1   = ipi(1,i)
c            write(*,*) 'potcalc_i mes',px,py,pz,gamma,masse,vx,vy,vz
            call linint1(x,y,z,deriv)
            j0    = deriv(0)
            j1    = deriv(1)
            j2    = deriv(2)
            j3    = deriv(3)

            if(j0 .gt. 1.0e-6) then
              betlrfx = j1/j0
              betlrfy = j2/j0
              betlrfz = j3/j0
            else
              betlrfx = 0.0
              betlrfy = 0.0
              betlrfz = 0.0
            end if

              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba potcalc11 lorentz, negative mass",
     &              j1,j2,j3,j0
c                 stop
              end if
            call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)
            icount = 0
            flag   = .true.
            do while(flag)
              icount = icount + 1

              pb1 = px
              pb2 = py
              pb3 = pz
              pb0 = sqrt((masse+scapo)**2+pb1**2+pb2**2+pb3**2)
              pbinv  = pb0**2 - pb1**2 - pb2**2 - pb3**2
          if(pb1**2+pb2**2+pb3**2.gt.pb0**2) then
            write(*,*) "hiba potcalc12 lorentz, negative mass",
     &                pb1,pb2,pb3,pb0
c            stop
          end if
            call lorentz(betlrfx, betlrfy, betlrfz, pb1, pb2, pb3, pb0)

              pabs   = sqrt(pb1**2 + pb2**2 + pb3**2)
              pin    = pabs
              vecpo  = potmes(  j0, pin , id1)
              scapo  = -masse + sqrt(masse**2 +
     +                 2.0*sqrt(pabs**2+masse**2)*vecpo + vecpo**2)

              psinv  = (masse+scapo)**2

              if(abs(pbinv-psinv).lt.1.0e-04) then
                flag = .false.
              end if
              if(icount.gt.100) then
                write(*,*)'hier laeuft was schief 4 in potcalc_i',
     1         icount, ibar, nparticle, x,y,z, px,py,pz
                stop
              end if
            end do

            pot= scapo
            betlrfbom(i,1) = betlrfx
            betlrfbom(i,2) = betlrfy
            betlrfbom(i,3) = betlrfz
            betlrfbom(i,4) = j0

      else
        write(*,*)'potcalc : something wrong with isplipi = ', isplipi
        stop
      ENDIF

c      write(*,*)'end of potcalc_i'

c      tout = secnds(t0)
      return
      end

