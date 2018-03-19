************************************************************************
*                                                                      *
      subroutine dens_4
*                                                                      *
*       purpose:     calculation of nuclear 4-density from spatial     *
*                    distribution of testparticles for momentum dep.   *
*                    forces for the baryons and determination of the   *
*                    charge current for 4-dim coulomb                  *
*                                                                      *
*                                                                      *
*       variables (all input, all integer)                             *
*         minnum  -  first testparticle treated in one run for density *
*         maxnum  -  last testparticle treated in one run for density  *
*         mass    -  total mass                                        *
*                                                                      *
*                                                                      *
*     for threefluid uncomment call sphere                             *
************************************************************************

      implicit none
      include"common"
      include"cominput"
      integer ismear, ismeasy
      parameter(ismear  = 1)
      parameter(ismeasy = 0)

      integer i, ix, iy, iz, ndim,   jx, jy, jz, nst
      integer kx, ky, kz, ic, ib, idtest
      real*8    pstore(0:3), veloc,meff , ehelp

****************************************************************
      real*8    pboo(0:3), pxstat, pystat, pzstat
      real*8    j0, j1, j2, j3, sd
      real*8    pxneu, pyneu, pzneu, pst, pabs
      real*8    p1x, p1y, p1z, potanal, epsi, deriv(0:4)

      parameter(epsi = 1.0e-05)

      real*8  pxtest,pytest, pztest, dummy, masse, scneu, scold, vecpot
      logical flagit, flag
      integer itcount, id1
      real*8    rho, rx, ry, rz
      real*8    plrf(1:3), rhap(1:3)


*****************************************************************
*
      real*8 tin, tout, t0
      t0 = 0.0
c      tin =secnds(t0)
*----------------------------------------------------------------------*
*       set arrays for 4-baryon-density to zero                        *
*
      do ix = -maxx,maxx
        do iy = -maxx,maxx
          do iz = -maxz,maxz
            do ndim = 0,5       !  hw
              rhob_4(ndim, ix, iy, iz) = 0.0
            end do
            rhob_4(6   , ix, iy, iz) = 1.0
          end do
        end do
      end do


*
*-----------------------------------------------------------------------
*       set arrays for charge 4-densities to zero
*
      if(ipou .eq. 1) then

*        do ndim = 0,3
*        do ix = -maxxc,maxxc
*        do iy = -maxxc,maxxc
*        do iz = -maxzc,maxzc
*           chrho(ndim, ix, iy, iz) = 0.0
*        end do
*        end do
*        end do
*        end do

      end if
*
*-----------------------------------------------------------------------
*
*      determination of densities for mom-dep forces
*       (implicit assumption that grid size = 1 fm )!!!!!!!!

c      if(imomdep .eq. 1) then

      if(ismear .eq. 1) then
*
*       this part takes the smearing weights which are set in entry
*       densin of the subroutine dens (dens.f)
*       parameters are stored in common block PQ

c        write(*,*) 'dens_4 smear',ismear
        do i = 1, maxpar
          if(id(1,i).gt.0) then
*
          ix     = nint( r(1,i) )
          iy     = nint( r(2,i) )
          iz     = nint( r(3,i) )
          idtest = id(1,i)
c          write(*,*) 'dens_4 0',i,idtest,r(1,i),r(2,i),r(3,i),
c     &     p(1,i),p(2,i),p(3,i)
*

          if(abs(ix).le.maxx.and.abs(iy).le.maxx.and.abs(iz).le.maxz)
     &       then

            kx=nint(float(2*ip+1)*(r(1,i)-float(ix)))
            if(abs(kx) .eq. ip+1) kx = kx/abs(kx) * ip
            ky=nint(float(2*ip+1)*(r(2,i)-float(iy)))
            if(abs(ky) .eq. ip+1) ky = ky/abs(ky) * ip
            kz=nint(float(2*ip+1)*(r(3,i)-float(iz)))
            if(abs(kz) .eq. ip+1) kz = kz/abs(kz) * ip
            ic=1+(kz+ip)+(ky+ip)*(2*ip+1)+(kx+ip)*(2*ip+1)**2
            ib=0
            do jx=ix-iq,ix+iq
              do jy=iy-iq,iy+iq
                do jz=iz-iq,iz+iq

                  ib=ib+1

                  if(cm(ic,ib).gt.0.0.and.
     &               abs(jx).le.maxx.and.abs(jy).le.maxx.and.
     &               abs(jz).le.maxz) then
*
                    pstore(1) = p(1,i)
                    pstore(2) = p(2,i)
                    pstore(3) = p(3,i)
                    meff      = e(i) + upot(i)
                    pstore(0) = sqrt(meff**2 + pstore(1)**2 +
     &                 pstore(2)**2 + pstore(3)**2)

c                    write(*,*)'dens_4 1',e(i),upot(i),idtest,pstore
                    do nst = 0,3
                      veloc = pstore(nst)/pstore(0)
                      rhob_4(nst,jx,jy,jz)=rhob_4(nst,jx,jy,jz)
     &                   + float(sign(1,idtest))*cm(ic,ib)* veloc
                    end do


                    if( (id(1,i).eq.1) .and. (id(2,i).eq.1)) then
                      rhob_4(4,jx,jy,jz)=rhob_4(4,jx,jy,jz)+cm(ic,ib)
                    end if

                    if( (id(1,i).eq.1) .and. (id(2,i).eq.0)) then
                      rhob_4(5,jx,jy,jz)=rhob_4(5,jx,jy,jz)+cm(ic,ib)
                    end if

                    rhob_4(6,jx,jy,jz)=1.0
                    if( rhob_4(0,jx,jy,jz) .gt. 1.e-4) then
                      rhob_4(6,jx,jy,jz)=
     &                   1.0/sqrt(1.-( (rhob_4(1,jx,jy,jz)**2+
     &                   rhob_4(2,jx,jy,jz)**2+rhob_4(3,jx,jy,jz)**2)/
     &                   rhob_4(0,jx,jy,jz)**2)  )
                    end if

                  end if
*
c         if (jx .eq. -4  .and. jy .eq. 2) write(*,*) ' dens_4 10',
c    1        rhob_4(1,jx,jy,jz), veloc, ix, iy, iz, jx, jy, jz
                end do
              end do
            end do

*
          end if
          end if
        end do

      else if(ismear .eq. 0) then


        do i = 1,maxpar
*
          if(id(1,i).ne.0) then
          ix = nint( r(1,i) )
          iy = nint( r(2,i) )
          iz = nint( r(3,i) )
          idtest = id(1,i)

*
          if(abs(ix).lt.maxx.and.abs(iy).lt.maxx.and.abs(iz).lt.maxz)
     &       then

            rx     = r(1,i)
            ry     = r(2,i)
            rz     = r(3,i)

            rhap(1) = rx
            rhap(2) = ry
            rhap(3) = rz

            meff = e(i) + upot(i)
            pstore(1) = p(1,i)
            pstore(2) = p(2,i)
            pstore(3) = p(3,i)
            pstore(0) = sqrt(meff**2 + pstore(1)**2 +
     &         pstore(2)**2 + pstore(3)**2)
*
            id1       = id(1,i)

****************************** iterate eff mass
c           write(*,*)'vor linint'

            call linint1(rx, ry, rz, deriv)

c           write(*,*)'nach linint'

            j0 = deriv(0)
            j1 = deriv(1)
            j2 = deriv(2)
            j3 = deriv(3)
            sd = deriv(4)

            pboo(1) = pstore(1)
            pboo(2) = pstore(2)
            pboo(3) = pstore(3)

            pxstat = pboo(1)
            pystat = pboo(2)
            pzstat = pboo(3)
            scneu = upot(i)
            masse = e(i)

               write(*,*)'1kkk',masse,scneu,deriv

            flagit = .false.
            itcount = 0
            do while(.not.flagit)

              itcount = itcount + 1

c               write(*,*)'2kkk'
              pboo(1) = pxstat
              pboo(2) = pystat
              pboo(3) = pzstat
              pboo(0) = sqrt( (masse+scneu)**2 + pboo(1)**2 +
     &           pboo(2)**2 + pboo(3)**2)
              scold   = scneu

              call lorlrf(j0,j1,j2,j3,pboo,rho,flag )

*      now the vector pboo contains the energy and momenta
*      of the particle in the LRF of the cell k.
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

c               write(*,*)'3kkk'
              p1x = pboo(1)
              p1y = pboo(2)
              p1z = pboo(3)
              pabs= sqrt(p1x**2+p1y**2+p1z**2)
              pst = pabs

              plrf(1) = p1x
              plrf(2) = p1y
              plrf(3) = p1z

c                rho = sd
              pst       = pabs
              id1 = 1
c               write(*,*)'4kkk'
              vecpot    = potanal(rho0,rho,plrf,id1,rhap)

              ehelp = sqrt(masse**2 + pabs**2)
              vecpot = - masse + sqrt((ehelp+vecpot)**2-pabs**2)

c               write(*,*)'5kkk'
              pboo(0) = sqrt( (masse+vecpot)**2 + pabs**2)


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*             boost back to calc. frame                              *

              call lorlrf(j0,-j1,-j2,-j3,pboo,dummy,flag )
c               write(*,*)'6kkk'
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
c               write(*,*)'7'

              if(pxtest.le.epsi .and. pytest.le.epsi .and.
     &           pztest.le.epsi) then

                flagit = .true.
              end if

              if(itcount.ge.100) then
                write(*,*)'in iteration gradu 1d1 ', itcount
                write(*,*)'old ', pxstat, pystat, pzstat, masse,scold
                write(*,*)'neu ', pxneu, pyneu, pzneu, masse,scneu
cc                 stop
                goto 777
              end if

            end do
 777        continue
            write(*,*)'dens_4 scneu',scneu

            pstore(0) = pboo(0)
            meff      = e(i) + scneu
            upot(i)   = scneu

c               write(*,*)'9'
            do nst = 0,3
              veloc = pstore(nst)/pstore(0)
              rhob_4(nst,ix,iy,iz)=rhob_4(nst,ix,iy,iz)
     &           + float(sign(1,idtest))*veloc/float(num)
            end do
c               write(*,*)'9'


            if( (id(1,i).eq.1) .and. (id(2,i).eq.1)) then
              rhob_4(4,ix,iy,iz)=rhob_4(4,ix,iy,iz)+1./float(num)
            end if

            if( (id(1,i).eq.1) .and. (id(2,i).eq.0)) then
              rhob_4(5,ix,iy,iz)=rhob_4(5,ix,iy,iz)+1./float(num)
            end if

          end if
          end if
*
        end do

        rhob_4(6,jx,jy,jz)=1.0
        if( rhob_4(0,ix,iy,iz) .gt. 1.e-4) then
          rhob_4(6,ix,iy,iz)=1.0/sqrt(1.-( (rhob_4(1,ix,iy,iz)**2+
     &       rhob_4(2,ix,iy,iz)**2+rhob_4(3,ix,iy,iz)**2)/
     &       rhob_4(0,ix,iy,iz)**2)  )
        end if

c               write(*,*)'10'
        if(ismeasy.eq.1) then
          call smear
        end if

      else
        write(*,*)"poblems in dens_4: parameter ismear wrong!"
      end if

c      stop
c      if(ithree.eq.1) then
c        write(*,*)'vor sphere'
c        call sphere(ismear)
c        write(*,*)'nach sphere'
c      end if

*                                                                     *
*                                                                     *

c      tout = secnds(t0)
c     tin = tout - tin
c     write(*,*)'time ellapsed in dens_4 = ',tin,'  sec.'

      return
      end
************************************************************************
*                                                                      *
      function vpot(i1)
*                                                                      *
*       purpose:     calculation of nuclear potential V from spatial   *
*                    distribution of testparticles for momentum dep.   *
*                    forces for the baryons                            *
*                                                                      *
*       variables (all input, all integer)                             *
*         minnum  -  first testparticle treated in one run for density *
*         maxnum  -  last testparticle treated in one run for density  *
*         mass    -  total mass                                        *
*                                                                      *
*                                                                      *
************************************************************************

      implicit none
      include"common"
      include"cominput"
      integer ismear, ismeasy
      parameter(ismear  = 1)
      parameter(ismeasy = 0)

      integer i1, ix, iy, iz, ndim,   jx, jy, jz, nst
      integer kx, ky, kz, ic, ib, idtest, ixx, iyy, izz, i
      real*8    vpot, pstore(0:3), veloc,meff , ehelp

************************************************************************
      real*8    pboo(0:3), pxstat, pystat, pzstat
      real*8    j0, j1, j2, j3, sd
      real*8    pxneu, pyneu, pzneu, pst, pabs
      real*8    p1x, p1y, p1z, potanal, epsi, deriv(0:4)

      parameter(epsi = 1.0e-05)

      real*8  pxtest,pytest, pztest, dummy, masse, scneu, scold, vecpot
      logical flagit, flag
      integer itcount, id1
      real*8    rho, rx, ry, rz
      real*8    plrf(1:3), rhap(1:3)

************************************************************************
*
      real*8 tin, tout, t0
      t0 = 0.0
c      tin =secnds(t0)
*----------------------------------------------------------------------*
      vpot = 0.0
      ixx  = nint( r(1,i1) )
      iyy  = nint( r(2,i1) )
      izz  = nint( r(3,i1) )

*-----------------------------------------------------------------------
*
*      determination of densities for mom-dep forces
*       (implicit assumption that grid size = 1 fm )!!!!!!!!

      if(ismear .eq. 1) then
*
*       this part takes the smearing weights which are set in entry
*       densin of the subroutine dens (dens.f)
*       parameters are stored in common block PQ

      do i = 1, maxpar
        if(id(1,i).ne.0) then
        ix     = nint( r(1,i) )
        iy     = nint( r(2,i) )
        iz     = nint( r(3,i) )
        idtest = id(1,i)
*
       if(abs(ix-ixx).le.iq.and.abs(iy-iyy).le.iq.and.abs(iz-izz).le.iq)
     &    then
          kx=nint(float(2*ip+1)*(r(1,i1)-float(ix)))
          if(abs(kx) .eq. ip+1) kx = kx/abs(kx) * ip
          ky=nint(float(2*ip+1)*(r(2,i1)-float(iy)))
          if(abs(ky) .eq. ip+1) ky = ky/abs(ky) * ip
          kz=nint(float(2*ip+1)*(r(3,i1)-float(iz)))
          if(abs(kz) .eq. ip+1) kz = kz/abs(kz) * ip
          ic=1+(kz+ip)+(ky+ip)*(2*ip+1)+(kx+ip)*(2*ip+1)**2
          ib=0
          do jx=ix-iq,ix+iq
          do jy=iy-iq,iy+iq
          do jz=iz-iq,iz+iq
            ib=ib+1
            if(jx.eq.ixx.and.jy.eq.iyy.and.jz.eq.izz) then
              if(cm(ic,ib).gt.0.0.and.
     &          abs(jx).le.maxx.and.abs(jy).le.maxx.and.
     &          abs(jz).le.maxz) then
*
                pstore(1) = p(1,i)
                pstore(2) = p(2,i)
                pstore(3) = p(3,i)
                meff      = e(i) + upot(i)
                pstore(0) = sqrt(meff**2 + pstore(1)**2 +
     &                          pstore(2)**2 + pstore(3)**2)


                do nst = 0,3
                  veloc = pstore(nst)/pstore(0)
                  rhob_4(nst,jx,jy,jz)=rhob_4(nst,jx,jy,jz)
     &                       + float(sign(1,idtest))*cm(ic,ib)* veloc
                end do


              end if
            end if
*
          end do
          end do
          end do

*
        end if
        end if
      end do


      else if(ismear .eq. 0) then


      do i = 1,maxpar
*
        if(id(1,i).ne.0) then
        ix = nint( r(1,i) )
        iy = nint( r(2,i) )
        iz = nint( r(3,i) )
        idtest = id(1,i)
*
        if(abs(ix).lt.maxx.and.abs(iy).lt.maxx.and.abs(iz).lt.maxz)then

        rx     = r(1,i)
        ry     = r(2,i)
        rz     = r(3,i)

        rhap(1) = rx
        rhap(2) = ry
        rhap(3) = rz

           meff = e(i) + upot(i)
           pstore(1) = p(1,i)
           pstore(2) = p(2,i)
           pstore(3) = p(3,i)
           pstore(0) = sqrt(meff**2 + pstore(1)**2 +
     &                      pstore(2)**2 + pstore(3)**2)
*
           id1       = id(1,i)

****************************** iterate eff mass
c           write(*,*)'vor linint'

               call linint1(rx, ry, rz, deriv)

c           write(*,*)'nach linint'

               j0 = deriv(0)
               j1 = deriv(1)
               j2 = deriv(2)
               j3 = deriv(3)
               sd = deriv(4)

               pboo(1) = pstore(1)
               pboo(2) = pstore(2)
               pboo(3) = pstore(3)


               pxstat = pboo(1)
               pystat = pboo(2)
               pzstat = pboo(3)
               scneu = upot(i)
               masse = e(i)

c               write(*,*)'1kkk'

               flagit = .false.
               itcount = 0
               do while(.not.flagit)

                itcount = itcount + 1

c               write(*,*)'2kkk'
                pboo(1) = pxstat
                pboo(2) = pystat
                pboo(3) = pzstat
                pboo(0) = sqrt( (masse+scneu)**2 + pboo(1)**2 +
     &                         pboo(2)**2 + pboo(3)**2)
                scold   = scneu

               call lorlrf(j0,j1,j2,j3,pboo,rho,flag )

*      now the vector pboo contains the energy and momenta
*      of the particle in the LRF of the cell k.
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

c               write(*,*)'3kkk'
                 p1x = pboo(1)
                 p1y = pboo(2)
                 p1z = pboo(3)
                 pabs= sqrt(p1x**2+p1y**2+p1z**2)
                 pst = pabs

                 plrf(1) = p1x
                 plrf(2) = p1y
                 plrf(3) = p1z


c                rho = sd
                pst       = pabs
                id1 = 1
c               write(*,*)'4kkk'
                vecpot    = potanal(rho0,rho,plrf,id1,rhap)

                ehelp = sqrt(masse**2 + pabs**2)
                vecpot = - masse + sqrt((ehelp+vecpot)**2-pabs**2)

c               write(*,*)'5kkk'
               pboo(0) = sqrt( (masse+vecpot)**2 + pabs**2)


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*             boost back to calc. frame                              *

               call lorlrf(j0,-j1,-j2,-j3,pboo,dummy,flag )
c               write(*,*)'6kkk'
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
c               write(*,*)'7'

               if(pxtest.le.epsi .and. pytest.le.epsi .and.
     &            pztest.le.epsi) then

                 flagit = .true.
               end if

               if(itcount.ge.100) then
                 write(*,*)'in iteration gradu 1d2 ', itcount
                 write(*,*)'old ', pxstat, pystat, pzstat, masse,scold
                 write(*,*)'neu ', pxneu, pyneu, pzneu, masse,scneu
cc                 stop
                goto 778
               end if

              end do
  778 continue
c               write(*,*)'8'
            write(*,*)'dens_4 2 scneu',scneu

              pstore(0) = pboo(0)
              meff      = e(i) + scneu
              upot(i)   = scneu

c               write(*,*)'9'
           do nst = 0,3
             veloc = pstore(nst)/pstore(0)
             rhob_4(nst,ix,iy,iz)=rhob_4(nst,ix,iy,iz)
     &                        + float(sign(1,idtest))*veloc/float(num)
           end do
c               write(*,*)'9'


              if( (id(1,i).eq.1) .and. (id(2,i).eq.1)) then
                rhob_4(4,ix,iy,iz)=rhob_4(4,ix,iy,iz)+1./float(num)
              end if

              if( (id(1,i).eq.1) .and. (id(2,i).eq.0)) then
                rhob_4(5,ix,iy,iz)=rhob_4(5,ix,iy,iz)+1./float(num)
              end if

              rhob_4(6,jx,jy,jz)=1.0
              if( rhob_4(0,ix,iy,iz) .gt. 1.e-4) then
                rhob_4(6,ix,iy,iz)=1.0/sqrt(1.-( (rhob_4(1,ix,iy,iz)**2+
     &                 rhob_4(2,ix,iy,iz)**2+rhob_4(3,ix,iy,iz)**2)/
     &                 rhob_4(0,ix,iy,iz)**2)  )
              end if





        end if
*
        end if
      end do

c               write(*,*)'10'
        if(ismeasy.eq.1) then
          call smear
        end if

      else
        write(*,*)"poblems in dens_4: parameter ismear wrong!"
      end if


c      if(ithree.eq.1) then
c        write(*,*)'vor sphere'
c        call sphere(ismear)
c        write(*,*)'nach sphere'
c      end if

*                                                                     *
*                                                                     *

c      tout = secnds(t0)
c     tin = tout - tin
c     write(*,*)'time ellapsed in dens_4 = ',tin,'  sec.'

      return
      end
