      subroutine propa(nt)

*-------------------------------------------------------------------
*
*       this routine does the propagatin of both
*       the baryons
*       and
*       the mesons
*
*-------------------------------------------------------------------
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"
      common /nthwhw/  nthw
      integer nthw
      real*8  chatot, rchar(3)
      common /couleft/  chatot, rchar !  part from coucom
      common /counthw/ ihw
      integer ihw

*     multiplication factor for pion-propa
      integer ndtfac
      parameter(ndtfac = 4)

      integer i, id1, id2, i_pert, ii
      integer ix, iy, iz
      integer nescpic, nescc, nt, idt

      real*8 rx, ry, rz, px, py, pz, velox, veloy, veloz
      real*8  etot, masslim,vrel,dens,rself,iself,sgamma,rn,JPsim0
      real*8 gradpx, gradpy, gradpz, gradm
      real*8 gradx, grady, gradz, tmass
      real*8 gradrx, gradry, gradrz
      real*8 rk(1:3,1:maxpar)
      real*8 pk(1:3,1:maxpar), uk(1:maxpar)
      real*8 pk1(1:3,1:maxpar)
      real*8 gradrx1(1:maxpar), gradry1(1:maxpar), gradrz1(1:maxpar)
      integer ipavp, mass, nlost, idn
      real*8    ttc
      real*8  velox1(1:maxpar), veloy1(1:maxpar), veloz1(1:maxpar)
      real*8    meff, dtprime, usc, j0,j1,j2,j3
c      real*8    rho_old(0:2, 0:maxpar)             !            hw

      integer idtp,ipropag
      real*8    dtfac,ddt,newepi,mesmass
*     potential parameters

      real*8  tin, tout

*-------------------------------------------------------------------
*       determine the 4-densities and 4-current-densities needed
*       for the momentum depedent forces acting on the baryons and
*       for evaluating the forces due to coulomb
*


      tin = 0.0
c      tin = secnds(tin)

      ipavp = 0
      mass   = massta + masspr
      nlost  = 0
      ttc = 0.0

      write(*,*)'vor dens'
      call dens_4
      write(*,*)'after dens_4'
c      call f77flush()

      do iz  = -maxz, maxz                          !  hw
      do iy  = -maxx, maxx
      do ix  = -maxx, maxx
        rhob_4(7,ix,iy,iz) = rhob_4(1,ix,iy,iz)
        rhob_4(8,ix,iy,iz) = rhob_4(2,ix,iy,iz)
        rhob_4(9,ix,iy,iz) = rhob_4(3,ix,iy,iz)
        rhob_4(10,ix,iy,iz) = rhob_4(0,ix,iy,iz)
c       write(*,*) ' propa bary 10 ', rhob_4(1,ix,iy,iz), ix,iy,iz
      enddo
      enddo
      enddo

************************************************************************

c        write(*,*)'  propa  1'
      if(isplipi.eq.1) then
        call spline3di1
        call spline3di
      end if
c        write(*,*)' propa 2'

      if(ipou .eq. 1 ) then
        call initcoul
        call cdens(nescc,nescpic)
c        if(iret .eq. 0) then
          call pois(ncont)
c        else if(iret .eq. 1) then
c          call reta1(nt, dt, ncont)
c        end if
      end if
c       write(*,*) ' upot(30/123)3', upot(30), upot(123), p(1,30),
c    1              r(1,30), p(2,123), r(2,123)
c        write(*,*)'  propa 3'
*
*-------------------------------------------------------------------
*        do baryon propagation
*

*     update momenta

      do i = 1, maxpar
        if(id(1,i).ne.0) then
          rk(1,i) = r(1,i)
          rk(2,i) = r(2,i)
          rk(3,i) = r(3,i)
          pk(1,i) = p(1,i)
          pk(2,i) = p(2,i)
          pk(3,i) = p(3,i)
          uk(i)   = upot(i)
        end if
      end do



***************************************************************
*     predictor step ******************************************
***************************************************************

      do i = 1,maxpar
        if(id(1,i).ne.0) then
        rx   =  rk(1,i)
        ry   =  rk(2,i)
        rz   =  rk(3,i)

        px   =  pk(1,i)
        py   =  pk(2,i)
        pz   =  pk(3,i)


        id1  =  id(1,i)
        id2  =  id(2,i)

        ix = nint(rx)
        iy = nint(ry)
        iz = nint(rz)

        idtp = i

        call gradu(rx, ry, rz, px, py, pz, etot, id1, id2,idtp,
     &             gradrx, gradry, gradrz,
     &             gradpx,gradpy,gradpz,nt,1)

        gradrx1(i) = gradrx
        gradry1(i) = gradry
        gradrz1(i) = gradrz

        p(1,i) = px - dt * gradrx1(i)
        p(2,i) = py - dt * gradry1(i)
        p(3,i) = pz - dt * gradrz1(i)
*
        pk1(1,i) = p(1,i)
        pk1(2,i) = p(2,i)
        pk1(3,i) = p(3,i)

        end if
      end do

***************************************************************
*     predictor step ******************************************
***************************************************************

      do i = 1,maxpar
        if(id(1,i).ne.0) then
        rx   =  rk(1,i)
        ry   =  rk(2,i)
        rz   =  rk(3,i)

        px   =  pk(1,i)
        py   =  pk(2,i)
        pz   =  pk(3,i)
        usc  =  uk(i)

        id1  =  id(1,i)
        id2  =  id(2,i)

        ix = nint(rx)
        iy = nint(ry)
        iz = nint(rz)
        idtp = i

        call gradu(rx, ry, rz, px, py, pz, etot, id1, id2,idtp,
     &             gradrx, gradry, gradrz,
     &             gradpx,gradpy,gradpz,nt,2)

*
*     update positions
*
        velox1(i) =  gradpx
        veloy1(i) =  gradpy
        veloz1(i) =  gradpz
*
        r(1,i) = rx + dt * velox1(i)
        r(2,i) = ry + dt * veloy1(i)
        r(3,i) = rz + dt * veloz1(i)
*
        end if
      end do

c       write(*,*)  '   propa 4  end of predictor '
**************** end of predictor step ******************************

      call dens_4

      if(isplipi.eq.1) then
        call spline3di1
        call spline3di
      end if

**************** do corrector step ***********************************

      do i = 1,maxpar
        if(id(1,i).ne.0) then

          id1  =  id(1,i)
          id2  =  id(2,i)
          idtp = i

          rx   =  r(1,i)
          ry   =  r(2,i)
          rz   =  r(3,i)

          px   =  p(1,i)
          py   =  p(2,i)
          pz   =  p(3,i)
          usc  =  upot(i)

c        write(*,*)  '  propa  5  - corrector '
          call gradu(rx, ry, rz, px, py, pz, etot, id1, id2,idtp,
     &             gradrx, gradry, gradrz,
     &             gradpx,gradpy,gradpz,nt,1 )

          p(1,i) = pk(1,i) - dt * (gradrx1(i) + gradrx)*0.5
          p(2,i) = pk(2,i) - dt * (gradry1(i) + gradry)*0.5
          p(3,i) = pk(3,i) - dt * (gradrz1(i) + gradrz)*0.5
        end if
      end do

      do i = 1,maxpar
        if(id(1,i).ne.0) then

          id1  =  id(1,i)
          id2  =  id(2,i)
          idtp = i

          rx   =  r(1,i)
          ry   =  r(2,i)
          rz   =  r(3,i)

          px   =  pk1(1,i)
          py   =  pk1(2,i)
          pz   =  pk1(3,i)

          call gradu(rx, ry, rz, px, py, pz, etot, id1, id2,idtp,
     &             gradrx, gradry, gradrz,
     &             gradpx,gradpy,gradpz,nt,2 )

          velox = gradpx
          veloy = gradpy
          veloz = gradpz

*
*     update positions
*
          r(1,i) = rk(1,i) + dt * 0.5*(velox1(i)+velox)
          r(2,i) = rk(2,i) + dt * 0.5*(veloy1(i)+veloy)
          r(3,i) = rk(3,i) + dt * 0.5*(veloz1(i)+veloz)
*
        end if
      end do

c 10   continue

*
      write(*,*)'nach baryon propagation'
      call dens_4
      write(*,*)'after  dens4 inside propa.f '
      call f77flush()



*       end of baryon propagation
*-------------------------------------------------------------------
*
*       do meson propagation

*-------------------------------------------------------------------
*
*  ipi(1,i)=particle type (pion, eta, ...)
*  ipi(2,i)=particle charge
*  ipi(3,i)=the serial number for the parent resonance of meson i
*  ipi(4,i)=how many times the meson was created + for pions - for e
*  ipi(5,i)=type of the parent resonance (id(1,ipi(i,3))))
*
*-------------------------------------------------------------------
*     update positions and momenta
*
      dtfac = 1.0/float(ndtfac)
      dtprime = dt * dtfac
      do idt = 1,ndtfac
        write(*,*) ' in propa mesons ',idt,ndtfac
        ddt = idt*dtprime/dt

       do  i = 1,maxppar
        if(ipi(1,i) .ne. 0) then
c       write(*,*) ' propa: cycle for mesons ',idt, i, ipi(1,i)
c       call f77flush()
           rx  =  rpi(1,i)
           ry  =  rpi(2,i)
           rz  =  rpi(3,i)

           id1 =  ipi(1,i)
           id2 =  ipi(2,i)

           px  =  ppi(1,i)
           py  =  ppi(2,i)
           pz  =  ppi(3,i)

           meff =  epi(i) + mpot(i)
           mesmass = epi(i)
           etot = sqrt(meff**2+px**2+py**2+pz**2)
           idn  = i

c          write(*,*) 'before gradupi: ',nthw,i,idt, id1, id2,idn,
c    &        rx, ry, rz,
c    &        px, py, pz, etot, epi(idn),
c    &        ddt,dtprime
c       call f77flush()
c---------
           ipropag = 1
           call gradupi(rx, ry, rz, px, py, pz, etot, id1, id2,mesmass,
     &              gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
     &              ddt,dtprime,ipropag)
c-------------------------------------
c          if (id1 .eq. 3)
c    1     write(*,*) 'after gradupi: ',nthw,i, idt,
c    &        rx, ry, rz,
c    &        px, py, pz, etot, id1, id2,idn,
c    &        gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
c    &        ddt,dtprime
c-------------------------------------------

           meff     = epi(i) + mpot(i)
           etot     = sqrt(meff**2+ppi(1,i)**2+ppi(2,i)**2+ppi(3,i)**2)
c
           if (id1.eq.3 .or. id1.eq.5) then
c -------------- no limit  hw-----------------
             masslim = float(id1+1) / 2.
c           that's because of the new spectral funtion : no masslimit
             if (id1.eq.3 .and. icbro .gt. 1)  masslim = .04/pmass
c -------------- no limit  hw-----------------

             newepi = epi(i) + dtprime * gradm * etot/meff  ! problem with
                                                            ! mpot(i)  ??!
             if ( newepi .gt. masslim*pmass+.01
     1          .and. newepi .lt. 1.5
     2          .and. abs(newepi-epi(i)).lt. 1.2)  then
               ppi(1,i) = px - dtprime * gradrx
               ppi(2,i) = py - dtprime * gradry
               ppi(3,i) = pz - dtprime * gradrz
               etot = sqrt(meff**2+ppi(1,i)**2+ppi(2,i)**2+ppi(3,i)**2)
               epi(i) = newepi
c             if (id1.eq.5) write(52,*) ' mass ',id1, i,gradm,newepi
             else
               write(*,*) 'propa: too small rho mass',i,id1,newepi,
     1            epi(i),dtprime,gradm,etot,
     2            meff,ppi(1,i),ppi(2,i),ppi(3,i)
               if (newepi .gt. 1.2) epi(i) =  1.23456 ! masslim*pmass + .01
               meff     = epi(i) + mpot(i)
               ppi(1,i) = px
               ppi(2,i) = py
               ppi(3,i) = pz
               etot = sqrt(meff**2+ppi(1,i)**2+ppi(2,i)**2+ppi(3,i)**2)
               gradpx = px/etot
               gradpy = py/etot
               gradpz = pz/etot
c                     what else can be done ???  HW
c             write(*,*) 'mass change: ',id1,i,nthw,
c    1           sqrt(ppi(1,i)**2+ppi(2,i)**2+ppi(3,i)**2),rz,epi(i)
             end if

             velox = gradpx
             veloy = gradpy
             veloz = gradpz
           else

             ppi(1,i) = px - dtprime * gradrx
             ppi(2,i) = py - dtprime * gradry
             ppi(3,i) = pz - dtprime * gradrz
             etot = sqrt(meff**2+ppi(1,i)**2+ppi(2,i)**2+ppi(3,i)**2)

             gradpx   = gradpx*(meff/etot)
             gradpy   = gradpy*(meff/etot)
             gradpz   = gradpz*(meff/etot)

             velox    = ppi(1,i) / etot + gradpx
             veloy    = ppi(2,i) / etot + gradpy
             veloz    = ppi(3,i) / etot + gradpz
           end if

           rpi(1,i) = rpi(1,i) + dtprime * velox
           rpi(2,i) = rpi(2,i) + dtprime * veloy
           rpi(3,i) = rpi(3,i) + dtprime * veloz

        end if
       end do
      end do

       write(*,*)'nach mesons in propa'
*-----------------------------------------------------------------------
*
*       hyperon propagation
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*     update positions and momenta
*
      do  i = 1,max_kminu
        if(nx_hyp(0,i) .eq. 1) then

            rx =  r_hyp(1,i)
            ry =  r_hyp(2,i)
            rz =  r_hyp(3,i)

c           id2=  1

            px =  p_hyp(1,i)
            py =  p_hyp(2,i)
            pz =  p_hyp(3,i)
            tmass = p_hyp(0,i)
            etot = sqrt(tmass**2+px**2+py**2+pz**2)
c           call gradu_hyp(nt, rx, ry, rz, px, py, pz, etot, id2,
c    &                  gradx, grady, gradz, velox, veloy, veloz)


           velox    = p_hyp(1,i) / etot                 !   hw
           veloy    = p_hyp(2,i) / etot
           veloz    = p_hyp(3,i) / etot

          r_hyp(1,i) = rx + dt * velox
          r_hyp(2,i) = ry + dt * veloy
          r_hyp(3,i) = rz + dt * veloz

        end if
      end do
*-----------------------------------------------------------------------
*
*       do kaon propagation
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*     update positions and momenta
*
      do  i = 1,maxkaon
        if(ika(1,i) .ne. 0) then

            rx =  rkao(1,i)
            ry =  rkao(2,i)
            rz =  rkao(3,i)

            id2= ika(1,i)

            px =  pkao(1,i)
            py =  pkao(2,i)
            pz =  pkao(3,i)
            tmass = xkmas
            etot = sqrt(tmass**2+px**2+py**2+pz**2)
          if (ikaonpot .gt. 0)  then
c         write(*,*)   ' propa, vor gradu kaon ',i
            call gradukao2(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                  gradx, grady, gradz, velox, veloy, veloz)

c            if ( i.le.5 ) then
c             write(*,7712)   i, nt, pkao(4,i),px, py, pz, gradz,
c    1        chatot, rchar, rx,ry,rz
c7712         format ( ' propa kaon i t ',2i4, e12.4, 3f8.4, e12.4,
c    1                    f8.1, ' r ',6f6.1)
c            endif
            pkao(1,i) = px - dt * gradx
            pkao(2,i) = py - dt * grady
            pkao(3,i) = pz - dt * gradz
          else
            velox    = pkao(1,i) / etot                 !   hw
            veloy    = pkao(2,i) / etot
            veloz    = pkao(3,i) / etot
          endif
            rkao(1,i) = rkao(1,i) + dt * velox
            rkao(2,i) = rkao(2,i) + dt * veloy
            rkao(3,i) = rkao(3,i) + dt * veloz
        end if  !  if kaon
      end do
*-----------------------------------------------------------------------
*
*       do kminus propagation
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*     update positions and momenta
*
      do  i = 1,max_kminu
        if(nx_kminu(0,i) .ne. 0) then

            rx =  r_kminu(1,i)
            ry =  r_kminu(2,i)
            rz =  r_kminu(3,i)

            id2=  -nx_kminu(0,i)

            px =  p_kminu(1,i)
            py =  p_kminu(2,i)
            pz =  p_kminu(3,i)
            tmass = xkmas
            etot = sqrt(tmass**2+px**2+py**2+pz**2)
c - - - - -
          if (i_kminu_pot .gt. 0)  then
          ihw = i
c         if (i .eq. 3197)
c    1    write(*,*)   ' propa, vor gradu kminus ',i, px, py,pz
            call gradukao2(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                  gradx, grady, gradz, velox, veloy, veloz)

c      if ( i.lt.4 ) then
c         write(*,7712)   i, nt, px, gradx, r_kminu(1,i), velox
c7712  format ( ' i t ',2i4, 4e12.4)
c      endif
             p_kminu(1,i) = px - dt * gradx
             p_kminu(2,i) = py - dt * grady
             p_kminu(3,i) = pz - dt * gradz
          else
             velox    = p_kminu(1,i) / etot                 !   hw
             veloy    = p_kminu(2,i) / etot
             veloz    = p_kminu(3,i) / etot
          endif
          r_kminu(1,i) = r_kminu(1,i) + dt * velox
          r_kminu(2,i) = r_kminu(2,i) + dt * veloy
          r_kminu(3,i) = r_kminu(3,i) + dt * veloz
        end if                ! end kminus
      end do

*-----------------------------------------------------------------------
*
*       do phi propagation (zm)
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*     update positions and momenta
*
      do  i = 1,max_pert
        if(nx_pert(id_phi,0,i) .ne. 0) then

            rx =  r_pert(id_phi,1,i)
            ry =  r_pert(id_phi,2,i)
            rz =  r_pert(id_phi,3,i)

c            id2=  0

            px =  p_pert(id_phi,1,i)
            py =  p_pert(id_phi,2,i)
            pz =  p_pert(id_phi,3,i)
            tmass = p_pert(id_phi,0,i)
            etot = sqrt(tmass**2+px**2+py**2+pz**2)
            call graduphi(nt, rx, ry, rz, px, py, pz, etot,
     &                  gradx, grady, gradz, velox, veloy, veloz)

c            velox    = px / etot
c            veloy    = py / etot
c            veloz    = pz / etot

            p_pert(id_phi,0,i) = tmass
            p_pert(id_phi,1,i) = px - dt * gradx
            p_pert(id_phi,2,i) = py - dt * grady
            p_pert(id_phi,3,i) = pz - dt * gradz

            r_pert(id_phi,1,i) = r_pert(id_phi,1,i) + dt * velox
            r_pert(id_phi,2,i) = r_pert(id_phi,2,i) + dt * veloy
            r_pert(id_phi,3,i) = r_pert(id_phi,3,i) + dt * veloz
          if(p_pert(id_phi,0,i).lt.0.05) then
            write(*,*) 'propa phimass',p_pert(id_phi,0,i)
            stop
          end if

        end if
      end do

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*       do J/Psi propagation
*
*     update positions and momenta
*
*-------------------------------------------------------------------
*
*  ipi(1,i)=particle type (pion, eta, ...)
*  ipi(2,i)=particle charge
*  ipi(3,i)=the serial number for the parent resonance of meson i
*  ipi(4,i)=how many times the meson was created + for pions - for e
*  ipi(5,i)=type of the parent resonance (id(1,ipi(i,3))))
*
*-------------------------------------------------------------------
*     update positions and momenta
*
      if(i_JPsi.ge.5 .and. i_JPsi.le.8) then
        dtfac = 1.0/float(ndtfac)
        dtprime = dt * dtfac
        write(*,*) ' in propa JPsi ',ndtfac
        do idt = 1,ndtfac
          ddt = idt*dtprime/dt

          do  i_pert = 1,max_pert
            do ii =1,3
              if(nx_pert(id_JPsi(ii),0,i_pert) .ne. 0) then
c       write(*,*) ' propa: cycle for mesons ',idt, i, ipi(1,i)
c       call f77flush()
                id1  = 100+id_JPsi(ii)
                id2 =  0
                rx = r_pert(id_JPsi(ii),1,i_pert)
                ry = r_pert(id_JPsi(ii),2,i_pert)
                rz = r_pert(id_JPsi(ii),3,i_pert)

                px = p_pert(id_JPsi(ii),1,i_pert)
                py = p_pert(id_JPsi(ii),2,i_pert)
                pz = p_pert(id_JPsi(ii),3,i_pert)
                mesmass = p_pert(id_JPsi(ii),0,i_pert)
                if(abs(mesmass-JPsi_prop(ii,1)).gt.0.5) write(*,*)
     &            'propa1 JPsimass',mesmass,JPsi_prop(ii,1),ii,i_pert
     &             ,px,py,pz,p_pert(id_JPsi(ii),4,i_pert),
     &            nx_pert(id_JPsi(ii),0,i_pert),
     &            nx_pert(id_JPsi(ii),2,i_pert)
                meff =  mesmass
                etot = sqrt(meff**2+px**2+py**2+pz**2)


c          write(*,*) 'JPsi before gradupi: ',id1,i_pert,mesmass,etot,
c     &        rx, ry, rz,
c     &        px, py, pz
                call f77flush()
c---------
                ipropag = 1
                call gradupi(rx,ry,rz, px, py, pz, etot,id1,id2,mesmass,
     &              gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
     &              ddt,dtprime,ipropag)

                newepi = mesmass + dtprime * gradm * etot/meff
c                if(abs(newepi-JPsi_prop(ii,1)).gt.0.8 .and. idt.eq.1)
c     &      write(*,*)'propa3 JPsimass',newepi,JPsi_prop(ii,1),ii,i_pert
              
                ix = nint(r_pert(id_JPsi(ii),1,i_pert))
                iy = nint(r_pert(id_JPsi(ii),2,i_pert))
                iz = nint(r_pert(id_JPsi(ii),3,i_pert))
                dens = 0.0
                j0 = 0.0
                if(iabs(ix).le.maxx.and.iabs(iy).le.maxx.and.
     &             iabs(iz).le.maxz)then
                  j0=rhob_4(0,ix,iy,iz)  
                  j1=rhob_4(1,ix,iy,iz)  
                  j2=rhob_4(2,ix,iy,iz)  
                  j3=rhob_4(3,ix,iy,iz)  
                  dens = sqrt(j0**2-j1**2-j2**2-j3**2)
c                  write(*,*) 'JPpropaens:', j0,dens,rhb(ix,iy,iz)
                  dens = rhb(ix,iy,iz)
                end if
                if (dens.lt.1.0e-3) then
                  dens=0.
                  vrel = 0.
                else
                  call f77flush()
                  if(j1**2+j2**2+j3**2.gt.j0**2) then
                    write(*,*) "hiba propa, mass<0",j0,j1,j2,j3
                    stop
                  end if
                  call lorentz(px/etot,py/Etot,pz/Etot,
     &                j1,j2,j3,j0)
                  vrel = sqrt(j1**2+j2**2+j3**2)/j0
                end if

                call self_JPsi(ii,vrel,dens,rself,iself,sgamma)
                JPsim0= JPsi_prop(ii,1)+rself/(2.0*JPsi_prop(ii,1))
                if (abs(newepi - JPsi_prop(ii,1)) .lt. 1.0) then
                  p_pert(id_JPsi(ii),0,i_pert) = newepi
                else
  20              meff = JPsim0 - (rn(iseed)-0.5)*iself/JPsi_prop(ii,1)
                  if(rn(iseed).gt. iself**2
     &                /((meff**2-JPsim0**2)**2+iself**2))      goto 20
c                  write(*,*)'propa: too small JPsi m',i_pert,id1,newepi,
c     1          p_pert(id_JPsi(ii),0,i_pert),meff,gradm,etot,iself,rself

                  p_pert(id_JPsi(ii),0,i_pert) = meff
                end if
                p_pert(id_JPsi(ii),1,i_pert) = px - dtprime * gradrx
                p_pert(id_JPsi(ii),2,i_pert) = py - dtprime * gradry
                p_pert(id_JPsi(ii),3,i_pert) = pz - dtprime * gradrz
c        if(i_pert.eq.4012) write(*,*) 'JP_4012propa1 p',
c     &   p_pert(id_JPsi(ii),1,i_pert),p_pert(id_JPsi(ii),2,i_pert),
c     &   p_pert(id_JPsi(ii),3,i_pert),px,py,pz,gradrx,gradry,gradrz
                etot = sqrt(p_pert(id_JPsi(ii),0,i_pert)**2
     &                +p_pert(id_JPsi(ii),1,i_pert)**2
     &                +p_pert(id_JPsi(ii),2,i_pert)**2
     &                +p_pert(id_JPsi(ii),3,i_pert)**2)

                velox = gradpx
                veloy = gradpy
                veloz = gradpz
                r_pert(id_JPsi(ii),1,i_pert) =
     &            r_pert(id_JPsi(ii),1,i_pert) + dtprime * velox
                r_pert(id_JPsi(ii),2,i_pert) =
     &            r_pert(id_JPsi(ii),2,i_pert) + dtprime * veloy
                r_pert(id_JPsi(ii),3,i_pert) =
     &            r_pert(id_JPsi(ii),3,i_pert) + dtprime * veloz
              end if
            end do
          end do
        end do
      end if
      write(*,*) 'end JPsi propa'
cc          call dens_4         !  used after baryon propagation hw
c      call thermodyn

c      write(*,*)'vor thermo2 '
c      call thermo2
c      write(*,*)'nach thremo2'
      tout = 0.0
c      tout = secnds(tout)
c     tout = tout - tin
c     write(*,*)'time ellapsed in propa : ',tout,' sec '
c     write(*,*) ' NACH PROPA - ', epi(6003),ppi(1,6003),
c    1                             ppi(2,6003), ppi(3,6003)


      return
      end

