      subroutine mesfin

*---------------------------------------------------------------------*
*     this routine propagtes the charged mesons in the Coulomb-field  *
*     that is left after the dynamical calculation is finished for    *
*     ntaft (parameter) timesteps.                                    *
*---------------------------------------------------------------------*
      implicit none
      include"common"
      include"coucom"
      include"cominput"
      include"com_kminu"

      integer ntaft, ndtfac
      parameter(ntaft = 5000)
      parameter(ndtfac = 5)

      integer i, k, ik, id2, nescc, nescpic, id1, idn
      real*8    rx, ry, rz, px, py, pz, etot, gradxr, gradyr, gradzr
      real*8    gradxp, gradyp, gradzp
      real*8    velox, veloy, veloz, mesmass, mepot, meff
      real*8    dilfac, dtprime, dummy1,dummy2,dummy3
      integer l,ipropag

      ipropag = 0

      if(ipou .eq. 1) then
        call cdens(nescc,nescpic)
        write(*,*)'in mesfin nach cdens'
*---------------------------------------------------------------------*
*     do the propagation of the mesons                                *
*     update positions and momenta                                    *
*
        dilfac = 1./float(ndtfac)
        dtprime = dt*dilfac
      do k = 1, ntaft
       do l = 1, ndtfac
       do  i = 1,maxppar
           if(ipi(1,i).ne.0 .and. ipi(2,i).ne.0) then


           rx =  rpi(1,i)
           ry =  rpi(2,i)
           rz =  rpi(3,i)

           id1=  ipi(1,i)
           id2=  ipi(2,i)
           idn=  i

           px =  ppi(1,i)
           py =  ppi(2,i)
           pz =  ppi(3,i)
           mesmass = epi(i)
           mepot   = mpot(i)
           meff    = mesmass + mepot
           etot = sqrt(meff**2+px**2+py**2+pz**2)


           call gradupi(rx, ry, rz, px, py, pz, etot, id1, id2, mesmass,
     &                gradxr,gradyr,gradzr, gradxp, gradyp,gradzp,
     &                dummy1,dummy2,dummy3,ipropag)


            ppi(1,i) = px - dtprime * gradxr
            ppi(2,i) = py - dtprime * gradyr
            ppi(3,i) = pz - dtprime * gradzr
            mepot   = mpot(i)
            meff    = mesmass + mepot
            etot = sqrt(meff**2+ppi(1,i)**2+ppi(2,i)**2+ppi(3,i)**2)
          velox    = ppi(1,i) / etot
          veloy    = ppi(2,i) / etot
          veloz    = ppi(3,i) / etot
          rpi(1,i) = rpi(1,i) + dtprime * velox
          rpi(2,i) = rpi(2,i) + dtprime * veloy
          rpi(3,i) = rpi(3,i) + dtprime * veloz

        end if
       end do
       end do
      end do            !     pions
c------------------------------------------  kaons
      if (ikaon .eq. 0)  goto 2100
      do k = 1, ntaft
       do l = 1, ndtfac
       do  2000  ik = 1, maxkaon
       if (ika(1,ik) .ne. 1) goto  2002

           rx =  rkao(1,ik)
           ry =  rkao(2,ik)
           rz =  rkao(3,ik)

           id1=  1
           id2=  1
           idn=  ik

           px =  pkao(1,ik)
           py =  pkao(2,ik)
           pz =  pkao(3,ik)
           mesmass =  xkmas
           meff    =  xkmas
           etot = sqrt(meff**2+px**2+py**2+pz**2)


           call gradupi(rx, ry, rz, px, py, pz, etot, id1, id2, mesmass,
     &                gradxr,gradyr,gradzr, gradxp, gradyp,gradzp,
     &                dummy1,dummy2,dummy3,ipropag)


            pkao(1,ik) = px - dtprime * gradxr
            pkao(2,ik) = py - dtprime * gradyr
            pkao(3,ik) = pz - dtprime * gradzr
            etot = sqrt(meff**2+
     1                  pkao(1,ik)**2+pkao(2,ik)**2+pkao(3,ik)**2)
          velox    = pkao(1,ik) / etot
          veloy    = pkao(2,ik) / etot
          veloz    = pkao(3,ik) / etot
          rkao(1,ik) = rkao(1,ik) + dtprime * velox
          rkao(2,ik) = rkao(2,ik) + dtprime * veloy
          rkao(3,ik) = rkao(3,ik) + dtprime * veloz

 2002 continue
 2000 continue   ! ik          end do
       end do
      end do
c
 2100 continue
c                     expansion         -        kminus
      if (i_kminu .eq. 0)  goto 3100
      do k = 1, ntaft
       do l = 1, ndtfac
       do  3000  ik = 1, max_kminu
       if (nx_kminu(0,ik) .ne. 1) goto  3002

           rx =  r_kminu(1,ik)
           ry =  r_kminu(2,ik)
           rz =  r_kminu(3,ik)

           id1=   1
           id2=  -1
           idn=  ik

           px =  p_kminu(1,ik)
           py =  p_kminu(2,ik)
           pz =  p_kminu(3,ik)
           mesmass =  xkmas
           meff    =  xkmas
           etot = sqrt(meff**2+px**2+py**2+pz**2)


           call gradupi(rx, ry, rz, px, py, pz, etot, id1, id2, mesmass,
     &                gradxr,gradyr,gradzr, gradxp, gradyp,gradzp,
     &                dummy1,dummy2,dummy3,ipropag)


            p_kminu(1,ik) = px - dtprime * gradxr
            p_kminu(2,ik) = py - dtprime * gradyr
            p_kminu(3,ik) = pz - dtprime * gradzr
            etot = sqrt(meff**2+
     1             p_kminu(1,ik)**2+p_kminu(2,ik)**2+p_kminu(3,ik)**2)
          velox    = p_kminu(1,ik) / etot
          veloy    = p_kminu(2,ik) / etot
          veloz    = p_kminu(3,ik) / etot
          r_kminu(1,ik) = r_kminu(1,ik) + dtprime * velox
          r_kminu(2,ik) = r_kminu(2,ik) + dtprime * veloy
          r_kminu(3,ik) = r_kminu(3,ik) + dtprime * veloz

 3002 continue
 3000 continue   ! ik          end do
       end do
      end do
c -----
 3100 continue
c-------
      end if
       write(*,*)'nach mesons in mesfin'

      return
      end
