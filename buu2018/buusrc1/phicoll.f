************************************************************************
*                                                                      *
      subroutine phicoll(isu,wref,nt)
*                                                                      *
*     purpose: phi baryon elastic and absorption scattering  (zm)      *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"

      real*8 ppel, ppp,wref,plab,p0,pzlab,thdeg,xtrav
      real*8 s,pcm,xx,yy,zz,r2,rr,betax,betay,betaz,gamma,pphiab
      real*8 x1, y1, z1, px1, py1, pz1, em1, em12, e1,e1cm,p1beta,transf
      real*8 dx, dy, dz, rsqare, px2, py2, pz2, em2, em22, e2, x2, y2,z2
      real*8 p12, p1dr, p2dr, b12, a12, c12, brel, t1, t2, ddlt, b21,dxm
      real*8 rn,phase,p2beta,e2cm,qx1,qy1,qz1,qx2,qy2,qz2,dxp,denst
      real*8 xkaolam, xradkl, alp, eps_abs, rad_abs,radmax,sigmax,decide
      real*8 a3,a2,ax0,adx,distr,pirphi,ttt
      integer ntag, numcoll,isu,nt
      integer collphi(2*maxkaon),collphid(2*maxkaon),u1,u2,u3,u4,u5
      integer maxk,irun,ink,ini,i1,i2,ii,ixn,iyn,izn
c      integer iz,jj,numabss,ix,iy

c      parameter (pirphi   = 0.134) ! -> 1000 mb (test)
c      parameter (pirphi   = 0.977)      !-> 30 mb
c      parameter (pirphi   = 0.399)  !-> 5 mb
c      parameter (pirphi   = 0.560)  !-> 10 mb
c      parameter (pirphi   = 0.736)  !-> 17 mb

c      parameter (pirphi = 0.6)   !elast phiN: -> 11.3 mb (20fache)
c      parameter (pirphi = 0.736)   !elast phiN: -> 5600.0 mb (100fache)

!      parameter (xkaolam = 1.611, xradkl = 0.425, alp = 1.66)
      parameter (xkaolam = 1.611, xradkl = 0.3989, alp = 1.66)

c      save numabss
      save u1,u2,u3,u4,u5,numcoll,ttt
*----------------------------------------------------------------------*
*    pirphi = 0.134 fm                corresponds  to 0.56 mb          *
*    xkaolam = kaon+lambda, xrad :: absorption-radius,                 *
*    alp::   sigma = pi*xkl^2(5.7mb)* exp(-2*alp*(sqrt(s)-threshold)   *
*----------------------------------------------------------------------*

!      return


      pirphi = 0.564       !manuell auf 10mb gesetzt

      write(*,*) 'call phicoll'
c      crosphi = 10.0 * pi * pirphi**2 ! not used
      maxk = max_pert/num
      kanum = 0
!      numcoll = 0
!      ttt = 0
      u1 = 0
      u2 = 0
      u3 = 0
      u4 = 0
      u5 = 0
!       numabss = 0
*     loop over all parallel runs
      do 1000 irun = 1,num
c        write(*,*) "phicoll, irun"
        ink = (irun-1) * maxk
        ini = (irun-1) * maxb
        do 800 ii  = 1,maxk
          i1  = ii + ink
          if(nx_pert(id_phi,0,i1) .eq. 0 )                     goto 800

          ttt = ttt + 1

          kanum = kanum + 1
          x1  = r_pert(id_phi,1,i1)
          y1  = r_pert(id_phi,2,i1)
          z1  = r_pert(id_phi,3,i1)
          px1 = p_pert(id_phi,1,i1)
          py1 = p_pert(id_phi,2,i1)
          pz1 = p_pert(id_phi,3,i1)
          em1 = p_pert(id_phi,0,i1)
          em12= em1**2
          e1  = sqrt( em12 + px1**2 + py1**2 + pz1**2 )
          if(em1.lt.0.1) then
             write(*,*) "phicoll, maxk",i1,em1,nx_pert(id_phi,0,i1),
     &            nx_pert(id_phi,1,i1), nx_pert(id_phi,2,i1)
          end if
*     look for a scattering pseudonucleon in the same run
          i2  = ini
  600     i2  = i2 + 1
          if(i2 .gt. ini+maxb)                                 goto 800
          if(id(1,i2).le.0)                                    goto 600
          if(i2.eq.nx_pert(id_phi,3,i1).or.i2.eq.nx_pert(id_phi,4,i1))
     &                                                         goto 600
          x2 = r(1,i2)
          dx     = x1 - x2
            if(nbound.eq.1) then
              dxp = dx+2.0*boxx
              dxm = dx-2.0*boxx
              if(abs(dx) .gt. abs(dxp)) dx=dxp
              if(abs(dx) .gt. abs(dxm)) dx=dxm
            end if
          if (abs(dx) .gt. 1.0*delpi)                          goto 600
          y2 = r(2,i2)
          dy     = y1 - y2
            if(nbound.eq.1) then
              dxp = dy+2.0*boxx
              dxm = dy-2.0*boxx
              if(abs(dy) .gt. abs(dxp)) dy=dxp
              if(abs(dy) .gt. abs(dxm)) dy=dxm
            end if
          if (abs(dy) .gt. 1.0*delpi)                          goto 600
          z2 = r(3,i2)
          dz     = z1 - z2
            if(nbound.eq.1) then
              dxp = dz+2.0*boxz
              dxm = dz-2.0*boxz
              if(abs(dz) .gt. abs(dxp)) dz=dxp
              if(abs(dz) .gt. abs(dxm)) dz=dxm
            end if
          if (abs(dz) .gt. 1.0*delpi)                          goto 600
          rsqare = dx**2 + dy**2 + dz**2
          if (rsqare .gt. 1.0*delpi**2)                        goto 600
*         now particles are close enough to each other !

c          write(*,*) "phicoll, particles are close"
          px2    = p(1,i2)
          py2    = p(2,i2)
          pz2    = p(3,i2)
          em2    = e(i2)
          em22   = e(i2)**2
          e2     = sqrt ( em22 + px2**2 + py2**2 + pz2**2 )
          s      = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                           - (pz1+pz2)**2
          eps_abs = sqrt(s) - em1 - em2

          u1 = u1 + 1
          if (eps_abs .gt. .0) then
!            rad_abs = xradkl * exp(-alp*eps_abs)
!            rad_abs = xradkl 					!ANKE fester absorptionsquerschnitt von 5mb
            pphiab = sqrt(px1**2+py1**2+pz1**2)

!            pirphi = sqrt((10./(1.+pphiab))/(pi*10.))           !NPA625,832
!           rad_abs = sqrt((15.*pphiab + 1.)/(pi*10.))        !from BUU extracted:f(p)=10mb/GeV/c + 4mb
!          rad_abs=sqrt((-25.*pphiabs**2+80.*pphiabs-39.)/(pi*10.))


            a3  = 0.453
            a2  = 26.32
            adx = 0.198
            ax0 = 0.948

c            Write(*,*) 'phicoll iso elott',u1
*********************** isospin asymmetry dependence ********************
            if     (id(2,i2) .eq. 0)  then          !neutron
              a2  = 26.32
              adx = 0.26
              ax0 = 0.91
            elseif (id(2,i2) .eq. 1)  then          !proton
              a2  = 25.00
              adx = 0.11
              ax0 = 0.98

            endif
*****************************************************************************

          rad_abs=sqrt((a2+(a3-a2)/(1.+exp((pphiab-ax0)/adx)))/(pi*10.))
c            Write(*,*) 'phicoll rad utan',rad_abs,em1,em2

 !            if(sqrt(s).gt.1.958 .and.sqrt(s).lt.2.086) then
 !               rad_abs = 0.
 !               else
 !                  rad_abs = 100000.
 !            endif
!        rad_abs = sqrt(17.0/(pi*10.))            !17mb = const


          else
            rad_abs = 0.0
          endif

          if (id(1,i2) .ne. 1) rad_abs = 0.0    !only nucleons

          u2 = u2 + 1
          sigmax = rad_abs**2 + pirphi**2      !ANKE Faktor * abs
          radmax = sqrt(sigmax)

************************** isospin dependent absorption( kapt/kampf) *************
!          if     (id(1,i2).eq.1 .and. id(2,i2).eq.0) then       !phi+n
!            sigmax = 1.0*(4.0)*rad_abs**2 + pirphi**2
!          elseif (id(1,i2).eq.1 .and. id(2,i2).eq.1) then       !phi+p
!            sigmax = 0.25*(4.0)*rad_abs**2 + pirphi**2
!          else                                                  !phi+X
!            sigmax = 1.0*(4.0)*rad_abs**2 + pirphi**2
!          endif
**********************************************************************************

*   is their impact parameter small enough?
          p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          p1dr   = px1 * dx + py1 * dy + pz1 * dz
          p2dr   = px2 * dx + py2 * dy + pz2 * dz
          a12    = 1.0 - ( em1 * em2 / p12 ) ** 2
          b12    = p1dr / em1 - p2dr * em1 / p12
          c12    = rsqare + ( p1dr / em1 )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
c_hw          if (brel .gt. pirphi)                            goto 600
          if (brel .gt. radmax)                                goto 600
          u3 = u3 +1
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / em1 - b12 / a12 ) * e1 / em1
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
          ddlt   = t1 - t2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt )                            goto 600
*   now  the phi will collide or be absorbed in this time step

c          write(*,*) "phicoll, collide"
           colphi(i1) = colphi(i1) + 1

c           density dependence
!           ixn = nint(x1)
!           iyn = nint(y1)
!           izn = nint(z1)
!       if(iabs(ixn).le.maxx.and.iabs(iyn).le.maxx.and.iabs(izn).le.maxz)
!      &      denst = rhb(ixn,iyn,izn) / rho0
!         collphi(i1)    = int(time*2.)
!         collphid(i1)   = int(denst*20.)
************************************************************************
c_hw----------------      absorption   --------
              decide = (pirphi/radmax)**2
             if(rn(iseed).gt. decide)  then
                   p_pert(id_phi,4,i1) = .0
                   nx_pert(id_phi,0,i1) = 0

c         write(*,*) "phicoll abs  :",pi*sigmax*10.
c           density dependence
          ixn = nint(x1)
          iyn = nint(y1)
          izn = nint(z1)
      if(iabs(ixn).le.maxx.and.iabs(iyn).le.maxx.and.iabs(izn).le.maxz)
     &      denst = rhb(ixn,iyn,izn) / rho0
          collphi(i1)    = int(time*2.)
          collphid(i1)   = int(denst*20.)
         distr = sqrt(x1**2 + y1**2 + z1**2)
c       write(*,*)'phicoll ANKE:',collphi(i1)/2.,collphid(i1)/20.

             goto 600
             endif

c_hw----------------  end of absorption   --------
          u4 = u4 + 1
          pcm = sqrt(0.25*(s-em1**2+em2**2)**2/s -em2**2)

  100     xx = 1.0 - 2.*rn(iseed)
          yy = 1.0 - 2.*rn(iseed)
          zz = 1.0 - 2.*rn(iseed)
          r2 = xx**2 + yy**2 + zz**2
          if(r2 .gt. 1.0 .or. r2 .lt. 0.000001) goto 100
          rr = sqrt(r2)

          betax = (px1+px2) / (e1+e2)
          betay = (py1+py2) / (e1+e2)
          betaz = (pz1+pz2) / (e1+e2)
          gamma  = 1.0 / sqrt(1.0-betax**2-betay**2-betaz**2)

          qx1   = pcm * xx/rr
          qy1   = pcm * yy/rr
          qz1   = pcm * zz/rr
          qx2   = -pcm * xx/rr
          qy2   = -pcm * yy/rr
          qz2   = -pcm * zz/rr

*             lorentz-transformation into lab frame
          e2cm  = sqrt (em2**2 + qx2**2 + qy2**2 + qz2**2)
          p2beta = qx2*betax + qy2*betay + qz2*betaz
          transf = gamma * ( gamma * p2beta / (gamma + 1.0) + e2cm )
          qx2 = betax * transf + qx2
          qy2 = betay * transf + qy2
          qz2 = betaz * transf + qz2

          ntag = 0
          if(id(1,i2).eq.1 .and. ipauli.eq.1)
     &      call pauli(i2,ntag,iseed,phase,r(1,i2),r(2,i2),r(3,i2),
     &                         qx2,qy2,qz2)
*
c          write(*,*) "phicoll, pauli után",ntag

          if(ntag .eq. -1) go to 600

          e1cm  = sqrt (em1**2 + qx1**2 + qy1**2 + qz1**2)
          p1beta = qx1*betax + qy1*betay + qz1*betaz
          transf = gamma * ( gamma * p1beta / (gamma + 1.0) + e1cm )
          p_pert(id_phi,1,i1) = betax * transf + qx1
          p_pert(id_phi,2,i1) = betay * transf + qy1
          p_pert(id_phi,3,i1) = betaz * transf + qz1
cc     note that the phi mass does not change by elastic collision
*************************************************************************

 !             write(20,*)  p_phi(1,i1), p_phi(2,i1) , p_phi(3,i1)

c          write(*,*) "phicoll, új imp"
          xtrav= sqrt(em1**2+p_pert(id_phi,1,i1)**2
     &                      +p_pert(id_phi,2,i1)**2)
          p0   = sqrt(xtrav**2 + p_pert(id_phi,3,i1)**2)
          pzlab= p0 * sinh(wref) + p_pert(id_phi,3,i1) * cosh(wref)
!          pzlab= p_phi(3,i1)
          plab = sqrt(pzlab**2 + p_pert(id_phi,1,i1)**2
     &                         + p_pert(id_phi,2,i1)**2)
          thdeg= 180./pi * acos(pzlab / plab)

          ppp =sqrt(px1**2+py1**2+pz1**2)
          ppel=sqrt(p_pert(id_phi,1,i1)**2+p_pert(id_phi,2,i1)**2
     &                                     +p_pert(id_phi,3,i1)**2)


          numcoll = numcoll + 1
          u5 = u5 +1

c          write(*,*) 'phi baryon scattering, nx_pert(id_phi,1,i1)=',
c     &     nx_pert(id_phi,1,i1)
c     &             'id(1,i2)=',id(1,i2)
c              write(*,*) 'phi probability:',p_pert(id_phi,4,i1)
c              write(*,*) 'old momenta (px,py,pz):',px1,py1,pz1
c              write(*,*) 'new momenta (px,py,pz):',
c     &             p_pert(id_phi,1,i1),p_pert(id_phi,2,i1),p_pert(id_phi,3,i1)
          if(abs(p_pert(id_phi,1,i1)+qx2-px1-px2).gt.1.e-3 .or.
     &           abs(p_pert(id_phi,2,i1)+qy2-py1-py2).gt.1.e-3 .or.
     &           abs(p_pert(id_phi,3,i1)+qz2-pz1-pz2).gt.1.e-3)
     &        write(*,*) 'hiba in phicoll, momentum nonconservation',
     &        p_pert(id_phi,1,i1),qx2,px1,px2,p_pert(id_phi,2,i1),qy2,
     &        py1,py2,p_pert(id_phi,3,i1),qz2,pz1,pz2
c          if(p_pert(id_phi,1,i1)**2+p_pert(id_phi,2,i1)**2+p_pert(id_phi,3,i1)**2.gt.2.5)
c     &      write(*,*) 'warning in phicoll, momentum is large',
c     &        p_pert(id_phi,1,i1),qx2,px1,px2,p_pert(id_phi,2,i1),qy2,py1,py2,
c     &        p_pert(id_phi,3,i1),qz2,pz1,pz2

          icollk = icollk + 1
!            collphi(int(time*2.)) = collphi(int(time*2.)) + 1
c                   write(54,*) ' phi elastic ',icollk,time
c             ix=nint(r_pert(id_phi,1,i1))
c             iy=nint(r_pert(id_phi,2,i1))
c             iz=nint(r_pert(id_phi,3,i1))
c             if(iabs(ix).gt.maxx) ix=ix/iabs(ix)*maxx
c             if(iabs(iy).gt.maxx) iy=iy/iabs(iy)*maxx
c             if(iabs(iz).gt.maxz) iz=iz/iabs(iz)*maxz
cc              ithermvar = ithermvar + 1
cc              thermok(2,i1) = ithermvar
cc              if(ithermo.eq.1) call tmunu(ix,iy,iz)
c              if(ithermo.eq.1) write(mterpri,'(i8,f6.2,3i4,14e11.3)')
c     &         ithermvar,time,ix,iy,iz,(avp(jj,ix,iy,iz),jj=1,14)
c-----------------------------------------------------------------------
*
************************************************************************
*                                                                      *

  800   continue
 1000 continue

c      write(*,*) "phicoll, coll over"
      if (nt.eq.ntmax .and. isu.eq.isubs) then
        write(20,*)
        write(20,*) "phi,coll",ttt,numcoll
        write(20,*)
      endif

c       write(54,*) 'number of phi+N coll,absorption:',numcoll,numabss
      write(*,*) 'total no. of fictive phi-s:',kanum


      return
      end
