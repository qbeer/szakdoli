************************************************************************
*
      subroutine phi_dec(wref)
*
*     slightly changed from zetenyi miklos's version  by  hw
*     this version is only used if kaons cannot appear via other channels
*     phi decays to dileptons and K+K-, 3pi decays, but pi's not generated
************************************************************************
      implicit none
      include"common"
      include"com_kminu"
      include"com_pert"
      include"cominput"

      integer irun,i1,maxk,maxkm,maxph, inem,maxem, ie_min
      integer ink, inkm, inkp ,iepair, idestroy
      integer  bin_dens,bin_time,ix,iy,iz
      integer ixi,iyi,izi,count_all,count_vac,count_med,cex
      real*8 x,y,z,px,py,pz,mas_ph,mas_kp,mas_km, checkdens
      real*8 gradx,grady,gradz,gradpx, gradpy, gradpz  ! not used
      real*8 ff,kk,kk0,Gtot0,Gtot,tau,tau0,gamma, br_epair,BRkao
      real*8 amkcha, amkneu, br_3pi, BRk0, denst
      real*8 gam_kao, gam_k, gam_fix, br_k, rck
      real*8 Ephi,decprob_tot,decprob_kao,prob_kao, prob_epair
      real*8 Ekp,Ekm,px_kp,py_kp,pz_kp,px_km,py_km,pz_km
      real*8 xx,yy,zz,rr,betax,betay,betaz
      real*8 rn, pxp,pyp,pzp,pkm,pka,sphiE
      integer kaonum,ikao_min,ikp,ikm

      real*8 chlab,shlab,ptra2,pcms,ptra,xtrav,p0,pzlab,E0lab
      real*8 plab,tlab,thc,thl,thldeg,wref

!       parameter(Gtot0=.004458)           ! GeV
!       parameter(BRkao=0.492, br_epair=0.000297,BRk0=0.340)
! !       parameter(BRkao=1.0, br_epair=0.000297,BRk0=0.340)
!       parameter(amkcha=.4937, amkneu=.4977, br_3pi = .17)
      parameter(Gtot0=.00426)           ! GeV
!       parameter(Gtot0=.04)           ! GeV
      parameter(BRkao=0.4920, br_epair=0.000297,BRk0=0.340)
      parameter(amkcha=.4937, amkneu=.4977, br_3pi = .153)
      save count_all,count_vac,count_med,cex,iepair,idestroy
c
      write(*,*) 'In phi_dec'
c
      gam_fix = br_3pi * Gtot0
      gam_kao = Gtot0 - gam_fix
      rck     = BRkao / (1.-br_3pi)
c
      maxk = maxkaon/num
      maxph = max_pert/num
      maxkm = max_kminu/num
      do 1000 irun = 1,num
        ink = (irun-1) * maxk
        inkp = (irun-1) * maxph
        inkm = (irun-1) * maxkm
        do 800 i1 = inkp+1,inkp+maxph
          if(nx_pert(id_phi,0,i1) .eq. 0)        goto 800
          if(nx_pert(id_phi,5,i1) .eq. 0)        goto 800
          x  = r_pert(id_phi,1,i1)
          y  = r_pert(id_phi,2,i1)
          z  = r_pert(id_phi,3,i1)
          px = p_pert(id_phi,1,i1)
          py = p_pert(id_phi,2,i1)
          pz = p_pert(id_phi,3,i1)
c.
         chlab = cosh(wref)
         shlab = sinh(wref)
         ptra2 = px**2 + py**2
         pcms  = sqrt(ptra2+pz**2)
         ptra  = sqrt(ptra2)
         xtrav = sqrt(p_pert(id_phi,0,i1)**2+ptra2)
         p0    = sqrt(xtrav**2 + pz**2)
         pzlab =  p0 * shlab + pz * chlab
         e0lab  =  p0 * chlab + pz * shlab
         plab  =  sqrt(pzlab**2 + ptra2)
         tlab  =  sqrt(p_pert(id_phi,0,i1)**2 + pzlab**2 + ptra2)
     &               - p_pert(id_phi,0,i1)
         thc   = acos(pz / pcms)
         thl   = acos(pzlab / plab)
         thldeg = 180./pi * thl
!        write(54,*) "phi ", thldeg

c calculate the effective masses:
          if (iphi_pot .ge. 1) then
            call graduphi(-1, x,y,z, px,py,pz,mas_ph,
     1         gradx,grady,gradz,gradpx, gradpy, gradpz)
          else
            mas_ph = xphimas
          endif
          ff = xkmas/xphimas
          if (ikaonpot .gt. 0) then
            call gradukaon(-1, x,y,z, px*ff,py*ff,pz*ff,mas_kp,1,
     1         gradx,grady,gradz,gradpx, gradpy, gradpz)
          else
            mas_kp = xkmas
          endif
          if (i_kminu_pot  .gt. 0) then
            call gradukao2(-1, x,y,z, px*ff,py*ff,pz*ff,mas_km,-1,
     1         gradx,grady,gradz,gradpx, gradpy, gradpz)
          else
            mas_km = xkmas
          endif
          if(mas_ph .lt. mas_kp+mas_km+0.00001) then
          write(*,*) 'phi below K+ K- threshold!: ',
     &         mas_ph,mas_kp,mas_km
           gam_k = .0
          else
            kk = sqrt((mas_ph**2 - mas_kp**2 - mas_km**2)**2
     &           - (2.*mas_km*mas_kp)**2) / (2.*mas_ph)
            kk0 = sqrt(xphimas**2/4. - xkmas**2)
            gam_k =  gam_kao * (xphimas/mas_ph)**2 * (kk/kk0)**3
          end if
          Gtot = gam_fix + gam_k
          br_k = rck * gam_k / Gtot
          Ephi = sqrt(mas_ph**2 + px**2 + py**2 + pz**2)
          gamma = Ephi / mas_ph
          tau0 = hbc/Gtot0 * gamma
          tau  = hbc/Gtot  * gamma
          decprob_tot = 1. - exp(-dt/tau)
          decprob_kao = p_pert(id_phi,4,i1)* br_k
!         write(54,*) ' tot, kao ',decprob_tot,decprob_kao,br_k,dt,tau
c         write(*,*)  '  before decay of phi no,',i1, nx_pert(id_phi,5,i1),
c    1                   p_pert(id_phi,4,i1)
c------------
      if (i_epair .eq. 1)  then           !  electron pair
      if (mas_ph .lt. .975) goto 590
      prob_epair  = p_pert(id_phi,4,i1)*dt*br_epair/tau0
          maxem = max_epair / num
          inem = (irun-1) * maxem
          do 530 ikm = inem+1,inem+maxem
            if(nx_epair(0,ikm).eq.0)                goto 530
            if (abs(p_epair(1,ikm) - px) .gt. .003) goto 530
            if (abs(p_epair(2,ikm) - py) .gt. .003) goto 530
            if (abs(p_epair(3,ikm) - pz) .gt. .003) goto 530
            p_epair(4,ikm) = p_epair(4,ikm) + prob_epair
            nx_epair(1,ikm) = nx_epair(1,ikm) + 1
!           write(54,*) '  phi - at time ', ikm, nx_epair(1,ikm),
!      1                  prob_epair 
            goto 590
  530     continue
          prob_kao = 1000.
          ie_min = 0
          do ikm = inem+1,inem+maxem
            if(nx_epair(0,ikm).eq.0) goto 520
            if(p_epair(4,ikm).lt.prob_kao) then
              ie_min = ikm
              prob_kao = p_epair(4,ikm)
            endif
          end do
          if(prob_epair .lt. prob_kao) write(54,*) 'no place for e-pair'
          if(prob_epair .lt. prob_kao) goto 590
          ikm = ie_min
 520      continue

          iepair = iepair + 1
          nx_epair(0,ikm) = 1
          nx_epair(1,ikm) = 1
          p_epair(1,ikm) = px
          p_epair(2,ikm) = py
          p_epair(3,ikm) = pz
          p_epair(4,ikm) = prob_epair

          write(54,*) 'iepair: ',iepair

c
      endif        !      electron pair
c-------------
 590  continue

      if(time.lt.dt*ntmax) then
        if(rn(iseed).gt.decprob_tot) goto 800  !  phi doesn't decay
      else
        idestroy = idestroy + 1
        goto 345
      endif

 345  continue
          nx_pert(id_phi,5,i1) = 0              ! phi destroyed
c         write(*,*) ' phi_dec prob ',i1, decprob_tot, Gtot,
c    1    p_pert(id_phi,4,i1), Ephi, tau, nx_pert(id_phi,5,i1)
c-----------------------
c   find a place for K+:
          prob_kao = 1000.
          kaonum = 0
          ikao_min = 0
          do ikp = ink+1,ink+maxk
            if(ika(1,ikp).eq.0) goto 700
              kaonum = kaonum+1
              if(kaonum.gt.maxk) goto 600
              if(pkao(4,ikp).lt.prob_kao) then
                ikao_min = ikp
                prob_kao = pkao(4,ikp)
              end if
          end do
 600      if(decprob_kao.lt.prob_kao) write(54,*) 'no place for K+'
          if(decprob_kao.lt.prob_kao) goto 800
          write(54,*) 'low prob. K+ deleted'
          ikp = ikao_min
!           if(ika(5,ikp).ne.0) nx_kminu(4,nx_pert(id_phi,5,ikp))=0    !!fehler??
          if(ika(5,ikp).ne.0) nx_kminu(4,ika(5,ikp))=0
 700      continue
c   find a place for K-:
          prob_kao = 1000.
          ikao_min = 0
          do ikm = inkm+1,inkm+maxkm
            if(nx_kminu(0,ikm).eq.0) goto 500
            if(p_kminu(4,ikm).lt.prob_kao) then
              ikao_min = ikm
              prob_kao = p_kminu(4,ikm)
            endif
          end do
          if(decprob_kao.lt.prob_kao) write(54,*) 'no place for K-'
          if(decprob_kao.lt.prob_kao) goto 800
          write(54,*) 'low prob. K- deleted'
          ikm = ikao_min
          if(nx_kminu(4,ikm).ne.0) ika(5,nx_kminu(4,ikm))=0
 500      continue
c   create K+ and K-:
          Ekp = sqrt(mas_kp**2 + kk**2)
          Ekm = sqrt(mas_km**2 + kk**2)
 60       continue
          xx       = 1. - 2. * rn(iseed)
          yy       = 1. - 2. * rn(iseed)
          zz       = 1. - 2. * rn(iseed)
          rr       = sqrt( xx**2 + yy**2 + zz**2 )
          if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 60
          px_kp = kk * xx/rr
          py_kp = kk * yy/rr
          pz_kp = kk * zz/rr
          px_km = -kk * xx/rr
          py_km = -kk * yy/rr
          pz_km = -kk * zz/rr
          betax = - px/Ephi
          betay = - py/Ephi
          betaz = - pz/Ephi
          
              if(px_kp**2+py_kp**2+pz_kp**2.gt.Ekp**2) then
                 write(*,*) "hiba phi_dec lorentz, negative mass",
     &              px_kp,py_kp,pz_kp,Ekp
c                 stop
              end if
          call lorentz(betax,betay,betaz,px_kp,py_kp,pz_kp,Ekp)
              if(px_km**2+py_km**2+pz_km**2.gt.Ekm**2) then
                 write(*,*) "hiba phi_dec2 lorentz, negative mass",
     &              px_km,py_km,pz_km,Ekm
c                 stop
              end if
          call lorentz(betax,betay,betaz,px_km,py_km,pz_km,Ekm)

c          write(*,*) 'fourmom. of phi:       ',px,py,pz,Ephi
c          write(*,*) 'fourmom. of kaon pair: ',
c     &       px_kp+px_km, py_kp+py_km, pz_kp+pz_km, Ekp+Ekm
c    K+:
          ika(1,ikp) = 1
          ika(2,ikp) = 77      ! ireac
          ika(3,ikp) = 0
          ika(4,ikp) = 0
          ika(5,ikp) = ikm    ! index of the assoc. K-
!           ika(6,ikp) = ikm    ! index of the assoc. K-
          rkao(1,ikp) = x
          rkao(2,ikp) = y
          rkao(3,ikp) = z
          pkao(1,ikp) = px_kp
          pkao(2,ikp) = py_kp
          pkao(3,ikp) = pz_kp
          pkao(4,ikp) = decprob_kao
c-----------------------------------------------------------------------
c           density dependence
          ix = nint(x)
          iy = nint(y)
          iz = nint(z)
          denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &      denst = rhb(ix,iy,iz) / rho0
          bin_dens = int(denst*10.)
          bin_time = int(time*1.)
c          kpl_phi(bin_dens,bin_time)=kpl_phi(bin_dens,bin_time)+1
c-----------------------------------------------------------------------
c    K-:
          nx_kminu(0,ikm) = 1
          nx_kminu(1,ikm) = 77
          nx_kminu(2,ikm) = 0
          nx_kminu(3,ikm) = 0
          nx_kminu(4,ikm) = ikp  ! index of the assoc. K+
!           nx_kminu(6,ikm) = ikp  ! index of the assoc. K+
          r_kminu(1,ikm) = x
          r_kminu(2,ikm) = y
          r_kminu(3,ikm) = z
          p_kminu(1,ikm) = px_km
          p_kminu(2,ikm) = py_km
          p_kminu(3,ikm) = pz_km
          p_kminu(4,ikm) = decprob_kao
c          write(*,*) 'decprob_kao: ',decprob_kao
c-----------------------------------------------------------------------
c           density dependence
          kmi_phi(bin_dens,bin_time)=kmi_phi(bin_dens,bin_time)+
     1    p_kminu(4,ikm)/real(num*isubs)
c-----------------------------------------------------------------------
          pxp = px_km + px_kp
          pyp = py_km + py_kp
          pzp = pz_km + pz_kp

          pkm = sqrt(xkmas**2 + px_km**2 + py_km**2 + pz_km**2)
          pka = sqrt(xkmas**2 + px_kp**2 + py_kp**2 + pz_kp**2)
          sphiE = abs(sqrt((pkm+pka)**2-pxp**2-pyp**2-pzp**2)-xphimas)
c-----------------------------------------------------------------------
c           density dependence
  233 format(2f10.6,4i10.6)
        checkdens = 0.0
        count_all = count_all + 1
        ixi = nint(x)
        iyi = nint(y)
        izi = nint(z)
        if(abs(ixi).lt.maxx.and.abs(iyi).lt.maxx.and.abs(izi).lt.
     &   maxz) then
        checkdens = rhb(ixi,iyi,izi)/rho0
!       pkao(9, ikm)  = checkdens    	!dens at K- index
        pkao(9, ikm)  = sphiE           !ex dens at K- index
        pkao(10,ikm)  = time            !time at K- index
!       pkao(11,ikp)  = checkdens       !dens at K+ index
        pkao(11,ikp)  = sphiE           !ex dens at K+ index
        pkao(12,ikp)  = time            !time at K+ index

! 	  write(54,*) 'K+K- direct ',tau,time,nx_kminu(4,ikm),ika(5,ikp)

        if(sphiE.lt.0.025) cex = cex + 1
        if(checkdens.lt.0.3333) then
          count_vac = count_vac + 1
        else if(checkdens.ge.0.3333) then
          count_med = count_med + 1
!        write(54,*) ' dec in med',tau,time,checkdens,count_med,sphiE
        endif
        write(53,233)time,checkdens,count_all,count_med,count_vac,cex        !all created K+K- pairs
        else
!         write(54,*) ' phi_dec raus ',tau,time,checkdens rhb(ix,iy,iz) / rho0

        endif
!         write(54,*) ' counter:all,vac,med',count_all,count_vac,count_med


!        write(54,*) time,rhb(0,0,0) / rho0

c-----------------
 800    continue
 1000 continue
      if (idestroy.gt.0 .and. irun.eq.num) then
        write(54,*) 'idestroy, b : ', idestroy,b
      endif
      return
      end

