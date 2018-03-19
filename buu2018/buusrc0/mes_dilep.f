************************************************************************
*
      subroutine mes_dilep
*
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      real*8 x1,y1,z1,px1,py1,pz1,em1,e1, branch, branch_ee
      real*8 gamma_kin,tau,decprob,prob_epair, eps, zmass, dec_tot
      real*8 rho,jrho1,jrho2,jrho3,density,vmlife,vmgam,rself,iself,sgam
      real*8 ymass,trmass, gam0, vrel0, rn, rnxx, decprob0
      integer iix,iiy,iiz
      integer irun,inp,ii,i1,id1,maxem,ikm,inem,ie_min
      integer ichannel

      write(*,*) 'In mes_dilep'
c      call f77flush()

*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        do 800 ii  = 1,maxp-1
          i1  = ii + inp
          if(ipi(1,i1) .eq. 0)                                 goto 800
          if(.not.(ipi(1,i1) .eq. 3 .or. ipi(1,i1) .eq. 5))    goto 800
          if(ipi(1,i1) .eq. 3 .and. ipi(2,i1) .ne. 0)          goto 800
c         write(*,*) 'In mes_dilep inside ',ipi(1,i1)
c          call f77flush()
          if (ipi(1,i1) .eq. 3) then ! rho -> e+e-
            gam0    = 0.149
            branch  = 7.02e-6/gam0 * romas**3
            ichannel = 1
          else if (ipi(1,i1) .eq. 5) then ! omega -> e+e-
            gam0    = 0.00844
            branch  = 0.60e-6/gam0 * omass**3
            ichannel = 2
          else
            write(*,*) 'strange meson in mes_dilep, ipi(,1)=',ipi(1,i1)
            stop
          end if
          id1 = ipi(1,i1)
          x1  = rpi(1,i1)
          y1  = rpi(2,i1)
          z1  = rpi(3,i1)
          px1 = ppi(1,i1)
          py1 = ppi(2,i1)
          pz1 = ppi(3,i1)
          em1 = epi(i1)
          e1  = sqrt( em1**2 + px1**2 + py1**2 + pz1**2 )
          gamma_kin = e1/em1

c---------------------------------we need density  for discrimination
          iix = nint(x1)
          iiy = nint(y1)
          iiz = nint(z1)
          if(abs(iix).lt.maxx .and. abs(iiy).lt.maxx .and.
     &        abs(iiz).lt.maxz) then
            rho=rhob_4(0,iix,iiy,iiz)
            jrho1 = rhob_4(1,iix,iiy,iiz)
            jrho2 = rhob_4(2,iix,iiy,iiz)
            jrho3 = rhob_4(3,iix,iiy,iiz)
          else
            rho = 0.0
            jrho1 = 0.0
            jrho2 = 0.0
            jrho3 = 0.0
          end if
          density = (rho**2 - jrho1**2 - jrho2**2 - jrho3**2)
          if (density .lt. .0) then
            density = .0
          else
            density = sqrt(density)
          endif
          IF (icbro .eq. 3)  THEN
            eps = g_col_null*(density/0.16)*(dt/gamma_kin)/hbc
            if (abs(eps) .lt. 1.e-4) then
              decprob0 = eps
            else
              decprob0 = 1.0 - exp(-eps)
            endif
            rnxx   = rn(iseed)
            write(*,*)  ' decprob0 ', rnxx, decprob0, dt, gamma_kin
     1                    ,gam_null
            if (rnxx .lt.  decprob0)   then
c                        switch off in case of counting rho's
              ipi(1,i1) = 0
              goto 800
            endif
          ENDIF
          if (id1.eq.3) then
            vrel0 = .0
            call self_rho(vrel0,density,em1,rself,iself,sgam)
            ymass=romas+0.5*rself/romas
            trmass=1.99*pmass
            if (icbro .eq. 2) trmass = .0
c----------------------------
c         that's because the new spectral fct has an e independent Gamma

            zmass = 1.5 * (em1 - ymass)/(em1-trmass)
            if (zmass .gt. .8)  zmass = .8
            if (icbro.eq.2 .and. em1.le.mrho_lim) zmass = .0
          endif
          if(id1.eq.5)then
            vrel0 = .0
            call self_omega(vrel0,density,em1,rself,iself,sgam)
            ymass=omass+0.5*rself/omass
            trmass=2.99*pmass
            if (em1.lt. trmass) then
              ipi(1,i1) = 0
              goto 800
            endif
            zmass = .0
          endif
          vmgam = sgam
          if(ireslife.eq.1 )
ccc             Danielewicz description
     1      vmgam = (sgam/2. + 2.0*(em1-ymass)**2/sgam)
     2          /(1. - zmass)
          if(ireslife.gt.1) write(*,*) 'problem with ireslife'

c           vmlife=gamma_kin*hbc/vmgam
c           if(vmlife.gt.dt) vmlife=dt
c         end if

          branch_ee = branch / em1**3
          tau = hbc*gamma_kin / vmgam
          eps = dt    /tau
          if (abs(eps) .gt. 1.e-3) then
             dec_tot = 1. - exp(-eps)
          else
             dec_tot = eps
          endif
          decprob = branch_ee * dec_tot
c
c           if (id1 .eq. 3  .and.  em1 .lt. 2.*pmass+.001) then
c
          rnxx = rn(iseed)
          if (rnxx .lt. dec_tot) then
c                     forbids meson  decay
            ipi(1,i1) = 0     !           meson decays
          endif
c           endif
c        write(*,*) 'ipi,mass,vmgam,tau,dt,decprob',id1,i1,ipi(1,i1),
c     1    em1, vmgam,tau,dt,dec_tot, decprob
          maxem = max_epair / num
          inem = (irun-1) * maxem
          do 530 ikm = inem+1,inem+maxem
            if (nx_epair(0,ikm).eq.0)                goto 530
            if (nx_epair(2,ikm) .ne. ichannel)       goto 530
            if (abs(p_epair(0,ikm) - em1) .gt. .003) goto 530
            if (abs(p_epair(1,ikm) - px1) .gt. .003) goto 530
            if (abs(p_epair(2,ikm) - py1) .gt. .003) goto 530
            if (abs(p_epair(3,ikm) - pz1) .gt. .003) goto 530
            if (abs(p_epair(5,ikm) - density) .gt. .005) goto 530
            p_epair(4,ikm) = p_epair(4,ikm) + decprob
            nx_epair(1,ikm) = nx_epair(1,ikm) + 1
            goto 590
  530     continue
          prob_epair = 1000.
          ie_min = 0
          do ikm = inem+1,inem+maxem
            if(nx_epair(0,ikm).eq.0) goto 520
            if(p_epair(4,ikm).lt.prob_epair) then
              ie_min = ikm
              prob_epair = p_epair(4,ikm)
            endif
          end do
          if(prob_epair .lt. decprob) write(*,*) 'no place for e-pair'
          if(prob_epair .lt. decprob) goto 590
          ikm = ie_min
 520      continue
          nx_epair(0,ikm) = 1
          nx_epair(1,ikm) = 1
          nx_epair(2,ikm) = ichannel
          p_epair(0,ikm) = em1
          p_epair(1,ikm) = px1
          p_epair(2,ikm) = py1
          p_epair(3,ikm) = pz1
          p_epair(4,ikm) = decprob
          p_epair(5,ikm) = density
c
c-------------
 590      continue
c          write(*,*) 'decprob: ',decprob,ikm,nx_epair(1,ikm),i1
c-----------------
 800    continue
 1000 continue
      write(*,*) 'end of mes_dilep'
      call f77flush()
      return
      end

