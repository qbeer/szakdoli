************************************************************************
      subroutine init_dilepton(yta,ypr,ipiNbrems)
*     initialize variables for continuous channels of dilepton prod.
*
*     let j the index of the integration point in the dilepton momentum*
*          space, then qq and qy give the appopriate momenta in the    *
*                                                                      *
*         qq(i,j) /i=1-3/   -  i. coord. of the j. point in mom.space  *
*         qq(4,j)           -  the energy of the j. point              *
*         qy(0,j)           -  the invariant mass belonging to j.      *
*         qy(1,j)           -  the transverse mom belonging to j.      *
*         qy(2,j)           -  the rapidity       belonging to j.      *
*         qy(3,j)           -  angle in trans.dir.belonging to j.      *
*       variables:                                                     *
*         yref    - -1 * target rapidity in the frame                  *
*         dlmasl  - minimal dilepton invariant mass       (real,input) *
*         ddilmas - step in dilepton invariant mass       (real,input) *
*         ndlmas  - number of dilepton invariant mass  (integer,input) *
*         ny      - number of rapidities               (integer,input) *
*         nqt     - number of transverse momenta       (integer,input) *
*         nf      - number of angles                   (integer,input) *
*         ifilt   - 0-> no filter, 1-> exp. filter used(integer,input) *
*                                                                      *
c-----------------------------------------------------------------------
      implicit none
      include 'common'
c      include 'cominput'
      include 'com_cont_epair'
      integer ich,ibin,ii,jj,kk,ntotq,ntoty,JP_totq,JP_toty,iqt,iy,iph
      integer ipiNbrems
      real*8 yta,ypr,yrap,dilmass,qtra,tramass,phi,dphi,szig,nul
      integer nDYim,IR_totq,IR_toty

      real*8 DrellYan_Cross(3,nDYmax*nDYsrtmax)
      real*8 xxa(nDYmax),yya(nDYmax),yy2(nDYmax)

      write(*,*) 'in dilep_in'
      nul=0.d0
      iproba = 0
      if(ipiNbrems.eq.1)      call piNdilep_read
c      write(*,*) 'after piNdilep_read',ivmesdil,
c     &  ibrems,ipiNbrems,iresdalitz,imesdalitz
      ndlmas = max0(1,ndlmas)
      ny     = max0(1,ny    )
      nqt    = max0(1,nqt   )
      nf     = max0(1,nf    )
      ntotq  = ny*nqt*nf*ndlmas
      ntoty  = ny*nqt*ndlmas
      if(ntoty .gt. maxy .or. ntotq.gt.maxq) then
        write(isum,'(''c: hiba too many q-s are to be calculated'')')
        write(*,'(''hiba dilepin too many q-s are to be calculated'')')
        stop
      end if
      ymin = 0.5*((yta+ypr)-yminscal * (ypr-yta))
      ymax = 0.5*((yta+ypr)+ymaxscal * (ypr-yta))
c      write(*,*) 'dilep init rap', yta,ypr,ymin,ymax,yminscal
      dy     = (ymax-ymin) / float(ny)
      dqt    = qtmaxi / float(nqt)
c      if(idiltra.eq.1) dqt    = (qtmaxi-dlmasl) / float(nqt)
      dphi   = 2. * pi / float(nf)
      do ich = 1,n_cont_channel
        do ibin = 1,ndlmas
          prob_cont_epair(ich,ibin) = 0.
c          prob_cont_epair_hades(ich,ibin) = 0.
c          prob_cont_epair_dls(ich,ibin) = 0.
        end do
        do ibin = 1,maxy
        do ii = 1,maxde
          sig(ich,ibin,ii) = 0.0
c         iqq(ich,ibin,ii) = 0
        end do
        end do
      end do
      ii = 0
      jj = 0
      dilmass = dlmasl-ddlmaslow
      do 1500 ibin = 1,ndlmas
        if(ibin.le.ndlmaslow) then
          dilmass = dilmass + ddlmaslow
        else
          dilmass = dilmass + ddlmashigh
        end if
        do 1400 iqt= 1,nqt
          qtra = (float(iqt) - 0.5) * dqt
          tramass  = sqrt(dilmass**2+qtra**2)
          if(idiltra.eq.1) then
            tramass = dilmass + (float(iqt) - 0.5) * dqt
            qtra    = tramass**2-dilmass**2
            if(qtra .le. 0.0) then
              qtra = 0.0
            else
              qtra    = sqrt(qtra)
            end if
          end if
          do 1300 iy = 1,ny
            yrap = ymin + (float(iy) - 0.5) * dy
c            iacc = 1
c            if(ifilt.eq.1) call dlsnacc(dilmass,qtra,yrap+yref,iacc)
            jj = jj + 1
            do 1200 iph = 0,nf-1
              ii = ii + 1
              phi = - pi + (float(iph)+0.5) * dphi
              qy(0,ii) = dilmass
              qy(1,ii) = qtra
c              if(idiltra.eq.1) qy(1,ii) = tramass
              qy(2,ii) = yrap
              qy(3,ii) = phi
              qq(1,ii) = qtra * cos(phi)
              qq(2,ii) = qtra * sin(phi)
              qq(3,ii) = tramass * sinh(yrap)
              qq(4,ii) = tramass * cosh(yrap)
 1200       continue
 1300     continue
 1400   continue
 1500 continue
cccccccc             charmonium      ccccccccccccccccccccccccccccccccccc
      JP_dlmas = max0(1,JP_dlmas)
      JP_y     = max0(1,JP_y    )
      JP_qt    = max0(1,JP_qt   )
      JP_f     = max0(1,JP_f    )
      JP_totq  = JP_y*JP_qt*JP_f*JP_dlmas
      JP_toty  = JP_y*JP_qt*JP_dlmas
      if(JP_toty .gt. JP_maxy .or. JP_totq.gt.JP_maxq) then
        write(isum,'(''c: hiba too many JP_q-s are to be calculated'')')
        write(*,'(''hiba dilepin too many JP_q are to be calculated'')')
        stop
      end if
      JP_ymin = 0.5*((yta+ypr)-JP_yminscal * (ypr-yta))
      JP_ymax = 0.5*((yta+ypr)+JP_ymaxscal * (ypr-yta))
c      write(*,*) 'dilep init rap', yta,ypr,ymin,ymax,yminscal
      JP_dy     = (JP_ymax-JP_ymin) / float(JP_y)
      JP_dqt    = JP_qtmaxi / float(JP_qt)
      dphi   = 2. * pi / float(JP_f)
      do ich = 1,JP_channel
        do ibin = 1,JP_dlmas
          JP_prob_epair(ich,ibin) = 0.
        end do
        do ibin = 1,JP_maxy
        do ii = 1,maxde
          JP_sig(ich,ibin,ii) = 0.0
c         iqq(ich,ibin,ii) = 0
        end do
        end do
      end do
      ii = 0
      jj = 0
      do 2500 ibin = 1,JP_dlmas
        dilmass = JP_dlmasl + float(ibin-1)*JP_ddlmas
        do 2400 iqt= 1,JP_qt
          qtra = (float(iqt) - 0.5) * JP_dqt
          tramass  = sqrt(dilmass**2+qtra**2)
          do 2300 iy = 1,JP_y
            yrap = JP_ymin + (float(iy) - 0.5) * JP_dy
            jj = jj + 1
            do 2200 iph = 0,JP_f-1
              ii = ii + 1
              phi = - pi + (float(iph)+0.5) * dphi
              JP_qy(0,ii) = dilmass
              JP_qy(1,ii) = qtra
              JP_qy(2,ii) = yrap
              JP_qy(3,ii) = phi
              JP_qq(1,ii) = qtra * cos(phi)
              JP_qq(2,ii) = qtra * sin(phi)
              JP_qq(3,ii) = tramass * sinh(yrap)
              JP_qq(4,ii) = tramass * cosh(yrap)
 2200       continue
 2300     continue
 2400   continue
 2500 continue
cccccccc        intermediate region     cccccccccccccccccccccccccccccccc
      IR_dlmas = max0(1,IR_dlmas)
      IR_y     = max0(1,IR_y    )
      IR_qt    = max0(1,IR_qt   )
      IR_f     = max0(1,IR_f    )
      IR_totq  = IR_y*IR_qt*IR_f*IR_dlmas
      IR_toty  = IR_y*IR_qt*IR_dlmas
      if(IR_toty .gt. IR_maxy .or. IR_totq.gt.IR_maxq) then
        write(isum,'(''#Dil hiba too many IR_q-s'')')
        write(*,'(''hiba dilepin too many IR_q-s'')')
        stop
      end if
      IR_ymin = 0.5*((yta+ypr)-IR_yminscal * (ypr-yta))
      IR_ymax = 0.5*((yta+ypr)+IR_ymaxscal * (ypr-yta))
c      write(*,*) 'dilep init rap', yta,ypr,ymin,ymax,yminscal
      IR_dy     = (IR_ymax-IR_ymin) / float(IR_y)
      IR_dqt    = IR_qtmaxi / float(IR_qt)
      dphi   = 2. * pi / float(IR_f)
      do ich = 1,IR_channel
        do ibin = 1,IR_dlmas
          IR_prob_epair(ich,ibin) = 0.
        end do
        do ibin = 1,IR_maxy
          do ii = 1,maxde
            IR_sig(ich,ibin,ii) = 0.0
c         iqq(ich,ibin,ii) = 0
          end do
        end do
      end do
      ii = 0
      jj = 0
      do 3500 ibin = 1,IR_dlmas
        dilmass = IR_dlmasl + float(ibin-1)*IR_ddlmas
        do 3400 iqt= 1,IR_qt
          qtra = (float(iqt) - 0.5) * IR_dqt
          tramass  = sqrt(dilmass**2+qtra**2)
          do 3300 iy = 1,IR_y
            yrap = IR_ymin + (float(iy) - 0.5) * IR_dy
            jj = jj + 1
            do 3200 iph = 0,IR_f-1
              ii = ii + 1
              phi = - pi + (float(iph)+0.5) * dphi
              IR_qy(0,ii) = dilmass
              IR_qy(1,ii) = qtra
              IR_qy(2,ii) = yrap
              IR_qy(3,ii) = phi
              IR_qq(1,ii) = qtra * cos(phi)
              IR_qq(2,ii) = qtra * sin(phi)
              IR_qq(3,ii) = tramass * sinh(yrap)
              IR_qq(4,ii) = tramass * cosh(yrap)
 3200       continue
 3300     continue
 3400   continue
 3500 continue
      open(41,file="buuinput/Drell-Yan_cross_sect.dat",status='old')
      read(41,'(//)')
      read(41,*) nDYsrt, nDYim
      do ii=1,nDYsrt
        do jj=1,nDYim
          ibin=(ii-1)*nDYim+jj 
          read(41,*) (DrellYan_Cross(kk,ibin),kk=1,3)
          xxa(jj) = DrellYan_Cross(2,ibin)
          yya(jj) = DrellYan_Cross(3,ibin)
c          write(*,*) 'DY 0',ii,jj,ibin,xxa(jj),yya(jj)
        end do
        call splined(xxa,yya,nDYim,nul,nul,yy2)
c        write(*,*) 'DY 1',yy2
        do jj=1,IR_dlmas
          dilmass = IR_qy(0,(jj-1)*IR_y*IR_qt*IR_f+1)
          if(dilmass.ge.xxa(1) .and. dilmass.le.xxa(nDYim)) then
            call splintd(xxa,yya,yy2,nDYim,dilmass,szig)
          else
            szig=0.0
          end if
          DY_InvMcross(1,ii,jj) = DrellYan_Cross(1,(ii-1)*nDYim+1)
          DY_InvMcross(2,ii,jj) = dilmass
          DY_InvMcross(3,ii,jj) = 2.0*dilmass*szig/1.d6
c          write(*,*)'DY_InvMcross',
c     &         DrellYan_Cross(1,(ii-1)*nDYim+1),dilmass,szig
        end do
      end do
      close(41)

      write(*,*) 'Drell-Yan in',nDYsrt, nDYim
c      do ii=1,nDYsrt
c        do jj=1,IR_dlmas
c          write(*,*) ii,jj,(DY_InvMcross(kk,ii,jj),kk=1,3)
c        end do
c      end do

      return
      end

************************************************************************
      subroutine npbrems(srt,szig,beta,gamma,xxx,yyy,zzz)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
c      include 'cominput'
      include 'com_cont_epair'
      real*8 szig,srt,beta(3),gamma,xxx,yyy,zzz
      real*8 dilmass,sig_cugnon,sigma,prob,qmax,q0max,qbeta
      real*8 eprim, spr, phasespace
      integer ibin,iacc,ii,jj,iqt,iy,iph,nde,dens_bin
c     real*8 dilmass_cont

***                                                                  ***
***         preparing    for the pauli blocking                      ***
***                                                                  ***
c      call paulpro0(xxx,yyy,zzz)
***                                                                  ***
c      write(*,*) 'In npbrems (srt,sig) - ',srt,sig
      if(srt.le.2.*rmass+0.0001) return
      sig_cugnon = 35./(1. + (srt-2.*rmass)*100.0) + 20.
      prob = sig_cugnon/szig
      ii    = 0
      jj    = 0

c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      

      do ibin = 1,ndlmas
        dilmass = qy(0,(ibin-1)*nqt*ny*nf+1)
cc The one below was used by ZM, reason????
c        sig_cugnon = 35./(1. + (srt-2.*dilmass)) + 20.
c        dsdm = 4.*alfa**2/(3.*pi**2) * 1./dilmass *
c     &     (srt**2/(4.*rmass**2) - 1.) *
c     &     log((srt-2.*rmass)/dilmass) *
c     &     sig_cugnon
c        sigma = ddlmas * dsdm
c        prob = sigma/szig
c        write(*,*) '->>',ibin,dilmass,sig_cugnon,dsdm,sigma,prob
c        if (prob.gt.1.) write(*,*) "prob>1 in npbrems ",sigma,sig,prob
c        write(*,*) '--+++-->',prob_cont_epair(ch_brems,ibin)
        q0max= 0.5 * (srt**2 + dilmass**2 - 4.0*rmass**2) / srt**2
        if(q0max.le.dilmass+0.001) return
        qmax = sqrt(max(0.0, q0max**2-dilmass**2))
        sigma = alfa**2/(3.*pi**2) *1./dilmass * 
     &            (4.0-16.*rmass**2/srt**2) *
     &            (log((qmax+q0max)/dilmass)-qmax/q0max)*prob
        if(log((qmax+q0max)/dilmass)-qmax/q0max.lt.0.0) then
          write(*,*) 'npbrems ',qmax,q0max,dilmass
          sigma = 0.0
        end if
        prob_cont_epair(ch_brems,ibin) =
     &     prob_cont_epair(ch_brems,ibin)
     &     + sigma
c        write(*,*) 'npbrems --->',ch_brems,num,isubs,
c     &     prob_cont_epair(ch_brems,ibin)
        do 4302 iqt=1,nqt
        do 4301 iy =1,ny
          jj = jj + 1
        do 4300 iph=1,nf
          ii = ii + 1
          qbeta = beta(1)*qq(1,ii) + beta(2)*qq(2,ii) + beta(3)*qq(3,ii)
          eprim = gamma * (qq(4,ii) - qbeta)
          spr   = srt * (srt - 2. * eprim) + qy(0,ii)**2
          if(spr.le.4.*rmass**2)    goto 4300
          phasespace= sqrt((1.-4.*rmass**2/spr)/(1.-4.*rmass**2/srt**2))
          sigma = alfa**2/(3.*pi**2)/dilmass / eprim**2 *
     &            0.5*(srt**2-4.0*rmass**2)/rmass**2 * prob * phasespace
          sig(ch_brems,jj,nde)= sig(ch_brems,jj,nde) 
     &                       + sigma
 4300   continue
 4301   continue
 4302   continue
      end do
      return
      end

************************************************************************
      subroutine baryon_dalitz(idres,charge,width,srt,dt0,beta,gamm,
     &       xxx,yyy,zzz)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
c      include 'cominput'
      include 'com_cont_epair'
      integer idres             ! id of resonance type
      integer charge            ! charge of the resonance
      real*8 srt                  ! resonance mass = sqrt(s)
      real*8 dt0                  ! length of time step in rest frame of res.
      real*8 width                ! total gamma - needed to check
                                !   that BR(gamma+N) is small
                                !   (-> perturbative method is justified)
      real*8 gamm                 ! Lorentz gamma
      integer ibin,iacc,ii,jj,nde,iqt,iy,iph,dens_bin,bin_num_cont
      real*8 dilmass
      real*8 gamma,prob, beta(3),xxx,yyy,zzz,pt_dil,rap_dil
      real*8 dgdm_bary_dalitz
      real*8 rn,getHadesPairAcceptance,acc
      real*8 rapydel,qbtra,trams,xyz,rap1,rap2,dif1,dif2,qzpr,q0pr,factu
      real*8 rap, sigm,q0,qqabs,qq2

c        write(*,*) 'baryon_dalitz', idres,charge,width,srt,dt0,beta,gamm
c     &     ,xxx,yyy,zzz
c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      
      do ibin = 1,ndlmas
        dilmass = qy(0,(ibin-1)*nqt*ny*nf+1)
        if ((srt-dilmass-rmass).lt.0.001) go to 11
        gamma = dgdm_bary_dalitz(idres,charge,srt,dilmass)
c        prob = 1. - exp(-dt0 * gamma *ddlmas/ hbc)
        prob = dt0 * gamma / hbc
        if(dt0.gt.10.) prob = gamma/width
        prob_cont_epair(idres+bary_chan_offset,ibin) =
     &     prob_cont_epair(idres+bary_chan_offset,ibin)
     &     + prob
        ii  = (ibin-1) * nqt * ny * nf
        jj  = (ibin-1) * nqt * ny
        q0      = 0.5 * (srt**2 + dilmass**2 - rmass**2) / srt
        qq2     = q0**2 - dilmass**2
        if(qq2.le.0.0)                   return
        qqabs   = sqrt(qq2)
************************************************************************
*                                                                      *
*         to fulfill the energy conservation - q=q0 in the delta frame *
*           we look for the appropriate rapidity at a given qt, y(qt)  *
*           function. it may happen that the domain of qt, where a     *
*           solution exist, is inside of the qt steps, so there will   *
*           be no contribution. when the statistic is good enough, this*
*           will not couse big problem, the error remain in 20 %.      *
*                                                                      *
*           when there is no solution, then xyz**2 < 1                 *
*                                                                      *
************************************************************************
        rapydel = 0.5 * log((1.+beta(3))/(1.-beta(3)))
        do 3300 iqt=1,nqt
          do 3200 iy=1,ny
            jj = jj + 1
            do 3100 iph=1,nf
              ii = ii + 1
              qbtra = qq(1,ii) * beta(1) + qq(2,ii) * beta(2)
              trams = sqrt(qy(0,ii)**2 + qy(1,ii)**2)
              if(idiltra.eq.1) trams = qy(1,ii)
              if(trams.le.qy(0,ii)) goto 3300
              xyz   = cosh(rapydel) * (q0/gamm+qbtra)/trams
              if(abs(xyz) .le. 1.0)                           goto 3100
              rap1= rapydel - log(xyz + sqrt(xyz**2 - 1.0))
              rap2= rapydel + log(xyz + sqrt(xyz**2 - 1.0))
              dif1= abs(rap1-qy(2,ii))
              dif2= abs(rap2-qy(2,ii))
              if((dif1.gt.dy/2.).and.(dif2.gt.dy/2.))         goto 3100
              if(dif1.le.dy/2.) rap = rap1
              if(dif2.le.dy/2.) rap = rap2
              qzpr = trams * sinh(rap)
              q0pr = trams * cosh(rap)
              factu = 4.*pi * qqabs * gamm * dy * abs(qzpr-beta(3)*q0pr)
              sigm = prob*2.*pi/factu
              sig(idres+bary_chan_offset,jj,nde) = 
     &           sig(idres+bary_chan_offset,jj,nde)+sigm
c HADES acceptance:
              pt_dil = qy(1,ii)
              rap_dil = qy(2,ii)
              if(ihades.eq.1) then
                acc =getHadesPairAcceptance(dilmass,pt_dil,rap_dil,0)
                prob_cont_epair_hades(idres+bary_chan_offset,ibin) =
     &              prob_cont_epair_hades(idres+bary_chan_offset,ibin)
     &              + acc * prob *dy * dqt * qy(1,ii)
              end if
c--------------
              iacc = 1
              if(ifilt.eq.1)call dlsnacc(dilmass,qy(1,ii),qy(2,ii),iacc)
              prob_cont_epair_dls(idres+bary_chan_offset,ibin) =
     &            prob_cont_epair_dls(idres+bary_chan_offset,ibin)
     &        + float(iacc) * prob *dy * dqt * qy(1,ii)
 3100       continue
 3200     continue
 3300   continue
      end do
 11   continue
      return
      end

************************************************************************
      subroutine meson_dalitz(idm,srt,dt0,beta,gamm,xxx,yyy,zzz)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
c      include 'cominput'
c      include 'commonthreef'
      include 'com_cont_epair'
      integer idm ! meson identification 1: pion, 2 eta, 5 omega
      real*8 srt         ! meson mass = sqrt(s)
      real*8 dt0    ! length of time step in rest frame of res. >10:final
      real*8 gamm                 ! Lorentz gamma
      integer ibin,iacc,ii,jj,nde,iqt,iy,iph,dens_bin,ch_mes
      real*8 dilmass,meswidth,restmass
      real*8 gamma,prob, beta(3),xxx,yyy,zzz,pt_dil,rap_dil
      real*8 dgdm_omega_dalitz,dgdm_eta_dalitz,dgdm_pion_dalitz
      real*8 rapydel,qbtra,trams,xyz,rap1,rap2,dif1,dif2,qzpr,q0pr,factu
      real*8 rap, sigm,q0,qqabs,qq2
      real*8 rn,getHadesPairAcceptance,acc

c      write(*,*) 'START meson_dalitz',srt,beta
c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      
      if(idm.eq.1) then
        ch_mes   = ch_pi
        meswidth = piwidth
        restmass = 0.0
      else if(idm.eq.2) then 
        ch_mes   = ch_eta
        meswidth = ewidth
        restmass = pmass
      else if(idm.eq.5) then 
        ch_mes   = ch_ome
        meswidth = owidth
        restmass = pmass
      else
        write(*,*) 'hiba in mesdalitz', idm
      end if

      do ibin = 1,ndlmas
        dilmass = qy(0,(ibin-1)*nqt*ny*nf+1)
        if ((srt-dilmass-restmass).lt.0.001) go to 11
        if(idm.eq.1) then
          gamma = dgdm_pion_dalitz(srt,dilmass)
        else if(idm.eq.2) then 
          gamma = dgdm_eta_dalitz(srt,dilmass)
        else if(idm.eq.5) then 
          gamma = dgdm_omega_dalitz(srt,dilmass)
        end if
c        prob = 1. - exp(-dt0 * gamma / hbc)
        prob = dt0 * gamma / hbc
        if(dt0.gt.10.0) prob=gamma/meswidth
        prob_cont_epair(ch_mes,ibin) =
     &     prob_cont_epair(ch_mes,ibin)
     &     + prob

        ii  = (ibin-1) * nqt * ny * nf
        jj  = (ibin-1) * nqt * ny
        q0      = 0.5 * (srt**2 + dilmass**2 - pmass**2) / srt
        qq2     = q0**2 - dilmass**2
        if(qq2.le.0.0)                return
        qqabs   = sqrt(qq2)
************************************************************************
*                                                                      *
*         to fulfill the energy conservation - q=q0 in the delta frame *
*           we look for the appropriate rapidity at a given qt, y(qt)  *
*           function. it may happen that the domain of qt, where a     *
*           solution exist, is inside of the qt steps, so there will   *
*           be no contribution. when the statistic is good enough, this*
*           will not couse big problem, the error remain in 20 %.      *
*                                                                      *
*           when there is no solution, then xyz**2 < 1                 *
*                                                                      *
************************************************************************
        rapydel = 0.5 * log((1.+beta(3))/(1.-beta(3)))
        do 3300 iqt=1,nqt
          do 3200 iy=1,ny
            jj = jj + 1
            do 3100 iph=1,nf
              ii = ii + 1
              qbtra = qq(1,ii) * beta(1) + qq(2,ii) * beta(2)
              trams = sqrt(qy(0,ii)**2 + qy(1,ii)**2)
c              if(idiltra.eq.1) trams = qy(1,ii)
c              if(trams.le.qy(0,ii)) goto 3300
              xyz   = cosh(rapydel) * (q0/gamm+qbtra)/trams
              if(abs(xyz) .le. 1.0)                           goto 3100
              rap1= rapydel - log(xyz + sqrt(xyz**2 - 1.0))
              rap2= rapydel + log(xyz + sqrt(xyz**2 - 1.0))
              dif1= abs(rap1-qy(2,ii))
              dif2= abs(rap2-qy(2,ii))
              if((dif1.gt.dy/2.).and.(dif2.gt.dy/2.))         goto 3100
              if(dif1.le.dy/2.) rap = rap1
              if(dif2.le.dy/2.) rap = rap2
              qzpr = trams * sinh(rap)
              q0pr = trams * cosh(rap)
              factu = 4.*pi * qqabs * gamm * dy * abs(qzpr-beta(3)*q0pr)
              sigm = prob*2.*pi/factu
              sig(ch_mes,jj,nde) = 
     &           sig(ch_mes,jj,nde)+sigm
c HADES acceptance:
              pt_dil = qy(1,ii)
              rap_dil = qy(2,ii)
              if(ihades.eq.1) then
                acc = getHadesPairAcceptance(dilmass,pt_dil,rap_dil,0)
c          smear according to detector resolution:
                call smearHadesPair(dilmass,pt_dil,rap_dil,1) ! 3 -> high resolution
c          store result:
                prob_cont_epair_hades(ch_mes,ibin) =
     &              prob_cont_epair_hades(ch_mes,ibin)
     &              + prob * acc *dy * dqt * qy(1,ii)
              end if
              iacc = 1
              if(ifilt.eq.1) call dlsnacc(dilmass,pt_dil,rap_dil,iacc)
              prob_cont_epair_dls(ch_mes,ibin) =
     &            prob_cont_epair_dls(ch_mes,ibin)
     &            + prob*float(iacc) *dy *dqt *qy(1,ii)
 3100       continue
 3200     continue
 3300   continue
      end do
 11   continue
      return
      end

************************************************************************
      subroutine piNdilep_read
*     read the piN-dilep data files
*
*                                                                      *
c-----------------------------------------------------------------------
      implicit none
      integer ibin,ii,isrt,nsrt,numdmass,kk,numdmax
      integer maxpiNsrt,maxpiNmass
      real*8 factrho,factome,factrhome,dilmass,srt
      parameter(maxpiNsrt=150, maxpiNmass=210)
      real*8 piNdilep(4,maxpiNsrt,maxpiNmass,5)
      common/piNdilep/ piNdilep,nsrt,numdmax
      open(41,file="buuinput/piNdilep_n_pi0_fine.dat",status='old')
      open(42,file="buuinput/piNdilep_n_pi+_fine.dat",status='old')
      open(43,file="buuinput/piNdilep_p_pi0_fine.dat",status='old')
      open(44,file="buuinput/piNdilep_p_pi-_fine.dat",status='old')
      numdmax = 0
      do ibin=1,4
        read(40+ibin,'(/////)')
        read(40+ibin,*) nsrt
        write(*,*) 'piNdilep readin',nsrt
        if(maxpiNsrt.lt.nsrt) then
          write(*,*) 'hiba:piNdilep_read, nsrt is too big'
        end if
        nsrt = min(nsrt,maxpiNsrt)
        do isrt=1,nsrt
          read(40+ibin,*) numdmass
          numdmax = max(numdmax,numdmass)
          if(maxpiNmass.lt.numdmass) then
            write(*,*) 'hiba:piNdilep_read, numdmass is too big'
            stop
          end if
          write(*,*) numdmass
          do ii=1,numdmass
            read(40+ibin,*) (piNdilep(ibin,isrt,ii,kk),kk=1,5)
c        from mikrobarn to mb
            piNdilep(ibin,isrt,ii,3) = 1.e-3 * piNdilep(ibin,isrt,ii,3)
            piNdilep(ibin,isrt,ii,4) = 1.e-3 * piNdilep(ibin,isrt,ii,4)
            piNdilep(ibin,isrt,ii,5) = 1.e-3 * piNdilep(ibin,isrt,ii,5)
          end do
c          write(*,*) ibin,isrt,numdmass,
c     &                  (piNdilep(ibin,isrt,1,kk),kk=1,5)
c          write(*,*) ibin,isrt,numdmass,
c     &                  (piNdilep(ibin,isrt,numdmass,kk),kk=1,5)
          read(40+ibin,'(//)')
        end do
        close(40+ibin)
      end do
      write(*,*) 'piNdilep kiiras nsrt',nsrt,numdmax
      do ibin=1,4
      do isrt=1,nsrt
      do ii=1,numdmax
        piNdilep(ibin,isrt,ii,1) = piNdilep(ibin,isrt,1,1)
        piNdilep(ibin,isrt,ii,2) = piNdilep(ibin,nsrt,ii,2)
      end do
      end do
      end do

      return
      end
************************************************************************
      subroutine piNdilepCreate(charge_pi,charge_N,ipara,srt,szig,beta,
     &       gamm,xxx,yyy,zzz)
c    momdep potential energy change at creation is not taken into account
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      include 'com_cont_epair'
      include 'cominput'
      integer charge_pi,charge_N ! charge of the resonance
      integer ipara              ! index of the paralel runs
      real*8 srt                   ! total energy sqrt(s)
      real*8 szig                  ! xsection for piN collision
      real*8 gamm                  ! Lorentz gamma
      integer ichannel,isrt,imass,ibin,nde,ix,iy,iz,dens_bin
      integer ii,iqt,iph,jj
      real*8 beta(3),xxx,yyy,zzz,xx,yy,zz,rr,p_rhome,e_rhome,dilmass
      real*8 rn,px_rhome,py_rhome,pz_rhome,density
      real*8 qzpr,rap,rap1,rap2,rapydel,sigm,trams,xyz
      real*8 q0pr,qbtra,dif1,dif2,factu,factrho,factome,factrhome
      integer maxrhomega
      parameter(maxrhomega=500000)
      integer irhomega(3,maxrhomega)
      real*8 rhomega(30,maxrhomega)
      common/rhomega_dat/ rhomega,irhomega
      integer maxpiNsrt,maxpiNmass,nsrt,numdmax
      parameter(maxpiNsrt=150, maxpiNmass=210)
      real*8 piNdilep(4,maxpiNsrt,maxpiNmass,5)
      common/piNdilep/ piNdilep,nsrt,numdmax

      if(rn(iseed).gt.1./float(npiNbrems)) return

      iproba = iproba + 1
      if(charge_pi+charge_N.lt.0 .or. charge_pi+charge_N.gt.1) return
      if(charge_pi.eq.0 .and. charge_N.eq.0) ichannel = 1
      if(charge_pi.eq.1 .and. charge_N.eq.0) ichannel = 2
      if(charge_pi.eq.0 .and. charge_N.eq.1) ichannel = 3
      if(charge_pi.eq.-1.and. charge_N.eq.1) ichannel = 4

      density = 0.0
      if(abs(xxx).lt.maxx .and. abs(yyy).lt.maxx .and.
     &            abs(zzz).lt.maxz) then
        ix = nint(xxx)
        iy = nint(yyy)
        iz = nint(zzz)
        density = sqrt(rhob_4(0,ix,iy,iz)**2-rhob_4(1,ix,iy,iz)**2-
     &                 rhob_4(2,ix,iy,iz)**2-rhob_4(3,ix,iy,iz)**2)
      end if

c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      

c   determination isrt for the given srt
c      write(*,*) 'piNCreate dens', srt, density, nde,ix,iy,iz,beta

      do isrt = 1,nsrt
        if(srt.lt.piNdilep(ichannel,isrt,1,1)) goto 11
      end do
      isrt = nsrt
 11   continue
      if(isrt.gt.1) then
        if(abs(srt-piNdilep(ichannel,isrt,1,1)).gt.
     &     abs(srt-piNdilep(ichannel,isrt-1,1,1))) isrt=isrt-1
      end if

c      if(iproba.eq.1) then
c      write(*,'("# piNdilep specfunc in matter in origo",f10.4)') srt 
c      write(*,'("# mass    rho    omega   interf    sum frho fome")')
c      do ii=1,numdmax
c        dilmass = piNdilep(2,isrt,ii,2)
c      call piNcrossmodfact(0.0,0.0,0.0,0.0,0.0,0.0,
c     &     dilmass,factrho,factome,factrhome)
c      write(*,'(f10.3,7e10.2)')dilmass,piNdilep(2,isrt,ii,3)*factrho,
c     & piNdilep(2,isrt,ii,4)*factome,piNdilep(2,isrt,ii,5)*factrhome,
c     & piNdilep(2,isrt,ii,3)*factrho+piNdilep(2,isrt,ii,4)*factome+
c     & piNdilep(2,isrt,ii,5)*factrhome,factrho,factome,factrhome
c      end do
c      end if
c      stop
      do ibin = 1,ndlmas
        dilmass = qy(0,(ibin-1)*nqt*ny*nf+1)
        if(dilmass+rmass+0.001 .gt. srt) return
cc  determination of imass to dilmass 
        do imass = 1,numdmax
          if(dilmass.lt.piNdilep(ichannel,isrt,imass,2)) goto 12
        end do
        imass = numdmax
 12     continue
        if(imass.gt.1) then
          if(abs(dilmass-piNdilep(ichannel,isrt,imass,2)).gt.
     &     abs(dilmass-piNdilep(ichannel,isrt,imass-1,2))) imass=imass-1
        end if

        e_rhome = (srt**2+dilmass**2 -rmass**2)/(2.*srt)
        if(e_rhome.le.dilmass) return 
        p_rhome = sqrt(e_rhome**2-dilmass**2)

c        write(*,*) 'piNdilepCreate1',e_rhome, p_rhome, dilmass
 13     continue
        xx = rn(iseed) - 0.5
        yy = rn(iseed) - 0.5
        zz = rn(iseed) - 0.5
        rr = sqrt(xx**2 + yy**2 + zz**2)
        if( (rr.lt.0.0001) .or. (rr.gt.0.25) )  goto 13
        px_rhome = p_rhome * xx/rr
        py_rhome = p_rhome * yy/rr
        pz_rhome = p_rhome * zz/rr
          if(px_rhome**2+py_rhome**2+pz_rhome**2.gt.e_rhome**2) then
            write(*,*) "hiba dilepton lorentz, negative mass",
     &                px_rhome,py_rhome,pz_rhome,e_rhome
c            stop
          end if
        call lorentz(-beta(1),-beta(2),-beta(3),
     1                 px_rhome,py_rhome,pz_rhome,e_rhome)
        factrho   = 1.0
        factome   = 1.0
        factrhome = 1.0
c        write(*,*) 'piNdilepCreate2',xxx,yyy,zzz,px_rhome,py_rhome,
c     &     pz_rhome,dilmass,factrho,factome,factrhome
        if(icbro.ge.1 .or. ivecmatt.ge.1) 
     &     call piNcrossmodfact(xxx,yyy,zzz,px_rhome,py_rhome,pz_rhome,
     &     dilmass,factrho,factome,factrhome)

cc     store the piN cross sections
c        if(iproba.eq.7)
c     &   write(*,*) 'piNcross',prob_cont_epair(ch_rho_piN,ibin),
c     &    prob_cont_epair(ch_ome_piN,ibin),
c     &    prob_cont_epair(ch_rhome_piN,ibin),
c     &    prob_cont_epair(ch_piN,ibin),ibin,imass,isrt
        prob_cont_epair(ch_rho_piN,ibin) =
     &     prob_cont_epair(ch_rho_piN,ibin)
     &     + piNdilep(ichannel,isrt,imass,3)/szig*factrho
     &     * float(npiNbrems)
        prob_cont_epair(ch_ome_piN,ibin) =
     &     prob_cont_epair(ch_ome_piN,ibin)
     &     + piNdilep(ichannel,isrt,imass,4)/szig*factome
     &     * float(npiNbrems)
        prob_cont_epair(ch_rhome_piN,ibin) =
     &     prob_cont_epair(ch_rhome_piN,ibin)
     &     + piNdilep(ichannel,isrt,imass,5)/szig*factrhome
     &     * float(npiNbrems)
        prob_cont_epair(ch_piN,ibin) =
     &     prob_cont_epair(ch_piN,ibin)
     &     + piNdilep(ichannel,isrt,imass,3)/szig*factrho
     &     * float(npiNbrems)
     &     + piNdilep(ichannel,isrt,imass,4)/szig*factome
     &     * float(npiNbrems)
     &     + piNdilep(ichannel,isrt,imass,5)/szig*factrhome
     &     * float(npiNbrems)
c        if(iproba.eq.7)
c     &   write(*,*) 'piNdilepcross2',prob_cont_epair(ch_rho_piN,ibin),
c     &    prob_cont_epair(ch_ome_piN,ibin),
c     &    prob_cont_epair(ch_rhome_piN,ibin),
c     &    prob_cont_epair(ch_piN,ibin),ibin,imass,isrt
ccc   piN contribution to the 3 dim diff cross section
        ii  = (ibin-1) * nqt * ny * nf
        jj  = (ibin-1) * nqt * ny
************************************************************************
*                                                                      *
*         to fulfill the energy conservation - q=q0 in the delta frame *
*           we look for the appropriate rapidity at a given qt, y(qt)  *
*           function. it may happen that the domain of qt, where a     *
*           solution exist, is inside of the qt steps, so there will   *
*           be no contribution. when the statistic is good enough, this*
*           will not couse big problem, the error remain in 20 %.      *
*                                                                      *
*           when there is no solution, then xyz**2 < 1                 *
*                                                                      *
************************************************************************
        rapydel = 0.5 * log((1.+beta(3))/(1.-beta(3)))
        do 3300 iqt=1,nqt
          do 3200 iy=1,ny
            jj = jj + 1
            do 3100 iph=1,nf
              ii = ii + 1
              qbtra = qq(1,ii) * beta(1) + qq(2,ii) * beta(2)
              trams = sqrt(qy(0,ii)**2 + qy(1,ii)**2)
c              if(idiltra.eq.1) trams = qy(1,ii)
c              if(trams.le.qy(0,ii)) goto 3300
              xyz   = cosh(rapydel) * (e_rhome/gamm+qbtra)/trams
              if(abs(xyz) .le. 1.0)                           goto 3100
              rap1= rapydel - log(xyz + sqrt(xyz**2 - 1.0))
              rap2= rapydel + log(xyz + sqrt(xyz**2 - 1.0))
              dif1= abs(rap1-qy(2,ii))
              dif2= abs(rap2-qy(2,ii))
              if((dif1.gt.dy/2.).and.(dif2.gt.dy/2.))         goto 3100
              if(dif1.le.dy/2.) rap = rap1
              if(dif2.le.dy/2.) rap = rap2
              qzpr = trams * sinh(rap)
              q0pr = trams * cosh(rap)
              factu = 4.*pi*p_rhome * gamm * dy * abs(qzpr-beta(3)*q0pr)

              sigm = piNdilep(ichannel,isrt,imass,3)/szig *2.*pi/factu
              sig(ch_rho_piN,jj,nde) = 
     &           sig(ch_rho_piN,jj,nde)+sigm * float(npiNbrems)

              sigm = piNdilep(ichannel,isrt,imass,4)/szig *2.*pi/factu
              sig(ch_ome_piN,jj,nde) = 
     &           sig(ch_ome_piN,jj,nde)+sigm * float(npiNbrems)
              sigm = piNdilep(ichannel,isrt,imass,5)/szig *2.*pi/factu
              sig(ch_rhome_piN,jj,nde) = 
     &           sig(ch_rhome_piN,jj,nde)+sigm * float(npiNbrems)
              sigm = (piNdilep(ichannel,isrt,imass,3)
     &               +piNdilep(ichannel,isrt,imass,4)
     &               +piNdilep(ichannel,isrt,imass,5))/szig *2.*pi/factu
              sig(ch_piN,jj,nde) = 
     &           sig(ch_piN,jj,nde)+sigm * float(npiNbrems)
c HADES acceptance: not finished yet
c              pt_dil = qy(1,ii)
c              rap_dil = qy(2,ii)
c              if(ihades.eq.1) then
c                acc =getHadesPairAcceptance(dilmass,pt_dil,rap_dil,0)
c                prob_cont_epair_hades(idres+bary_chan_offset,ibin) =
c     &              prob_cont_epair_hades(idres+bary_chan_offset,ibin)
c     &              + acc * prob *dy * dqt * qy(1,ii)/ eventnum
c              end if
c--------------
 3100       continue
 3200     continue
 3300   continue

c***********************************************************************
c
c rhomega(i,j) the jth particle
c i:1    rho_mass
c   2-4  rho_momentum
c   5-7  rho_position
c   8    creation probability
c   9    summed probability of the decay
c   10   potential energy
c   11-20 the same for omega
c   21-29 the same for interference terms (mass, position, mom. dont evolve) 
c   30   local density of creation
c
c***********************************************************************

        do ii=1,maxrhomega
          if(irhomega(1,ii).eq.0.and.irhomega(2,ii).eq.0) goto 14
        end do
        write(*,*) 'hiba: piNdilep, too many rhomegas'
        ii=maxrhomega
        stop
c        write(*,*) 'piNdilep creation ',ipara,dilmass,px_rhome,py_rhome
c     &    ,pz_rhome,xxx,yyy,zzz,szig,factrho,factome,factrhome
c        if(iproba.eq.7)
c     &  write(*,*) 'piNdilep 1',iproba,ii,dilmass,imass,srt,isrt,
c     &      piNdilep(ichannel,isrt,imass,3)/szig
 14     continue
        irhomega(1,ii) = 1
        irhomega(2,ii) = 1
        irhomega(3,ii) = ipara
        rhomega(1,ii)  = dilmass    ! rho part
        rhomega(2,ii)  = px_rhome
        rhomega(3,ii)  = py_rhome
        rhomega(4,ii)  = pz_rhome
        rhomega(5,ii)  = xxx
        rhomega(6,ii)  = yyy
        rhomega(7,ii)  = zzz
        rhomega(8,ii)  = piNdilep(ichannel,isrt,imass,3)/szig*factrho
     &                   *float(npiNbrems)
        rhomega(9,ii)  = 0.0
        rhomega(10,ii) = 0.0  ! rho potential energy
        rhomega(11,ii) = dilmass  ! starts the omega part
        rhomega(12,ii) = px_rhome
        rhomega(13,ii) = py_rhome
        rhomega(14,ii) = pz_rhome
        rhomega(15,ii) = xxx
        rhomega(16,ii) = yyy
        rhomega(17,ii) = zzz
        rhomega(18,ii) = piNdilep(ichannel,isrt,imass,4)/szig*factome
     &                   *float(npiNbrems)
        rhomega(19,ii) = 0.0
        rhomega(20,ii) = 0.0  ! omega potential energy
        rhomega(21,ii) = dilmass
        rhomega(22,ii) = px_rhome
        rhomega(23,ii) = py_rhome
        rhomega(24,ii) = pz_rhome
        rhomega(25,ii) = xxx
        rhomega(26,ii) = yyy
        rhomega(27,ii) = zzz
        rhomega(28,ii) = piNdilep(ichannel,isrt,imass,5)/szig*factrhome
     &                   *float(npiNbrems)
        rhomega(29,ii) = 0.0
        rhomega(30,ii) = density
c        write(*,*) 'piNdileprhomegacreate', ii,rhomega(1,ii),
c     &    rhomega(8,ii),rhomega(18,ii),rhomega(28,ii),isrt,imass    
c 1000   continue
      end do
      return
      end

************************************************************************
      subroutine piNdilepDecay(dt0)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      include 'cominput'
      include 'com_cont_epair'
      real*8 dt0                   ! length of time step in rest frame of res.

      integer ibin,idm
      real*8 vmass,gamma,bwmes,prob,prob_rho,prob_ome,en_rho,en_ome
      real*8 px,py,pz,xxx,yyy,zzz,theta
      integer maxrhomega
      parameter(maxrhomega=500000)
      integer irhomega(3,maxrhomega)
      real*8 rhomega(30,maxrhomega)
      common/rhomega_dat/ rhomega,irhomega
c
c      write(*,*) 'piNdilepdecay1',dt0
      do ibin=1,maxrhomega
cc        if(ibin.eq.7) write(*,*) 'piNdilepdecay7',dt0,irhomega(1,ibin),
c     & irhomega(2,ibin),irhomega(3,ibin),rhomega(1,ibin),rhomega(8,ibin)
c     & ,rhomega(9,ibin),prob_cont_epair(ch_rho_piN_prop,1),
c     &    prob_cont_epair(ch_ome_piN_prop,1),
c     &    prob_cont_epair(ch_rhome_piN_prop,1),
c     &    prob_cont_epair(ch_piN_prop,1)
        if(irhomega(1,ibin).eq.1.or.irhomega(2,ibin).eq.1) then
          prob_rho= 0.
          vmass=rhomega(1,ibin)
          px  = rhomega(2,ibin)
          py  = rhomega(3,ibin)
          pz  = rhomega(4,ibin)
          en_rho = sqrt(vmass**2+px**2+py**2+pz**2)
          if(irhomega(1,ibin).eq.1) then
            xxx = rhomega(5,ibin)
            yyy = rhomega(6,ibin)
            zzz = rhomega(7,ibin)
            gamma = 0.0
            idm = 1
            if(vmass.gt.2.*pmass)
     &        gamma=bwmes(vmass**2,idm,idec2pi,iresmode,0,0,iwidth,0)
            prob_rho= dt0 * gamma / hbc
            if(gamma.gt.1.e-2 .and. dt0.lt.10.) 
     &              prob_rho=1.-exp(-dt0 * gamma / hbc)
            if(dt0.ge.10.0) prob_rho=1.0
            prob_rho= prob_rho * (1.-rhomega(9,ibin))
            prob = prob_rho * rhomega(8,ibin)
            idm = 3
c            write(*,*) 'piNdecay cal vect',idm,vmass,px,py,pz,
c     &       xxx,yyy,zzz,prob
            call vectmes_dilep(idm,vmass,px,py,pz,xxx,yyy,zzz,prob,2)
c            write(*,*) 'piNdecay-rho',vmass,prob_rho,prob
c            if(ibin.eq.7) write(*,*) 'piNdilep 1D',dt0,vmass,ibin,
c     &      rhomega(8,ibin),prob_rho,rhomega(9,ibin),prob,
c     &      prob_cont_epair(ch_rho_piN_prop,1)
          end if
          prob_ome= 0.
          vmass=rhomega(11,ibin)
          px  = rhomega(12,ibin)
          py  = rhomega(13,ibin)
          pz  = rhomega(14,ibin)
          en_ome = sqrt(vmass**2+px**2+py**2+pz**2)
          if(irhomega(2,ibin).eq.1) then
            xxx = rhomega(15,ibin)
            yyy = rhomega(16,ibin)
            zzz = rhomega(17,ibin)
            gamma = 0.0
            idm = 3
            if(vmass.gt.3.*pmass)
     &        gamma=bwmes(vmass**2,idm,idec2pi,iresmode,0,0,iwidth,0)
            prob_ome= dt0 * gamma / hbc
            if(gamma.gt.1.e-2 .and. dt0.lt.10.) 
     &              prob_ome=1.-exp(-dt0 * gamma / hbc)
            if(dt0.ge.10.0) prob_ome=1.0
            prob_ome= prob_ome * (1.-rhomega(19,ibin))
            prob = prob_ome * rhomega(18,ibin)
            idm = 5
            call vectmes_dilep(idm,vmass,px,py,pz,xxx,yyy,zzz,prob,2)
c            write(*,*) 'piNdecay-rho',vmass,prob_ome,prob,gamma,dt0
          end if
          vmass=rhomega(21,ibin)
          px  = rhomega(22,ibin)
          py  = rhomega(23,ibin)
          pz  = rhomega(24,ibin)
          xxx = rhomega(25,ibin)
          yyy = rhomega(26,ibin)
          zzz = rhomega(27,ibin)
          prob = prob_rho* prob_ome
     &      + prob_rho*rhomega(19,ibin) + prob_ome*rhomega(9,ibin)
          rhomega(9,ibin)  = rhomega(9,ibin)  + prob_rho
          rhomega(19,ibin) = rhomega(19,ibin) + prob_ome
          rhomega(29,ibin) = rhomega(29,ibin) + prob
          theta = abs(en_ome-en_rho) * dt0/hbc
          theta = min(theta,0.5*pi)
c          theta = 0.0
          prob = rhomega(28,ibin) * cos(theta) * prob
          idm = -1
          call vectmes_dilep(idm,vmass,px,py,pz,xxx,yyy,zzz,prob,2)
c           write(*,*) 'piNdecay', rhomega(9,ibin),rhomega(19,ibin),
c     & rhomega(29,ibin)
        if(rhomega(9,ibin).ge.1.001 .or. rhomega(19,ibin).ge.1.001
     &    .or. rhomega(19,ibin).ge.1.001) stop
        end if
        if(dt0.gt.10.0) then
          irhomega(1,ibin) = 0
          irhomega(2,ibin) = 0
          irhomega(3,ibin) = 0
        end if
      end do
      return
      end

************************************************************************
      subroutine piNdilepAbs(dt0)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      include 'cominput'
      include 'com_cont_epair'
      real*8 dt0                   ! length of time step in rest frame of res.
      real*8 pirr
      parameter(pirr = 1.3)
      integer ibin,i2,idt,ii,idm
      real*8 dtfac,dtprime,ddt,etot,meff
      real*8 gradrx,gradry,gradrz,gradpx,gradpy,gradpz,gradm
      real*8 vmass,gamma,bwmes,prob,prob_rho,prob_ome,en_rho,en_ome
      real*8 px,py,pz,xxx,yyy,zzz,theta,dummyf(1:9),bwdist
      real*8 dxxx,dyyy,dzzz,rsqare,px2,py2,pz2,em2,e2,srt2,p12,p1dr,p2dr
      real*8 a12,b12,c12,brel,b21,t1,t2,qq2,sigt,isofac,help,rn
      integer maxrhomega
      parameter(maxrhomega=500000)
      integer irhomega(3,maxrhomega)
      real*8 rhomega(30,maxrhomega)
      common/rhomega_dat/ rhomega,irhomega

      write(*,*) 'pindilepabs',nres
      if(nres.ne.24) stop
c                rho
      do ibin=1,maxrhomega
        if(irhomega(1,ibin).eq.1 .or. irhomega(2,ibin).eq.1) then
        prob_rho = 0.0
        vmass=rhomega(1,ibin)
c       if(vmass.le.2.*pmass) return
        px  = rhomega(2,ibin)
        py  = rhomega(3,ibin)
        pz  = rhomega(4,ibin)
        en_rho  = sqrt(vmass**2+px**2+py**2+pz**2)
        if(irhomega(1,ibin).eq.1) then
          xxx = rhomega(5,ibin)
          yyy = rhomega(6,ibin)
          zzz = rhomega(7,ibin)

c--- selection of baryons:
          i2  = (irhomega(3,ibin)-1) * totmass
c          write(*,*)'piNdilepabs', ibin,vmass,i2,totmass,delpi
  600     i2  = i2 + 1
          if(i2 .gt. irhomega(3,ibin)*totmass)               goto 700
          if(id(1,i2).ne.1)                                  goto 600
          dxxx     = xxx - r(1,i2)
          if (abs(dxxx) .gt. delpi)                          goto 600
          dyyy     = yyy - r(2,i2)
          if (abs(dyyy) .gt. delpi)                          goto 600
          dzzz     = zzz - r(3,i2)
          if (abs(dzzz) .gt. delpi)                          goto 600
          rsqare = dxxx**2 + dyyy**2 + dzzz**2
          if (rsqare .gt. delpi**2)                          goto 600
*         now particles are close enough to each other !
c          if(ipi(1,i1).eq.3) write(*,*) 'zm: rho+bar. are close enough'

          px2    = p(1,i2)
          py2    = p(2,i2)
          pz2    = p(3,i2)
          em2    = e(i2)
          e2     = sqrt (em2**2+px2**2+py2**2+pz2**2 )
          srt2   = (en_rho+e2)**2 - (px+px2)**2 - (py+py2)**2
     &                                           - (pz+pz2)**2
*   is their impact parameter small enough?
          p12    = en_rho * e2 - px * px2 - py * py2 - pz * pz2
          p1dr   = px * dxxx + py * dyyy + pz * dzzz
          p2dr   = px2 * dxxx + py2 * dyyy + pz2 * dzzz
          a12    = 1.0 - ( vmass * em2 / p12 )**2
          b12    = p1dr / vmass - p2dr * vmass / p12
          c12    = rsqare + ( p1dr / vmass )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if(brel .gt. pirr)                                  goto 600
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / vmass - b12 / a12 ) * en_rho / vmass
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
c   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt0 )                          goto 600
*   now  the pion may be absorbed in this time step
*
          qq2= 0.25*(srt2-em2**2+vmass**2)**2/srt2-vmass**2
          sigt = 0.0
          do ii = 1, nres
*------- calculate the isospin factors for the reactions ------------*
**       isospin 1/2 x 1
            if(resprop2(ii,1).eq.3) then
              isofac=2./3.
            else if(resprop2(ii,1).eq.1) then
              isofac=1./3.
            end if

** 3/2 x 1 delta rho is not included
c            write(*,*) 'before bwrho',srt2,ii,idec2pi,iresmode,4,
c     &                     0,iwidth,1,dummyf
            help = bwdist(srt2,ii,idec2pi,iresmode,4,
     &                     0,iwidth,1,dummyf)
            sigt = sigt + 40.0*pi/qq2 * isofac * (resprop2(ii,3)+1.)/6.0
     &                 * help * hbc**2
c            write(*,*) 'after bwrho',help,sigt
          end do
          if(rn(iseed).le. sigt/(10.*pi*pirr**2)) then
            irhomega(1,ibin)=0
            prob_rho = 1.0 - rhomega(9,ibin)
c            write(*,*) 'piNdilep ome absorbed', sigt,srt2
          end if
        end if
 700    continue

        prob_ome = 0.0
        vmass=rhomega(11,ibin)
c       if(vmass.le.3.*pmass) return
        px  = rhomega(12,ibin)
        py  = rhomega(13,ibin)
        pz  = rhomega(14,ibin)
        en_ome  = sqrt(vmass**2+px**2+py**2+pz**2)
        if(irhomega(2,ibin).eq.1) then
          xxx = rhomega(15,ibin)
          yyy = rhomega(16,ibin)
          zzz = rhomega(17,ibin)

c--- selection of baryons:
          i2  = (irhomega(3,ibin)-1) * totmass
c          write(*,*)'t1 ', irun, ii
  800     i2  = i2 + 1
          if(i2 .gt. irhomega(3,ibin)*totmass)                 goto 900
          if(id(1,i2).ne.1)                                    goto 800
          dxxx     = xxx - r(1,i2)
          if (abs(dxxx) .gt. delpi)                            goto 800
          dyyy     = yyy - r(2,i2)
          if (abs(dyyy) .gt. delpi)                            goto 800
          dzzz     = zzz - r(3,i2)
          if (abs(dzzz) .gt. delpi)                            goto 800
          rsqare = dxxx**2 + dyyy**2 + dzzz**2
          if (rsqare .gt. delpi**2)                            goto 800
*        now particles are close enough to each other !
c          if(ipi(1,i1).eq.3) write(*,*) 'zm: rho+bar. are close enough'

          px2        = p(1,i2)
          py2        = p(2,i2)
          pz2        = p(3,i2)
          em2        = e(i2)
          e2         = sqrt (em2**2+px2**2+py2**2+pz2**2 )
          srt2       = (en_ome+e2)**2 - (px+px2)**2 - (py+py2)**2
     &                                           - (pz+pz2)**2
*   is their impact parameter small enough?
          p12    = en_ome * e2 - px * px2 - py * py2 - pz * pz2
          p1dr   = px * dxxx + py * dyyy + pz * dzzz
          p2dr   = px2 * dxxx + py2 * dyyy + pz2 * dzzz
          a12    = 1.0 - ( vmass * em2 / p12 )**2
          b12    = p1dr / vmass - p2dr * vmass / p12
          c12    = rsqare + ( p1dr / vmass )**2
          brel   = sqrt( abs(c12 - b12**2/a12) )
*   is their impact parameter small enough?
          if(brel .gt. pirr)                                  goto 800
          b21    = - p2dr / em2 + p1dr * em2 / p12
          t1     = ( p1dr / vmass - b12 / a12 ) * en_ome / vmass
          t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2
*   will particles get closest point in this time interval ?
          if ( abs(t1+t2) .gt. dt0 )                            goto 800
*   now  the pion may be absorbed in this time step
*
          qq2= 0.25*(srt2-em2**2+vmass**2)**2/srt2-vmass**2
          sigt = 0.0
          do ii = 1, nres

*------- calculate the isospin factors for the reactions ------------*
**       isospin 1/2 x 0
            isofac=1.0
c            write(*,*) 'before bwome',srt2,ii,idec2pi,iresmode,4,
c     &                     0,iwidth,1,dummyf
            help = bwdist(srt2,ii,idec2pi,iresmode,9,
     &                     0,iwidth,1,dummyf)
            sigt = sigt + 40.0*pi/qq2 * isofac * (resprop2(ii,3)+1.)/6.
     &                 * help * hbc**2
c            write(*,*) 'after bwome',help,sigt
          end do
          if(rn(iseed).le. sigt/(10.*pi*pirr**2)) then
            irhomega(2,ibin)=0
            prob_ome = 1.0 - rhomega(19,ibin)
c            write(*,*) 'piNdilep ome absorbed', sigt,srt2
          end if
        end if
 900    continue
        if(icbro.ge.1) then
          theta = abs(en_ome-en_rho) * dt0/hbc
          theta = min(theta,0.5*pi)
          prob = prob_rho* prob_ome
     &      + prob_rho*rhomega(19,ibin) + prob_ome*rhomega(9,ibin)
          prob = rhomega(28,ibin) * cos(theta) * prob
          idm = -1
          call vectmes_dilep(idm,vmass,px,py,pz,xxx,yyy,zzz,prob,2)
        end if
      end if
      end do
      return
      end
************************************************************************
      subroutine piNdilepProp(dt0)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
      include 'cominput'
      include 'com_cont_epair'
      real*8 dt0                   ! length of time step in rest frame of res.

      integer ndtfac
      parameter(ndtfac = 2)
      integer ibin,idt,ipropag,id1,id2
      real*8 dtfac,dtprime,ddt,etot,meff,newepi,masslim
      real*8 gradrx,gradry,gradrz,gradpx,gradpy,gradpz,gradm
      real*8 vmass,gamma,bwmes
      real*8 px,py,pz,xxx,yyy,zzz
c      real*8 dxxx,dyyy,dzzz,rsquare,px2,py2,pz2,em2,e2,srt2,p12,p1dr,p2dr
      integer maxrhomega
      parameter(maxrhomega=500000)
      integer irhomega(3,maxrhomega)
      real*8 rhomega(30,maxrhomega)
      common/rhomega_dat/ rhomega,irhomega
c
c rhomega(i,j) the jth particle
c i:1    rho_mass
c   2-4  rho_momentum
c   5-7  rho_position
c   8    creation probability
c   9    summed probability of the decay
c   10   potential energy
c   11-20 the same for omega
c   21   original mass (in case of spect. evolution relevant)
c   22   "probability" of the interference
c   23   summed "probability" of the interference
c   30   local density of creation
c
c***********************************************************************
c
c             propagation
      dtfac = 1.0/float(ndtfac)
      dtprime = dt0 * dtfac
      id2=0
      do idt = 1,ndtfac
c        write(*,*) ' in propa piNdil',idt,ndtfac
        ddt = idt*dtfac

        do ibin = 1,maxrhomega
          if(irhomega(1,ibin).eq.1) then
            vmass=rhomega(1,ibin)

            xxx  =  rhomega(5,ibin)
            yyy  =  rhomega(6,ibin)
            zzz  =  rhomega(7,ibin)

            px  =  rhomega(2,ibin)
            py  =  rhomega(3,ibin)
            pz  =  rhomega(4,ibin)
 
            meff =  vmass + rhomega(10,ibin)
            etot = sqrt(meff**2+px**2+py**2+pz**2)
c---------
            ipropag = 1
            id1=3
c            write(*,*) 'piNDilep1 ', xxx, yyy, zzz, px, py, pz, etot, 3, 
c     &      0, vmass,
c     &              gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
c     &              ddt,dtprime,ipropag
            call gradupi(xxx, yyy, zzz, px, py, pz, etot,id1,id2,vmass,
     &              gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
     &              ddt,dtprime,ipropag)

            etot = sqrt(meff**2+rhomega(2,ibin)**2+rhomega(3,ibin)**2
     &                     +rhomega(4,ibin)**2)
c
            newepi = vmass + dtprime * gradm * etot/meff  ! problem with
                                                            ! mpot(i)  ??!
            if (newepi .gt. .01 .and. newepi .lt. 1.5) then
c              write(*,*) 'piNdilepPropa ',px,py,pz,gradrx,gradry,gradrz
c     &               ,dtprime
              rhomega(2,ibin) = px - dtprime * gradrx
              rhomega(3,ibin) = py - dtprime * gradry
              rhomega(4,ibin) = pz - dtprime * gradrz
              etot = sqrt(meff**2+rhomega(2,ibin)**2+rhomega(3,ibin)**2
     &                     +rhomega(4,ibin)**2)
              rhomega(1,ibin) = newepi
c             if (id1.eq.5) write(52,*) ' mass ',id1, i,gradm,newepi
            else
c              write(*,*) 'propa: too small rho mass',ibin,newepi,
c     1          dtprime,gradm,etot,meff
            end if

            rhomega(5,ibin) = rhomega(5,ibin) + dtprime * gradpx
            rhomega(6,ibin) = rhomega(6,ibin) + dtprime * gradpy
            rhomega(7,ibin) = rhomega(7,ibin) + dtprime * gradpz
          end if
          if(irhomega(2,ibin).eq.1) then
            vmass=rhomega(11,ibin)

            xxx  =  rhomega(15,ibin)
            yyy  =  rhomega(16,ibin)
            zzz  =  rhomega(17,ibin)

            px  =  rhomega(12,ibin)
            py  =  rhomega(13,ibin)
            pz  =  rhomega(14,ibin)
 
            meff =  vmass + rhomega(20,ibin)
            etot = sqrt(meff**2+px**2+py**2+pz**2)
c---------
            ipropag = 1
            id1=5
c            write(*,*) 'piNDilep2 ', xxx, yyy, zzz, px, py, pz, etot, 5, 
c     &      0, vmass,
c     &              gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
c     &              ddt,dtprime,ipropag
            call gradupi(xxx, yyy, zzz, px, py, pz, etot,id1,id2,vmass,
     &              gradrx, gradry, gradrz,gradpx,gradpy,gradpz,gradm,
     &              ddt,dtprime,ipropag)

            etot = sqrt(meff**2+rhomega(12,ibin)**2+rhomega(13,ibin)**2
     &                     +rhomega(14,ibin)**2)
c
            masslim = 3.0 * pmass
            if (icbro .gt. 1)  masslim = .04
            newepi = vmass + dtprime * gradm * etot/meff  ! problem with
                                                            ! mpot(i)  ??!
            if (newepi .gt.01 .and. newepi .lt. 1.5) then
c              write(*,*) 'piNdilepPropa2',px,py,pz,gradrx,gradry,gradrz
c     &               ,dtprime
              rhomega(12,ibin) = px - dtprime * gradrx
              rhomega(13,ibin) = py - dtprime * gradry
              rhomega(14,ibin) = pz - dtprime * gradrz
              etot =sqrt(meff**2+rhomega(12,ibin)**2+rhomega(13,ibin)**2
     &                     +rhomega(14,ibin)**2)
              rhomega(11,ibin) = newepi
            else
c              write(*,*) 'propa hiba: too small ome mass',ibin,newepi,
c     1          dtprime,gradm,etot,meff
            end if

            rhomega(15,ibin) = rhomega(15,ibin) + dtprime * gradpx
            rhomega(16,ibin) = rhomega(16,ibin) + dtprime * gradpy
            rhomega(17,ibin) = rhomega(17,ibin) + dtprime * gradpz
          end if
        end do
      end do
      return
      end

************************************************************************
      subroutine piNdilepRun(dt0)
c        
c-----------------------------------------------------------------------
      implicit none
      real*8 dt0                   ! length of time step in rest frame of res.
      call f77flush()
      write(*,*) 'piNdilepRun', dt0
c             decay
      call piNdilepDecay(dt0)

      if(dt0 .gt.5.0) return
c             propagation
      call piNdilepProp(dt0)

cccc          absorption
      call piNdilepAbs(dt0)
*-----------------------------------------------------------------------
      return
      end

************************************************************************
      subroutine vectmes_dilep(id,vmass,px,py,pz,xxx,yyy,zzz,prob,ich)
c-----------------------------------------------------------------------
      implicit none
c      include 'common'
c      include 'cominput'
      include 'com_cont_epair'
      integer id                ! id=1: rho, 3: omega, -1 interference term
      integer i_reg !1 low mass, 2: intermediate
      real*8 vmass       ! mass of the vector meson
      real*8 px,py,pz,xxx,yyy,zzz ! monemtum and position of vector mes
      integer ich      ! 1 direct vector meson, 2 piNdilep
      integer nde, dens_bin,ichannel,ibin,imass,iqt,iy
      real*8 rap,en,qt,dvmass,volume,prob
      real*8 br_rho, br_ome
      parameter(br_rho=4.7e-5, br_ome=7.3e-5)

      integer mass_bin_num,qt_bin_num,y_bin_num
      real*8 mass_bin_size,qt_bin_size,y_bin_size
c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      

      i_reg = 1
      if(ich.eq.1 .and. id.eq.3)   ichannel = ch_rho_dir
      if(ich.eq.1 .and. id.eq.5)   ichannel = ch_ome_dir
      if(ich.eq.2 .and. id.eq.3)   ichannel = ch_rho_piN_prop
      if(ich.eq.2 .and. id.eq.5)   ichannel = ch_ome_piN_prop
      if(ich.eq.2 .and. id.eq.-1)  ichannel = ch_rhome_piN_prop
c      write(*,*) 'vectmesdecay', ichannel, id, vmass,px,py,pz,xxx,yyy,
c     &   zzz,prob,ich,nde
c      call f77flush()

c      write(*,*) 'piNdilep Vectmes', ichannel,vmass
      imass = mass_bin_num(vmass,i_reg)
      if(imass.le.0 .or. imass.gt.ndlmas) return
      dvmass = mass_bin_size(imass)
      if(ich.eq.1 .and. id.eq.3) prob=prob*br_rho/dvmass
      if(ich.eq.1 .and. id.eq.5) prob=prob*br_ome/dvmass
c      write(*,*) 'piNdilep Vectmes3', ichannel,vmass,imass,iqt,iy,prob
      prob_cont_epair(ichannel,imass)=prob_cont_epair(ichannel,imass)
     &                 + prob
      
      en = sqrt(vmass**2+px**2+py**2+pz**2)
      qt = sqrt(px**2+py**2)
      rap = 0.5 * log((en+pz)/(en-pz))
      iqt = qt_bin_num(qt)
      iy = y_bin_num(rap)
c      write(*,*) 'piNdilep Vectmes', ichannel,vmass,imass,iqt,iy,prob
      if(iqt.lt.1 .or. iqt.gt.nqt
     &  .or. iy.lt.1 .or. iy.gt.ny) return
      ibin = (imass-1)*nqt*ny+(iqt-1)*ny+(iy-1)+1
      if(ibin.gt.maxy) then
        write(*,*) 'hiba vec-dec ', ibin,ndlmas,nqt,ny,nf,
     &   imass,iqt,iy
      stop
      end if
c      dqt = qt_bin_size(iqt)
c      dy = y_bin_size(iy)
      volume = qt*dqt*dy
      sig(ichannel,ibin,nde) = sig(ichannel,ibin,nde) 
     &                     + prob/volume
      return
      end

************************************************************************
      subroutine JPsi_dilep(id,vmass,px,py,pz,xxx,yyy,zzz,prob)
c-----------------------------------------------------------------------
      implicit none
      include 'com_pert'
      include 'com_cont_epair'
      integer id       ! id=1: rho, 3: omega, -1 interference term
      real*8 vmass       ! mass of the vector meson
      real*8 px,py,pz,xxx,yyy,zzz ! momentum and position of vector mes
      integer nde, dens_bin,ibin,imass,iqt,iy
      real*8 rap,en,qt,dvmass,volume,prob

      integer JP_mass_bin_num,JP_qt_bin_num,JP_y_bin_num
c      real*8 mass_bin_size,qt_bin_size,y_bin_size
c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      

c      call f77flush()
      if(id.gt.3) then
           write(*,*) 'hiba in JP_dilep, id>3',id
           stop
      end if
c      if(abs(vmass-JPsi_prop(id,1)).gt.0.2) write(*,*)
c     & 'JPsidil JPsimass',vmass,JPsi_prop(id,1),id
               
      imass = JP_mass_bin_num(vmass)
c      write(*,*) 'JP_dilep',vmass,id,imass,prob
      if(imass.le.0 .or. imass.gt.JP_dlmas) return
      dvmass = JP_ddlmas
      prob=prob*JPsi_prop(id,4)/dvmass
      JP_prob_epair(id,imass)=JP_prob_epair(id,imass)
     &                 + prob
      
      en = sqrt(vmass**2+px**2+py**2+pz**2)
      qt = sqrt(px**2+py**2)
      rap = 0.5 * log((en+pz)/(en-pz))
      iqt = JP_qt_bin_num(qt)
      iy = JP_y_bin_num(rap)
c      write(*,*) 'piNdilep Vectmes', ichannel,vmass,imass,iqt,iy,prob
      if(iqt.lt.1 .or. iqt.gt.JP_qt
     &  .or. iy.lt.1 .or. iy.gt.JP_y) return
      ibin = (imass-1)*JP_qt*JP_y+(iqt-1)*JP_y+(iy-1)+1
      if(ibin.gt.JP_maxy) then
        write(*,*) 'hiba JP vec-dec ', ibin,JP_dlmas,JP_qt,JP_y,JP_f,
     &   imass,iqt,iy
      stop
      end if
c      dqt = qt_bin_size(iqt)
c      dy = y_bin_size(iy)
      volume = qt*JP_dqt*JP_dy
      JP_sig(id,ibin,nde) = JP_sig(id,ibin,nde) 
     &                     + prob/volume
      return
      end
************************************************************************
      subroutine Drell_Yan(srt,szig,beta,gamma,xxx,yyy,zzz)
c-----------------------------------------------------------------------
      implicit none
      include 'common'
c      include 'cominput'
      include 'com_cont_epair'
      real*8 szig,srt,beta(3),gamma,xxx,yyy,zzz
      real*8 dilmass,sigma,prob,qmax,q0max,qbeta
      real*8 eprim, spr
      integer ibin,ii,jj,kk,iqt,iy,iph,nde,dens_bin

c      write(*,*) 'In Drell_Yan (srt,sig) - ',srt,szig
c      if(srt.le.2.*rmass+0.0001) return
      ii    = 0
      jj    = 0

c***********************************************************************
c           density dependence
      nde = dens_bin(xxx,yyy,zzz)

c***********************************************************************      

      kk = 1
      do while (kk.le.nDYsrt .and. srt.gt.DY_InvMcross(1,kk,1))
        kk = kk + 1
      end do
c      if(kk.gt.1 .and. kk.le.nDYsrt) write(*,*)  'Drell Yan interp', kk,
c     &  DY_InvMcross(1,kk-1,1),srt,DY_InvMcross(1,kk,1)
      q0max= srt
      do ibin = 1,IR_dlmas
        dilmass = IR_qy(0,(ibin-1)*IR_qt*IR_y*IR_f+1)

        if(q0max.le.dilmass+0.001) return
        qmax = sqrt(max(0.0, q0max**2-dilmass**2))
        if(kk.gt.1 .and. kk.le.nDYsrt) then
           sigma = DY_InvMcross(3,kk-1,ibin)+
     &            (srt-DY_InvMcross(1,kk-1,ibin))/
     &     (DY_InvMcross(1,kk,ibin)-DY_InvMcross(1,kk-1,ibin)) * 
     &     (DY_InvMcross(3,kk,ibin)-DY_InvMcross(3,kk-1,ibin))
        else if (kk.eq.1) then
           sigma = srt/DY_InvMcross(1,1,ibin) * DY_InvMcross(3,1,ibin)
        else if(kk.gt.nDYsrt) then
           sigma=DY_InvMcross(3,nDYsrt,ibin)+
     &           (srt-DY_InvMcross(1,nDYsrt,ibin))/
     &     (DY_InvMcross(1,nDYsrt,ibin)-DY_InvMcross(1,nDYsrt-1,ibin)) * 
     &     (DY_InvMcross(3,nDYsrt,ibin)-DY_InvMcross(3,nDYsrt-1,ibin))
        end if

        IR_prob_epair(ch_DY,ibin) =
     &     IR_prob_epair(ch_DY,ibin)
     &     + sigma/szig
c        write(*,*) 'Drell-Yan ->',ch_DY,kk,dilmass,sigma,
c     &     IR_prob_epair(ch_DY,ibin)
        do 4302 iqt=1,nqt
        do 4301 iy =1,ny
          jj = jj + 1
        do 4300 iph=1,nf
          ii = ii + 1
          qbeta = beta(1)*qq(1,ii) + beta(2)*qq(2,ii) + beta(3)*qq(3,ii)
          eprim = gamma * (qq(4,ii) - qbeta)
          spr   = srt * (srt - 2. * eprim) + qy(0,ii)**2
          if(spr.le.0.0)    goto 4300
          IR_sig(ch_DY,jj,nde)= IR_sig(ch_DY,jj,nde) 
     &                       + 0.0
 4300   continue
 4301   continue
 4302   continue
      end do
      return
      end
************************************************************************
      subroutine OpenCharm_dilep(id,vmass,px,py,pz,prob)
c-----------------------------------------------------------------------
      implicit none
c      include 'common'
      include 'com_pert'
      include 'com_cont_epair'
      integer id       ! id=1: charged 2: neutral
      real*8 vmass     ! mass of the vector meson
      real*8 px,py,pz ! monemtum and position of vector mes
      integer ich      ! 1 nonpert, 2 pert
      integer nde, dens_bin,ichannel,ibin,imass,iqt,iy
      real*8 rap,en,qt,dvmass,volume,prob

      integer IR_mass_bin_num,IR_qt_bin_num,IR_y_bin_num,i_reg
      real*8 mass_bin_size,qt_bin_size,y_bin_size

      i_reg=2
c***********************************************************************
c           density dependence
      nde=1
c***********************************************************************      

      imass = IR_mass_bin_num(vmass)
      if(imass.le.0 .or. imass.gt.IR_dlmas) return
      dvmass = IR_ddlmas
      prob=prob*Dmes_prop(id,3)**2/dvmass
c      write(*,*) 'IR_dilep',vmass,id,imass,prob,Dmes_prop(id,3),
c     &     IR_prob_epair(ch_OC,imass)
      IR_prob_epair(ch_OC,imass)=IR_prob_epair(ch_OC,imass)
     &                 + prob
      en = sqrt(vmass**2+px**2+py**2+pz**2)
      qt = sqrt(px**2+py**2)
      rap = 0.5 * log((en+pz)/(en-pz))
      iqt = IR_qt_bin_num(qt)
      iy = IR_y_bin_num(rap)
c      write(*,*) 'OpenCharm Vectmes', ichannel,vmass,imass,iqt,iy,prob
      if(iqt.lt.1 .or. iqt.gt.IR_qt
     &  .or. iy.lt.1 .or. iy.gt.IR_y) return
      ibin = (imass-1)*IR_qt*IR_y+(iqt-1)*IR_y+(iy-1)+1
      if(ibin.gt.IR_maxy) then
        write(*,*) 'hiba IR vec-dec ',ibin,IR_dlmas,IR_qt,IR_y,IR_f,
     &   imass,iqt,iy
      stop
      end if
c      dqt = qt_bin_size(iqt)
c      dy = y_bin_size(iy)
      volume = qt*IR_dqt*IR_dy
      IR_sig(ch_OC,ibin,nde) = IR_sig(ch_OC,ibin,nde) 
     &                     + prob/volume

      return
      end

************************************************************************
      subroutine dilepout(scala,yref)
c-----------------------------------------------------------------------
      implicit none
      include 'com_cont_epair'
      include 'cominput'
      include 'common'
      include 'com_pert'
      real*8 scala, scal,yref
      real*8 dsdm(n_cont_channel,n_mass_bins),syi(n_cont_channel)
      real*8 JP_dsdm(JP_channel,JP_mass_bins)
      integer ii,jj,idm,ichannel,nde,iqt,iy,ij

*** convert the results from mb to mikrob ***
      scal   = scala * 1000.0
c      scal = 1.0
      write(*,*) 'Dilepout 1', scal, dqt, dy
      ii = 0
      jj = 0
      do 5900 idm = 1,ndlmas
        do ichannel = 1,n_cont_channel
          dsdm(ichannel,idm) = 0.0
        end do
        do 5700 iqt= 1,nqt
          do ichannel=1,n_cont_channel
            syi(ichannel) = 0.
          end do
          do 5500 iy = 1,ny
            jj = jj + 1
            ii = ii + nf
            do 5301 nde=1,maxde
c              if(iqq(0,jj,nde) .eq. 0)                        goto 5301
              do 5300 ichannel =1,n_cont_channel
                syi(ichannel) = syi(ichannel) + sig(ichannel,jj,nde)
 5300         continue
 5301       continue
 5500     continue
          do ichannel = 1,n_cont_channel
            dsdm(ichannel,idm)  =dsdm(ichannel,idm) 
     &       +syi(ichannel)  *qy(1,ii)*dy*dqt
c            write(*,*) 'qt: ', iqt, ii, qy(1,ii)
          end do
 5700   continue
        prob_cont_epair(ch_piN_prop,idm)=
     &           prob_cont_epair(ch_rho_piN_prop,idm)+
     &           prob_cont_epair(ch_ome_piN_prop,idm)+
     &           prob_cont_epair(ch_rhome_piN_prop,idm)
        dsdm(ch_piN_prop,idm)=
     &           dsdm(ch_rho_piN_prop,idm)+
     &           dsdm(ch_ome_piN_prop,idm)+
     &           dsdm(ch_rhome_piN_prop,idm)
c        write(*,*) 'dilepout piNout', ch_piN_prop,ch_rho_piN_prop,
c     &    ch_ome_piN_prop,ch_rhome_piN_prop,idm,
c     &    prob_cont_epair(ch_piN_prop,idm),
c     &    prob_cont_epair(ch_rho_piN_prop,idm),
c     &    prob_cont_epair(ch_ome_piN_prop,idm),
c     &    prob_cont_epair(ch_rhome_piN_prop,idm)
c      write(*,*) 'dilepout',nres
        do ichannel = bary_chan_offset+1,bary_chan_offset+nres
          prob_cont_epair(bary_chan_offset,idm)=
     &           prob_cont_epair(bary_chan_offset,idm)
     &          +prob_cont_epair(ichannel,idm)
          dsdm(bary_chan_offset,idm)=
     &           dsdm(bary_chan_offset,idm) + dsdm(ichannel,idm)
        end do
 5900 continue
c        stop
*
      write(isum,'(//79(''*''))')
      write(isum,'(''#c:      dilepton production''/''#''/
     &            ''#c:bombarding energy:'',f8.2,'' GeV''/
     &            ''#c:impact parameter: '',f8.2,'' fm'')') elab,b
      write(isum,'(''#c: ny:'',i5,5x,''nqt:'',i5,5x,''nphi:'',i5)')
     &                                  ny,nqt,nf
      write(isum,'(''#c:rapidity interval:'',f8.4,1x,''-'',f8.4,3x,
     &             ''lab. rapidity:'',f8.4,3x,''qtmax:'',f7.2)')
     &                 ymin,ymax,yref,qtmaxi
*
      if(idiltra.eq.1) write(isum,'(''#c: with transverse mass'')')
      write(isum,100)
      do idm = 1,ndlmas
        ii = (idm-1) * nqt * ny * nf + 1
        write(isum,101)qy(0,ii),(prob_cont_epair(ij,idm)*scal,ij=1,6),
     &    prob_cont_epair(ch_piN_prop,idm)*scal,
     &    prob_cont_epair(bary_chan_offset,idm)*scal
      end do
c      write(isum,'(/''#c: integrated over phase space'')')
c      write(isum,100)
c      do idm = 1,ndlmas
c        ii = (idm-1) * nqt * ny * nf + 1
c        write(isum,101)qy(0,ii),(dsdm(ij,idm)*scal,ij=1,6),
c     &    dsdm(ch_piN_prop,idm)*scal,
c     &    dsdm(bary_chan_offset,idm)*scal
c      end do

      write(isum,103)
      do idm = 1,ndlmas
        ii = (idm-1) * nqt * ny * nf + 1
        write(isum,101)qy(0,ii),(prob_cont_epair(ij,idm)*scal,ij=7,9),
     &    (prob_cont_epair(ij,idm)*scal,ij=11,13)
      end do
      write(isum,'(/''#c: integrated over phase space'')')
      write(isum,103)
      do idm = 1,ndlmas
        ii = (idm-1) * nqt * ny * nf + 1
        write(isum,101)qy(0,ii),(dsdm(ij,idm)*scal,ij=7,9),
     &    (dsdm(ch_piN_prop,idm)*scal,ij=11,13)
      end do
 100  format(/'#n:      ',
     &  ' pi0dec  etadec    omegadec   pn        rh0      omega    ',
     &  ' piN       ResDalitz')
 101  format(f7.3,10e10.3)
 103  format(//'#n:  mass',
     &  ' piN     piN-rho   piN-omega piN-prop  piN-rhoPr piN-omePr')
*
      ii = 0
      jj = 0
      do 6900 idm = 1,JP_dlmas
        do ichannel = 1,JP_channel
          JP_dsdm(ichannel,idm) = 0.0
        end do
        do 6700 iqt= 1,JP_qt
          do ichannel=1,JP_channel
            syi(ichannel) = 0.
          end do
          do 6500 iy = 1,JP_y
            jj = jj + 1
            ii = ii + JP_f
            do 6301 nde=1,maxde
c              if(iqq(0,jj,nde) .eq. 0)                        goto 5301
              do 6300 ichannel =1,JP_channel
                syi(ichannel) = syi(ichannel) + JP_sig(ichannel,jj,nde)
 6300         continue
 6301       continue
 6500     continue
          do ichannel = 1,JP_channel
            JP_dsdm(ichannel,idm)  =JP_dsdm(ichannel,idm) 
     &       +syi(ichannel)  *JP_qy(1,ii)*JP_dy*JP_dqt
c            write(*,*) 'qt: ', iqt, ii, qy(1,ii)
          end do
 6700   continue
 6900 continue
c        stop
*
      write(isum,'(//79(''*''))')
      write(isum,'(''#     charmonium''/''#''/
     &            ''# bombarding energy:'',f8.2,'' GeV''/
     &            ''# impact parameter: '',f8.2,'' fm''/
     &            ''# coll. broadening: '',i3,'' em branching'',i3)')
     & elab,b,icbro,i_charm_matt_dec
      write(isum,'(''#c: ny:'',i5,5x,''nqt:'',i5,5x,''nphi:'',i5)')
     &                                  JP_y,JP_qt,JP_f
      write(isum,'(''#c:rapidity interval:'',f8.4,1x,''-'',f8.4,3x,
     &             ''lab. rapidity:'',f8.4,3x,''qtmax:'',f7.2)')
     &                 JP_ymin,JP_ymax,yref,JP_qtmaxi
*
      write(*,'(/''#     charmonium''/''#''/
     &            ''# bombarding energy:'',f8.2,'' GeV''/
     &            ''# impact parameter: '',f8.2,'' fm''/
     &            ''# coll. broadening: '',i3,'' em branching'',i3)')
     & elab,b,icbro,i_charm_matt_dec
      write(*,'(''#c: ny:'',i5,5x,''nqt:'',i5,5x,''nphi:'',i5)')
     &                                  JP_y,JP_qt,JP_f
      write(*,'(''#c:rapidity interval:'',f8.4,1x,''-'',f8.4,3x,
     &             ''lab. rapidity:'',f8.4,3x,''qtmax:'',f7.2)')
     &                 JP_ymin,JP_ymax,yref,JP_qtmaxi
*
      write(isum,104)
      write(*,104)
      do idm = 1,JP_dlmas
        ii = (idm-1) * JP_qt * JP_y * JP_f + 1
        write(isum,101)JP_qy(0,ii),(JP_prob_epair(ij,idm)*scal,ij=1,3)
        write(*,101)JP_qy(0,ii),(JP_prob_epair(ij,idm)*scal,ij=1,3)
      end do
      call f77flush()

      write(isum,'(//''#c: integrated over phase space'')')
      write(*,'(//''#c: integrated over phase space'')')
      write(*,104)
      write(isum,104)
      do idm = 1,JP_dlmas
        ii = (idm-1) * JP_qt * JP_y * JP_f + 1
        write(isum,101)JP_qy(0,ii),(JP_dsdm(ij,idm)*scal,ij=1,3)
        write(*,101)JP_qy(0,ii),(JP_dsdm(ij,idm)*scal,ij=1,3)
      end do

 104  format('#  mass  J/Psi     Psi1      Psi2')
*
*
      write(isum,'(//"#",78(''*''))')
      write(isum,'(''#c: Intermediate mass region:Drell-Yan''/''#''/
     &            ''#c:bombarding energy:'',f8.2,'' GeV''/
     &            ''#c:impact parameter: '',f8.2,'' fm'')') elab,b
      write(*,'(//''#c: Intermediate mass region:Drell-Yan''/
     &            ''#c:bombarding energy:'',f8.2,'' GeV''/
     &            ''#c:impact parameter: '',f8.2,'' fm'')') elab,b
c      write(isum,'(''#c: ny:'',i5,5x,''nqt:'',i5,5x,''nphi:'',i5)')
c     &                                  JP_y,JP_qt,JP_f
c      write(isum,'(''#c:rapidity interval:'',f8.4,1x,''-'',f8.4,3x,
c     &             ''lab. rapidity:'',f8.4,3x,''qtmax:'',f7.2)')
c     &                 JP_ymin,JP_ymax,yref,JP_qtmaxi
*
      write(isum,105)
      write(*,105)
      do idm = 1,IR_dlmas
        ii = (idm-1) * IR_qt * IR_y * IR_f + 1
        write(isum,101)IR_qy(0,ii),(IR_prob_epair(ij,idm)*scal,ij=1,2)
        write(*,101)IR_qy(0,ii),(IR_prob_epair(ij,idm)*scal,ij=1,2)
      end do
      call f77flush()
 105  format('# mass  Drell-Yan OpenCharm',f8.3)

      return
      end
c----------------------------------------------------------------------c
      SUBROUTINE SPLINED(X,Y,N,YP1,YPN,Y2)
      implicit none
      integer n,nDYmax,i,k
      real*8 x,y,y2,u,yp1,ypn,sig,p,qn,un
      PARAMETER (nDYmax=200)
      DIMENSION X(nDYmax),Y(nDYmax),Y2(nDYmax),U(nDYmax)
      IF (YP1.GT..99d30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99d30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

      SUBROUTINE SPLINTD(XA,YA,Y2A,N,X,Y)
      implicit none
      integer n,klo,khi,k,nmax
      real*8 xa,ya,y2a,x,y,h,a,b
      PARAMETER (NMAX=200)
      DIMENSION XA(Nmax),YA(Nmax),Y2A(Nmax)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
