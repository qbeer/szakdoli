************************************************************************
*                                                                      *
      subroutine JPsi_pert_init
*                                                                      *
*       purpose:    calculating JPsi perturbative production from:     *
*                           N+N  collision      (entry JPsi_pert_NN)   *
*                           N+pi collision      (entry JPsi_pert_Npi)  *
*--------------------------------------------------------------------  *
*                   initialization by subroutine JPsi_pert_init        *
*--------------------------------------------------------------------  *
*         ireac: 10->  N_A + N_B  ->  N_A' + N_B' + JPsi  + X          *
*         ireac: 20->  N + pi  ->  N + JPsi                            *
*         ireac: 30->  anti_p + p  ->  anti_p' + p' + JPsi             *
*         ireac: 31->  anti_p + p  ->  pi_0 + JPsi                     *
*         ireac: 32->  anti_p + p  ->  rho_0 + JPsi                    *
*                                                                      *
*                propag     "nonpert"       abs                        *
*         i_JPsi    1           1            1        +1               *
*                   0           0            0                         *
*                propag:  5-8                                          *
*                nonpert: 3,4,7,8                                      *
*                abs:     2,4,6,8                                      *
************************************************************************

      implicit none
      include 'common'
      include 'cominput'
      include 'com_pert'
*----------------------------------------------------------------------*

c ---------------------------------------------------------------------
c ----------- variables used in the N+N collision ---------------------
c ---------------------------------------------------------------------

      real*8 prob_JPsi_scale
      real*8 mmass ! J/Psi meson mass
      real*8 pmax,pmax2 ! maximal momentum and momentum square of the J/Psi (NN --> J/Psi + NN)
      integer i1, i2 ! id numbers of nucleons
      integer id1,id2 ! type of the colliding nucleons
      integer id21,id22 ! charges of the nucleons
      integer ireac ! id of the reaction
      integer ib ! number of times the i1 nucleon has collided
      integer irun              ! index of parallel runs
      integer ii
      real*8 beta(3) ! velocity of the two nucleons in the observable frame
      real*8 gamma ! 1/sqrt(1-beta^2)
      real*8 srt,s ! cms energy of the two nucleons, s=srt**2
      real*8 xxx,yyy,zzz ! position of the collision in the observable frame
      real*8 sig0 ! cross section corresponding to the maximal impact parameter for calling this routine
      real*8 sig1, sig2, sig3, sig4 ! auxiliary cross section variables
      real*8 sigma ! cross section for NN --> J/Psi + X  
      integer srt_MeV_int ! closest integer of 'srt' in MeV
      integer srt1_l, srt1_h, srt2_l, srt2_h ! auxiliary 'srt' variables
      integer srt1_id, srt2_id ! auxiliary 'srt' variables
      real*8 alfa1, beta1, f_JPsi, a_fac ! parameters of the N + N --> J/Psi + X cross section (Eq. 1 of O. Linnyk et al. Nucl. Phys. A 786 (2007) 183)
      parameter (alfa1=10.0, beta1=0.775, f_JPsi=0.581, a_fac=0.2) ! a_fac = 0.2 microbarn, sig0 in milibarn


      integer size_of_run, i_pert_min, i_pert, nt, iz0
      real*8 probab, fact_pauli ! probability of JPsi creation
      real*8 prob_min, zz0

      real*8 rn ! random number generation
      real*8 ranp ! random variable for JPsi momentum generation
      real*8 p_abs_JPsi, E_JPsi ! momentum and energy of JPsi in  the 3 particle c.m.s.
      real*8 xx,yy,zz,rr ! for random unit vector generation
      real*8 xxn,yyn,zzn,rrn, pnuc2 ! for random unit vector generation for nucleons
      real*8 p_x_JPsi, p_y_JPsi, p_z_JPsi ! JPsi momentum components
      real*8 s_prime ! see later
      real*8 E_nucl_prime,p_nucl_prime ! momentum and energy of the nucleon in the 2 nucleon subsystem
      real*8 factn1,factn2 ! transformation factors for the N_A' momentum calculation
      real*8 srtmin, JPsiMass
      real*8 vrel, rself,iself,sgamma,MassJPsi,mass2,j0,j1,j2,j3
      real*8  p3(3),p4(3)             ! momentum components of N_A'
      real*8  e3, e4, E_mes, pbeta ! energy of N_A', pbeta = p_vec*beta_vec
      real*8 traf,transf ! variables for Lorentz transformation
      real*8 phase ! from Pauli exclusion principle
      integer ntag ! flag which tells if phase-space is pauli-blocked (-1: blocked)
      integer ix,iy,iz ! integer coordinates of the collision  
      real*8 denst ! density
      integer bin_dens,bin_time ! bins for the histogram
      integer p_dbb(0:999,0:999) ! for density histogram 

c ---------------------------------------------------------------------
c ------ additional variable used in N+pi collisions ------------------
c ---------------------------------------------------------------------      

      integer iz1,iz2 ! charges of the meson and nuleon

      real*8 E_Npi, px_Npi, py_Npi, pz_Npi ! Total energy and momentum of the N+pi system
      real*8 p_abs_JPsi2
      real*8 betax, betay, betaz

      real*8 gamma1, b_fac ! parameters of the N + pi --> J/Psi + N cross section (Eq. 1 of O. Linnyk et al. Nucl. Phys. A 786 (2007) 183)
      parameter (gamma1=7.3, b_fac=1.24) ! b_fac = 1.24 microbarn, sig0 in milibarn

      return

*----------------------------------------------------------------------*
      entry JPsi_pert_pbarN(i1,i2,id1,id2,beta,srt,xxx,yyy,zzz,sig0,
     &  irun)

*       variables:                                                     *
*  CHANGED!!     ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*     cross-sections are taken from sibirtsev   or  Chung PL B401(97)1 *
*----------------------------------------------------------------------*

      srtmin = JPsi_prop(1,1)-0.5
      if(srt.le.srtmin) return

c      write(*,*) 'in JPsi_pert_cr',i1,i2,id1,id2,srt
      prob_JPsi_scale = 1./JPsi_scale_factor
      if (rn(iseed) .ge. prob_JPsi_scale)                         return
cc ----------------------------------------------------------------------
      size_of_run = (max_pert/num) ! num: number of parallel runs, size_of_run: array size for one parallel run (copy)

      ix = nint(xxx)            ! nint: closest integer
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx.and.iabs(iy).le.maxx.and.iabs(iz).le.maxz)then
        j0=rhob_4(0,ix,iy,iz)  
        j1=rhob_4(1,ix,iy,iz)  
        j2=rhob_4(2,ix,iy,iz)  
        j3=rhob_4(3,ix,iy,iz)  
        denst = sqrt(j0**2-j1**2-j2**2-j3**2)
c        write(*,*) 'JPcreate dens:', j0,denst,rhb(ix,iy,iz)
        denst = rhb(ix,iy,iz)
      else
        denst = 0.0
        j0 = 0.0
      end if

      if(denst/rho0.gt.1.5) write(*,*) 'JPsi_create dens',
     &   denst/rho0,ix,iy,iz

      s     = srt**2
      srt_MeV_int = nint(srt*1000.) ! closest integer of sqrt(s) in MeV
ccc   loop over the 3 charmonium states
      vrel = 0.0
      gamma = 1.0/sqrt(1.0-(beta(1)**2+beta(2)**2+beta(3)**2)) 
      traf  = gamma / (gamma+1.0)

      do ii=1,3                 ! do fo JPsi types
         
c           direction of the charmonium mom. xx/rr, yy/rr, zz/rr
 51     xx       = 1. - 2. * rn(iseed) ! choosing a unit vector for the momentum direction   
        yy       = 1. - 2. * rn(iseed)
        zz       = 1. - 2. * rn(iseed)
        rr       = sqrt( xx**2 + yy**2 + zz**2 )
        if((rr .lt. 0.001) .or. (rr .gt. 1.) ) goto 51

        pmax2=.25*(s-(JPsi_prop(ii,1)+pmass)**2)*
     &       (s-(JPsi_prop(ii,1)-pmass)**2)/s
        write(*,*) 'JPsicreate pmax2', pmax2
        if(pmax2.lt.0.0001) then
          pmax = 0.0 
          vrel = 0.0
        else
          pmax = sqrt(pmax2)
          p_x_JPsi = pmax * xx / rr ! JPsi momentum components for vrel       
          p_y_JPsi = pmax * yy / rr
          p_z_JPsi = pmax * zz / rr
          E_JPsi   = sqrt(pmax2 + JPsi_prop(ii,1)**2) ! for vrel JPsi energy in the ppbar c.m. system
          call lorentz(-beta(1),-beta(2),-beta(3),
     &       p_x_JPsi,p_y_JPsi,p_z_JPsi,E_JPsi)
        
          if (denst.lt.1.0e-3) then
            denst=0.
            vrel = 0.
          else
            call f77flush()
c         write(*,*) 'lorentz call gradupi',betacm(1),betacm(2),betacm(3),
c     &     j1,j2,j3,j0
            if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba gradupi lorentz, mass<0",j0,j1,j2,j3
              stop
            end if
            call lorentz(p_x_JPsi/E_JPsi,p_y_JPsi/E_JPsi,p_z_JPsi/E_JPsi
     &          ,j1,j2,j3,j0)
            vrel = sqrt(j1**2+j2**2+j3**2)/j0
          end if
        end if

        call self_JPsi(ii,vrel, denst,rself,iself,sgamma)

c     size of the charmonium mom. ranp*pmax  
 50     continue
        ranp     = rn(iseed)
        xx = ranp**2 * sqrt(1.00001-ranp**2) / 0.385 ! p^2*sqrt(1-p^2) is the momentum distribution 
c          0.385 is the normalization factor, the maximal value of xx is 1.
        if(xx .lt. rn(iseed)) goto 50


        sig1 = 0.0

        goto 95
c
ccc                 pbarp -> pbarp charmonium
c
        ireac = 30
        mass2 = 2.*rmass
c        write(*,*) 'JP create1',srt,mass2+JPsi_prop(ii,1)+
c     &       0.5*rself/JPsi_prop(ii,1)+2.*iself/JPsi_prop(ii,1)
        if(srt.le.mass2+JPsi_prop(ii,1)+rself/(2.0*JPsi_prop(ii,1))+
     &       2.*iself/JPsi_prop(ii,1)) goto 95
        
        mmass = JPsiMass(ii,srt,mass2,rself,iself,iseed)
        pmax2=.25*(s-(mass2+mmass)**2)*(s-(mass2-mmass)**2)/s
        if(pmax2 .le. 0.0)                      goto 95
c     pmax  = the maximal J/Psi momentum
        pmax  = sqrt(pmax2)

c--------------------------------------------------
c------ Creation probability ----------------------
c--------------------------------------------------
        fact_pauli = 1.0

c---------------------------------------------
c---------------JPsi momentum ----------------
c---------------------------------------------

        p_abs_JPsi = pmax*ranp ! outgoing JPsi momentum absolute value  
        E_JPsi   = sqrt(p_abs_JPsi**2 + mmass**2) ! JPsi energy in the 3 particle c.m. system
        s_prime  = s - 2.0*srt*E_JPsi+ mmass**2 ! s'=(p_N_A' + p_N_B')^2 : two nucleon s in the 3 particle c. m. system
        E_nucl_prime = 0.5* sqrt(s_prime) ! energy of the nucleons in the two nucleon c. m. system
        p_nucl_prime = sqrt(E_nucl_prime**2-rmass**2) ! absolute value of the nucleon momentum in the two nucleon c. m. system 
        p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
        p_y_JPsi = p_abs_JPsi * yy / rr
        p_z_JPsi = p_abs_JPsi * zz / rr

c----------------------------------------------
c--------------- N_A' momentum ----------------
c----------------------------------------------

 52     xxn       = 1. - 2. * rn(iseed)
        yyn       = 1. - 2. * rn(iseed)
        zzn       = 1. - 2. * rn(iseed)
        rrn       = sqrt( xxn**2 + yyn**2 + zzn**2 )
        if((rrn .lt. 0.001) .or. (rrn .gt. 1.) ) goto 52
        factn1   = s_prime + sqrt(s_prime)*(srt-E_JPsi)
        factn2   = p_nucl_prime*(p_x_JPsi*xxn+p_y_JPsi*yyn+p_z_JPsi*zzn)
     $        /factn1/rrn-E_nucl_prime/sqrt(s_prime)
*     p3:                    nucleon momentum in N_A - N_B c.m. system (it is the same as the 3 particle c.m.s.)
c            to the 3 particle c.m.s (JPsi + N_A' + N_B')
        p3(1)    = p_nucl_prime*xxn/rrn + factn2*p_x_JPsi ! Lorentz transformation from the two nucleon (A' + B') system
        p3(2)    = p_nucl_prime*yyn/rrn + factn2*p_y_JPsi
        p3(3)    = p_nucl_prime*zzn/rrn + factn2*p_z_JPsi
        e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*     p3:                    nucleon momentum in observable system
         
        pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
        transf = gamma * (traf * pbeta + e3)
        p3(1)  = p3(1) + beta(1) * transf
        p3(2)  = p3(2) + beta(2) * transf
        p3(3)  = p3(3) + beta(3) * transf
        if(ipauli.eq.1)
     &     call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))

        fact_pauli = (1.0-phase) ! modification of the probability according to the Pauli exclusion principle of N_A'

        sig1 = sig1 * fact_pauli

 95     continue
        sig2 = 0.0
        sig3 = 0.0
        sig4 = 0.0

c     p antip ->  J/Psi + pi_0
        if(id1*id2.eq.-1 .and. id1+id2.eq.0) then
          if(ii.eq.1) then
            srt2_l = nint(sig_ppbar_JPsi_pi0(1,1)*1000) ! smallest integer 'srt' value for the process ireac=30 in MeV   
            srt2_h = nint(sig_ppbar_JPsi_pi0(NoL_2,1)*1000) ! largest integer 'srt' value for the process ireac=30 in MeV   
            if((srt_MeV_int .ge. srt2_l) .and. 
     $        (srt_MeV_int .le. srt2_h)) then 
              srt2_id = srt_MeV_int - srt2_l + 1
              sig2 = sig_ppbar_JPsi_pi0(srt2_id,2)
            elseif(srt_MeV_int .lt. srt2_l) then
c              write(*,*) 'srt below threshold of J/Psi + pi_0', srt
              sig2 = 0
            elseif(srt_MeV_int .gt. srt2_h) then
              write(*,*) 'srt above largest existing value in table', 
     $           '--> sigma(srt_max) value is used', srt
              sig2 = sig_ppbar_JPsi_pi0(NoL_2,2)
            endif
          else if(ii.eq.2) then
            srt2_l = nint(sig_ppbar_Psi1_pi0(1,1)*1000) ! smallest integer 'srt' value for the process ireac=30 in MeV   
            srt2_h = nint(sig_ppbar_Psi1_pi0(NoL_4,1)*1000) ! largest integer 'srt' value for the process ireac=30 in MeV   
            if((srt_MeV_int .ge. srt2_l) .and. 
     $        (srt_MeV_int .le. srt2_h)) then 
              srt2_id = srt_MeV_int - srt2_l + 1
              sig2 = sig_ppbar_Psi1_pi0(srt2_id,2)
            else if(srt_MeV_int .lt. srt2_l) then
c              write(*,*) 'srt below threshold of Psi1 + pi_0', srt
              sig2 = 0
            else if(srt_MeV_int .gt. srt2_h) then
              write(*,*) 'srt above largest existing value in table', 
     $           '--> sigma(srt_max) value is used', srt
              sig2 = sig_ppbar_Psi1_pi0(NoL_4,2)
            end if
          else if(ii.eq.3) then
            srt2_l = nint(sig_ppbar_Psi2_pi0(1,1)*1000) ! smallest integer 'srt' value for the process ireac=30 in MeV   
            srt2_h = nint(sig_ppbar_Psi2_pi0(NoL_6,1)*1000) ! largest integer 'srt' value for the process ireac=30 in MeV
            if((srt_MeV_int .ge. srt2_l) .and. 
     $        (srt_MeV_int .le. srt2_h)) then 
              srt2_id = srt_MeV_int - srt2_l + 1
              sig2 = sig_ppbar_Psi2_pi0(srt2_id,2)
            else if(srt_MeV_int .lt. srt2_l) then
c              write(*,*) 'srt below threshold of Psi2 + pi_0', srt
              sig2 = 0
            else if(srt_MeV_int .gt. srt2_h) then
              write(*,*) 'srt above largest existing value in table', 
     $           '--> sigma(srt_max) value is used', srt
              sig2 = sig_ppbar_Psi2_pi0(NoL_6,2)
            end if
          end if

c     p antip ->  J/Psi + rho_0
          if(ii.eq.1) then
            srt2_l = nint(sig_ppbar_JPsi_rho0(1,1)*1000) ! smallest integer 'srt' value for the process ireac=30 in MeV   
            srt2_h = nint(sig_ppbar_JPsi_rho0(NoL_3,1)*1000) ! largest integer 'srt' value for the process ireac=30 in MeV   
            if((srt_MeV_int .ge. srt2_l) .and. 
     $        (srt_MeV_int .le. srt2_h)) then 
              srt2_id = srt_MeV_int - srt2_l + 1
              sig3 = sig_ppbar_JPsi_rho0(srt2_id,2)
            elseif(srt_MeV_int .lt. srt2_l) then
c              write(*,*) 'srt below threshold of J/Psi + rho_0', srt,ii
              sig3 = 0
            elseif(srt_MeV_int .gt. srt2_h) then
              write(*,*) 'srt above largest existing value in table', 
     $           '--> sigma(srt_max) value is used', srt
              sig3 = sig_ppbar_JPsi_rho0(NoL_3,2)
            endif
          else if(ii.eq.2) then
            srt2_l = nint(sig_ppbar_Psi1_rho0(1,1)*1000) ! smallest integer 'srt' value for the process ireac=30 in MeV   
            srt2_h = nint(sig_ppbar_Psi1_rho0(NoL_5,1)*1000) ! largest integer 'srt' value for the process ireac=30 in MeV   
            if((srt_MeV_int .ge. srt2_l) .and. 
     $        (srt_MeV_int .le. srt2_h)) then 
              srt2_id = srt_MeV_int - srt2_l + 1
              sig3 = sig_ppbar_Psi1_rho0(srt2_id,2)
            else if(srt_MeV_int .lt. srt2_l) then
c              write(*,*) 'srt below threshold of Psi1 + rho_0', srt
              sig3 = 0
            else if(srt_MeV_int .gt. srt2_h) then
              write(*,*) 'srt above largest existing value in table', 
     $           '--> sigma(srt_max) value is used', srt
              sig3 = sig_ppbar_Psi1_rho0(NoL_5,2)
            end if
          else if(ii.eq.3) then
            srt2_l = nint(sig_ppbar_Psi2_rho0(1,1)*1000) ! smallest integer 'srt' value for the process ireac=30 in MeV   
            srt2_h = nint(sig_ppbar_Psi2_rho0(NoL_7,1)*1000) ! largest integer 'srt' value for the process ireac=30 in MeV
            if((srt_MeV_int .ge. srt2_l) .and. 
     $        (srt_MeV_int .le. srt2_h)) then 
              srt2_id = srt_MeV_int - srt2_l + 1
              sig3 = sig_ppbar_Psi2_rho0(srt2_id,2)
            else if(srt_MeV_int .lt. srt2_l) then
c              write(*,*) 'srt below threshold of Psi2 + rho_0', srt
              sig3 = 0
            else if(srt_MeV_int .gt. srt2_h) then
              write(*,*) 'srt above largest existing value in table', 
     $           '--> sigma(srt_max) value is used', srt
              sig3 = sig_ppbar_Psi2_rho0(NoL_7,2)
            end if
          end if
c             p + antip -> J/Psi
          if(abs(srt-(JPsi_prop(ii,1)+rself/(2.0*JPsi_prop(ii,1))))
     &     .lt.-iself/JPsi_prop(ii,1)) then
            pnuc2 = 0.25*s-rmass**2 
          sig4=120.*pi/pnuc2*hbc**2*s*JPsi_prop(ii,2)**2*JPsi_prop(ii,3)
     &          /((s-(JPsi_prop(ii,1)**2+rself)**2)+iself**2)
            if(abs(srt-JPsi_prop(ii,1)).gt.0.5)write(*,*)'JPsi_cr pbarp'
     &         ,ii,srt,rself/(2.0*JPsi_prop(ii,1)),iself/JPsi_prop(ii,1)
          end if
        end if  ! end of the if(ireac.eq.30)
        sig1=0.0
        sigma = sig1 + sig2 + sig3 + sig4
        if(sigma.lt.1.e-9) return
c        write(*,*) 'JPcreate4',sig1,sig2,sig3,sig4,ii
        xxn = rn(iseed)
        if(xxn .le. sig1/sigma) return
        if(xxn .gt. sig1/sigma .and. xxn .lt. (sig1+sig2)/sigma) then
c     p+pbar -> JPsi + pi0
          ireac = 31
          mass2 = pmass
          mmass=JPsiMass(ii,srt,mass2,rself,iself,iseed)
c          if(abs(mmass-JPsi_prop(ii,1)).gt.0.5) write(*,*)
c     &      'JPsiCr JPsimass',mmass,JPsi_prop(ii,1),ii,rself,iself,srt
          p_abs_JPsi2=.25*(s-(mmass+mass2)**2)*(s-(mmass-mass2)**2)/s
          if(p_abs_JPsi2 .le. 0.0)    p_abs_JPsi2 = 0.0                       
c     p_abs_JPsi : momentum of the J/Psi (process: p + pbar -> J/Psi + pi_0)
          p_abs_JPsi = sqrt(p_abs_JPsi2)
          E_JPsi = sqrt(p_abs_JPsi2 + mmass**2)
          p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
          p_y_JPsi = p_abs_JPsi * yy / rr
          p_z_JPsi = p_abs_JPsi * zz / rr
          if(i_JPsi.eq.3.or.i_JPsi.eq.4.or.i_JPsi.eq.7.or.i_JPsi.eq.8)
     &      then
c            irup = (i1-1)/maxb
            ib = (irun-1)*maxp
  10        ib=ib+1
            if(ipi(1,ib).ne.0 .and. ib .lt.(irun+1)*maxp  ) goto 10
            rpi(1,ib)= r(1,i1)
            rpi(2,ib)= r(2,i1)
            rpi(3,ib)= r(3,i1)
            pbeta=-(beta(1)*p_x_JPsi+beta(2)*p_y_JPsi+beta(3)*p_z_JPsi)
            E_mes = sqrt(p_abs_JPsi2 + mass2**2)
            transf = gamma * (traf * pbeta + E_mes)

            ppi(1,ib)= -p_x_JPsi + beta(1) * transf
            ppi(2,ib)= -p_y_JPsi + beta(2) * transf
            ppi(3,ib)= -p_z_JPsi + beta(3) * transf
            epi(ib)  = mass2
            rpie(4,ib)= denst
            rpie(5,ib)= mass2
            rpie(6,ib)= time
            mpot(ib) = 0.0
            ipi(1,ib)= 1
            ipi(2,ib)= 0
            ipi(3,ib)= i1
            ipi(5,ib)= -10
            ipi(6,ib)= id(6,i1)
            ipi(7,ib)= id(4,i1)
            ipi(8,ib)= 0

            id(1,i1) = 0
            id(1,i2) = 0
          end if
        else if(xxn.gt.(sig1+sig2)/sigma .and. xxn.le.(sig1+sig2+sig3)
     &     /sigma) then
c     p+pbar -> JPsi + rho0
          ireac = 32
          mass2 = romas
          mmass= JPsiMass(ii,srt,mass2,rself,iself,iseed)
c          if(abs(mmass-JPsi_prop(ii,1)).gt.0.3) write(*,*)
c     &      'JPsiC2 JPsimass',mmass,JPsi_prop(ii,1),ii,rself,iself,srt
          p_abs_JPsi2=.25*(s-(mmass+mass2)**2)*(s-(mmass-mass2)**2)/s
          if(p_abs_JPsi2 .le. 0.0)   p_abs_JPsi2 = 0.0                
c     p_abs_JPsi : momentum of the J/Psi (process: p + pbar -> J/Psi + rho_0)
          p_abs_JPsi = sqrt(p_abs_JPsi2)
          E_JPsi = sqrt(p_abs_JPsi2 + mmass**2)
          p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
          p_y_JPsi = p_abs_JPsi * yy / rr
          p_z_JPsi = p_abs_JPsi * zz / rr
          if(i_JPsi.eq.3.or.i_JPsi.eq.4.or.i_JPsi.eq.7.or.i_JPsi.eq.8)
     &     then
            ib = (irun-1)*maxp
  20        ib=ib+1
            if(ipi(1,ib).ne.0 .and. ib .lt.(irun+1)*maxp  ) goto 20
            rpi(1,ib)= r(1,i1)
            rpi(2,ib)= r(2,i1)
            rpi(3,ib)= r(3,i1)
            pbeta=-(beta(1)*p_x_JPsi+beta(2)*p_y_JPsi+beta(3)*p_z_JPsi)
            E_mes = sqrt(p_abs_JPsi2 + mass2**2)
            transf = gamma * (traf * pbeta + E_mes)

            ppi(1,ib)= -p_x_JPsi + beta(1) * transf
            ppi(2,ib)= -p_y_JPsi + beta(2) * transf
            ppi(3,ib)= -p_z_JPsi + beta(3) * transf
            epi(ib)  = mass2
            rpie(4,ib)= denst
            rpie(5,ib)= mass2
            rpie(6,ib)= time
            mpot(ib) = 0.0
            ipi(1,ib)= 3
            ipi(2,ib)= 0
            ipi(3,ib)= i1
            ipi(5,ib)= -10
            ipi(6,ib)= id(6,i1)
            ipi(7,ib)= id(4,i1)
            ipi(8,ib)= 0

            id(1,i1) = 0
            id(1,i2) = 0
          end if
        else if(xxn .gt. (sig1+sig2+sig3)/sigma) then

c     p+pbar -> JPsi
          ireac = 35
          mmass= srt
          if(abs(mmass-JPsi_prop(ii,1)).gt.0.3) write(*,*)
     &      'JPsiC3 JPsimass',mmass,JPsi_prop(ii,1),ii,rself,iself,srt
c     p_abs_JPsi : momentum of the J/Psi (process: p + pbar -> J/Psi
          p_abs_JPsi = 0.0
          E_JPsi = mmass
          p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
          p_y_JPsi = p_abs_JPsi * yy / rr
          p_z_JPsi = p_abs_JPsi * zz / rr
          if(i_JPsi.eq.3.or.i_JPsi.eq.4.or.i_JPsi.eq.7.or.i_JPsi.eq.8)
     &     then
            id(1,i1) = 0
            id(1,i2) = 0
          end if
        else if(xxn .le. sig1/sigma) then
          if(i_JPsi.eq.3.or.i_JPsi.eq.4.or.i_JPsi.eq.7.or.i_JPsi.eq.8)
     &     then
            p(1,i1) = p3(1)
            p(2,i1) = p3(2)
            p(3,i1) = p3(3)
            p(1,i2) = p4(1)
            p(2,i2) = p4(2)
            p(3,i2) = p4(3)
          end if
        end if
        probab = sigma/sig0*JPsi_scale_factor ! probability of producing a J/Psi
        
c----------------------------------------------
c-- Find a place for the created particle -----
c----------------------------------------------
        
        i_pert_min = 0
        prob_min  = 1000.0
        do i_pert = (irun-1) * size_of_run + 1, irun * size_of_run ! i_pert: running index of nx_pert array containing perturbativ particles  
c                                                                   ! (irun-1) * size_of_run + 1: beginning index of current run in nx_pert array 
c          write(*,*)'JPsicreate2',ii,i_pert,id_JPsi(ii),irun,size_of_run
          if (nx_pert(id_JPsi(ii),0,i_pert) .eq. 0) goto 54 ! true: found an empty slot for the just created J/Psi
          if (prob_min .gt. p_pert(id_JPsi(ii),4,i_pert)) then
            i_pert_min = i_pert
            prob_min  = p_pert(id_JPsi(ii),4,i_pert)
          endif
        enddo
        i_pert = i_pert_min
        if (probab .lt. prob_min) goto 98
 54     continue

*   J/Psi momentum in observable system
        pbeta  = beta(1)*p_x_JPsi + beta(2)*p_y_JPsi + beta(3)*p_z_JPsi
        transf = gamma * (traf * pbeta + E_JPsi)

        p_pert(id_JPsi(ii),0,i_pert)= mmass
        p_pert(id_JPsi(ii),1,i_pert)= p_x_JPsi + beta(1) * transf
        p_pert(id_JPsi(ii),2,i_pert)= p_y_JPsi + beta(2) * transf
        p_pert(id_JPsi(ii),3,i_pert)= p_z_JPsi + beta(3) * transf
        p_pert(id_JPsi(ii),4,i_pert)= probab
c        if(i_pert.eq.4012) write(*,*) 'JP_4012create1 p',ireac,
c     &   p_pert(id_JPsi(ii),1,i_pert),p_pert(id_JPsi(ii),2,i_pert),
c     &   p_pert(id_JPsi(ii),3,i_pert),beta(1),beta(2),beta(3),E_JPsi,
c     &   p_x_JPsi,p_y_JPsi,p_z_JPsi,transf,p_abs_JPsi,gamma,traf,pbeta   
        r_pert(id_JPsi(ii),1,i_pert)= xxx
        r_pert(id_JPsi(ii),2,i_pert)= yyy
        r_pert(id_JPsi(ii),3,i_pert)= zzz
        r_pert(id_JPsi(ii),4,i_pert)= denst/rho0 ! baryon density at the point of the J/Psi creation
        r_pert(id_JPsi(ii),5,i_pert)= time ! time of the J/Psi creation

        nx_pert(id_JPsi(ii),0,i_pert) = ii ! iitype JPsi at the i_pert position
        nx_pert(id_JPsi(ii),1,i_pert) = 1 ! JPsi exists at the i_pert position
        nx_pert(id_JPsi(ii),2,i_pert) = ireac ! reaction type
        nx_pert(id_JPsi(ii),3,i_pert) = i1 ! which nucleon, N_A
        nx_pert(id_JPsi(ii),4,i_pert) = i2 ! which nucleon, N_B
        nx_pert(id_JPsi(ii),5,i_pert) = 1 ! ????????????????????????
        nx_pert(id_JPsi(ii),6,i_pert) = id1 ! type of N_A
        nx_pert(id_JPsi(ii),7,i_pert) = id2 ! type of N_B
        nx_pert(id_JPsi(ii),8,i_pert) = id21 ! charge of N_A
        nx_pert(id_JPsi(ii),9,i_pert) = id22 ! charge of N_B
        mass_evol(1,1,id_JPsi(ii),i_pert) = mmass
        mass_evol(2,1,id_JPsi(ii),i_pert) = sqrt(mmass**2+
     &      p_pert(id_JPsi(ii),1,i_pert)**2+
     &      p_pert(id_JPsi(ii),2,i_pert)**2+
     &      p_pert(id_JPsi(ii),3,i_pert)**2)
        mass_evol(3,1,id_JPsi(ii),i_pert) = denst/rho0
c        if(denst.lt.0.01) write(*,*) 'JPsi_create2 dens',mmass,time,ii
c     &        ,i_pert,denst,ix,iy,iz
c     numprodd = numprodd + 1
        zz0 = sqrt(max(radius**2-b**2,0.0))
        iz0 = nint((zzz+zz0)/0.1)+10
        iz0 = max(0,iz0)
        iz0 = min(iz0,200)
c        iz0 = min(max(0,nint((zzz+zz0)/0.1)+10.),200)
        JPsi_init_pos(1,iz0)= JPsi_init_pos(1,iz0) + 1.0
        iz0 = nint(denst/rho0/0.1)
        iz0 = min(iz0,200)
c        iz0 = min(nint(denst/rho0/0.1),200)
        JPsi_init_pos(2,iz0)= JPsi_init_pos(2,iz0) + 1.0
        
 98     continue
c      write(20,*) "#J/Psi prod: ",i_pert,i1,i2,irun,numprodd
c-----------------------------------------------------------------------
c----------------------------------------------
c-- Find a place for the created particle -----
c----------------------------------------------

c-----------------------------------------------------------------------
      end do ! end of ii loop
      return

*************************************************************************

      entry JPsi_pert_Npi(i1,i2,id1,id2,iz1,iz2,E_Npi,px_Npi,py_Npi,
     $     pz_Npi,srt,xxx,yyy,zzz,sig0,irun)
*    ireac = 20: N + pi -> N + J/Psi                                    *
*    variables:    1 = pion       2 = baryon                            *
*    id1: meson type (1 -> pion), id2: baryon type, i2: baryon id       *
*    iz1: charge of the pion, iz2: charge of the baryon                 *
*    E_Npi, px_Npi, py_Npi, pz_Npi- energy, momentum of the pi+baryon sys.(real,input) *
*    srt     - cms energy (mass of the pi+baryon sys.(real,input)       *
*    xxx...zzz - coordinates of the event                               *
*-----------------------------------------------------------------------*

c      write(*,*) 'in JPsi_pert_Npi',srt

      if(iz1+iz2 .gt. 1)                                          return ! Q_pi + Q_N = {0,1}
      if(iz1+iz2 .lt. 0)                                          return ! Q_pi + Q_N = {0,1}
      if(id1.ne.1)                                                return ! id1 = 1 refers to pion
      if(id2.gt.nres+1)                                           return ! in this way baryon can be a nucleon or a nucleon resonance

      prob_JPsi_scale = 1./JPsi_scale_factor
      if (rn(iseed) .ge. prob_JPsi_scale)                         return

      if(srt-rmass .lt. JPsi_prop(1,1)-0.2)                       return
      betax = px_Npi / E_Npi
      betay = py_Npi / E_Npi
      betaz = pz_Npi / E_Npi
      gamma = E_Npi / srt

      ix = nint(xxx)            ! nint: closest integer
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx.and.iabs(iy).le.maxx.and.iabs(iz).le.maxz)then
        j0=rhob_4(0,ix,iy,iz)  
        j1=rhob_4(1,ix,iy,iz)  
        j2=rhob_4(2,ix,iy,iz)  
        j3=rhob_4(3,ix,iy,iz)  
        denst = sqrt(j0**2-j1**2-j2**2-j3**2)
c        write(*,*) 'JPcreate dens:', j0,denst,rhb(ix,iy,iz)
        denst = rhb(ix,iy,iz)
      else
        denst = 0.0
        j0 = 0.0
      end if

c           direction of the charmonium mom. xx/rr, yy/rr, zz/rr
 72   continue
      xx       = 1. - 2. * rn(iseed) ! unit vector for the momentum direction
      yy       = 1. - 2. * rn(iseed)
      zz       = 1. - 2. * rn(iseed)
      rr       = sqrt( xx**2 + yy**2 + zz**2 )
      if((rr .lt. 0.001) .or. (rr .gt. 1.) ) goto 72

      s     = srt**2
      mmass = JPsi_prop(1,1)

      p_abs_JPsi2=.25*(s-(rmass+mmass)**2)*(s-(rmass-mmass)**2)/s
c     p_abs_JPsi : momentum of the J/Psi (process: N + pi -> R -> N + J/Psi)
      p_abs_JPsi = sqrt(max(0.0000001,p_abs_JPsi2))
      E_JPsi = sqrt(p_abs_JPsi2 + mmass**2)

      p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
      p_y_JPsi = p_abs_JPsi * yy / rr
      p_z_JPsi = p_abs_JPsi * zz / rr

      call lorentz(-betax,-betay,-betaz,
     &       p_x_JPsi,p_y_JPsi,p_z_JPsi,E_JPsi)
        
      if (denst.lt.1.0e-3) then
        denst=0.
        vrel = 0.
      else
        call f77flush()
c         write(*,*) 'lorentz call gradupi',betacm(1),betacm(2),betacm(3),
c     &     j1,j2,j3,j0
        if(j1**2+j2**2+j3**2.gt.j0**2) then
           write(*,*) "hiba JPsi piN, mass<0",j0,j1,j2,j3
            stop
        end if
        call lorentz(p_x_JPsi/E_JPsi,p_y_JPsi/E_JPsi,p_z_JPsi/E_JPsi,
     &      j1,j2,j3,j0)
          vrel = sqrt(j1**2+j2**2+j3**2)/j0
        end if
      vrel = 0.0
      
      call self_JPsi(1,vrel, denst,rself,iself,sgamma) ! calculates for 3 ccbar
      mmass = JPsiMass(1,srt,rmass,rself,iself,iseed)
      p_abs_JPsi2=.25*(s-(rmass+mmass)**2)*(s-(rmass-mmass)**2)/s
      if(p_abs_JPsi2 .le. 0.0)                                    return
c     p_abs_JPsi : momentum of the J/Psi (process: N + pi -> R -> N + J/Psi)
      p_abs_JPsi = sqrt(max(0.0000001,p_abs_JPsi2))
      E_JPsi = sqrt(p_abs_JPsi2 + mmass**2)

      p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
      p_y_JPsi = p_abs_JPsi * yy / rr
      p_z_JPsi = p_abs_JPsi * zz / rr
c      if(abs(mmass-JPsi_prop(1,1)).gt.0.5) write(*,*)
c     & 'JPsiC4 JPsimass',srt,JPsi_prop(1,1),rself,iself

c--------------------------------------------------
c------ Creation probability ----------------------
c--------------------------------------------------

      if((id2 .eq. 1))          ! meson was already checked
     $     ireac = 20           ! N + pi

      if(ireac .eq. 20) then
         sigma = f_JPsi*b_fac*(1-mmass/srt)**gamma1/1000. ! sigma in milibarn; Eq. 1 of Nucl. Phys. A 786 (2007) 183
      endif

      probab = sigma/sig0*JPsi_scale_factor

c----------------------------------------------
c-- Find a place for the created particle -----
c----------------------------------------------


      size_of_run = (max_pert/num) ! num: number of parallel runs, size_of_run: array size for one parallel run (copy)

      i_pert_min = 0
      prob_min  = 1000.0
      do    i_pert = (irun-1) * size_of_run + 1, irun * size_of_run ! i_pert: running index of nx_pert array containing perturbativ particles  
c                                                                   ! (irun-1) * size_of_run + 1: beginning index of current run in nx_pert array 
        if (nx_pert(id_JPsi(1),0,i_pert) .eq. 0) goto 71 ! true: found an empty slot for the just created J/Psi
        if (prob_min .gt. p_pert(id_JPsi(1),4,i_pert)) then
           i_pert_min = i_pert
           prob_min  = p_pert(id_JPsi(1),4,i_pert)
        endif
      enddo
      i_pert = i_pert_min
      if (probab .lt. prob_min) goto 101
 71   continue

c----------------------------------------------
c--------------- JPsi momentum ----------------
c----------------------------------------------

c     J/Psi momentum in observable system

      pbeta  = betax*p_x_JPsi + betay*p_y_JPsi + betaz*p_z_JPsi
      transf = gamma * (gamma / (gamma + 1.0) * pbeta + E_JPsi) 

      p_pert(id_JPsi(1),0,i_pert)= mmass
      p_pert(id_JPsi(1),1,i_pert)= p_x_JPsi + betax * transf
      p_pert(id_JPsi(1),2,i_pert)= p_y_JPsi + betay * transf
      p_pert(id_JPsi(1),3,i_pert)= p_z_JPsi + betaz * transf
c        if(i_pert.eq.4012) write(*,*) 'JP_4012create2 p',
c     &   p_pert(id_JPsi(1),1,i_pert),p_pert(id_JPsi(1),2,i_pert),
c     &   p_pert(id_JPsi(1),3,i_pert),p_x_JPsi,betax,transf,E_JPsi

c----------------------------------------------
c--------------- Nucleon  momentum ------------
c----------------------------------------------

c p_N = p_Npi - p_JSPi in the observable system

      p3(1)  = px_Npi - p_pert(id_JPsi(1),1,i_pert)
      p3(2)  = py_Npi - p_pert(id_JPsi(1),2,i_pert)
      p3(3)  = pz_Npi - p_pert(id_JPsi(1),3,i_pert)

      if(ipauli.eq.1)
     &  call pauli(i2,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
      probab = probab * (1.0-phase) ! modifiacation of the probability according to the Pauli exclusion principle of the nucleon

c----------------------------------------------

      p_pert(id_JPsi(1),4,i_pert)= probab ! JPsi creation probability

      r_pert(id_JPsi(1),1,i_pert)= xxx
      r_pert(id_JPsi(1),2,i_pert)= yyy
      r_pert(id_JPsi(1),3,i_pert)= zzz
      r_pert(id_JPsi(1),4,i_pert)= denst/rho0 ! baryon density at the point of the J/Psi creation
      r_pert(id_JPsi(1),5,i_pert)= time ! time of the J/Psi creation

      nx_pert(id_JPsi(1),0,i_pert) = 1 ! A JPsi exists at the i_pert position
      nx_pert(id_JPsi(1),1,i_pert) = 1 ! A JPsi exists at the i_pert position
      nx_pert(id_JPsi(1),2,i_pert) = ireac ! reaction type
      nx_pert(id_JPsi(1),3,i_pert) = 0 ! which meson -> no data?
      nx_pert(id_JPsi(1),4,i_pert) = i2 ! which nucleon
      nx_pert(id_JPsi(1),5,i_pert) = 1 ! ????????????????????????
      nx_pert(id_JPsi(1),6,i_pert) = id1 ! type of meson
      nx_pert(id_JPsi(1),7,i_pert) = id2 ! type of nucleon
      nx_pert(id_JPsi(1),8,i_pert) = iz1 ! charge of meson
      nx_pert(id_JPsi(1),9,i_pert) = iz2 ! charge of nucleon
      mass_evol(1,1,id_JPsi(1),i_pert) = mmass
      mass_evol(2,1,id_JPsi(1),i_pert) = sqrt(mmass**2+
     &      p_pert(id_JPsi(1),1,i_pert)**2+
     &      p_pert(id_JPsi(1),2,i_pert)**2+
     &      p_pert(id_JPsi(1),3,i_pert)**2)
      mass_evol(3,1,id_JPsi(1),i_pert) = denst/rho0
        if(denst.lt.0.01) write(*,*) 'JPsi_createpi dens',mmass,time,ii
     &        ,i_pert,denst,ix,iy,iz
c      numprodd = numprodd + 1
      zz0 = sqrt(max(radius**2-b**2,0.0))
      iz0 = nint((zzz+zz0)/0.1)+10
      iz0 = min(max(0,iz0),200)
c      iz0 = min(max(0,nint((zzz+zz0)/0.1)+10.),200)
      JPsi_init_pos(1,iz0)= JPsi_init_pos(1,iz0) + 1.0
      iz0 = nint(denst/rho0/0.1)
      iz0 = min(iz0,200)
c      iz0 = min(nint(denst/rho0/0.1),200)
      JPsi_init_pos(2,iz0)= JPsi_init_pos(2,iz0) + 1.0
        
c      write(20,*) "#J/Psi prod: ",i_pert,0,i2,irun,numprodd
      if(i_JPsi.eq.3.or.i_JPsi.eq.4.or.i_JPsi.eq.7.or.i_JPsi.eq.8)
     &     then
        ipi(1,i1) = 0
        p(1,i2) = p3(1)
        p(2,i2) = p3(2)
        p(3,i2) = p3(3)
      end if
c-----------------------------------------------------------------------

 101   continue
c      write(*,*) 'end of JPsi_pert_Npi'
*----------------------------------------------------------------------*
      return
      end

c***********************************************************************
      subroutine Nbarn_DD_pert(i1,i2,id1,id2,beta,srt,xxx,yyy,zzz,sig0,
     &    irun)

*       variables:                                                     *
*  CHANGED!!     ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*     cross-sections are taken from sibirtsev   or  Chung PL B401(97)1 *
*----------------------------------------------------------------------*
      implicit none
      include 'common'
      include 'cominput'
      include 'com_pert'
*----------------------------------------------------------------------*

c ----------- variables used in the N+N collision ---------------------
c ---------------------------------------------------------------------

      real*8 prob_JPsi_scale
      real*8 mmass ! Dmes meson mass
      real*8 pmax,pmax2 ! maximal momentum and momentum square
      integer i1, i2 ! id numbers of nucleons
      integer id1,id2 ! type of the colliding nucleons
      integer id21,id22 ! charges of the nucleons
      integer ireac ! id of the reaction
      integer irun              ! index of parallel runs
      integer ii
      real*8 beta(3) ! velocity of the two nucleons in the observable frame
      real*8 gamma ! 1/sqrt(1-beta^2)
      real*8 srt,s, srtmin ! cms energy of the two nucleons, s=srt**2
      real*8 xxx,yyy,zzz ! position of the collision in the observable frame
      real*8 sig0 ! cross section corresponding to the maximal impact parameter for calling this routine
      real*8 sigma ! cross section for NN --> DD + X  

      integer size_of_run, i_pert_min, i_pert, nt, i_srt
      real*8 probab
      real*8 prob_min

      real*8 rn ! random number generation
      real*8 ranp ! random variable for JPsi momentum generation
      real*8 p_abs_Dmes, E_Dmes, E_Dmes2 ! momentum and energy of D
      real*8 xx,yy,zz,rr ! for random unit vector generation
      real*8 p_x_Dmes, p_y_Dmes, p_z_Dmes ! JPsi momentum components
      real*8 vrel, rself,iself,sgamma,MassJPsi,mass2,pbeta
      real*8 traf,transf ! variables for Lorentz transformation

      integer ix,iy,iz ! integer coordinates of the collision  
      real*8 denst ! density
cc ----------------------------------------------------------------------

c      write(*,*) 'in NbarN_DD_pert',i1,i2,id1,id2,srt
      if((id1+id2 .ne. 0) .or. (id1*id2 .ne. -1)) return
      prob_JPsi_scale = 1./JPsi_scale_factor
      if (rn(iseed) .ge. prob_JPsi_scale)                         return

      size_of_run = (max_pert/num) ! num: number of parallel runs, size_of_run: array size for one parallel run (copy)

      ix = nint(xxx)            ! nint: closest integer
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)
     &     denst = rhb(ix,iy,iz)

      srtmin=10.0
      do ii=1,4                 ! do fo D meson types
        rself=0.0
        iself=0.0
c        call self_Dmes(ii,vrel, denst,rself,iself,sgamma)
        if(Dmes_prop((ii+1)/2,1)+rself/(2.0*Dmes_prop((ii+1)/2,1))
     &       +iself/Dmes_prop((ii+1)/2,1).lt.srtmin)
     &       srtmin = Dmes_prop((ii+1)/2,1)+rself/
     &         (2.*Dmes_prop((ii+1)/2,1))+2.*iself/Dmes_prop((ii+1)/2,1)
      end do
      if(srt.le.srtmin) return

      s     = srt**2
      vrel = 0.0
      gamma = 1.0/sqrt(1.0-(beta(1)**2+beta(2)**2+beta(3)**2)) 
      traf  = gamma / (gamma+1.0)

      i_srt = 1
      do while(i_srt.le.NoL_DD .and. srt.gt.sig_ppbar_DpDm(i_srt,1))
         i_srt=i_srt+1
c         write(*,*) 'DD srt',srt,i_srt,sig_ppbar_DpDm(i_srt,1)
      end do

      sigma = 0.0

      do ii=1,2                 ! do for D meson channels: D+D-, D0D0

c        call self_Dmes(2*ii-1,vrel, denst,rself,iself,sgamma)
c        call self_Dmes(2*ii  ,vrel, denst,rself,iself,sgamma)

        mmass = Dmes_prop(ii,1)
        mass2 = Dmes_prop(ii,1)
c        mmass = JPsiMass(ii,srt,mass2,rself,iself,iseed)
        pmax2=.25*(s-(mass2+mmass)**2)*(s-(mass2-mmass)**2)/s
        if(pmax2 .le. 0.0)                      goto 95
c     pmax  = the maximal D meson momentum
        pmax  = sqrt(pmax2)

c--------------------------------------------------
c------ Creation probability ----------------------
c--------------------------------------------------
        if(ii.eq.1) then
          if(i_srt .gt.1 .and. i_srt.le. NoL_DD+1) then
            sigma =  sig_ppbar_DpDm(i_srt-1,2)
          else if(i_srt .eq.1) then
            sigma =  sig_ppbar_DpDm(i_srt,2)
          end if
        else if(ii.eq.2) then
          if(i_srt .gt.1 .and. i_srt.le. NoL_DD+1) then
            sigma =  sig_ppbar_D0D0(i_srt-1,2)
          else if(i_srt .eq.1) then
            sigma =  sig_ppbar_D0D0(i_srt,2)
          end if
        endif
c---------------------------------------------
c---------------Dmes momentum ----------------
c---------------------------------------------

c        size of the D meson mom. ranp*pmax  

c           direction of the D meson mom. xx/rr, yy/rr, zz/rr
 51     xx       = 1. - 2. * rn(iseed) ! choosing a unit vector for the momentum direction   
        yy       = 1. - 2. * rn(iseed)
        zz       = 1. - 2. * rn(iseed)
        rr       = sqrt( xx**2 + yy**2 + zz**2 )
        if((rr .lt. 0.001) .or. (rr .gt. 1.) ) goto 51

        p_abs_Dmes = pmax
        E_Dmes  = sqrt(p_abs_Dmes**2 + mmass**2)
        E_Dmes2 = sqrt(p_abs_Dmes**2 + mass2**2)

        p_x_Dmes = p_abs_Dmes * xx / rr 
        p_y_Dmes = p_abs_Dmes * yy / rr
        p_z_Dmes = p_abs_Dmes * zz / rr

        probab = sigma/sig0*JPsi_scale_factor ! probability of producing a DD
c        write(*,*)'DDprod',sigma,i_srt,srt,JPsi_scale_factor,sig0,probab
c----------------------------------------------
c-- Find a place for the created particle -----
c----------------------------------------------
        
        i_pert_min = 0
        prob_min  = 1000.0
        do i_pert = (irun-1) * size_of_run + 1, irun * size_of_run ! i_pert: running index of nx_pert array containing perturbativ particles  
c                                                                   ! (irun-1) * size_of_run + 1: beginning index of current run in nx_pert array 

          if (nx_pert(id_Dmes(2*ii-1),0,i_pert) .eq. 0) goto 54 ! true: found an empty slot for the just created J/Psi
          if (prob_min .gt. p_pert(id_Dmes(2*ii-1),4,i_pert)) then
            i_pert_min = i_pert
            prob_min  = p_pert(id_Dmes(2*ii-1),4,i_pert)
          endif
        enddo
        i_pert = i_pert_min
        if (probab .lt. prob_min) goto 95
 54     continue

*   Dmes1 momentum in observable system
        pbeta  = beta(1)*p_x_Dmes + beta(2)*p_y_Dmes + beta(3)*p_z_Dmes
        transf = gamma * (traf * pbeta + E_Dmes)

        p_pert(id_Dmes(2*ii-1),0,i_pert)= mmass
        p_pert(id_Dmes(2*ii-1),1,i_pert)= p_x_Dmes + beta(1) * transf
        p_pert(id_Dmes(2*ii-1),2,i_pert)= p_y_Dmes + beta(2) * transf
        p_pert(id_Dmes(2*ii-1),3,i_pert)= p_z_Dmes + beta(3) * transf
        p_pert(id_Dmes(2*ii-1),4,i_pert)= probab
        r_pert(id_Dmes(2*ii-1),1,i_pert)= xxx
        r_pert(id_Dmes(2*ii-1),2,i_pert)= yyy
        r_pert(id_Dmes(2*ii-1),3,i_pert)= zzz
        r_pert(id_Dmes(2*ii-1),4,i_pert)= denst/rho0 ! baryon density at the point of the J/Psi creation
        r_pert(id_Dmes(2*ii-1),5,i_pert)= time ! time of the J/Psi creation

        nx_pert(id_Dmes(2*ii-1),0,i_pert) = ii ! iitype Dmes at the i_pert pos
        nx_pert(id_Dmes(2*ii-1),1,i_pert) = 1 ! Dmes exists at the i_pert pos
        nx_pert(id_Dmes(2*ii-1),2,i_pert) = ireac ! reaction type
        nx_pert(id_Dmes(2*ii-1),3,i_pert) = i1 ! which nucleon, N_A
        nx_pert(id_Dmes(2*ii-1),4,i_pert) = i2 ! which nucleon, N_B
        nx_pert(id_Dmes(2*ii-1),5,i_pert) = 1 ! ????????????????????????
        nx_pert(id_Dmes(2*ii-1),6,i_pert) = id1 ! type of N_A
        nx_pert(id_Dmes(2*ii-1),7,i_pert) = id2 ! type of N_B
        nx_pert(id_Dmes(2*ii-1),8,i_pert) = id21 ! charge of N_A
        nx_pert(id_Dmes(2*ii-1),9,i_pert) = id22 ! charge of N_B

*   Dmes1 momentum in observable system
        pbeta  = beta(1)*p_x_Dmes + beta(2)*p_y_Dmes + beta(3)*p_z_Dmes
        transf = gamma * (-traf * pbeta + E_Dmes2)

        p_pert(id_Dmes(2*ii),0,i_pert)= mass2
        p_pert(id_Dmes(2*ii),1,i_pert)= -p_x_Dmes + beta(1) * transf
        p_pert(id_Dmes(2*ii),2,i_pert)= -p_y_Dmes + beta(2) * transf
        p_pert(id_Dmes(2*ii),3,i_pert)= -p_z_Dmes + beta(3) * transf
        p_pert(id_Dmes(2*ii),4,i_pert)= probab
        r_pert(id_Dmes(2*ii),1,i_pert)= xxx
        r_pert(id_Dmes(2*ii),2,i_pert)= yyy
        r_pert(id_Dmes(2*ii),3,i_pert)= zzz
        r_pert(id_Dmes(2*ii),4,i_pert)= denst/rho0 ! baryon density at the point of the J/Psi creation
        r_pert(id_Dmes(2*ii),5,i_pert)= time ! time of creation

        nx_pert(id_Dmes(2*ii),0,i_pert) = ii ! iitype Dmes at the i_pert pos
        nx_pert(id_Dmes(2*ii),1,i_pert) = 1 ! Dmes exists at the i_pert pos
        nx_pert(id_Dmes(2*ii),2,i_pert) = ireac ! reaction type
        nx_pert(id_Dmes(2*ii),3,i_pert) = i1 ! which nucleon, N_A
        nx_pert(id_Dmes(2*ii),4,i_pert) = i2 ! which nucleon, N_B
        nx_pert(id_Dmes(2*ii),5,i_pert) = 1 ! ????????????????????????
        nx_pert(id_Dmes(2*ii),6,i_pert) = id1 ! type of N_A
        nx_pert(id_Dmes(2*ii),7,i_pert) = id2 ! type of N_B
        nx_pert(id_Dmes(2*ii),8,i_pert) = id21 ! charge of N_A
        nx_pert(id_Dmes(2*ii),9,i_pert) = id22 ! charge of N_B

 95     continue
      end do ! end of ii loop
      return
      end
*----------------------------------------------------------------------*
      subroutine JPsi_pert_NN(i1,i2,id1,id2,beta,srt,xxx,yyy,zzz,sig0,
     &       irun)

*       variables:                                                     *
*  CHANGED!!     ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*     cross-sections are taken from sibirtsev   or  Chung PL B401(97)1 *
*----------------------------------------------------------------------*
      implicit none
      include 'common'
      include 'cominput'
      include 'com_pert'
*----------------------------------------------------------------------*

c ---------------------------------------------------------------------
c ----------- variables used in the N+N collision ---------------------
c ---------------------------------------------------------------------

      real*8 prob_JPsi_scale
      real*8 mmass ! J/Psi meson mass
      real*8 pmax,pmax2 ! maximal momentum and momentum square of the J/Psi (NN --> J/Psi + NN)
      integer i1, i2 ! id numbers of nucleons
      integer id1,id2 ! type of the colliding nucleons
      integer id21,id22 ! charges of the nucleons
      integer ireac ! id of the reaction
      integer ib ! number of times the i1 nucleon has collided
      integer irun              ! index of parallel runs
      integer ii
      real*8 beta(3) ! velocity of the two nucleons in the observable frame
      real*8 gamma ! 1/sqrt(1-beta^2)
      real*8 srt,s ! cms energy of the two nucleons, s=srt**2
      real*8 xxx,yyy,zzz ! position of the collision in the observable frame
      real*8 sig0 ! cross section corresponding to the maximal impact parameter for calling this routine
      real*8 sig1, sig2, sig3, sig4 ! auxiliary cross section variables
      real*8 sigma ! cross section for NN --> J/Psi + X  
      integer srt_MeV_int ! closest integer of 'srt' in MeV
      integer srt1_l, srt1_h, srt2_l, srt2_h ! auxiliary 'srt' variables
      integer srt1_id, srt2_id ! auxiliary 'srt' variables

      real*8 alfa1, beta1, f_JPsi, a_fac ! parameters of the N + N --> J/Psi + X cross section (Eq. 1 of O. Linnyk et al. Nucl. Phys. A 786 (2007) 183)
      parameter (alfa1=10.0, beta1=0.775, f_JPsi=0.581, a_fac=0.2) ! a_fac = 0.2 microbarn, sig0 in milibarn

      integer size_of_run, i_pert_min, i_pert, nt, iz0
      real*8 probab, fact_pauli ! probability of JPsi creation
      real*8 prob_min, zz0

      real*8 rn ! random number generation
      real*8 ranp ! random variable for JPsi momentum generation
      real*8 p_abs_JPsi, E_JPsi ! momentum and energy of JPsi in  the 3 particle c.m.s.
      real*8 xx,yy,zz,rr ! for random unit vector generation
      real*8 xxn,yyn,zzn,rrn, pnuc2 ! for random unit vector generation for nucleons
      real*8 p_x_JPsi, p_y_JPsi, p_z_JPsi ! JPsi momentum components
      real*8 s_prime ! see later
      real*8 E_nucl_prime,p_nucl_prime ! momentum and energy of the nucleon in the 2 nucleon subsystem
      real*8 factn1,factn2 ! transformation factors for the N_A' momentum calculation
      real*8 srtmin, JPsiMass
      real*8 vrel, rself,iself,sgamma,MassJPsi,mass2,j0,j1,j2,j3
      real*8  p3(3),p4(3)             ! momentum components of N_A'
      real*8  e3, e4, E_mes, pbeta ! energy of N_A', pbeta = p_vec*beta_vec
      real*8 traf,transf ! variables for Lorentz transformation
      real*8 phase ! from Pauli exclusion principle
      integer ntag ! flag which tells if phase-space is pauli-blocked (-1: blocked)
      integer ix,iy,iz ! integer coordinates of the collision  
      real*8 denst ! density
      integer bin_dens,bin_time ! bins for the histogram
      integer p_dbb(0:999,0:999) ! for density histogram 


c      write(*,*) 'in JPsi_pert_cr',i1,i2,id1,id2,srt

      if(id1*id2.lt.0) return
c
ccc                 pp -> pp charmonium
c

      srtmin = JPsi_prop(1,1)+2.0*rmass - 0.500
      if(srt.le.srtmin) return

      prob_JPsi_scale = 1./JPsi_scale_factor
      if (rn(iseed) .ge. prob_JPsi_scale)                         return
cc ----------------------------------------------------------------------
      size_of_run = (max_pert/num) ! num: number of parallel runs, size_of_run: array size for one parallel run (copy)

      ix = nint(xxx)            ! nint: closest integer
      iy = nint(yyy)
      iz = nint(zzz)
      denst = 0.0
      if(iabs(ix).le.maxx.and.iabs(iy).le.maxx.and.iabs(iz).le.maxz)then
        j0=rhob_4(0,ix,iy,iz)  
        j1=rhob_4(1,ix,iy,iz)  
        j2=rhob_4(2,ix,iy,iz)  
        j3=rhob_4(3,ix,iy,iz)  
        denst = sqrt(j0**2-j1**2-j2**2-j3**2)
c        write(*,*) 'JPcreate dens:', j0,denst,rhb(ix,iy,iz)
        denst = rhb(ix,iy,iz)
      else
        denst = 0.0
        j0 = 0.0
      end if

c     if(denst.lt.0.01) write(*,*) 'JPsi_create dens',
c     &   denst,ix,iy,iz

      s     = srt**2
      srt_MeV_int = nint(srt*1000.) ! closest integer of sqrt(s) in MeV
ccc   loop over the 3 charmonium states
      vrel = 0.0
      gamma = 1.0/sqrt(1.0-(beta(1)**2+beta(2)**2+beta(3)**2)) 
      traf  = gamma / (gamma+1.0)
      mass2 = 2.*rmass

      do ii=1,3                 ! do fo JPsi types
         
c           direction of the charmonium mom. xx/rr, yy/rr, zz/rr
 51     xx       = 1. - 2. * rn(iseed) ! choosing a unit vector for the momentum direction   
        yy       = 1. - 2. * rn(iseed)
        zz       = 1. - 2. * rn(iseed)
        rr       = sqrt( xx**2 + yy**2 + zz**2 )
        if((rr .lt. 0.001) .or. (rr .gt. 1.) ) goto 51

        pmax2=.25*(s-(JPsi_prop(ii,1)+mass2)**2)*
     &       (s-(JPsi_prop(ii,1)-mass2)**2)/s
        if(pmax2.gt.0.0001) then
          pmax = sqrt(pmax2)
          p_x_JPsi = pmax * xx / rr ! JPsi momentum components for vrel       
          p_y_JPsi = pmax * yy / rr
          p_z_JPsi = pmax * zz / rr
          E_JPsi   = sqrt(pmax**2 + JPsi_prop(ii,1)**2) ! for vrel JPsi energy in the pp c.m. system
          call lorentz(-beta(1),-beta(2),-beta(3),
     &       p_x_JPsi,p_y_JPsi,p_z_JPsi,E_JPsi)
        
          if (denst.lt.1.0e-3) then
            denst=0.
            vrel = 0.
          else
            call f77flush()
c         write(*,*) 'lorentz call gradupi',betacm(1),betacm(2),betacm(3),
c     &     j1,j2,j3,j0
            if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba JPsi_cr pp, mass<0",j0,j1,j2,j3
              stop
            end if
            call lorentz(p_x_JPsi/E_JPsi,p_y_JPsi/E_JPsi,p_z_JPsi/E_JPsi
     &        ,j1,j2,j3,j0)
            vrel = sqrt(j1**2+j2**2+j3**2)/j0
          end if
        else
          pmax = 0.0
          vrel = 0.0 
        end if

        call self_JPsi(ii,vrel, denst,rself,iself,sgamma)

c     size of the charmonium mom. ranp*pmax  
 50     continue
        ranp     = rn(iseed)
        xx = ranp**2 * sqrt(1.00001-ranp**2) / 0.385 ! p^2*sqrt(1-p^2) is the momentum distribution 
c          0.385 is the normalization factor, the maximal value of xx is 1.
        if(xx .lt. rn(iseed)) goto 50

        sig1 = 0.0

c        write(*,*) 'JP create1',srt,mass2+JPsi_prop(ii,1)+
c     &       0.5*rself/JPsi_prop(ii,1)+2.*iself/JPsi_prop(ii,1)
        if(srt.le.mass2+JPsi_prop(ii,1)+rself/(2.0*JPsi_prop(ii,1))+
     &       2.*iself/JPsi_prop(ii,1)) goto 95

        mmass = JPsiMass(ii,srt,mass2,rself,iself,iseed)
        pmax2=.25*(s-(mass2+mmass)**2)*(s-(mass2-mmass)**2)/s
        if(pmax2 .le. 0.0)                      goto 95
c     pmax  = the maximal J/Psi momentum
        pmax  = sqrt(pmax2)

c--------------------------------------------------
c------ Creation probability ----------------------
c--------------------------------------------------
        fact_pauli = 1.0
        if((id1 .eq. 1) .and. (id2 .eq. 1)) then  ! N + N
          ireac = 10
        else if((id1 .gt. 1) .or. (id2 .gt. 1)) then  ! N + R
          ireac = 11
        endif
        if(ii.eq.1) sig1 = f_JPsi*a_fac * (1-mmass/srt)**alfa1 *
     $      (srt/mmass)**beta1/1000. ! sigma in milibarn; Eq. 1 of Nucl. Phys. A+ 786 (2007) 183

c---------------------------------------------
c---------------JPsi momentum ----------------
c---------------------------------------------

        p_abs_JPsi = pmax*ranp ! outgoing JPsi momentum absolute value  
        E_JPsi   = sqrt(p_abs_JPsi**2 + mmass**2) ! JPsi energy in the 3 particle c.m. system
        s_prime  = s - 2.0*srt*E_JPsi+ mmass**2 ! s'=(p_N_A' + p_N_B')^2 : two nucleon s in the 3 particle c. m. system
        E_nucl_prime = 0.5* sqrt(s_prime) ! energy of the nucleons in the two nucleon c. m. system
        p_nucl_prime = sqrt(E_nucl_prime**2-rmass**2) ! absolute value of the nucleon momentum in the two nucleon c. m. system 
        p_x_JPsi = p_abs_JPsi * xx / rr ! JPsi momentum components
        p_y_JPsi = p_abs_JPsi * yy / rr
        p_z_JPsi = p_abs_JPsi * zz / rr

c----------------------------------------------
c--------------- N_A' momentum ----------------
c----------------------------------------------

 52     xxn       = 1. - 2. * rn(iseed)
        yyn       = 1. - 2. * rn(iseed)
        zzn       = 1. - 2. * rn(iseed)
        rrn       = sqrt( xxn**2 + yyn**2 + zzn**2 )
        if((rrn .lt. 0.001) .or. (rrn .gt. 1.) ) goto 52
        factn1   = s_prime + sqrt(s_prime)*(srt-E_JPsi)
        factn2   = p_nucl_prime*(p_x_JPsi*xxn+p_y_JPsi*yyn+p_z_JPsi*zzn)
     $        /factn1/rrn-E_nucl_prime/sqrt(s_prime)
*     p3:                    nucleon momentum in N_A - N_B c.m. system (it is the same as the 3 particle c.m.s.)
c            to the 3 particle c.m.s (JPsi + N_A' + N_B')
        p3(1)    = p_nucl_prime*xxn/rrn + factn2*p_x_JPsi ! Lorentz transformation from the two nucleon (A' + B') system
        p3(2)    = p_nucl_prime*yyn/rrn + factn2*p_y_JPsi
        p3(3)    = p_nucl_prime*zzn/rrn + factn2*p_z_JPsi
        e3       = sqrt(rmass**2+p3(1)**2+p3(2)**2+p3(3)**2)
*     p3:                    nucleon momentum in observable system
         
        pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
        transf = gamma * (traf * pbeta + e3)
        p3(1)  = p3(1) + beta(1) * transf
        p3(2)  = p3(2) + beta(2) * transf
        p3(3)  = p3(3) + beta(3) * transf
        if(ipauli.eq.1)
     &     call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))

        fact_pauli = (1.0-phase) ! modification of the probability according to the Pauli exclusion principle of N_A'

c----------------------------------------------
c--------------- N_B' momentum ----------------
c----------------------------------------------
        if(ireac.eq.10) then
          factn2= -p_nucl_prime*(p_x_JPsi*xxn+p_y_JPsi*yyn+p_z_JPsi*zzn)
     $        /factn1/rrn-E_nucl_prime/sqrt(s_prime) ! Eddig hibas volt; itt a factn2-ben is
c         meg kell cserelni p_nucl_prime elojelet !
*   p4:                    nucleon momentum in i1-i2-c.m. system
          p4(1)    = -p_nucl_prime*xxn/rrn + factn2*p_x_JPsi 
          p4(2)    = -p_nucl_prime*yyn/rrn + factn2*p_y_JPsi 
          p4(3)    = -p_nucl_prime*zzn/rrn + factn2*p_z_JPsi 
          e4       = sqrt(rmass**2+p4(1)**2+p4(2)**2+p4(3)**2)
*   p4:                   nucleon momentum in observable system
          pbeta  = beta(1)*p4(1) + beta(2)*p4(2) + beta(3)*p4(3)
          transf = gamma * (traf * pbeta + e4)
          p4(1)  = p4(1) + beta(1) * transf
          p4(2)  = p4(2) + beta(2) * transf
          p4(3)  = p4(3) + beta(3) * transf
          if(ipauli.eq.1)
     &     call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,p3(1),p3(2),p3(3))
          fact_pauli = fact_pauli * (1.0-phase) ! modifiacation of the probability according to the Pauli exclusion principle of N_B'
        end if
c-------------------------------------------------------
c--------End of the momentum determination -------------
c-------------------------------------------------------

        sig1 = sig1 * fact_pauli

 95     continue
        sigma = sig1
        if(sigma.lt.1.e-9) return
c        write(*,*) 'JPcreate4',sig1,sig2,sig3,sig4,ii
        probab = sigma/sig0*JPsi_scale_factor ! probability of producing a J/Psi
        
c----------------------------------------------
c-- Find a place for the created particle -----
c----------------------------------------------
        
        i_pert_min = 0
        prob_min  = 1000.0
        do i_pert = (irun-1) * size_of_run + 1, irun * size_of_run ! i_pert: running index of nx_pert array containing perturbativ particles  
c                                                                   ! (irun-1) * size_of_run + 1: beginning index of current run in nx_pert array 
c          write(*,*)'JPsicreate2',ii,i_pert,id_JPsi(ii),irun,size_of_run
          if (nx_pert(id_JPsi(ii),0,i_pert) .eq. 0) goto 54 ! true: found an empty slot for the just created J/Psi
          if (prob_min .gt. p_pert(id_JPsi(ii),4,i_pert)) then
            i_pert_min = i_pert
            prob_min  = p_pert(id_JPsi(ii),4,i_pert)
          endif
        enddo
        i_pert = i_pert_min
        if (probab .lt. prob_min) goto 98
 54     continue

*   J/Psi momentum in observable system
        pbeta  = beta(1)*p_x_JPsi + beta(2)*p_y_JPsi + beta(3)*p_z_JPsi
        transf = gamma * (traf * pbeta + E_JPsi)

        p_pert(id_JPsi(ii),0,i_pert)= mmass
        p_pert(id_JPsi(ii),1,i_pert)= p_x_JPsi + beta(1) * transf
        p_pert(id_JPsi(ii),2,i_pert)= p_y_JPsi + beta(2) * transf
        p_pert(id_JPsi(ii),3,i_pert)= p_z_JPsi + beta(3) * transf
        p_pert(id_JPsi(ii),4,i_pert)= probab
c        if(i_pert.eq.4012) write(*,*) 'JP_4012create1 p',ireac,
c     &   p_pert(id_JPsi(ii),1,i_pert),p_pert(id_JPsi(ii),2,i_pert),
c     &   p_pert(id_JPsi(ii),3,i_pert),beta(1),beta(2),beta(3),E_JPsi,
c     &   p_x_JPsi,p_y_JPsi,p_z_JPsi,transf,p_abs_JPsi,gamma,traf,pbeta   
        r_pert(id_JPsi(ii),1,i_pert)= xxx
        r_pert(id_JPsi(ii),2,i_pert)= yyy
        r_pert(id_JPsi(ii),3,i_pert)= zzz
        r_pert(id_JPsi(ii),4,i_pert)= denst/rho0 ! baryon density at the point of the J/Psi creation
        r_pert(id_JPsi(ii),5,i_pert)= time ! time of the J/Psi creation

        nx_pert(id_JPsi(ii),0,i_pert) = ii ! iitype JPsi at the i_pert position
        nx_pert(id_JPsi(ii),1,i_pert) = 1 ! JPsi exists at the i_pert position
        nx_pert(id_JPsi(ii),2,i_pert) = ireac ! reaction type
        nx_pert(id_JPsi(ii),3,i_pert) = i1 ! which nucleon, N_A
        nx_pert(id_JPsi(ii),4,i_pert) = i2 ! which nucleon, N_B
        nx_pert(id_JPsi(ii),5,i_pert) = 1 ! ????????????????????????
        nx_pert(id_JPsi(ii),6,i_pert) = id1 ! type of N_A
        nx_pert(id_JPsi(ii),7,i_pert) = id2 ! type of N_B
        nx_pert(id_JPsi(ii),8,i_pert) = id21 ! charge of N_A
        nx_pert(id_JPsi(ii),9,i_pert) = id22 ! charge of N_B
        mass_evol(1,1,id_JPsi(ii),i_pert) = mmass
        mass_evol(2,1,id_JPsi(ii),i_pert) = sqrt(mmass**2+
     &      p_pert(id_JPsi(ii),1,i_pert)**2+
     &      p_pert(id_JPsi(ii),2,i_pert)**2+
     &      p_pert(id_JPsi(ii),3,i_pert)**2)
        mass_evol(3,1,id_JPsi(ii),i_pert) = denst/rho0
c        if(denst.lt.0.01) write(*,*) 'JPsi_create2 dens',mmass,time,ii
c     &        ,i_pert,denst,ix,iy,iz
c     numprodd = numprodd + 1
        zz0 = sqrt(max(radius**2-b**2,0.0))
        iz0 = nint((zzz+zz0)/0.1)+10
        iz0 = min(max(0,iz0),200)
c        iz0 = min(max(0,nint((zzz+zz0)/0.1)+10.0),200)
        JPsi_init_pos(1,iz0)= JPsi_init_pos(1,iz0) + 1.0
        iz0 = nint(denst/rho0/0.1)
        iz0 = min(iz0,200)
        JPsi_init_pos(2,iz0)= JPsi_init_pos(2,iz0) + 1.0
        
 98     continue
c      write(20,*) "#J/Psi prod: ",i_pert,i1,i2,irun,numprodd
c-----------------------------------------------------------------------
c----------------------------------------------
c-- Find a place for the created particle -----
c----------------------------------------------

c-----------------------------------------------------------------------
      end do ! end of ii loop
      return
      end
