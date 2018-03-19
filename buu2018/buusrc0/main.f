************************************************************************
*                 #####    Full BUU        #####                       *

*                                                                      *
*  last modified in  2016                                              *
*                                                                      *
*  ##    direct pions                                                  *
*         possible change of the width of resonances                   *
*         nstar(1440) resonance is included                            *
*                                                                      *
*  ##    momentum dependent force gbdmom                               *
*         effective mass and effective momentum in it                  *
*                                                                      *
*        this code includes the finit e range force and coulomb         *
*        potential(new version)                                        *
*                                                                      *
*        this version includes  125*125  gaussian smearing for density *
*        with width = 1.0 fm, cutoff = sqrt(5.0) fm                    *
*        and woods-saxon initial distribution with                     *
*            saa = 0.02444 * a**1./3. + 0.2864   for with yukawa       *
*            saa = 0.13                          for others            *
*        and local thomas fermi with average density for the initial   *
*        distribution of momenta                                       *
*                                                                      *
*        new version of the boltzmann-uehling-uhlenbeck-program        *
*        made by       stefan teis  1994                               *
*                                                                      *
*                                                                      *
************************************************************************
*  output units:                                                       *
*     15 (=isum)   -  normal control output; buu + particle production *
*----------------------------------------------------------------------*
*  parameters:                                                         *
*     maxpar  -  maximum number of testparticles program can handle    *
*     maxx    -  number of meshpoints in x and y direction = 2 maxx + 1*
*     maxz    -  number of meshpoints in z direction       = 2 maxz + 1*
*     rmass   -  1 atomic mass unit "gev/c**2"                         *
*----------------------------------------------------------------------*
      implicit  none
      real*8  alp, bet
c      integer k2,kk,ncputim
c      real*8 abababa
      integer   ncatas, isu, isubsold, ni1, nj1, irun, j1,
     1          nii, jn, jj, ix,iy,iz, k6, k7, is, js, ij, i0,
     2          j0, jseed, i, ii, ik, jcoll,
     3          masspec, nlost
      integer mass
      integer ictime, nrspl, ncqn, mcelpri, icell,
     1        nnuclpr, mdel1, mdel0,mdel2, mprot,
     1        mneut, mppan, mdelm, mrspl, mrsze, mcopt, mcodd,
     1        mqspl, mqsze, lcnu, lcde, lcrn, lcqn,
     2        mstanu, msprnu, nproton, isubsqt,npbar, nneutron
      integer napart, npionpr, nod, non, no1,no2, npi,
     1        npidi,
     1        neta, ipavp, ilast, iendel, nprot, nneut, ndel1,
     2        ndel2, ndel0, mcnne, mcond, mcopn, mcbre, merho,
     3        meeta, mesig, mpidi, neeta, nerho, nesig,
     3        ncnu, ncde, ncrn, ncopn, ncodd, ncond,
     4        ndelm, nrsze, nqspl, nqsze, ncbre, ncopt, nppan,
     5        jctime, nt, lcond, lcopn, lcopt, lppan, lcodd, lcbre
      integer  collkplus(2*40000),collkplusd(2*40000),cotkpl(200)
      integer ikairun
      real*8    mdpartini
      real*8    eventnum, scala, scalpro, cpuisuloop, cputime,
     1        tramb, tram,dt0
      real*8  en, bmass, ylab, ylabb, eganz, etot, dendel, wmass,
     1      punt, ydel, ypr, yta, ymim, ymam, yv,
     1      yref, gammta, gammpr, totcros,
     2      aaa, bbb, ccc, rkk, u00, vzero,
     3      radta, radpr, rmax, rxta, rzta, rxpr, rzpr,
     4      betata, betapr, surfta, surfpr, p0
c      real*8 tin,tout
      real*8  t0, vol, dpmin, enerc, dummy
      parameter     (alp=0.33,  bet=5.7)
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"
      include"commonthreef"
      include"com_cont_epair"

      integer numomega
      common /omega/  numomega
      common /nthwhw/  nthw
      integer nthw
      common /rancheckhw/ rancheck
      integer rancheck
      common /counthw/ ihw
      integer  ihw
      common /correlkaon/ ikairun(0:maxkaon)
      COMMON/self_stored/ totom,totro
      real*8 totom(7,2000),totro(7,2000)

      integer dimlcoll
      parameter(dimlcoll = 3*nres+6+nres**2)
      integer dimlmesc
      parameter(dimlmesc = nres+9)
*----------------------------------------------------------------------*
*
      real*8         dengamro(10000)
      real*8     ptxyz(3,0:200),ppxyz(3,0:200),pixyz(3,0:200),
     &         rtxyz(3,0:200),rpxyz(3,0:200),rixyz(3,0:200),
     &         etta(0:200),etpr(0:200),ekta(0:200),ekpr(0:200),
     &         ecta(0:200),ecpr(0:200),emta(0:200),empr(0:200),
     &         eipr(0:200),eita(0:200),eket(0:200),eztt(0:200)
      real*8     edif(0:200),eptt(0:200),epct(0:200),eftt(0:200),
     &         ettt(0:200),ektt(0:200),ectt(0:200),eitt(0:200),
     &         etpi(0:200),ekpi(0:200),ecpi(0:200),eipi(0:200),
     &         emtt(0:200),eflo(0:200),eflpi(0:200),tij(3,3,0:200),
     &         tijpi(3,3,0:200),fptp(6,20),time1(0:200),rmeson(0:200)
      real*8     rhmax(0:200),rhavr(0:200),reovr(0:200),distp(0:200),
     &         rmc(2), rmy(2)
      real*8     suru(0:200),ener(0:200),rati(0:200),time2(0:200)
      real*8     resdens(-2:nres,0:200,50),rootmsq(0:200)
      integer  icptimep
      integer  ndis(0:20,0:200), nlst(0:200)
      integer  numes(18),lumes(18),mumes(18,0:200)
      integer  lcoll(-6:dimlcoll),ncoll(-1:dimlcoll),
     &             mcoll(-1:dimlcoll)
      integer  lmesc(dimlmesc),nmesc(dimlmesc),mmesc(dimlmesc)
      integer  lmesa(dimlmesc),nmesa(dimlmesc),mmesa(dimlmesc),lmesa2
      integer  mcol(-1:dimlcoll,0:200), mesc(dimlmesc,0:200),
     &             mesa(dimlmesc,0:200)
      integer  mend(0:200),medd(0:200),mepn(0:200),mept(0:200),
     &         mdip(0:200),meta(0:200),mron(0:200),msin(0:200),
     &         mpan(0:200)
      integer  mprt(0:200),mntr(0:200),mrpl(0:200),mrze(0:200),
     &         mde2(0:200),mde1(0:200),mde0(0:200),mdem(0:200),
     &         mqpl(0:200),mqze(0:200),
     &         mcnu(0:200),mcde(0:200),mcrn(0:200),mcqn(0:200)
      integer  cres(1:nres,1:9)
      integer nelse
*----------------------------------------------------------------------*
*  input-section:                                                      *
*                                                                      *
*  1) target-related quantities                                        *
*       massta   -  target mass                             (integer)  *
*       mstapr   -  proton number in target                 (integer)  *
*                                                                      *
*  2) projectile-related quantities                                    *
*       masspr   -  projectile mass                         (integer)  *
*       msprpr   -  proton number in projectile             (integer)  *
*       mpion    -  0-> nucleus - nucleus collision         (integer)  *
*                -  1-> pi+ - nucleus; -1-> pi- -nucleus; 2->pi0-nucl. *
*       elab     -  beam energy in "gev/nucleon"               (real)  *
*       b        -  impact parameter "fm"                      (real)  *
*       rdist    -  rmax =radta+radpr+rdist; initial distance  (real)  *
*                                                                      *
*  3) program-control parameters                                       *
*       iseed    -  seed for random number generator        (integer)  *
*       dt       -  time-step-size "fm/c"                      (real)  *
*       ntmax    -  total number of time steps              (integer)  *
*       icoll    -  (= 1 -> mean field only,                           *
*                -   =-1 -> cascade only, = 0 -> full buu)  (integer)  *
*       num      -  number of testparticles per nucleon     (integer)  *
*       insys    -  (=0 -> lab-system, =1-> c.m. system)    (integer)  *
*       iavoid   -  (=1 -> avoid first coll. within same nucl.         *
*                    =0 -> allow them)                      (integer)  *
*       iangc    -  0-> isotropic; 1-> elastic; 2-> experimental;      *
*                    angular distribution for delta         (integer)  *
*       ggam     -  change the resonance width; 0-> width=width*(1-f), *
*                    else-> width=width*ggam                   (real)  *
*  4) choice of potential                                              *
*       ipot     -  0-> t0 and t3 terms  1-> t0,t3 and yukawa          *
*                   2-> t0,t3 and gbd momentum dep. force              *
*       ipou     -  0-> no coulomb, 1->with coulomb                    *
*       rpot     -  sigma=rpot+1 power of density depedent potential   *
*                   in mean field potential                    (real)  *
*       smu      -  range of yukawa potential (1/fm)                   *
*                                                                      *
*  5) test of coulomb and yukawa potential                             *
*      (alp)     -  minmum eigenvalue  ; fixed  0.33                   *
*      (bet)     -  maxmum eigenvalue  ; fixed  5.70                   *
*       initit   -  number of iteration for the initial condition      *
*       nity     -  number of iteration for yukawa                     *
*                                                                      *
*  6) control-printout options                                         *
*       isumm    -  0-> no summry, 1-> summry                          *
*       nfreq    -  number of time steps after which printout          *
*                   is required for time evolution for coll.(integer)  *
*       nfrep    -  number of time steps after which printout          *
*                   is required for figures                 (integer)  *
*                                                                      *
*  7) particle production                                              *
*                                                                      *
*  8) direct pions                                                     *
*       idipi    -  0=> no, 1 => pion emission at the end,             *
*                   2=> direct pions;                       (integer)  *
*                                                                      *
*======================================================================*
*                   ( following parameters are read in sub patini:)    *
*======================================================================*
      common /cputime/ cptimepr, cputim0, cputim, clockmax
      real*8 cptimepr,cputim0, cputim, clockmax, second
      external second
      numomega = 0
c
      ihw = 0
      rancheck = 0
      write(*,*)'BUU starts'
c TEST!!!:
c      call test_acc
c      stop
      clockmax = 0              !         old value 4294.97
c     call otime(ncputim)
c     cputim0=0.01*float(ncputim)
      cputim0  = second()
      cptimepr = cputim0
      icptimep = 0
      write(*,*)'before readin main'
      call readin
      write(*,*)'after readin main '
*-----------------------------------------------------------------------
*
      jseed = iseed/1000
cr    iseed = -abs(iseed - 1000*jseed)
c      write(46,*)  '   iseed - jseed ', iseed, jseed
c
      isubs=max(1,isubs)
      write(*,*)'nach rin',isubs,num
*
      mstanu=massta-mstapr
      msprnu=masspr-msprpr
      nproton=(mstapr+msprpr)*num
      maxp = maxppar / num
      maxb = maxpar / num
      mass   = abs(massta) + abs(masspr)
      totmass =  mass
      ntotal = num * mass
*
      if((masspr.eq.0) .and. (mpion.eq.0)) elab=0.0
      if(mpion.ne.0) masspr = 0
      if(mpion.ne.0) idipi  = 2
      if(mpion.ne.0) msprnu = 0
      if(mpion.eq.2) msprpr = 0
      if(mpion.eq.3) b = 0.0
      if(iabs(mpion).eq.1) msprpr = mpion
*
      if(i_JPsi.ge.1) idilepton = 1
      if(idilepton.eq.0) then
        ivmesdil   = 0
        ibrems     = 0
        ipiNbrems  = 0
        iresdalitz = 0
        imesdalitz = 0
        idilper = 0
      endif
*
*-----------------------------------------------------------------------
*   determine the coefficients of the potential
      write(*,*)'before potential',isubs

c      call potential(ipot,rpot,smu,tt0,tt03,yv,u00,aaa,bbb,rkk,idelpot,
c     &  vzero)

*
*----------------------------------------------------------------------*
      call front(isum,massta,masspr,elab,mstapr,mstanu,msprpr,msprnu,
     &           mpion,b)
      write(*,*)'after front'
      call inecho(isum  , massta, masspr, elab  , isubs ,
     &            iseed , dt    , ntmax , icoll , num   , insys ,
     &            iavoid, iwidth, ipou  , mstapr,mstanu,msprpr,msprnu,
     &   idipi,iabso,iempi,ipiab,ipico,irescol,ipauli,dlife)
      write(*,*)' inecho utan'
*----------------------------------------------------------------------*
*
c     call zero_pot   !      is that needed ?? or wrong??
*
      if (ntotal .gt. maxpar) then
        write(isum,'(//''c:'',10x,
     &  ''**** fatal error: too many test part. ****'')')
        stop
      end if
*
*-----------------------------------------------------------------------
*   initialization of kinematical variables                            *
*                                                                      *
      call kinema(massta,masspr,mstapr,msprpr,mpion,elab,b,ipou,icoll,
     &                insys,rdist,radta,radpr,rmax,rxta,rzta,rxpr,rzpr,
     &       yref,yta,ypr,pxta,pzta,pxpr,pzpr,gammta,gammpr,betata,
     &       betapr,surfta,surfpr)
      write(*,*)'nach kinema',radta,radpr
      radius = amax1(radta,radpr)
*-----------------------------------------------------------------------
*   some constants for summry
      p0 = abs(pzta)
      if (abs(pzpr) .gt. abs(pzta)) p0 = abs(pzpr)
      punt=(p0+.60)/hbc/24.0*20.0
      factor=float(nint(punt))/20.0
      fact = 1./factor/hbc
*-----------------------------------------------------------------------
*   initialization of subroutine dens, pauli and relcol
      call denin(num)

      write(*,*)'nach denin'
      call paulin(num,numpaul,isum)
      write(*,*)'nach paulin'
c      call fixbm
      call resdyn(iresdalitz)
      write(*,*)'nach resdyn'
*=======================================================================
*
*   initialization of coulomb and yukawa potential
*
*=======================================================================
c      call yupri(alp,bet,rmc,rmy)
c     call initcoulen                  !   put into  time loop


c     write(*,*)'nach cypri'

*-----------------------------------------------------------------------
      ymim = 1.2* yta
      ymam = 1.2* ypr
      if(imeson.eq.1) call mesonin(ymim,ymam)
      write(*,*)'nach mesonin'
*=======================================================================
*   initialization of mesout                             zm
      write(*,*) 'imestimpri: ', imestimpri
      if(imestimpri.eq.1) call mesout_in
      write(*,*)'after mesout_in'
*=======================================================================
c   initialization for dilepton production:
      if(idilepton.eq.1)
     &  call init_dilepton(yta,ypr,ipiNbrems)
*=======================================================================
c   initialization for perturbative meson, JPsi production:
      if(i_JPsi.ne.0)
     &  call pert_meson_init
*=======================================================================
      ydel = ypr - yta
*
*=======================================================================
*
c      call sumbain(ibaspe,iparti,binkba,bintba,binmba,iframeb,
c     &                   nangleb,anba,dangleb,yref,ydel,icoll)
c      write(*,*)'nach sumbain'
*
*=======================================================================
c      if(idipi .ge.1)
c     &  call sumpiin(ipispe,binkpi,bintpi,binmpi,iframep,
c     &                   nanglep,anpi,danglep,yref,ydel)
*
c      write(*,*)'nach sumpiin'
*=======================================================================
*
*   initialization of time-loop variables
*
      do 1111 i0 = 0,200
        do 1101 j0 =-1,dimlcoll
          mcol(j0,i0)= 0
 1101   continue
        do k6 =1,50
        do k7 =-2,nres
          resdens(k7,i0,k6) = 0.0
        enddo
        enddo
        mend(i0)= 0
        medd(i0)= 0
        mepn(i0)= 0
        mept(i0)= 0
        do  j0=1,dimlmesc
          mesc(j0,i0)=0
          mesa(j0,i0)=0
        end do
        mdip(i0)= 0
        meta(i0)= 0
        mron(i0)= 0
        msin(i0)= 0
        mpan(i0)= 0
        mprt(i0)= 0
        mntr(i0)= 0
        mde2(i0)= 0
        mde1(i0)= 0
        mde0(i0)= 0
        mdem(i0)= 0
        mrpl(i0)= 0
        mrze(i0)= 0
        mqpl(i0)= 0
        mqze(i0)= 0
        do j0=1,15
          mumes(j0,i0) = 0
        end do
        mcnu(i0)= 0
        mcde(i0)= 0
        mcrn(i0)= 0
        mcqn(i0)= 0
        nlst(i0)= 0
        rati(i0)=0.0
        ener(i0)=0.0
        suru(i0)=0.0
        rootmsq(i0)=0.0
        do js=0,20
          ndis(js,i0)=0
        end do
 1111 continue
      ndis(0,0)=ntotal*isubs
      do 4012 is=1,8
      do 4011 js=0,50
        mlife(is,js) = 0
        mdens(is,js) = 0
        if(is.le.2) mpath(is,js) = 0
 4011 continue
 4012 continue
      do 4013 js=1,100
        mbirt(1,js) = 0
        mbirt(2,js) = 0
        mbirt(3,js) = 0
        mbirt(4,js) = 0
 4013 continue
      do 4014 js=1,20
        mposi(1,js) = 0
        mposi(2,js) = 0
        mposi(3,js) = 0
        mposi(4,js) = 0
 4014 continue
      do 4015 ij=1,15
        lumes(ij) = 0
 4015 continue

      mestopi=0
      mestomes = 0
      pitomes=0
      mcnne  = 0
      do ij=-1,dimlcoll
        mcoll(ij) = 0
      end do
      mcond  = 0
      mcopn  = 0
      mcbre  = 0
      mcopt  = 0
      mcodd  = 0
      do ij = 1,dimlmesc
        mmesc(ij) = 0
        mmesa(ij) = 0
      end do
      mppan  = 0
      mprot = 0
      mneut = 0
      mdel2 = 0
      mdel1 = 0
      mdel0 = 0
      mdelm = 0
      mrspl = 0
      mrsze = 0
      mqspl = 0
      mqsze = 0
      lcnu  = 0
      lcde  = 0
      lcrn  = 0
      lcqn  = 0
      time2(0)=0.0
      rhmax(0)=0.0
      rhavr(0)=0.0
      reovr(0)=0.0
      if(ipou.eq.1) then
        distp(0)=rmax
      else
        distp(0)=sqrt(rmax**2+b**2)
      end if
*
*=======================================================================
*   initialization and initial output for thermo
*
      write(*,*)'vor thermo'
c      if(ithermo.eq.1) then
c        write(mterpri,'(4i5,2f8.3)') massta,masspr,num,isubs,elab,b
c.......................................................................
c    to avoid that only npartcel testparticle per cell gives higher
c             occupation (distribution) than one: dpmin is the minimal
c             momentum step length
c.......................................................................
c        vol  = 4./3. * pi * delr**3 * float(num)
c        dpmin = hbc * (float(npartcel)/vol)**0.333333333
c        if(dpmin.gt.delp)
c     &  write(mterpri,'(''too small step in mom. new dp='',f6.3)') dpmin
c        delp = amax1(dpmin,delp)
c       write(mterpri,'(''c: max.num. particle per cell '',i5)') npartcel
c        write(mterpri,'(''c: space cell radius (fm) '',f8.3)') delr
c        write(mterpri,'(''c: momentum step length (gev) '',f8.3)') delp
c      endif
      write(*,*)'nach thermo'
*-----------------------------------------------------------------------
*       ============== loop over subsequent runs ===============       *
*

      if(icoll.ne.1) then
      call bwini
      end if
      call croswread
      write(*,*)'vor subs do loop'
      write(*,*)'isubs = ',isubs
      do ij = 1,2
         m_birth(ij) = 0.
         t_birth(ij) = 0.
         m2_birth(ij) = 0.
         m_final(ij) = .0
         t_final(ij) = .0
         m2_final(ij) = .0
         n_birth(ij) = 0
         n_final(ij) = 0
      enddo
      isubsqt = 0
      do 50000 isu = 1,isubs
c          write(46,*) ' isu loop ',isu
cr        rancheck = 1
c         abababa  = rn(iseed)
c         rancheck = 0
c       if (isu .eq. 1) then
c         if (abs(jseed) .gt. 10) iseed = jseed
c       else
c         kk = abs(iseed)
c         k2 = kk / 1000000
c         iseed = -abs(kk - 1000000 * k2)
cr      endif
c        write(46,*)'in isubs do loop random  ', massta,isu,iseed
      call f77flush()
         write(*,*)'in isubs do loop random  ', massta,isu,iseed
*
*=======================================================================
*   initialization of coordinate space and momentum distribut2ion
*
      isubsqt = isubsqt+1
      do 10 ii = 1,maxpar
        id(1,ii) = 0
  10  continue
      do  ii = 1, max_kminu
        nx_hyp(0,ii)   = 0
        nx_kminu(0,ii) = 0
        nx_kminu(4,ii) = 0
        p_hyp(4,ii)   = 0.
        p_kminu(4,ii) = 0.
      enddo
      do  ii = 1, max_ksi
        nx_ksi(0,ii)   = 0
        p_ksi(4,ii)   = 0.
      enddo
      do  ii = 1, max_epair
        nx_epair(0,ii)   = 0
        nx_epair(1,ii)   = 0
        nx_epair(2,ii)   = 0
        p_epair (4,ii)   = 0.
      enddo
      ikmi_abs = 0
      icollk   = 0
      icollh   = 0
      icollk_mi= 0
*
      if(icomp.ne.3 .and. icomp.ne.4 .and. icomp.ne.6) then
       dummy = mdpartini(rho0,icomp)
      end if

      write(*,*)'vor inir',massta   ,num     ,radta,
     &          rxta    ,rzta     ,
     &          iseed    ,mass    ,mstapr,  massta, surfta
*   target
      if(massta.ne.0)
     & call inir(1       ,massta   ,num     ,radta,
     &          rxta    ,rzta     ,
     &          iseed    ,mass    ,mstapr,  massta, surfta)
*

      write(*,*)'vor dens '
      call dens(1,massta,num,nlost)
*
      write(*,*)'vor inip', massta, rzta, pzta, pxta, gammta, betata
      jcoll=icoll
      if(massta.eq.1) jcoll=-1
      call inip(1, massta, rzta, pzta, pxta, gammta, betata)
*
c      write(19,*)'after target',1,massta
c      do ii=1,maxpar
c         write(19,*) ii,(jj,id(jj,ii),jj=1,6)
c      end do
*   projectile
      write(*,*)'vor inirp',mass   ,num     ,radpr,
     &          rxpr    ,rzpr     ,
     &          iseed    ,mass    ,msprpr,  masspr, surfpr
      if(masspr.ne.0) then
        call inir(1+massta,mass     ,num     ,radpr,
     &          rxpr    ,rzpr     ,
     &          iseed    ,mass    ,msprpr,  masspr, surfpr)
*
        call dens(1+massta,mass,num,nlost)
*
        jcoll=icoll
        if(masspr.eq.1) jcoll=-1
        call inip(1+massta, mass, rzpr, pzpr, pxpr, gammpr,betapr)
      end if
      write(*,*)'after projectile',1,massta,1+massta,mass
c      do ij=1,num
c      do ii=(ij-1)*maxb+1,(ij-1)*maxb+mass
c         write(*,*) ii,(jj,id(jj,ii),jj=1,6),(r(jj,ii),jj=1,3),
c     &    (p(jj,ii),jj=1,3)
c      end do
c      end do
*
      write(*,*)'nach inip'
      nproton = 0
      nneutron = 0
      npbar = 0
      nelse=0
      do ij= 1, maxpar
        if(id(1,ij).eq.1 .and. id(2,ij).eq.1) nproton=nproton+1 
        if(id(1,ij).eq.1 .and. id(2,ij).eq.0) nneutron=nneutron+1 
        if(id(1,ij).eq.-1 .and. id(2,ij).eq.-1) npbar=npbar+1 
        if((id(1,ij).eq.1 .and. (id(2,ij).gt.1 .or. id(2,ij).lt.0)) .or.
     &       (id(1,ij).eq.-1.and.id(2,ij).ne.-1).or.abs(id(1,ij)).gt.1)
     &     nelse=nelse+1 
      end do
      write(*,*)'proton=',nproton,' neutron=',nneutron,' antipr=',npbar,
     &     '  else=',nelse
      call dens(1,mass,num,nlost)
c      if(idenspr.eq.1) call densir(0.0,0)
c      if(isu.eq.1) call densir(0.0,0)
*------------------------------------------
c        check  r-p
c       write(77,*)  '  gammpr, gammta ',gammpr, gammta
c       do irun = 1,num
c       init     = (irun-1)*maxb
c       do i =1 ,mass
c       i1  = i + init
c         e1  = sqrt(rmass**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
c    1          -rmass
c         write(77, 7717) irun, i, id(1,i1), id(2,i1),
c    1        (r(m,i1),m=1,3),(p(m,i1),m=1,3), e1
c7717   format(2i4,2i2, 3f6.2,2x,3f7.3,2x,f7.3)
c       enddo
c       enddo
c       stop
************** initialize cross sections
c    -----------------------------------
c        new  initialisation for coulomb   tried by hw
c      call zero_pot   !
      call initcoulen
      call initemrpot
c    -----------------------------------

      t0 = 0.0
c      tin = secnds(t0)
      call mesdecini
c      tout = secnds(t0)
c     tin = tout - tin
c     write(*,*)'time ellapsed in main cross read = ',tin,'  sec.'

*=======================================================================
*
*            direct pion part (initialization for piondi)              *
      if(idipi .ge.1)
     &      call pioini(iabso,iempi,ipico,idec2pi,imeson,ivmesdil
     &      ,ipiNbrems,imesdalitz)
      if(mpion.ne.0)
     &              call inipion(mpion,num,rxpr,rzpr,pxpr,pzpr)
      if(mpion.eq.3)
     &  call gamrhoin(num,iseed,massta,elab,dt,dengamro,masspec,.025)
*=======================================================================
c   initialization for JPsi production:
        if(i_JPsi.ne.0 .or. i_phi.ne.0)
     &    call pert_meson_init_isu
*=======================================================================
*   initialization of coulomb and yukawa potential
c      call yuint(vzero,rmc,rmy)
*
      write(*,*)'nach yuint'
*=======================================================================
*
*   output of density, potential and momentum distribution
*     --- ground state properties and pauli test ---
c      if(masspr.eq.0.and.mpion.eq.0) then
c       call denden(1,num,radta,isum,0.0,ipot,ipou,rpot,tt0,tt3)
c       call ptest(num,iseed,isum,radta)
c      end if
*
*-----------------------------------------------------------------------
*   control printout of initial configuration
*
c      if(isumm.eq.1) then
c      write(*,*) 'summry elott'
c      call summry(isum,masspr,massta,mstapr,msprpr,num,mpion,icoll,
c     &            0,0.,0,time1,idipi,
c     &            ipot,ipou,rpot,tt0,tt3,
c     &            ptxyz,ppxyz,pixyz,rtxyz,rpxyz,rixyz,rmeson,
c     &            edif,eptt,epct,ettt,ektt,ectt,eitt,eftt,emtt,
c     &            etpi,ekpi,ecpi,eipi,eflo,eflpi,tij,tijpi,fptp,
c     &            etta,etpr,ekta,ekpr,ecta,ecpr,emta,empr,
c     &            eipr,eita,eket,eztt)
c      write(*,'(/''c:energy deviation, kinetic, pion'',
c     &'' potential, coulomb, flow, mass/nucl.'')')
c      write(*,'(''c:time[fm/c]'')')
c      write(*,'(
c     &''c:time    edif/n    kin e/n   pion e/n '',
c     &  ''pot  e/n  coul e/n  flow e/n  mass e/n'')')
c      write(*,'(
c     &''n:  x     n(edif)    n(kin),  n(pio),  '',
c     &  ''n(pot),0  n(cou),0  n(flo),0  n(mass)'')')
c      write(*,'(f7.2,7e10.3)') time1(0),
c     &edif(0),ektt(0),etpi(0)*rmeson(0),eptt(0),epct(0),eflo(0),emtt(0)
*
c      write(*,*) 'summry utan'
c      end if
*-----------------------------------------------------------------------
      if(ntmax.le.0) stop
c     initialisation of iphapc(0:10) and iphapt(0:10)
c     phasespace for delta decay with local density aproximation
c     and with testparticle summation
      do 556 i=0,9
        iphapc(i)=0
        iphapt(i)=0
556   continue
*
      write(*,*) 'amplread elott'
      call amplread
      write(*,*) 'amplread utan'

      if (isu .eq. 1) call print_spect
      write(*,*) 'print_spect utan'


      write(*,*)'vor timestep loop', ' seed ',iseed
*=======================================================================
*
*   initialization of time-loop variables
*
      do 17 ij=-1,dimlcoll
        ncoll(ij) = 0
  17  continue
      ncond  = 0
      ncodd  = 0
      ncopn  = 0
      ncbre  = 0
      ncopt  = 0
      do ij = 1, dimlmesc
        nmesc(ij) = 0
        nmesa(ij) = 0
      end do
      nppan  = 0
*-----------------------------------------------------------------------
*       ============== loop over all time steps ================       *
*                                                                      *
         ictime = 0
         jctime = 0
*
*=======================================================================
         time = 0.0
*
          itip = 1
          itet = 1
c          call sumbarq(time,enerc)
c          call sumpioq(time)
*=======================================================================
*

      do 10000 nt = 1,ntmax
        nthw = nt
c     ---  hw  -------
        time=dt*float(nt)
c         write(46,*)'in time loop', isu, nt, iseed
         write(*,*)'in time loop', isu, nt, iseed,rhb(0,0,0)/rho0
     &    ,rhob_4(0,0,0,-5),rhob_4(0,0,0,0),rhob_4(0,0,0,5)    
c      do ij=1,num
c      do ii=(ij-1)*maxb+1,(ij-1)*maxb+mass
c         write(*,*) ii,(jj,id(jj,ii),jj=1,6),(r(jj,ii),jj=1,3),
c     &    (p(jj,ii),jj=1,3)
c      end do
c      end do
*
c        call dens_4    !    is called   in relcol
*=======================================================================
        if (idenspr.eq.1 .and.
     &    ((nt/nfrba)*nfrba.eq.nt) .and. (time.ge.timinsp))
     &    call densir(time,1)

*
*=======================================================================
*
*     collision term
*

        if (icoll .ne. 1) then
*
*     particle production part
*
*
*=======================================================================
        call checkl
        write(*,*)'vor relcol', ' seed ',iseed
        call relcol(lcoll,nt)
        lcond = lcoll(-3)
        lcopn = lcoll(-2)
        lcodd = lcoll(-4)
        lcbre = lcoll(-6)
        lcopt = lcoll(-5)
        write(*,*)'nach relcol', iseed
        call checkl
        write(*,*)'nach check'
*
*=======================================================================
c      call mes_dilep ! called from piondi instead

*            direct pion part

        write(*,*)'vor piondi', isu,nt, idipi, mass, ' seed ',iseed
        if(idipi .eq.2)
     &      call piondi(num,dt,lmesc,lmesa,lppan,numes,lmesa2,
     &                  cres)
*
      write(*,*)'nach piondi', ' seed ',iseed

*=======================================================================
c      call pert_mesons !

*           pert mesons

        write(*,*)'vor pert_mesons', isu,nt,i_JPsi, mass, ' seed ',iseed
        dt0=dt
        if(i_JPsi .ge.1.or.i_phi.ge.1) call pert_meson(nt,dt0,isu,yref)
*
      write(*,*)'nach pert_mesons', ' seed ',iseed

*=======================================================================
*            direct kaon part

c        write(*,*)'vor kaoncoll', kanum, icollk
c        if(ikaondi.ge.1.and.ikaondi.ne.5) then
c            call kaoncoll(collkplus,collkplusd,colkpl)
c        endif
c        write(44,*)'  vor  kminu_prod   '
!        if (i_kminu .eq. 1 .and. iphi_dec  .ne. 1)  call kminu_pi_prod
!        if (i_kminu .eq. 1 .and. iphi_dec  .ne. 1)  call hyp_coll
c        if (i_kminu .eq. 1)  call kminu_pi_prod(kmi_hpi)
*=======================================================================
*            Kminus-baryon elastic scattering + absorption
      write(*,*)  ' vor kminu_coll ', nt, mass
c      if( (ikaon.ge.1 .or. imeson.eq.1) .and. i_kminu.eq.1 .and.
c     & i_kminu_coll.eq.1 ) call kminu_coll(collkminus,collkminusd,
c     &                             colkmi,abskminus,abskminusd,abskmi)
*
      write(*,*) 'after kminu_coll', mass

*=======================================================================
*
c        if ( ((nt/nfrba)*nfrba .eq. nt) .and. (time .ge. timinsp)) then
c          itip = itip + 1
c          itip = min0(10,itip)
c          call sumbary(time,ntotal,nproton)
*
c          write(*,*)'vor idipi loop'
c      if(idipi .eq.2) then
c          call sumpion(time,numes(1),numes(2),numes(3))
*
c          write(*,*)'in idipi loop nfrpira = ',nfrpira
c        if (((nt/nfrpira)*nfrpira.eq. nt) .and. (time .ge.timinsp))
c     &  call pidensxz(time)
c      end if
c      end if

*=======================================================================
c      if ( (nt/nfreq)*nfreq .eq. nt ) then
c          itet = itet + 1
c          itet = min0(50,itet)
c          call sumbarq(time,enerc)
c          if(idipi .eq.2) then
c            call sumpioq(time)
c          end if
c      end if
      write(*,*)'vor deltadec', ' seed ',iseed
*=======================================================================
*
          do 11 ij = -1,dimlcoll
          ncoll(ij) = ncoll(ij) + lcoll(ij)
  11      continue
          ncond = ncond + lcond
          ncodd = ncodd + lcodd
          ncopn = ncopn + lcopn
          ncbre = ncbre + lcbre
          ncopt = ncopt + lcopt

          do  ij = 1,dimlmesc
            nmesc(ij) = nmesc(ij) + lmesc(ij)
            nmesa(ij) = nmesa(ij) + lmesa(ij)
          end do

          npidi = nmesc(1)+nmesc(2)+nmesc(3)+nmesc(4)+nmesc(5)+nmesc(6)
     &           -nmesa(1)-nmesa(2)-nmesa(3)-nmesa(4)-nmesa(5)-nmesa(6)
     &           +2*(nmesc(12)-nmesa(12)+nmesc(13)-nmesa(13))
          neeta = nmesc(7) - nmesa(7)
          nerho = nmesc(8)-nmesa(8)+nmesc(9)-nmesa(9)-nmesc(12)+
     &            nmesa(12)
          nesig = nmesc(10)-nmesa(10)+nmesc(11)-nmesa(11)-nmesc(13)+
     &            nmesa(13)
          nppan = nppan + lppan
*
          write(*,*)'main stop 1'

*
*=======================================================================
           icell = nt/ntcell
           if(icelpri.eq.1 .and. icell*ntcell.eq.nt .and. icell.le.10)
     &                    then
            mcelpri = ncelpri + (isu-1)*10 + icell
            write(6,*) mcelpri
            open(mcelpri)
            write(mcelpri,'(4i5,2f8.3)') massta,masspr,num,isubs,elab,b
c            write(mcelpri,'(3f8.3)') tt0,tt3,rpot
            write(mcelpri,'(f6.2,/)') time
*  id(6,i) = abs: number of baryon collision; sign: target or projectile
            do i =1,maxpar
               if(id(1,i).ne.0 .and.
     &            r(1,i)**2+r(2,i)**2+r(3,i)**2.lt.cell**2) then
            write(mcelpri,'(2i5,7f10.4)') id(1,i),id(6,i),
     &                         (r(ik,i),ik=1,3),e(i),(p(ik,i),ik=1,3)
             endif
             enddo
            close(mcelpri)
           endif
*=======================================================================
c------  store the counting variables       ----------------------------
*

           if ( (nt/nfreq)*nfreq .eq. nt ) then

           ictime=ictime+1
           time2(ictime)=time
           do 12 ij=-1,dimlcoll
             mcol(ij,ictime)=mcol(ij,ictime)+ncoll(ij)
  12       continue
           mend(ictime)=mend(ictime)+ncond
           medd(ictime)=medd(ictime)+ncodd
           mepn(ictime)=mepn(ictime)+ncopn
           mept(ictime)=mept(ictime)+ncopt

           do ij=1,dimlmesc
             mesc(ij,ictime) = mesc(ij,ictime)+nmesc(ij)
             mesa(ij,ictime) = mesa(ij,ictime)+nmesa(ij)
           end do
           mdip(ictime)=mdip(ictime)+npidi
           meta(ictime)=meta(ictime)+neeta
           mron(ictime)=mron(ictime)+nerho
           msin(ictime)=msin(ictime)+nesig
           mpan(ictime)=mpan(ictime)+nppan
           do  23 ij=1,15
             mumes(ij,ictime) = mumes(ij,ictime) + numes(ij)
 23        continue

           nprot = 0
           nneut = 0
           ndel2 = 0
           ndel1 = 0
           ndel0 = 0
           ndelm = 0
           nrspl = 0
           nrsze = 0
           nqspl = 0
           nqsze = 0
           ncnu  = 0
           ncde  = 0
           ncrn  = 0
           ncqn  = 0
           write(*,*)'ictime = ',ictime, ' seed ',iseed
           do 4999 i =1,maxpar
             if(id(1,i).eq.0) goto 4999
             if((id(1,i).eq.1).and.(id(2,i).eq. 1)) nprot = nprot + 1
             if((id(1,i).eq.1).and.(id(2,i).eq. 0)) nneut = nneut + 1
             if((id(1,i).eq.2).and.(id(2,i).eq. 2)) ndel2 = ndel2 + 1
             if((id(1,i).eq.2).and.(id(2,i).eq. 1)) ndel1 = ndel1 + 1
             if((id(1,i).eq.2).and.(id(2,i).eq. 0)) ndel0 = ndel0 + 1
             if((id(1,i).eq.2).and.(id(2,i).eq.-1)) ndelm = ndelm + 1
             if((id(1,i).eq.3).and.(id(2,i).eq. 1)) nrspl = nrspl + 1
             if((id(1,i).eq.3).and.(id(2,i).eq. 0)) nrsze = nrsze + 1
             if((id(1,i).eq.4).and.(id(2,i).eq. 1)) nqspl = nqspl + 1
             if((id(1,i).eq.4).and.(id(2,i).eq. 0)) nqsze = nqsze + 1
c************** if none of them
c             if((id(1,i).gt.4).or. (id(2,i).lt.-1) .or. (id(2,i).gt.2)
c     &      .or.(id(1,i).ne.2 .and.(id(2,i).eq.-1.or.id(2,i).eq.2)) )
c     &         write(6,'(''hiba:main'',5i6)') id(1,i),id(2,i),i,isu,nt
             rootmsq(ictime) = rootmsq(ictime)
     &                              +r(1,i)**2+r(2,i)**2+r(3,i)**2
             if(r(1,i)**2+r(2,i)**2+r(3,i)**2 .gt. 4)          goto 4999
             if(id(1,i).eq.1) ncnu = ncnu + 1
             if(id(1,i).eq.2) ncde = ncde + 1
             if(id(1,i).eq.3) ncrn = ncrn + 1
             if(id(1,i).eq.4) ncqn = ncqn + 1
 4999      continue
           do 4997 i=1,maxpar
chw                  --- not interested in id > nres  causing index overflow
             if (id(1,i).gt. nres) goto   4997
             if (id(1,i).eq. 0) goto   4997
c
             ix = nint( r(1,i) )
             iy = nint( r(2,i) )
             iz = nint( r(3,i) )
             dendel = 0.
             if(abs(ix).lt.maxx.and.abs(iy).lt.maxx.and.abs(iz).lt.maxz)
     &         dendel = rhb(ix,iy,iz)/rho0
             iendel = nint(5.0 * dendel) + 1
             if(iendel.gt.50) iendel = 50
             resdens(id(1,i)-1,0,iendel) =
     &             resdens(id(1,i)-1,0,iendel) + 1.0
             resdens(id(1,i)-1,ictime,iendel) =
     &             resdens(id(1,i)-1,ictime,iendel) + 1.0
 4997      continue
           mprt(ictime)=mprt(ictime)+nprot
           mntr(ictime)=mntr(ictime)+nneut
           mde2(ictime)=mde2(ictime)+ndel2
           mde1(ictime)=mde1(ictime)+ndel1
           mde0(ictime)=mde0(ictime)+ndel0
           mdem(ictime)=mdem(ictime)+ndelm
           mrpl(ictime)=mrpl(ictime)+nrspl
           mrze(ictime)=mrze(ictime)+nrsze
           mqpl(ictime)=mqpl(ictime)+nqspl
           mqze(ictime)=mqze(ictime)+nqsze
           mcnu(ictime)=mcnu(ictime)+ncnu
           mcde(ictime)=mcde(ictime)+ncde
           mcrn(ictime)=mcrn(ictime)+ncrn
           mcqn(ictime)=mcqn(ictime)+ncqn
           if(mcnu(ictime).gt.0)
     &       rati(ictime)=float(mcde(ictime)+mcrn(ictime)+mcqn(ictime))
     &                   /float(mcnu(ictime))
*
           do 4000 is = 1,maxpar
             if(id(1,is).gt.0) then
               js = abs(id(6,is))-1
               if (js .le. 19)
     &         ndis(js,ictime)=ndis(js,ictime) + 1
             endif
 4000      continue
           endif
        endif
*
*     update momenta
************************************************************************
        write(*,*)'vor propa', ' seed ',iseed
*     do complete propagation  *************************

        call propa( nt)

        write(*,*)'nach propa ',nthw, ' seed ',iseed


********************************************************
*
      if (icoll.ne.-1) then
*
*     update density
*
        ipavp=0
        if ( (nt/nfreq)*nfreq .eq. nt ) ipavp=1
      call dens(1,mass,num,nlost)
      write(*,*) 'after dens in main ',nthw,rhb(0,0,0),enerc,nfreq,rho0
      suru(ictime) = suru(ictime)+rhb(0,0,0)/rho0/float(isubs*nfreq)
      ener(ictime) =ener(ictime)+enerc/float(nfreq)
*
*     yukawa and coulomb
*
c      if(ipot.eq.1) then
c        call yukw(2*nity,rmy,vzero)
c      end if
      end if
*-----------------------------------------------------------------------
*
*       control-printout of configuration (if required)
*
c       if (((nt/nfreq)*nfreq .eq. nt ).and.(isumm.eq.1)) then
c         jctime=jctime+1
c         jctime=min(jctime,200-2)
c         nlst(jctime)=nlost-nlst(jctime+1)
c         nlst(jctime+2)=nlost
c         ilast=0
c         if(nt+nfreq.gt.ntmax) ilast=1
c         write(*,*) 'summry elott'
c         call summry(isum,masspr,massta,mstapr,msprpr,num,mpion,icoll,
c     &               ilast,time,jctime,time1,idipi,
c     &               ipot,ipou,rpot,tt0,tt3,
c     &               ptxyz,ppxyz,pixyz,rtxyz,rpxyz,rixyz,rmeson,
c     &               edif,eptt,epct,ettt,ektt,ectt,eitt,eftt,emtt,
c     &               etpi,ekpi,ecpi,eipi,eflo,eflpi,tij,tijpi,fptp,
c     &               etta,etpr,ekta,ekpr,ecta,ecpr,emta,empr,
c     &               eipr,eita,eket,eztt)
c      write(*,'(/''c:energy deviation, kinetic, pion'',
c     &'' potential, coulomb, flow, mass/nucl.'')')
c      write(*,'(''c:time[fm/c]'')')
c      write(*,'(
c     &''c:time    edif/n    kin e/n   pion e/n '',
c     &  ''pot  e/n  coul e/n  flow e/n  mass e/n'')')
c      write(*,'(
c     &''n:  x     n(edif)    n(kin),  n(pio),  '',
c     &  ''n(pot),0  n(cou),0  n(flo),0  n(mass)'')')
c      write(*,'(f7.2,7e10.3)') time1(jctime),
c     &edif(jctime),ektt(jctime),etpi(jctime)*rmeson(jctime),
c     &eptt(jctime),epct(jctime),eflo(jctime),emtt(jctime)
c         write(*,*) 'summry utan'
c      end if
*=======================================================================

c zm:  printout of mesons
       write(*,*) 'Before mesout', ' seed ',iseed
       call mesout

        eganz = 0.0
         do i = 1, maxpar
           if(id(1,i).ne.1) then
           etot = sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
           eganz = eganz + etot
           end if
         end do
*-----------------------------------------------------------------------
         write(*,*)'etotal = ',nt,eganz, eganz/float(num)
         cputim = second()
         cputim = cputim - cputim0
         if(cputim.lt.cptimepr) icptimep = icptimep + 1
         cptimepr = cputim
         cputim = cputim + float(icptimep) * clockmax

         write(*,*)' end of timestep ',isu,nt
         write(*,*)
         write(*,*)
         if(cputim.gt.timecpu) goto 60000

      nod = 0
      non = 0
      no1 = 0
      no2 = 0
      npi = 0
      neta = 0

      do i = 1, maxpar
        if(id(1,i).eq.1) non = non + 1
        if(id(1,i).eq.2) nod = nod + 1
        if(id(1,i).eq.3) no1 = no1 + 1
        if(id(1,i).eq.4) no2 = no2 + 1
      end do

      do i = 1, maxppar
        if(ipi(1,i).eq.1) npi = npi +1
        if(ipi(1,i).eq.2) neta = neta +1
      end do

*      write(117,117)nt, non, nod, no1, no2, npi, neta

c 117  format(7(1x,i5))

      write(*,*) ' end of time step loop', isu, nt, mass,
     1             ' seed ',iseed
      call f77flush()
10000 continue
cc----------------  new hw -----------------
      call prohwout (yref, yta, ypr,isu)
      call pihwout (yref, yta, ypr,isu)
cc----------------  new hw -----------------
c     write(*,*)  ' in main 8801 K1', ipi(1,8801)

*       ==============  end of time step loop   ================
*
      write(*,*)'vor piondif'
      call piondif(num,lmesc,numes,imesdalitz,igamma,
     &               lmesa2,cres)
      write(*,*)'nach piondif'
*=======================================================================
c      call pert_mesons !

        write(*,*)'vegso pert_mesons', isu,nt,i_JPsi, mass
        dt0=1000.0
        if(i_JPsi .ge.1 .or.i_phi.ge.1) call pert_meson(nt,dt0,isu,yref)
*
      write(*,*)'nach pert_mesons', ' seed ',iseed
*=======================================================================


      do ij=-1,dimlcoll
         mcoll(ij) = ncoll(ij) + mcoll(ij)
      enddo
          mcond = ncond + mcond
          mcopn = ncopn + mcopn
          mcbre = ncbre + mcbre
          mcopt = ncopt + mcopt
          mcodd = ncodd + mcodd

          do ij=1,dimlmesc
            mmesa(ij) = nmesa(ij)+mmesa(ij)
            mmesc(ij) = nmesc(ij)+mmesc(ij)
          enddo
          do ij=1,15
            lumes(ij) =numes(ij)+lumes(ij)
          enddo

          mpidi = npidi + mpidi
          meeta = neeta + meeta
          merho = nerho + merho
          mesig = nesig + mesig
          mppan = nppan + mppan
          mprot = mprot + nprot
          mneut = mneut + nneut
          mdel2 = mdel2 + ndel2
          mdel1 = mdel1 + ndel1
          mdel0 = mdel0 + ndel0
          mdelm = mdelm + ndelm
          mrspl = mrspl + nrspl
          mrsze = mrsze + nrsze
          mqspl = mqspl + nqspl
          mqsze = mqsze + nqsze
          lcnu  = lcnu  + ncnu
          lcde  = lcde  + ncde
          lcrn  = lcrn  + ncrn
          lcqn  = lcqn  + ncqn


      itip = itip + 1
      itip = min0(10,itip)
      itet = itet + 1
      itet = min0(50,itet)

       if(idenspr.eq.1) call densir(time,1)
*
*
      if(idenspr.eq.1) call densout
*
*=======================================================================
*
*=======================================================================
      if(inuclpri .eq. 1) then
        nnuclpr  = max0(nnuclpri-(isu-1)*num,0)
        nnuclpr  = min0(nnuclpr,num)
        write(mnucpri,'(4i8,f10.6)') massta,mstapr,masspr,msprpr,elab
        write(mnucpri,'(2i8,2f10.6)') isubs,nnuclpr,b,yref
        do 6112 ii = 1,nnuclpr
        write(mnucpri,'(i8)') mass
        do 6111 jj = (ii-1)*maxb+1,ii*maxb
          if(id(1,jj).ne.0) then
          tram = sqrt(e(jj)**2+p(1,jj)**2+p(2,jj)**2)
          etot = sqrt(tram**2+p(3,jj)**2)
          ylab = yref + 0.5 * log((etot+p(3,jj))/(etot-p(3,jj)))
c id 1 - nukleon, id 2 - töltés, e - tömeg, p1, p2, p3 impulzusok lab rendszerben, id 6 ütközési szám
         write(mnucpri,'(2i3,7f10.6,i4)') id(1,jj),id(2,jj),e(jj),
     &         r(1,jj), r(2,jj), r(3,jj),
     &         p(1,jj),p(2,jj),tram*sinh(ylab),id(6,jj)
c     r(1,jj), r(2,jj), r(3,jj) - helykoordináta
          end if
 6111   continue
 6112   continue
      end if
      write(*,*)'nach if loof'
*=======================================================================
      if(ipionpri .eq. 1.or.ipionpri.eq.3) then
        npionpr  = max0(npionpri-(isu-1)*num,0)
        npionpr  = min0(npionpr,num)
        write(mpiopri,'(4i8,f10.6)') massta,mstapr,masspr,msprpr,elab
        write(mpiopri,'(2i8,2f10.6)') isubs,npionpr,b,yref
        do 6115 ii = 1,npionpr
        nii = 0
        do 6113 jj = (ii-1)*maxp+1,ii*maxp
          if(ipi(1,jj) .ne. 0) nii = nii + 1
 6113   continue
        write(mpiopri,'(i8)') nii
        do 6114 jj = (ii-1)*maxp+1,ii*maxp
          if(ipi(1,jj) .ne.0) then
            wmass = epi(jj)
            tram = sqrt(wmass**2+ppi(1,jj)**2+ppi(2,jj)**2)
            etot = sqrt(tram**2+ppi(3,jj)**2)
            ylab = yref + 0.5 * log((etot+ppi(3,jj))/(etot-ppi(3,jj)))
            write(mpiopri,'(i3,5f10.6)') ipi(2,jj),wmass,ppi(1,jj),
     &                                      ppi(2,jj),tram*sinh(ylab)
          end if
 6114   continue
 6115   continue
        endif

*=====================================================================*
      if(ideltpri .eq. 1) then
        if(isu.eq.1) then
          write(mdelpri,'(4i8,f10.6)') massta,mstapr,masspr,msprpr,elab
          write(mdelpri,'(2i8,4f10.6)') isubs,npionpr,b,yref,gammta
        endif
        do ii = 1,num
          nii = 0
          do jj = (ii-1)*maxp+1,ii*maxp
             if(ipi(1,jj).ne.0) then
              if(ipi(3,jj).ne.0 .and.jj.eq.id(8,ipi(3,jj)).and.
     &           ipi(6,jj).eq.id(6,ipi(3,jj)) .and.
     &           ipi(7,jj).eq.id(5,ipi(3,jj)) ) nii = nii + 1
            endif
          enddo
          write(*,*) 'observable delta number: nii',nii
          write(mdelpri,'(i8)') nii
          do jj = (ii-1)*maxp+1,ii*maxp
            if(ipi(1,jj).ne.0) then
              if(ipi(3,jj).ne.0 .and.jj.eq.id(8,ipi(3,jj)).and.
     &           ipi(6,jj).eq.id(6,ipi(3,jj)) .and.
     &           ipi(7,jj).eq.id(5,ipi(3,jj)) ) then
                wmass = epi(jj)
                tram = sqrt(wmass**2+ppi(1,jj)**2+ppi(2,jj)**2)
                etot = sqrt(tram**2+ppi(3,jj)**2)
                ylab = yref+ 0.5*log((etot+ppi(3,jj))/(etot-ppi(3,jj)))
                jn = ipi(3,jj)
                bmass = e(jn)
                tramb = sqrt(bmass**2+p(1,jn)**2+p(2,jn)**2)
                etot = sqrt(tramb**2+p(3,jn)**2)
                ylabb = yref + 0.5*log((etot+p(3,jn))/(etot-p(3,jn)))
                write(mdelpri,'(2i3,4f10.6,2i3,4f10.6,i6,6f10.6)')
     &            ipi(1,jj),ipi(2,jj),wmass,ppi(1,jj),ppi(2,jj)
     &            ,tram*sinh(ylab),id(1,jn),id(2,jn),e(jn),p(1,jn),
     &           p(2,jn),tramb*sinh(ylabb),id(6,jn),(rpie(ij,jj),ij=1,6)
              endif
            end if
          enddo
        enddo
      endif

*=======================================================================
*
c      call pert_mesons !

        write(*,*)'vegso pert_meson out', isu,nt,i_JPsi
        if(i_JPsi .ge.1) call pert_meson_store
*
      write(*,*)'nach pert_meson out', ' seed ',iseed
*=====================================================================*
*     propagate the charged pions in the coulomb-field that is left   *
*                                                                     *
      call mesfin
*                                                                     *
*=====================================================================*
*                                                                     *
*    Do the fserat-output                                             *

        if( ifserat.eq.1.and. isu.eq.1 ) then

        if(isu .eq. 1 ) then
          write(mfserat,91) isubs, num, insys
          write(mfserat,92)massta, mstapr, masspr, msprpr
          write(mfserat,93) elab, b,-betata, gammta
          write(mfserat,95) ipou, ncont, icomp
          write(mfserat,96) ipipot, ireapl
        end if


         do irun = 1, num
            write(mfserat,*)irun
           do  j1 = 1,maxb
              jj  = j1 + (irun - 1) * maxb
              if(id(1,jj).ne.0) then
              write(mfserat,90)e(jj), p(1,jj), p(2,jj), p(3,jj)
              write(mfserat,89) upot(jj), id(1,jj), id(2,jj),id(6,jj)
              end if
           end do
         end do
         end if

 90      format(4e16.8,2i4)
 89      format(1e16.8,i8,i8,i8)

*     end for fserat-output
*=====================================================================*
*                                                                     *
*    Do the fserat-output                                             *

        if( ifserat.eq.2) then

        if(isu .eq. 1 ) then
          write(mfserat,91) isubs, num, insys
          write(mfserat,92)massta, mstapr, masspr, msprpr
          write(mfserat,93) elab, b,-betata, gammta
        end if

         do irun = 1, num
            write(mfserat,*) (isu-1)*num+irun
           do  j1 = 1,maxb
              jj  = j1 + (irun - 1) * maxb
              if(id(1,jj).ne.0)
     &        write(mfserat,90)e(jj), p(1,jj), p(2,jj), p(3,jj),id(2,jj)
           end do
           do  j1 = 1,maxp
             jj  = j1 + (irun - 1) * maxp
             if(ipi(1,jj).gt.0)
     &        write(mfserat,90)epi(jj), ppi(1,jj), ppi(2,jj), ppi(3,jj)
     &          ,ipi(2,jj), ipi(5,jj)
           end do
         end do
         end if

*     end for fserat-output
*======================================================================*
*                                                                     *
*    Do the pion-output                                               *

        if( icpipri.eq.1 ) then

        if(isu .eq. 1 ) then
          write(mcpipri,91) isubs, num, insys
          write(mcpipri,92) massta, mstapr, masspr, msprpr
          write(mcpipri,93) elab, b,-betata, gammta
          write(mcpipri,95) ipou, ncont, icomp
        end if
 91     format(3i8)
 92     format(4i8)
 93     format(4f10.5)
c 94     format(2i4,4e16.8,i4)
 95     format(3i8)
 96     format(2i8)
c 97     format(7e16.8)


        npionpr  = max0(npionpri-(isu-1)*num,0)
        npionpr  = min0(npionpr,num)

        do  ii = 1,num

          napart = 0
          do nj1 = 1,maxb
            ni1  = nj1 + (ii-1)*maxb
            if(id(1,ni1).ne.0 .and. abs(id(6,ni1)).gt.1) then
              napart = napart + 1
            end if
          end do

          nii = 0
          do jj = (ii-1)*maxp+1,ii*maxp
            if(ipi(1,jj) .ne. 0) nii = nii + 1
          end do
          write(mcpipri,96) nii, napart
          do jj = (ii-1)*maxp+1,ii*maxp
            if(ipi(1,jj).ne. 0) then
              wmass = epi(jj)
              en = ppi(1,jj)**2+ppi(2,jj)**2+ppi(3,jj)**2
              en = en + wmass**2
              en = sqrt(en)

            end if
          end do
        end do
      endif
*     end for fserat-output
*======================================================================*

         cputim = second()
         cputim=cputim - cputim0
         if(cputim.lt.cptimepr) icptimep = icptimep + 1
         cptimepr = cputim
         cputim = cputim + float(icptimep) * clockmax
         cpuisuloop = cputim/float(isu)
         write(*,'(''c: end of isubs loop, time for a loop:'',f15.2)')
     &      cpuisuloop
         write(*,'(''c: cputime spent'',2i6,f15.2,i20)')
     &     isu, nt-1, cputim

         write(*,*)
         ncatas  = 0
         if(isu.lt.isubs .and. cputim+cpuisuloop.gt.timecpu) then
           ncatas   = 1
           isubsold = isubs
           isubs    = isu
         endif
*
      if (ncatas .eq. 1)  then
        isubs = isubsold
        goto 60000
      endif
50000 continue  ! end of isubs loop

60000 continue
      
      if(isu.lt.isubs) write(*,*) 'hiba0 isu<isubs',isu,isubs
      if(b .lt.0.01) then
        totcros= radius**2*pi*10.
      else
        totcros= 2.*b*pi*10.
      end if
      eventnum = float(isubsqt*num)
      scala  = totcros
c      if(masspr.ne.0) scala  = totcros/float(massta*masspr)
      scalpro= totcros/eventnum

      write(*,*)'vor end of job',isu,isubs,isubsqt,num,eventnum

      if(idilepton.eq.1) then
        call dilepout(scalpro,yref)
c        call el_pos_pair_out(yref, yta, ypr,isu, 20)
c        call epair_out (yref, yta, ypr,isu)
      end if
*
*=======================================================================
*
c      call pert_mesons !

        write(*,*)'vegso pert_meson out', isu,nt,i_JPsi
        if(i_JPsi .ge.1.or.i_phi.ge.1) call pert_meson_out
*
      write(*,*)'nach pert_meson out', ' seed ',iseed
*=======================================================================
      if(imeson .eq. 1)
     &  call mesonout(scalpro,yref)
*=======================================================================

      call inputwri(0)

*       ===================  output  ==========================
      if(ikaopri .eq. 1 .and. ikaon .eq. 1) then
        close(mkaopri)
      end if
*=======================================================================
      if(imespri .eq. 1 .and. imeson .eq. 1) then
        write(mmespri,'(4i8,f10.6)') massta,mstapr,masspr,msprpr,elab
        write(mmespri,'(2f10.6,e12.2,i5)') b,yref,sigpiom,icomp
        call finalmes(eventnum)
        close(mmespri)
      end if
*=======================================================================
cc      if(ippipri*ipertpi .eq. 1) then
cc        write(mppipri,'(4i8,f10.6)') massta,mstapr,masspr,msprpr,elab
cc        write(mppipri,'(f10.6)') b
cc        call finalpip
cc        close(mppipri)
cc      end if
      write(isum,'(//''c:total cpu time: '')')
      cputime = second() - cputim0
      if(cputim.lt.cptimepr) icptimep = icptimep + 1
      cputim = cputim + float(icptimep) * clockmax
      write(isum,'(''c: cpu time spent: isu,nt'',2i6,f15.2)')
     &     isu-1, nt-1, cputim


*
*=======================================================================
*

      write(*,'(''c: cpu time spent'',2i6,f15.2)')
     &     isu-1, nt-1, cputim
      if(isubs .eq. isu-1) write(*,*) 'vegemvan'
      if(isubs .gt. isu-1) write(*,*) 'idoelottvegemvan'
      if(isubs .eq. isu-1) write(isum,*) 'vegemvan'
      if(isubs .gt. isu-1) write(isum,*) 'idoelottvegemvan'

      stop
      end
c===========================================================
      real*8 function second()
      real*8 t
      call cpu_time(t)
      second = t
      return
      end
c=====================================================
