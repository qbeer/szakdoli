************************************************************************
*                                                                      *
      subroutine piondi(num,dt,lmesc,lmesa,lppan,numes,lmesa2,
     &                  cres)
*                                                                      *
************************************************************************
      implicit none
      include"common"
      logical gridte, flag,final
      integer dimlmesc
      parameter(dimlmesc = nres+9)
      integer lmesc(1:dimlmesc), lmesa(1:dimlmesc), numes(18)
      integer i
      integer jj,dummy(1:100)
      integer iempi, iabso, ipico, imeson, idec2pi
      integer lpiro, lppan,lpisi, lropi, lprom, lompi, lmesa2, ile, idj
      integer ixx, iyy, izz, iendel, il, ibir, ipos,km, jd
      integer jabso, jdec2pi, jmeson, num, inp
      integer jpico, imesdalitz, jmesdalitz,lsipi, jempi, mpion
      integer jvmesdil, ivmesdil, igamma,jpiNbrems,ipiNbrems
      real*8    dt, pos, dendel, x0, z0, px0, pz0
      real*8    phi,rap,gamma,probpi0,dens, ee ,pt
      integer cres(1:nres,1:9)
*----------------------------------------------------------------------*
*
*  id(1,i) = particle type (nucleon,2:delta, n*...nres+2:Lambda,Sigma)
*  id(2,i) = particle charge
*  id(3,i) = last colliding partner of particle i
*  id(4,i) = for resonances: how many times the meson was created producing i
*  id(5,i) = number of pion absorption
*  id(6,i) = abs: number of baryon coll.; sign:target or proj.;
*  id(7,i) = type of parent meson (>0), or nucl-nucl coll. (0) for res.
*  id(8,i) = the serial number of the daughter meson
*
*----------------------------------------------------------------------*
*
*  ipi(1,i) = particle type (pion, eta, ...)
*  ipi(2,i) = particle charge
*  ipi(3,i) = the serial number for the parent resonance of meson i
*  ipi(4,i) = how many times the meson was created
*  ipi(5,i) = type of the parent resonance (id(1,ipi(i,3)))),sstate:1,mesondec-
*  ipi(6,i) = collision number of the parent resonance;
*  ipi(7,i) = pion absorption number of the parent resonance
*  ipi(8,i) = in case of state prod. the serial number of the 2. father
*
*   the anti p + p -> JPsi + meson is not completely correct: ipi(5,...)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      save  iabso,iempi,ipico,idec2pi,imeson,ivmesdil,ipiNbrems
     &      ,imesdalitz
*----------------------------------------------------------------------*
      write(*,*)'vor pionem', iempi
      call f77flush()
      call checkl
      flag =.false.
      if(iempi .eq. 1)     call pionem(dt,lmesc,flag,cres)

      write(*,*) 'in piondi calling mes_dilep'
      if(ivmesdil.eq.1) then
        call mes_dilep
      endif

      if(ivmesdil.eq.1 .or. imesdalitz.eq.1) then
        call mes_dalitz(flag)
      endif

      if(ipiNbrems.eq.1) call piNdilepRun(dt)

      write(*,*)'nach pionem-1', dt
      write(*,*)'nach pionem', iabso
      call f77flush()
      call checkl
      write(*,*)'nach checkl pionab', iabso
      if(iabso .eq. 1)     call pionab(lmesa)

      call checkl
      write(*,*)'vor pionco'
      if(ipico.ge.1 .or. imeson.eq.1 .or. idec2pi.ge.1)
     &                     call pionco(lppan,lpiro,lpisi,lprom)

      write(*,*)'vor mesdec'
      call f77flush()
      call checkl
      flag = .false.
      if(idec2pi.ge.1)   call mesdec(lropi,lsipi,lompi,flag,lmesa2)

      write(*,*)'nach mesdec'
      call f77flush()
*
      call checkl
      write(*,*)
      write(*,*)'writeout in piondi '
      call f77flush()

      lmesc(12) = lropi
      lmesc(13) = lsipi
      lmesc(14) = lompi
      lmesa(12) = lpiro
      lmesa(13) = lpisi
      lmesa(14) = lprom
      do 10 i = 1,18
        numes(i) = 0
  10  continue
*
*     update positions and momenta in propa
*
      do 100 i = 1,maxppar
        if(ipi(1,i) .eq. 0)                                     goto 100
*
        if(ipi(1,i) .eq. 1) then                         ! pion
          if(ipi(2,i).eq. 1) numes(1) = numes(1) + 1
          if(ipi(2,i).eq. 0) numes(2) = numes(2) + 1
          if(ipi(2,i).eq.-1) numes(3) = numes(3) + 1
        end if
        if(ipi(1,i) .eq. 2)  numes(4) = numes(4) + 1     ! eta
        if(ipi(1,i) .eq. 3) then                         ! rho
          if(ipi(2,i).eq. 1) numes(5) = numes(5) + 1
          if(ipi(2,i).eq. 0) numes(6) = numes(6) + 1
          if(ipi(2,i).eq.-1) numes(7) = numes(7) + 1
        end if
        if(ipi(1,i) .eq. 4)  numes(8) = numes(8) + 1     ! sigma
        if(ipi(1,i) .eq. 5)  numes(9) = numes(9) + 1     ! omega
        if(ipi(1,i) .eq. 6)  numes(10) = numes(10) + 1   ! kaon
        if(rpi(1,i)**2+rpi(2,i)**2+rpi(3,i)**2 .lt. 4) then
          if(ipi(1,i).eq. 1) numes(11)  = numes(11) + 1
          if(ipi(1,i).eq. 2) numes(12)  = numes(12) + 1
          if(ipi(1,i).eq. 3) numes(13)  = numes(13) + 1
          if(ipi(1,i).eq. 4) numes(14)  = numes(14) + 1
          if(ipi(1,i).eq. 5) numes(15)  = numes(15) + 1
          if(ipi(1,i).eq. 6) numes(16)  = numes(16) + 1
        end if
        if(ipi(1,i).gt.6 .or. ipi(1,i).lt.0 .or. abs(ipi(2,i)).gt.1 .or.
     &    ( (ipi(1,i).eq.2.or.ipi(1,i).eq.4.or.ipi(1,i).eq.5)
     &       .and.ipi(2,i).ne.0 ) )
     &  write(*,*) 'hiba: piondi ', ipi(1,i),ipi(2,i),ipi(5,i)
************************************************************************
*       meson dileptonic decay
        if(ivmesdil.eq.1 .and. (ipi(1,i).eq.3 .or. ipi(1,i).eq.5)) then
          pt  = sqrt(ppi(1,i)**2 + ppi(2,i)**2)
          ee  = sqrt(epi(i)**2 + ppi(1,i)**2 + ppi(2,i)**2 +ppi(3,i)**2)
          rap = 0.5 * log( (ee+ppi(3,i))/(ee-ppi(3,i)) )
          phi = atan2(ppi(2,i),ppi(1,i))
          gamma=ee/epi(i)
          ixx  = nint(rpi(1,i))
          iyy  = nint(rpi(2,i))
          izz  = nint(rpi(3,i))
          dens= 0.0
          if(iabs(ixx).le.maxx.and.iabs(iyy).le.maxx.and.iabs(izz)
     &        .le.maxz)
     &         dens = rhb(ixx,iyy,izz)
*
c         call vmesdil(epi(i),pt,rap,phi,dens,gamma,dt,ipi(1,i),ipi(5,i))
        end if
************************************************************************
*
  100 continue
      write(*,*) 'return from piondi'
      call f77flush()
      return
*
*----------------------------------------------------------------------*
*                                                                      *
      entry piondif(num,lmesc,numes,igamma,
     &               lmesa2,cres)
*                                                                      *
*----------------------------------------------------------------------*

      write(*,*)'in piondif'
      write(*,*)'vor pionem',num, lmesc,numes
      flag = .true.

       call pionem(30.,lmesc,flag,cres)
*
      write(*,*)'nach pionem in piondif '

      if(ivmesdil.eq.1 .or. imesdalitz.eq.1)
     &   call mes_dalitz(flag)
      if(ipiNbrems.eq.1) call piNdilepRun(50.0)

      flag = .true.
      if(idec2pi.ge.1) call mesdec(lropi,lsipi,lompi,flag,lmesa2)
      write(*,*)'nach mesdec'
      do 20 i = 1,18
        numes(i) = 0
  20  continue
*
      do 200 i = 1,maxppar
*         write(*,*)'maxppar = ', i
        if(ipi(1,i) .eq. 0)                                     goto 200
*
        if(ipi(1,i) .eq. 1) then
          if(ipi(2,i).eq. 1) numes(1) = numes(1) + 1
          if(ipi(2,i).eq. 0) numes(2) = numes(2) + 1
          if(ipi(2,i).eq.-1) numes(3) = numes(3) + 1
        end if
        if(ipi(1,i) .eq. 2)  numes(4) = numes(4) + 1
        if(ipi(1,i) .eq. 3) then
          if(ipi(2,i).eq. 1) numes(5) = numes(5) + 1
          if(ipi(2,i).eq. 0) numes(6) = numes(6) + 1
          if(ipi(2,i).eq.-1) numes(7) = numes(7) + 1
        end if
        if(ipi(1,i) .eq. 4)  numes(8) = numes(8) + 1
        if(ipi(1,i) .eq. 5)  numes(9) = numes(9) + 1
        if(ipi(1,i) .eq. 6)  numes(10) = numes(10) + 1
        if(rpi(1,i)**2+rpi(2,i)**2+rpi(3,i)**2 .lt. 4) then
          if(ipi(1,i).eq. 1) numes(11) = numes(11) + 1
          if(ipi(1,i).eq. 2) numes(12) = numes(12) + 1
          if(ipi(1,i).eq. 3) numes(13) = numes(13) + 1
          if(ipi(1,i).eq. 4) numes(14) = numes(14) + 1
          if(ipi(1,i).eq. 5) numes(15) = numes(15) + 1
          if(ipi(1,i).eq. 6) numes(16) = numes(16) + 1
        end if
*
        if(ipi(1,i).gt.6 .or. ipi(1,i).lt.0 .or. abs(ipi(2,i)).gt.1 .or.
     &    ( (ipi(1,i).eq.2.or.ipi(1,i).eq.4.or.ipi(1,i).eq.5)
     &       .and.ipi(2,i).ne.0 ) )
     &  write(*,*) 'hiba: piondif', ipi(1,i),ipi(2,i),ipi(5,i)
*
c        write(*,*) 'piondif old stat routine',i,maxpar
************************************************************************
*       meson collision chain
        ile = nres+1+ipi(1,i)
        idj = min(iabs(ipi(4,i)), 50)
        mlife(ile,idj) = mlife(ile,idj) + 1
*
************************************************************************
        ixx = nint(rpi(1,i))
        iyy = nint(rpi(2,i))
        izz = nint(rpi(3,i))
        if(abs(ixx).lt.maxx .and. abs(iyy).lt.maxx
     &                .and.abs(izz).lt.maxz) then
          gridte = .true.
        else
          gridte = .false.
        endif
        dendel = 0.0
        if(gridte) dendel = rpie(4,i)/rho0
        iendel = nint(5.0*dendel)
        if(iendel .gt. 50) iendel = 50
        il = ipi(1,i) + nres + 1
        mdens(il,iendel) = mdens(il,iendel) + 1
************************************************************************
        ibir = max(1, nint(rpie(6,i)))
        if (ibir .gt. 100) ibir = 100
        mbirt(il-4,ibir) = mbirt(il-4,ibir) + 1
        pos = sqrt(rpie(1,i)**2+rpie(2,i)**2+rpie(3,i)**2)
        ipos = nint(2.0*pos)
        ipos = min(ipos,20)
        ipos = max(ipos,1)
        mposi(il-nres -1,ipos) = mposi(il-nres-1,ipos) + 1
*
c       km = nint(40.0 * (rpie(5,i)-1.050) )
c        if(km.gt.0) then
c          jd = ipi(5,i) - 1
c          if(jd.ge.1) idelmas(jd+3,km)=idelmas(jd+3,km)+1
c        endif
c***********************************************************************
c        write(*,*) 'piondif: vor dilepton decay'
*       pion or eta electromagnetic decay
        if(ipi(1,i) .eq. 1 .and. imesdalitz.ge.1) then
          probpi0  = 0.0
          if((imesdalitz.eq.1) .and. (ipi(2,i).eq.0)) probpi0 = 1.0
c          if((ipi0d .eq. 1) .or. ((ipi0d.eq.2) .and. (ipi(2,i).eq.0)))
c     &       call pi0dec(ppi(1,i),ppi(2,i),ppi(3,i),0,probpi0)
c        elseif(ipi(1,i) .eq. 2 .and. ietad.ge.1) then
c          call etadec(ppi(1,i),ppi(2,i),ppi(3,i),1.0)
        end if

*
c***********************************************************************
c       write(*,*) 'piondif: vor gamma decay'
*       pion or eta gamma-decay
cc      if(igamma.eq.1 .and.
cc     &   ((ipi(1,i).eq.1 .and. ipi(2,i).eq.0) .or. ipi(1,i).eq.2)) then
cc       call pi0gam(ipi(1,i),epi(i),ppi(1,i),ppi(2,i),ppi(3,i),1.0,0)
cc        endif

  200 continue
      write(*,*)  '  end of piondif '
      return
*-----------------------------------------------------------------------
      entry pioini(jabso,jempi,jpico,jdec2pi,jmeson,jvmesdil,
     &    jpiNbrems,jmesdalitz)
*
*       this entry reads the range of the pion absorbtion cross section
*
*     this entry is used to initialize the common block ipi
*-----------------------------------------------------------------------
*
*
      imesdalitz= jmesdalitz
      ipiNbrems= jpiNbrems
      ivmesdil= jvmesdil
      ipico  = jpico
      iabso  = jabso
      iempi  = jempi
      idec2pi= jdec2pi
      imeson = jmeson
*
      do 30 i = 1, maxppar
        ipi(1,i) = 0
  30  continue
      return
*-----------------------------------------------------------------------
      entry inipion(mpion,num,x0,z0,px0,pz0)
      do 11 i = 1,num
        inp = (i-1) * maxp + 1
        epi(inp) = pmass
        mpot(inp) = 0.0
        rpi(1,inp) = x0
        rpi(2,inp) = 0.0
        rpi(3,inp) = z0
        ppi(1,inp) = px0
        ppi(2,inp) = 0.0
        ppi(3,inp) = pz0
        ipi(1,inp) = 1
        ipi(2,inp) = mpion
        ipi(3,inp) = 0
        ipi(4,inp) = 0
c-hw       rpie(5,inp)= rmass+pmass
        rpie(6,inp)= 1.0
        if(mpion.eq.2) ipi(2,inp) = 0
  11  continue
      return
      end
