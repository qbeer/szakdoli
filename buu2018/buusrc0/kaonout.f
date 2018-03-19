************************************************************************
*                                                                      *
      subroutine kaonout(wref,wmin,wmax,isu)
*                                                                      *
*         qq(i,j) /i=1-4/   -  i. coord. of the j. point in mom.space  *
*         qy(1,j)           -  the transverse mom belonging to j.      *
*         qy(2,j)           -  the rapidity       belonging to j.      *
*         qy(3,j)           -  angle in trans.dir.belonging to j.      *
*         iqq(0,j)          -  1 -> allowed by filter, 0-> not allowed;*
*         iqq(1,j)          -  number of n+n  with    pauli            *
*         iqq(2,j)          -  number of n+n  without pauli            *
*         iqq(3,j)          -  number of n+pi with    pauli            *
*         iqq(4,j)          -  number of n+pi without pauli            *
*         iqq(5,j)          -  number of resonance decay               *
*       variables:                                                     *
*         yref    - -1 * target rapidity in the frame                  *
*         nyk     - number of rapidities               (integer,input) *
*         nqtk    - number of transverse momenta       (integer,input) *
*         nfk     - number of angles                   (integer,input) *
*         isu     - loop number of runs                                *
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 thetal(0:9), wcm(0:9), phideg(0:9), dzeta(0:9), dangla(0:9)
      real*8 siginv(0:50,0:9), sigazim(0:20,0:9), pcm(0:50,0:9)
      real*8 v0(0:9), v1(0:9), v2(0:9), v3(0:9), v1pt(0:5), v0pt(0:5)
      real*8 sigpt(0:50,0:9), sigpt2(0:50,0:9)
      real*8 px_mean(0:9), multip(0:20)
      integer  i,j,ik, it, ith, iphik, nq2, nq3
      integer  ithetal, nykk, nkao, nphi, nrap, nplab
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0, ptra2, ptra
      real*8  xtrav, pzlab, thl, xmt, tlab, ph, rap, gew, sum, dangle
      real*8 dw, yta, ypr, dphi, plab, e0lab, thldeg
      real*8  prob, qq, dang1, dang2
      integer ntrav, kanums, kzero, icollks,  icollhs, nptra
c
      save  siginv, sigazim, sigpt, sigpt2, v0, v1, v2, v3, px_mean
      save  kanums, kzero, icollks, icollhs, multip
c-------------------------------------
      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*32
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*48
      thetal(5) = pi/180.*56
      thetal(6) = pi
      ithetal   = 5
      do  i = 1, ithetal
        dang1 = thetal(i) - thetal(i-1)
        dang2 = thetal(i+1) - thetal(i)
        dang1 = .5* min(dang1,dang2, dangle)
        dangla(i) = dang1
        dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
      if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
          dzeta(1) = 1.0 - cos(dangla(1))
      if (insys .eq. 0 )  then
       write(isum,*)' NO calculation in the lab system at the moment'
         return
      endif
      nykk = nyk
      if(nyk .gt. 6)    nykk = 6
      if(nyk .le. 2)    nykk = 3
      dw  = (wmax - wmin) / (nykk-2)
      wcm(0) = - dw
      do  i = 1,nykk
        wcm(i) = wcm(i-1) + dw
      enddo
      iphik  = 4
      dphi  = 180./(iphik)
      phideg(0) = .0
      do  i = 1,iphik
        phideg(i) = phideg(i-1) + dphi
      enddo
      dqtk = qtmaxk / nqtk
c
      IF (isu .eq. 1)   then
        kanums  = 0
        kzero   = 0
        icollks = 0
        icollhs = 0
        do  i = 0,20
          multip(i) = .0
        enddo
        do  nplab = 0, nqtk
          do  it    = 0, ithetal
            siginv(nplab,it) = .0
          enddo
        enddo
        do  nrap = 0, nykk
          v0     (nrap) = .0
          v1     (nrap) = .0
          v2     (nrap) = .0
          v3     (nrap) = .0
          px_mean(nrap) = .0
          do  nphi = 0, iphik
            sigazim(nphi, nrap) = .0
          enddo
          do  ntrav = 0, nqtk
            sigpt(ntrav, nrap) = .0
            sigpt2(ntrav, nrap) = .0
          enddo
        enddo
        do  nptra = 0, 5
          v0pt(nptra) = .0
          v1pt(nptra) = .0
        enddo
      ENDIF
c
      icollks = icollks + icollk
      icollhs = icollhs + icollh
      nkao = 0
      tmass = xkmas
      chlab = cosh(wref)
      shlab = sinh(wref)
      do  3000  ik = 1, max_kminu
          j =nx_hyp(0,ik)
          if (j .eq. 0) goto 3000
          prob = p_hyp(4,ik) / real(num*isubs)
          multip(3+j) = multip(3+j) + prob
 3000 continue
      do  2000  ik = 1, maxkaon
        if (ika(1,ik) .eq. 0) goto  2000
        if (ika(1,ik) .eq. 2) then
          kzero = kzero + 1
          multip(0) = multip(0) + pkao(4,ik) / real(num*isubs)
          goto  2000
        endif
c       if (ika(2,ik) .ne. 1) goto  2000
        kanums  = kanums  + 1
        nkao = nkao + 1
c     write(isum,*) '  kaons ',ik, kanums, ika(2,ik), pkao(4,ik)
        prob = pkao(4,ik) / real(num*isubs)
        i = ika(2,ik)
        if (i .eq. 1)  multip(1) = multip(1) + prob
        if (i .gt. 1 .and. i .le. 4)  multip(2) = multip(2) + prob
        if (i .gt. 4 .and. i .lt. 10)  multip(3) = multip(3) + prob
c
c      write(*,*)  '  kaonout check  ',isu, ik, ika(2,ik), prob
c
        px = pkao(1,ik)
        py = pkao(2,ik)
        pz = pkao(3,ik)
        ptra2 = px**2 + py**2
        ptra  = sqrt(ptra2)
        xtrav = sqrt(tmass**2+ptra2)
        p0    = sqrt(xtrav**2 + pz**2)
        pzlab =  p0 * shlab + pz * chlab
        e0lab  =  p0 * chlab + pz * shlab
        plab  =  sqrt(pzlab**2 + ptra2)
        tlab  =  sqrt(tmass**2 + pzlab**2 + ptra2) - tmass
        thl   = acos(pzlab / plab)
        thldeg = 180./pi * thl
        xmt   =  xtrav - tmass
        ph    =  px/ptra
        ph    =  min(1.,ph)
        ph    =  max(-1.,ph)
        ph    =  acos(ph)
        nphi =  nint(180./pi * ph/dphi)
        rap   =  pzlab / e0lab
        rap   =  0.5 * log((1.+rap)/(1.-rap))
        nrap  =   nint ((rap + 1.0*dw) / dw)
        nplab = nint(plab/dqtk)
        ntrav = nint(xmt/dqtk)
        nptra = nint(ptra/dqtk)
c----
c      if (nkao .lt. 10)  then
c      write(isum, *) ' kaon ', nkao, prob,px,py,pz
c      write(isum, *) nphi, plab, ph, rap, nrap
c      write(isum, *) nplab, pzlab, thldeg, dangle
c      endif
        ith = 0
        do  it = 1, ithetal
          if (abs(thl - thetal(it)) .lt. dangla(it))  ith = it
        enddo
        if (ith .ne. 0 .and. nplab .le. nqtk)
     1      siginv(nplab,ith) = siginv(nplab,ith) + prob
c--------
        if (nrap .le. nykk  .and. nrap .ge. 0)  then
          v0(nrap) = v0(nrap) + prob
          v1(nrap) = v1(nrap) + prob * cos(ph)
          v2(nrap) = v2(nrap) + prob * cos(2.*ph)
          v3(nrap) = v3(nrap) + prob * cos(3.*ph)
          px_mean(nrap) = px_mean(nrap) + px*prob
          if (nphi .le. iphik)
     1      sigazim(nphi, nrap) = sigazim(nphi, nrap) + prob
          if (ntrav .le. nqtk)
     1      sigpt(ntrav,nrap) = sigpt(ntrav,nrap) + prob
          if (nptra .le. nqtk)
     1      sigpt2(nptra,nrap) = sigpt2(nptra,nrap) + prob
        endif
        if (nptra .lt. 6 .and. nrap .lt. nykk/2 .and. nrap .ge. 0)
     1   then
	  v0pt(nptra) = v0pt(nptra) + prob
	  v1pt(nptra) = v1pt(nptra) + prob * cos(ph)
        endif
 2000 continue
      if (isu .lt. isubs)   return
c
      write(isum,852)   wmax-wmin, wref, yta, dqtk, dangleb
  852 format(/,' rapidity dw, y_tar, y_ref:',3f7.3,
     1                             ' dqtk, dangleb', 2f7.3)
      write(isum,850)  kanums, kzero, multip(0), icollks, icollhs
  850 format(/,'    kaon - production ',/, ' fictive kaons K+/K0:',2i9,
     1        '(',e11.4,')',' No of collisions: K+',i9,'  Hyp',i9)
      write(isum, 854) (multip(i), i=1,5)
  854 format(' multiplicities from  NN', e12.4, '  BB',e12.4,'  piB',
     1        e12.4, ' L/S ', 2e12.4)
c
      do  nplab = 0, nqtk
          pp = (nplab+0.0) * dqtk
          p0 = sqrt(tmass**2+pp*pp)
          if (nplab .eq. 0) pp = dqtk / sqrt(24.)
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1                          (pp*sin(thetal(it)))   **2  )
          siginv(nplab,it) = siginv(nplab,it) * sqrt(tmass**2+pp**2) /
     1                       (dzeta(it)*pp**2 * 2.*pi * dqtk)
        enddo
        do  nrap = 0,nykk
          sigpt(nplab,nrap) = sigpt(nplab,nrap) / (dw * qq * dqtk)
          sigpt2(nplab,nrap) = sigpt2(nplab,nrap) / (dw * pp * dqtk)
        enddo
      enddo
c
      do  nrap = 0, nykk
      do nphi = 0, iphik
      gew = 1.
      if (nphi .eq. 0)      gew = 2.
      if (nphi .eq. iphik)  gew = 2.
        sigazim(nphi, nrap) = sigazim(nphi, nrap) * gew
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
      do nphi = 0, iphik
        sum = sum + sigazim(nphi, nrap)
      enddo
      enddo
      do  nrap = 0, nykk
      do nphi = 0, iphik
        if (sum .gt. .0) then
           sigazim(nphi, nrap) = sigazim(nphi, nrap) / sum
        endif
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        sum = sum + v0(nrap)
        if (v0(nrap) .gt. .0) then
            px_mean(nrap) =     px_mean(nrap)/v0(nrap)
            v1     (nrap) = 2.* v1     (nrap)/v0(nrap)
            v2     (nrap) = 2.* v2     (nrap)/v0(nrap)
            v3     (nrap) = 2.* v3     (nrap)/v0(nrap)
        endif
      enddo
      do  nrap = 0, nykk
        if (sum .gt. .0) v0(nrap) = v0(nrap) / sum
      enddo
      do   nptra= 1,5
       if (v0pt(nptra) .gt. .0) v1pt(nptra) =v1pt(nptra)/v0pt(nptra)
      enddo
      write(isum,*) '  invariant kaon spectra '
      write(isum,*)
      nq2 = ithetal/2
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900 format(30x, ' theta_lab ',/,(10x, 9f22.1))
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=1,nq2)
  910 format(f10.3, 10(f10.3,e12.4))
      enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=nq3,ithetal)
      enddo
c....
      write(isum,*)
      write(isum,*) '  azimuthal  kaon distribution '
      write(isum,950)  (wcm(i), i = 0,nykk)
      write(isum,*)
  950 format(30x, ' rapidity_lap ',/,(10x, 9f10.3))
      do i = 0, iphik
      write(isum,960) phideg(i), (sigazim(i,j),j=0,nykk)
  960 format(f10.1,10f10.4)
      enddo
c....
      write(isum,*)
      write(isum,972) sum
  972 format(15x,' rapidity  distribution --','kaon number=',e12.4)
      write(isum,970)  (v0(j),j=0,nykk)
  970 format(10x,10f10.4)
      write(isum,*) '       p_transverse  distribution -- v1 - v2- v3 -'
      write(isum,970)  (px_mean(j),j=0,nykk)
      write(isum,970)  (v1     (j),j=0,nykk)
      write(isum,970)  (v2     (j),j=0,nykk)
      write(isum,970)  (v3     (j),j=0,nykk)
c..
      write(isum,*)  ' v1  as a function of p_t '
      write(isum,956)  (0.1*nptra, v1pt(nptra),nptra = 1,5)
  956 format(f12.3, f10.4)
      write(isum,*)
      write(isum,*) '                  transverse mass distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
  952 format( ' rapidity_lab',(11f12.3))
      do nplab = 0, nqtk
      pp  = nplab * dqtk
      write(isum,954) pp, (sigpt(nplab,nrap),nrap=0,nykk)
  954 format(f12.3,10e12.4)
      enddo
c----------
c     write(isum,*) '    transverse mass distribution  as fct of pt'
c     write(isum,*)
c     do nplab = 0, nqtk
c     pp  = nplab * dqtk
c     write(isum,954) pp, (sigpt2(nplab,nrap),nrap=0,nykk)
c     enddo
c----------
      end
************************************************************************
*                                                                      *
      subroutine kminusout(wref,wmin,wmax,isu)
*                                                                      *
*       variables:                                                     *
*         yref    - -1 * target rapidity in the frame                  *
*         nyk     - number of rapidities               (integer,input) *
*         nqtk    - number of transverse momenta       (integer,input) *
*         nfk     - number of angles                   (integer,input) *
*         isu     - loop number of runs                                *
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 thetal(0:9), wcm(0:9), phideg(0:9), dzeta(0:9), dangla(0:9)
      real*8 siginv(0:50,0:9), sigazim(0:20,0:9), pcm(0:50,0:9)
      real*8 v0(0:9), v1(0:9), v2(0:9), v3(0:9), v1pt(0:5), v0pt(0:5)
      real*8 sigpt(0:50,0:9), sigpt2(0:50,0:9)
      real*8 px_mean(0:9), multip(0:20)
      integer  i,j,ik, it, ith, iphik, nq2, nq3
      integer  ithetal, nykk, nkao, nphi, nrap, nplab
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0, ptra2, ptra
      real*8  xtrav, pzlab, thl, xmt, tlab, ph, rap, gew, sum, dangle
      real*8 dw, yta, ypr, dphi, plab, e0lab, thldeg
      real*8  prob, qq, dang1, dang2
      integer ntrav, kanums, icollks,  nptra, ikmi_abss
c
      save  siginv, sigazim, sigpt2, sigpt, v0, v1, v2, v3, px_mean
      save  kanums, icollks, multip, ikmi_abss
c-------------------------------------
      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*32
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*48
      thetal(5) = pi/180.*56
      thetal(6) = pi
      ithetal   = 5
      do  i = 1, ithetal
	  dang1 = thetal(i) - thetal(i-1)
	  dang2 = thetal(i+1) - thetal(i)
	  dang1 = .5* min(dang1,dang2, dangle)
	  dangla(i) = dang1
          dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
	  if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
          dzeta(1) = 1.0 - cos(dangla(1))
      if (insys .eq. 0 )  then
       write(isum,*)' NO calculation in the lab system at the moment'
         return
      endif
         nykk = nyk
         if(nyk .gt. 6)    nykk = 6
         if(nyk .lt. 2)    nykk = 3
         dw  = (wmax - wmin) / (nykk-2)
         wcm(0) =  - dw
         do  i = 1,nykk
         wcm(i) = wcm(i-1) + dw
         enddo
      iphik  = 4
      dphi  = 180./(iphik)
      phideg(0) = .0
      do  i = 1,iphik
      phideg(i) = phideg(i-1) + dphi
      enddo
      dqtk = qtmaxk / nqtk
c
      IF (isu .eq. 1)   then
      kanums  = 0
      ikmi_abss = 0
      icollks = 0
      do  i = 0,20
          multip(i) = .0
      enddo
      do  nplab = 0, nqtk
      do  it    = 0, ithetal
          siginv(nplab,it) = .0
      enddo
      enddo
      do  nrap = 0, nykk
          v0     (nrap) = .0
          v1     (nrap) = .0
          v2     (nrap) = .0
          v3     (nrap) = .0
          px_mean(nrap) = .0
      do  nphi = 0, iphik
          sigazim(nphi, nrap) = .0
      enddo
      do  ntrav = 0, nqtk
          sigpt(ntrav, nrap) = .0
          sigpt2(ntrav, nrap) = .0
      enddo
      enddo
      do  nptra = 0, 5
         v0pt(nptra) = .0
         v1pt(nptra) = .0
      enddo
      ENDIF
c
      icollks = icollks + icollk_mi
      ikmi_abss = ikmi_abss + ikmi_abs
      nkao = 0
      tmass = xkmas
      chlab = cosh(wref)
      shlab = sinh(wref)
      do  2000  ik = 1, max_kminu
      if (nx_kminu(0,ik) .ne. 1) goto  2000
      kanums  = kanums  + 1
      nkao = nkao + 1
c     write(isum,*) '  kminus ',ik, kanums, nx_kminu(1,ik),
c    1                           p_kminu(4,ik)
      prob = p_kminu(4,ik) / real(num*isubs)
      if (prob .lt. .0) then
          write(isum,*) ' kminusout negativ ', isu, ik, prob,
     1      nx_kminu(1,ik)
          stop
      endif
      i = nx_kminu(1,ik)
      if (i .eq. 1)                 multip(1) = multip(1) + prob
      if (i .gt. 1 .and. i .le. 4)  multip(2) = multip(2) + prob
      if (i .eq. 10)                multip(3) = multip(3) + prob
c     if (i .eq. 5 .or.  i .eq. 6)  multip(4) = multip(4) + prob
      if (i .eq. 5 )  multip(4) = multip(4) + prob
      if (i .eq. 6 )  multip(5) = multip(5) + prob
      if (i .eq. 7 )  multip(6) = multip(6) + prob
      px = p_kminu(1,ik)
      py = p_kminu(2,ik)
      pz = p_kminu(3,ik)
c
c     write(isum,*)  '  kaon-kminu-out check',isu, ik, nx_kminu(1,ik),
c    1              nx_kminu(2,ik),  px,py,pz,       prob
c
      ptra2 = px**2 + py**2
      ptra  = sqrt(ptra2)
      xtrav = sqrt(tmass**2+ptra2)
      p0    = sqrt(xtrav**2 + pz**2)
      pzlab =  p0 * shlab + pz * chlab
      e0lab  =  p0 * chlab + pz * shlab
      plab  =  sqrt(pzlab**2 + ptra2)
      tlab  =  sqrt(tmass**2 + pzlab**2 + ptra2) - tmass
      thl   = acos(pzlab / plab)
      thldeg = 180./pi * thl
      xmt   =  xtrav - tmass
      ph    =  px/ptra
      ph    =  min(1.,ph)
      ph    =  max(-1.,ph)
      ph    =  acos(ph)
      nphi =  nint(180./pi * ph/dphi)
      rap   =  pzlab / e0lab
      rap   =  0.5 * log((1.+rap)/(1.-rap))
      nrap  =   nint ((rap + 1.0*dw) / dw)
      nplab = nint(plab/dqtk)
      ntrav = nint(xmt/dqtk)
      nptra = nint(ptra/dqtk)
c----
c     write(isum, *) ' kaon ', nkao, prob,px,py,pz
c     write(isum, *) nphi, plab, ph, rap, nrap
c     write(isum, *) nplab, pzlab, thldeg, dangle
      ith = 0
      do  it = 1, ithetal
        if (abs(thl - thetal(it)) .lt. dangla(it))  ith = it
      enddo
      if (ith .ne. 0 .and. nplab .le. nqtk)
     1      siginv(nplab,ith) = siginv(nplab,ith) + prob
c--------
      if (nrap .le. nykk  .and. nrap .ge. 0)  then
          v0(nrap) = v0(nrap) + prob
          v1(nrap) = v1(nrap) + prob * cos(ph)
          v2(nrap) = v2(nrap) + prob * cos(2.*ph)
          v3(nrap) = v3(nrap) + prob * cos(3.*ph)
          px_mean(nrap) = px_mean(nrap) + px*prob
         if (nphi .le. iphik)
     1      sigazim(nphi, nrap) = sigazim(nphi, nrap) + prob
         if (ntrav .le. nqtk)
     1      sigpt(ntrav,nrap) = sigpt(ntrav,nrap) + prob
         if (nptra .le. nqtk)
     1      sigpt2(nptra,nrap) = sigpt2(nptra,nrap) + prob
      endif
      if (nptra .lt. 6 .and. nrap .lt. nykk/2 .and. nrap .ge. 0)
     1   then
	  v0pt(nptra) = v0pt(nptra) + prob
	  v1pt(nptra) = v1pt(nptra) + prob * cos(ph)
      endif
 2000 continue
      if (isu .lt. isubs)   return
c
      write(isum,852)   wmax-wmin, wref, dqtk, dangleb
  852 format(/,' rapidity difference:',2f7.3, ' dqtk, dangleb', 2f7.3)
      write(isum,850)  kanums, icollks, ikmi_abss
  850 format(//,'    kminus - production ',/, ' fictive kminus:', i9,
     1        ' No of collisions:',i9, ' No of absorptions:',i9)
      write(isum, 854) (multip(i), i=1,6)
  854 format(' multipl.:  NN', e11.3,
     1       '  BB',e11.3,' piB', e11.3,' piY',2e11.3,' NY',e11.3)
c
      do  nplab = 0, nqtk
          pp = (nplab+0.0) * dqtk
          p0 = sqrt(tmass**2+pp*pp)
          if (nplab .eq. 0) pp = dqtk / sqrt(24.)
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1                          (pp*sin(thetal(it)))   **2  )
          siginv(nplab,it) = siginv(nplab,it) * sqrt(tmass**2+pp**2) /
     1                       (dzeta(it)*pp**2 * 2.*pi * dqtk)
        enddo
        do  nrap = 0,nykk
          sigpt(nplab,nrap) = sigpt(nplab,nrap) / (dw * qq * dqtk)
          sigpt2(nplab,nrap) = sigpt2(nplab,nrap) / (dw * pp * dqtk)
        enddo
      enddo
c
      do  nrap = 0, nykk
      do nphi = 0, iphik
      gew = 1.
      if (nphi .eq. 0)      gew = 2.
      if (nphi .eq. iphik)  gew = 2.
        sigazim(nphi, nrap) = sigazim(nphi, nrap) * gew
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
      do nphi = 0, iphik
        sum = sum + sigazim(nphi, nrap)
      enddo
      enddo
      do  nrap = 0, nykk
      do nphi = 0, iphik
        if (sum .gt. .0) then
           sigazim(nphi, nrap) = sigazim(nphi, nrap) / sum
        endif
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        sum = sum + v0(nrap)
        if (v0(nrap) .gt. .0) then
            px_mean(nrap) =     px_mean(nrap)/v0(nrap)
            v1     (nrap) = 2.* v1     (nrap)/v0(nrap)
            v2     (nrap) = 2.* v2     (nrap)/v0(nrap)
            v3     (nrap) = 2.* v3     (nrap)/v0(nrap)
        endif
      enddo
      do  nrap = 0, nykk
        if (sum .gt. .0) v0(nrap) = v0(nrap) / sum
      enddo
      do   nptra= 1,5
       if (v0pt(nptra) .gt. .0) v1pt(nptra) =v1pt(nptra)/v0pt(nptra)
      enddo
      write(isum,*) '  invariant K- spectra '
      write(isum,*)
      nq2 = ithetal/2
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900 format(30x, ' theta_lab ',/,(10x, 9f22.1))
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=1,nq2)
  910 format(f10.3, 10(f10.3,e12.4))
      enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=nq3,ithetal)
      enddo
c....
      write(isum,*)
      write(isum,*) '  azimuthal  K- distribution '
      write(isum,950)  (wcm(i), i = 0,nykk)
      write(isum,*)
  950 format(30x, ' rapidity_lab  ',/,(10x, 9f10.3))
      do i = 0, iphik
      write(isum,960) phideg(i), (sigazim(i,j),j=0,nykk)
  960 format(f10.1,10f10.4)
      enddo
c....
      write(isum,*)
      write(isum,972) sum
  972 format(15x,' rapidity  distribution --','K- number=',e12.4)
      write(isum,970)  (v0(j),j=0,nykk)
  970 format(10x,10f10.4)
      write(isum,*) '       p_transverse  distribution -- v1 - v2- v3 -'
      write(isum,970)  (px_mean(j),j=0,nykk)
      write(isum,970)  (v1     (j),j=0,nykk)
      write(isum,970)  (v2     (j),j=0,nykk)
      write(isum,970)  (v3     (j),j=0,nykk)
c..
      write(isum,*)  ' v1  as a function of p_t '
      write(isum,956)  (0.1*nptra, v1pt(nptra),nptra = 1,5)
  956 format(f12.3, f10.4)
      write(isum,*)
      write(isum,*) '                  transverse mass distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
  952 format( ' rapidity_lab',(11f12.3))
      do nplab = 0, nqtk
      pp  = nplab * dqtk
      write(isum,954) pp, (sigpt(nplab,nrap),nrap=0,nykk)
  954 format(f12.3,10e12.4)
      enddo
c----------
c     write(isum,*) '    transverse mass distribution  as fct of pt'
c     write(isum,*)
c     do nplab = 0, nqtk
c     pp  = nplab * dqtk
c     write(isum,954) pp, (sigpt2(nplab,nrap),nrap=0,nykk)
c     enddo
c----------
      end
************************************************************************
************************************************************************
*                                                                      *
      subroutine phi_out(wref,wmin,wmax,isu)
*                                                                      *
*     variables:                                                       *
*     yref    - -1 * target rapidity in the frame                      *
*     nyk     - number of rapidities               (integer,input)     *
*     nqtk    - number of transverse momenta       (integer,input)     *
*     nfk     - number of angles                   (integer,input)     *
*     isu     - loop number of runs                                    *
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
      include 'com_pert'
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 thetal(0:9), wcm(0:22), phideg(0:9), dzeta(0:9),dangla(0:9)
      real*8 siginv(0:50,0:9), sigazim(0:20,0:22), pcm(0:50,0:9)
      real*8 siginv2(0:50,0:9), facsig, fac
      real*8 v0(0:22), v1(0:22), v2(0:22), v3(0:22)
      real*8 sigpt(0:50,0:22), sigpt_det(0:50,0:22)
      real*8 sigtramom(0:50,0:22), sigtramom2(0:50,0:22)
      real*8 sigtramom_det_Hel(0:50,0:22), sigtramom_det_CDC(0:50,0:22)
      real*8 px_mean(0:22)
      integer  i,j,k, ik, ikk, it, ith, iphik,
     1         nq2, nq3, nqtkk, ikm, supmax
      integer  ithetal, nykk, n_phi, nphi, nrap, nplab, ntramom
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0, ptra2, ptra
      real*8  xtrav, pzlab, thl, xmt, tlab, ph, rap, gew, sum, dangle
      real*8 dw, yta, ypr, dphi, plab, e0lab, thldeg, tramom, dtmk, pcms
      real*8  prob, effic_Hel, effic_CDC, qq, thc, thcl
      real*8 x_phi_piN, x_phi_piD, x_phi_bary, x_phi_nucl
      real*8  x_phi_rhoN, x_phi_rhoD, x_phi_pipi, x_phi_pirho
      real*8  x_phi_piN1440, x_phi_piN1520
      real*8  x_phi_cre, x_phi_kpkm, x_phi_surv
      real*8  x_phi_epair, x_phi_kpair, x_phi_kpkm_good
      real*8 x_phi_det_Hel, x_phi_det_CDC
      real*8  x_phi_old_det_Hel, x_phi_old_det_CDC, detprob, boxsize
      integer ntrav, kanums, icollks, ntmk
      real*8 densmin, densmax, ddens, dens, brkao, br_epair
      integer numdens, idens, nzeta, npcm, npcl
      real*8 zeta, dang1, dang2
      real*8  angular(0:9,0:5)
      parameter (brkao = 0.492)
      parameter (br_epair = .000291)
      real*8  pka, pkm, sphi, phi_mas(0:100)
      integer nsphi
      integer n_phi_piN, n_phi_piD, n_phi_bary, n_phi_nucl
      integer  n_phi_rhoN, n_phi_rhoD, n_phi_pipi, n_phi_pirho
      integer  n_phi_piN1440, n_phi_piN1520
      integer  n_phi_cre, n_phi_kpkm, n_phi_surv
      integer  n_phi_epair, n_phi_kpair, n_phi_kpkm_good
c
      save  siginv, siginv2, sigazim, v0, v1, v2, v3, px_mean
      save  sigpt, sigpt_det, sigtramom, sigtramom2, n_phi, icollks
      save  sigtramom_det_Hel, sigtramom_det_CDC
      save  x_phi_piN, x_phi_piD, x_phi_bary, x_phi_nucl
      save  x_phi_cre, x_phi_kpkm, x_phi_surv
      save  x_phi_epair, x_phi_kpair, x_phi_kpkm_good
      save  x_phi_rhoN, x_phi_rhoD, x_phi_piN1440, x_phi_piN1520
      save  x_phi_pipi, x_phi_pirho, x_phi_det_Hel, x_phi_det_CDC
      save  x_phi_old_det_Hel, x_phi_old_det_CDC
      save x_phi_dens, angular, phi_mas
      save  n_phi_piN, n_phi_piD, n_phi_bary, n_phi_nucl
      save  n_phi_cre, n_phi_kpkm, n_phi_surv
      save  n_phi_epair, n_phi_kpair, n_phi_kpkm_good
      save  n_phi_rhoN, n_phi_rhoD, n_phi_piN1440, n_phi_piN1520
      save  n_phi_pipi, n_phi_pirho
c
      parameter(densmin=0., densmax=3.5)
      parameter(numdens=14)
      real*8 x_phi_dens(1:numdens)
c
      character*4 reac, reac0
      reac = 'X X '
      reac0 = reac
c-------------------------------------
      if(massta.eq.58 .and. masspr.eq.58) reac='NiNi'
      if(massta.eq.96 .and. masspr.eq.96) reac='RuRu'

      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*32
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*48
      thetal(5) = pi/180.*56
      thetal(6) = pi/180.*90
      thetal(7) = pi
      ithetal   = 6
      do  i = 1, ithetal
	  dang1 = thetal(i) - thetal(i-1)
	  dang2 = thetal(i+1) - thetal(i)
	  dang1 = .5* min(dang1,dang2, dangle)
	  dangla(i) = dang1
          dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
          if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
          dzeta(1) = 1.0 - cos(dangla(1))
      if (insys .eq. 0 )  then
        write(isum,*)' NO calculation in the lab system at the moment'
        return
      endif
      nykk = nyk
c      if(nyk .gt. 22)    nykk = 22
      if(nyk .gt. 6)    nykk = 6
      if(nyk .lt. 3)    nykk = 3
c      write(*,*) 'nyk, nykk:', nyk, nykk
      dw  = (wmax - wmin) / (nykk-2)
      wcm(0) = - dw
      do  i = 1,nykk
        wcm(i) = wcm(i-1) + dw
      enddo
      iphik  = 4
      dphi  = 180./(iphik)
      phideg(0) = .0
      do  i = 1,iphik
        phideg(i) = phideg(i-1) + dphi
      enddo
      ntmk = nqtk
      nqtkk = 2*nqtk
      dtmk =  2. * qtmaxk / ntmk   !   for phi's
      dqtk = qtmaxk / nqtkk
      ddens = (densmax - densmin)/numdens
c
      IF (isu .eq. 1)   then
        n_phi = 0
        x_phi_surv = 0
        x_phi_cre = 0
        x_phi_kpkm = 0
        x_phi_epair = 0
        x_phi_kpair = 0
        x_phi_kpkm_good = 0
        x_phi_piN = 0
        x_phi_piD = 0
        x_phi_nucl = 0
        x_phi_bary = 0
        x_phi_rhoN = 0
        x_phi_rhoD = 0
        x_phi_pipi = 0
        x_phi_pirho = 0
        x_phi_piN1440 = 0
        x_phi_piN1520 = 0
        x_phi_det_Hel = 0
        x_phi_old_det_Hel = 0
        x_phi_det_CDC = 0
        x_phi_old_det_CDC = 0
        n_phi_surv = 0
        n_phi_cre = 0
        n_phi_kpkm = 0
        n_phi_piN = 0
        n_phi_piD = 0
        n_phi_nucl = 0
        n_phi_bary = 0
        n_phi_rhoN = 0
        n_phi_rhoD = 0
        n_phi_pipi = 0
        n_phi_pirho = 0
        n_phi_piN1440 = 0
        n_phi_piN1520 = 0
        do  nplab = 0, ntmk
          do  it    = 0, ithetal
            siginv(nplab,it) = .0
            siginv2(nplab,it) = .0
          enddo
        enddo
        do  nrap = 0, nykk
          v0     (nrap) = .0
          v1     (nrap) = .0
          v2     (nrap) = .0
          v3     (nrap) = .0
          px_mean(nrap) = .0
          do  nphi = 0, iphik
            sigazim(nphi, nrap) = .0
          enddo
          do  ntrav = 0, nqtkk
            sigpt(ntrav, nrap) = .0
            sigpt_det(ntrav, nrap) = .0
          enddo
          do ntramom = 0, ntmk
            sigtramom(ntramom,nrap) = 0.
            sigtramom2(ntramom,nrap) = 0.
            sigtramom_det_Hel(ntramom,nrap) = 0.
            sigtramom_det_CDC(ntramom,nrap) = 0.
          enddo
        enddo
        do i = 0,100
          phi_mas(i) = 0.
        end do
        do idens = 1,numdens
          x_phi_dens(idens) = 0.
        end do
        do  i = 0, 9
        do  j = 0, 5
            angular(i,j) = 0.0
        enddo
        enddo
      ENDIF
c
      tmass = xphimas
      chlab = cosh(wref)
      shlab = sinh(wref)
c---------------------------------------------- search
      supmax = maxkaon + max_pert
      prob = 0.0
      do  2000  ik = 1, supmax
        sphi = .0
        if (ik .gt. max_pert) then
            ikk = ik - max_pert
            if (ika(1,ikk) .eq. 1) goto  2500
            goto 2000
        endif
        if (nx_pert(id_phi,0,ik) .eq. 0) goto  2000
          prob = p_pert(id_phi,4,ik) / real(num*isubs)
          x_phi_cre = x_phi_cre + prob
          n_phi_cre = n_phi_cre + 1
        if (nx_pert(id_phi,5,ik) .eq. 0) goto 2000
          nsphi = nint(20.*xphimas)
          phi_mas(nsphi) = phi_mas(nsphi) + prob
          x_phi_surv = x_phi_surv + prob
          n_phi_surv = n_phi_surv + 1
          x_phi_epair = x_phi_epair + prob
          n_phi_epair = n_phi_epair + 1
          x_phi_kpair = x_phi_kpair + prob
          n_phi_kpair = n_phi_kpair + 1
        goto  2599
 2500 continue
      if (ika(2,ikk) .ne. 77) goto 2000
      if (ika(5,ikk) .eq. 0)  goto 2000
      ikm = ika(5,ikk)
c     write(isum, *) ' kminus  pair ',ik, ikm, p_pert(id_phi,4,ik)
c    1         ,nx_kminu(4, ikm)
      if (nx_kminu(4, ikm) .eq. ikk) then
          prob = pkao(4,ikk) / real(num*isubs) /  brkao
          x_phi_kpkm = x_phi_kpkm + prob
          n_phi_kpkm = n_phi_kpkm + 1
          px = p_kminu(1,ikm) + pkao(1,ikk)
          py = p_kminu(2,ikm) + pkao(2,ikk)
          pz = p_kminu(3,ikm) + pkao(3,ikk)
          pkm = sqrt(xkmas**2+p_kminu(1,ikm)**2 +
     1               p_kminu(2,ikm)**2 + p_kminu(3,ikm)**2)
          pka = sqrt(xkmas**2 + pkao(1,ikk)**2  +
     1               pkao(2,ikk)**2 + pkao(3,ikk)**2)
          sphi = sqrt((pkm+pka)**2 -px**2 - py**2 - pz**2)
          if (abs(sphi - xphimas) .lt. 0.05) then
            x_phi_kpkm_good = x_phi_kpkm_good + prob
            n_phi_kpkm_good = n_phi_kpkm_good + 1
            x_phi_kpair = x_phi_kpair + prob
            n_phi_kpair = n_phi_kpair + 1
          end if
          nsphi = nint (20.*sphi)
          phi_mas(nsphi) = phi_mas(nsphi) + prob
      endif
      if (i_epair .eq. 0)    goto 2000  !  only used if e+e- pairs are considered
      goto  2597
c
 2599 continue      !           here real phi's  are used
        if (nx_pert(id_phi,2,ik) .eq. 1) then
          x_phi_nucl = x_phi_nucl + prob
          n_phi_nucl = n_phi_nucl + 1
        end if
        if (nx_pert(id_phi,2,ik).ge.2.and.nx_pert(id_phi,2,ik).le.4)then
          x_phi_bary  = x_phi_bary + prob
          n_phi_bary  = n_phi_bary + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 5) then
          x_phi_piN = x_phi_piN + prob
          n_phi_piN = n_phi_piN + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 6) then
          x_phi_piD = x_phi_piD + prob
          n_phi_piD = n_phi_piD + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 8) then
          x_phi_rhoN = x_phi_rhoN + prob
          n_phi_rhoN = n_phi_rhoN + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 9) then
          x_phi_rhoD = x_phi_rhoD + prob
          n_phi_rhoD = n_phi_rhoD + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 10) then
          x_phi_pipi = x_phi_pipi + prob
          n_phi_pipi = n_phi_pipi + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 11) then
          x_phi_pirho = x_phi_pirho + prob
          n_phi_pirho = n_phi_pirho + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 12) then
          x_phi_piN1440 = x_phi_piN1440 + prob
          n_phi_piN1440 = n_phi_piN1440 + 1
        end if
        if (nx_pert(id_phi,2,ik) .eq. 13) then
          x_phi_piN1520 = x_phi_piN1520 + prob
          n_phi_piN1520 = n_phi_piN1520 + 1
        end if

        dens = pkao(8,ik)
        idens = nint(dens/ddens+0.5)
        if(idens.ge.1 .and. idens.le.numdens) then
          x_phi_dens(idens) = x_phi_dens(idens) + prob
        else
          write(*,*) 'Density in phi_out outside the interval.',
     &       dens, idens
        end if
c
c     write(*,*)  '  phi_out check  ',isu, ik, nx_pert(id_phi,2,ik), prob
c
        px = p_pert(id_phi,1,ik)
        py = p_pert(id_phi,2,ik)
        pz = p_pert(id_phi,3,ik)
 2597 continue      !   now we use K+K-  pairs and real phi-s
        n_phi = n_phi + 1
        ptra2 = px**2 + py**2
	pcms  = sqrt(ptra2+pz**2)
        ptra  = sqrt(ptra2)
        xtrav = sqrt(tmass**2+ptra2)
        p0    = sqrt(xtrav**2 + pz**2)
        pzlab =  p0 * shlab + pz * chlab
        e0lab  =  p0 * chlab + pz * shlab
        plab  =  sqrt(pzlab**2 + ptra2)
        tlab  =  sqrt(tmass**2 + pzlab**2 + ptra2) - tmass
        thc   = acos(pz / pcms)
        thl   = acos(pzlab / plab)
        thldeg = 180./pi * thl
        xmt   =  xtrav - tmass
        ph    =  px/ptra
        ph    =  min(1.,ph)
        ph    =  max(-1.,ph)
        ph    =  acos(ph)
        nphi =  nint(180./pi * ph/dphi)
        rap   =  pzlab / e0lab
        rap   =  0.5 * log((1.+rap)/(1.-rap))
        nrap  =  nint ((rap + 1.0*dw) / dw)
        nplab = nint(plab/dtmk)
	npcm  = nint(pcms/dtmk)
        ntrav = nint(xmt/dqtk)
        ntramom = nint(ptra/dtmk)
c     write(isum,*)  '  phi direct ', px,py,pz, pzlab, rap, shlab,
c    1                                ntrav,  nrap
        effic_Hel = .0
        effic_CDC = .0
        if (ik .le. max_pert) then
          if (reac .ne. reac0) then
            effic_Hel = detprob(ptra,rap,'Hel',reac) ! detector efficiency
            effic_CDC = detprob(ptra,rap,'CDC',reac) ! detector efficiency

c        write(*,*) 'acc - ptra,rap,effic_Hel',ptra,rap,effic_Hel
c        write(*,*) 'acc - ptra,rap,effic_CDC',ptra,rap,effic_CDC

            x_phi_det_Hel = x_phi_det_Hel + prob*effic_Hel
            x_phi_det_CDC = x_phi_det_CDC + prob*effic_CDC
            if(nx_pert(id_phi,2,ik).ge.1 .and.nx_pert(id_phi,2,ik).le.7)
     &               then
              x_phi_old_det_Hel = x_phi_old_det_Hel + prob*effic_Hel
              x_phi_old_det_CDC = x_phi_old_det_CDC + prob*effic_CDC
            end if
          endif                 !   reac0 check
        endif
c----
c     if (n_phi .lt. 10)  then
c     write(isum, *) ' phi_ ', n_phi, prob,px,py,pz
c     write(isum, *) nphi, plab, ph, rap, nrap
c     write(isum, *) nplab, pzlab, thldeg, dangle
c     endif
        ith = 0
        npcl = nplab
        do  it = 1, ithetal
        if (it .lt. 1      ) then
           thcl = thc
	   npcl = npcm
	else
           thcl = thl
           npcl = nplab
	endif
         if (abs(thcl - thetal(it)) .lt. dangla(it))  ith = it
        enddo
c       write(isum,3377) isu, ik, ith,prob, dangla(ith), thetal(ith)
c    1                   ,npcm, dtmk, px, py, pz
c3377   format(' phi_out ',3i6,e12.4,2f10.5,i5,4f10.5)
        if (ith .ne. 0 .and. npcl .le. nqtk) then
           siginv(npcl,ith) = siginv(npcl,ith) + prob
           siginv2(npcl,ith) = siginv2(npcl,ith) + prob**2
        endif
c--------
        if (nrap .le. nykk  .and. nrap .ge. 0)  then
          v0(nrap) = v0(nrap) + prob
          v1(nrap) = v1(nrap) + prob * cos(ph)
          v2(nrap) = v2(nrap) + prob * cos(2.*ph)
          v3(nrap) = v3(nrap) + prob * cos(3.*ph)
          px_mean(nrap) = px_mean(nrap) + px*prob
          if (nphi .le. iphik)
     1       sigazim(nphi, nrap) = sigazim(nphi, nrap) + prob
          if (ntrav .le. nqtkk) then
chw -------------------------------------------
chw  if (nx_pert(id_phi,2,ik) .ge. 2 .and. nx_pert(id_phi,2,ik) .le. 4)
chw     1      sigpt(ntrav,nrap) = sigpt(ntrav,nrap) + prob
            sigpt(ntrav,nrap) = sigpt(ntrav,nrap) + prob
c           sigpt_det(ntrav,nrap) = sigpt_det(ntrav,nrap)
c    &         + prob*effic_Hel
          endif
         if (ntramom .le. nqtkk) then
         sigtramom(ntramom,nrap) = sigtramom(ntramom,nrap) + prob
         sigtramom2(ntramom,nrap) = sigtramom2(ntramom,nrap) + prob**2
         sigtramom_det_Hel(ntramom,nrap) =
     &         sigtramom_det_Hel(ntramom,nrap) + prob*effic_Hel
            sigtramom_det_CDC(ntramom,nrap) =
     &         sigtramom_det_CDC(ntramom,nrap) + prob*effic_CDC
         endif
        endif
c   angular distribution
      if (ik .le. max_pert) then
      zeta = abs(pz) / sqrt(pz**2 + ptra2)
      nzeta = int(zeta / 0.2)
      angular(0,nzeta) = angular(0,nzeta)+prob
      if(nx_pert(id_phi,2,ik).le. 3) angular(1,nzeta) = angular(1,nzeta)
     &                                                 +prob
      if(nx_pert(id_phi,2,ik).eq. 5) angular(2,nzeta) = angular(2,nzeta)
     &                                                 +prob
      if(nx_pert(id_phi,2,ik).eq. 6) angular(3,nzeta) = angular(3,nzeta)
     &                                                 +prob
      if(nx_pert(id_phi,2,ik).eq. 8) angular(4,nzeta) = angular(4,nzeta)
     &                                                 +prob
      if(nx_pert(id_phi,2,ik).eq. 9) angular(5,nzeta) = angular(5,nzeta)
     &                                                 +prob
      if(nx_pert(id_phi,2,ik).eq.13) angular(6,nzeta) = angular(6,nzeta)
     &                                                 +prob
      endif
c--------------------------------------------
c                look  at K+
c--------------------------------------------
 2000 continue
c-------------- e+e- pairs: -----------------
      do ik = 1,max_epair
        if (nx_epair(0,ik).eq.1) then
          prob = p_epair(4,ik) / real(num*isubs) / br_epair
          x_phi_epair = x_phi_epair + prob
          n_phi_epair = n_phi_epair + 1
        end if
      end do

      if (isu .lt. isubs)   return

c---------- printout: -------------------
      write(isum,852)   wmax-wmin, yta, wref, dqtk, dangleb, tmass
 852  format(//' rapidity dw, y_tar, y_ref:',3f7.3, ' dqtk, dangleb',
     1           2f7.3, ' m_phi=',f7.2 )
      write(isum,'(///,''phi_ production:'',/)')
      write(isum,'(''fictive phi_s: '',i9)') n_phi
      write(isum,'(''phi numbers:'')')
      write(isum,'('' created:                  '',e14.4,i9)')
     &   x_phi_cre,n_phi_cre
      write(isum,'('' reconstructed from e+e-:  '',e14.4,i9)')
     &   x_phi_epair,n_phi_epair
      write(isum,'('' reconstructed from K+K-:  '',e14.4,i9)')
     &   x_phi_kpair,n_phi_kpair
      write(isum,'('' survived phi:             '',e14.4,i9)')
     &   x_phi_surv,n_phi_surv
      write(isum,'('' phi undestroyed K+K-:     '',e14.4,i9)')
     &   x_phi_kpkm,n_phi_kpkm
      write(isum,'('' ... with |m-mphi|<50 MeV  '',e14.4,i9)')
     &   x_phi_kpkm_good,n_phi_kpkm_good
      write(isum,*)
      write(isum,'(''number of reconstructed phi-s: '')')
      write(isum,'(''   in Helitron: '',e14.4)') x_phi_det_Hel
      write(isum,'(''   in CDC:      '',e14.4)') x_phi_det_CDC
      write(isum,*) 'phi-s produced by the various channels:'
      write(isum,'(''NN:          '',e14.4,i9)') x_phi_nucl,n_phi_nucl
      write(isum,'(''ND, DD, ...: '',e14.4,i9)') x_phi_bary,n_phi_bary
      write(isum,'(''piN:         '',e14.4,i9)') x_phi_piN,n_phi_piN
      write(isum,'(''piD:         '',e14.4,i9)') x_phi_piD,n_phi_piD
      write(isum,'(''rhoN:        '',e14.4,i9)') x_phi_rhoN,n_phi_rhoN
      write(isum,'(''rhoD:        '',e14.4,i9)') x_phi_rhoD,n_phi_rhoD
      write(isum,'(''pi+pi:       '',e14.4,i9)') x_phi_pipi,n_phi_pipi
      write(isum,'(''pi+rho:      '',e14.4,i9)') x_phi_pirho,n_phi_pirho
      write(isum,'(''N(1440)pi:   '',e14.4,i9)')
     &   x_phi_piN1440,n_phi_piN1440
      write(isum,'(''N(1520)pi:   '',e14.4,i9)')
     &   x_phi_piN1520,n_phi_piN1520
      write(isum,'(''number of reconstructed phi-s'',/,
     &   ''in the NN, ND, DD, ..., piN, piD channels: '')')
      write(isum,'(''   in Helitron: '',e14.4)') x_phi_old_det_Hel
      write(isum,'(''   in CDC:      '',e14.4)') x_phi_old_det_CDC
      write(isum,*)
c
      write(isum,'(''phi-mass  spectrum:'',/, (f6.2,e13.3))')
     &   (i*0.05, phi_mas(i),i=16,30)
      write(isum,*)
      do  nplab = 0, ntmk
        pp = (nplab+0.0) * dtmk
        p0 = sqrt(tmass**2+pp*pp)
        if (nplab .eq. 0) pp = dtmk / sqrt(24.)
        qq = pp + tmass
        if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1       (pp*sin(thetal(it)))   **2  )
          facsig = sqrt(tmass**2+pp**2)/(dzeta(it)*pp**2 * 2.*pi*dtmk)
          fac = siginv2(nplab,it) - siginv(nplab,it)**2/(isubs*num)
          fac = sqrt(abs(fac)) * facsig
          siginv(nplab,it) = siginv(nplab,it) * facsig
          if (siginv(nplab,it) .gt. 0) then
          siginv2(nplab,it) = 100. * fac / siginv(nplab,it)
	  else
	  siginv2(nplab,it) = .0
	  endif
        enddo
      enddo
      do nplab = 0, nqtkk
          pp = (nplab+0.0) * dqtk
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  nrap = 0,nykk
          sigpt(nplab,nrap) = sigpt(nplab,nrap) / (dw * qq * dqtk)
          sigpt_det(nplab,nrap) =
     &       sigpt_det(nplab,nrap) / (dw * qq * dqtk)
        enddo
      enddo
c
      do ntramom = 0, ntmk
        tramom = ntramom*dtmk
        boxsize = tramom * dtmk * dw
        if (ntramom.eq.0) boxsize = dtmk**2/8.*dw
        do nrap = 0, nykk
	pp = sigtramom2(ntramom,nrap)
     1         - sigtramom(ntramom,nrap)**2/(isubs*num)
        sigtramom(ntramom,nrap) = sigtramom(ntramom,nrap)/boxsize
	sigtramom2(ntramom,nrap)  = sqrt(abs(pp)) / boxsize
          sigtramom_det_Hel(ntramom,nrap) =
     &       sigtramom_det_Hel(ntramom,nrap)/boxsize
          sigtramom_det_CDC(ntramom,nrap) =
     &       sigtramom_det_CDC(ntramom,nrap)/boxsize
        end do
      end do
c
      do  nrap = 0, nykk
        do nphi = 0, iphik
          gew = 1.
          if (nphi .eq. 0)      gew = 2.
          if (nphi .eq. iphik)  gew = 2.
          sigazim(nphi, nrap) = sigazim(nphi, nrap) * gew
        enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        do nphi = 0, iphik
          sum = sum + sigazim(nphi, nrap)
        enddo
      enddo
      do  nrap = 0, nykk
        do nphi = 0, iphik
          if (sum .gt. .0) then
            sigazim(nphi, nrap) = sigazim(nphi, nrap) / sum
          endif
        enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        sum = sum + v0(nrap)
        if (v0(nrap) .gt. .0) then
            px_mean(nrap) =     px_mean(nrap)/v0(nrap)
            v1     (nrap) = 2.* v1     (nrap)/v0(nrap)
            v2     (nrap) = 2.* v2     (nrap)/v0(nrap)
            v3     (nrap) = 2.* v3     (nrap)/v0(nrap)
        endif
      enddo
      do  nrap = 0, nykk
        if (sum .gt. .0) v0(nrap) = v0(nrap) / sum
      enddo
      write(isum,*) '  invariant phi_ spectra '
      write(isum,*)
      nq2 = 3
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900  format(30x, ' theta_cm  ',/,(7x, 9f18.1))
       do  nplab = 0, ntmk
      write(isum,910) nplab*dtmk,
     1       (siginv(nplab,it),siginv2(nplab,it),it=1,nq2)
  910    format(f10.3,5(e12.4,'(',f6.1,')'))
       enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
       do  nplab = 0, ntmk
      write(isum,910) nplab*dtmk,
     1       (siginv(nplab,it),siginv2(nplab,it),it=nq3,ithetal)
c  911    format(f10.3, 10(f10.3,e12.4))
       enddo
c....
      write(isum,*)
      write(isum,*) '  azimuthal  phi_ distribution '
      write(isum,950)  (wcm(i), i = 0,nykk)
      write(isum,*)
 950  format(30x, ' rapidity_lab  ',/,(10x, 23f10.3))
      do i = 0, iphik
        write(isum,960) phideg(i), (sigazim(i,j),j=0,nykk)
 960    format(f10.1,23f10.4)
      enddo
c....
      write(isum,*)
      write(isum,972) sum
 972  format(15x,' rapidity  distribution --','phi_ number=',e12.4)
      write(isum,970)  (v0(j),j=0,nykk)
 970  format(10x,22f10.4)
      write(isum,*) '       p_transverse  distribution -- v1 - v2- v3 -'
      write(isum,970)  (px_mean(j),j=0,nykk)
      write(isum,970)  (v1     (j),j=0,nykk)
      write(isum,970)  (v2     (j),j=0,nykk)
      write(isum,970)  (v3     (j),j=0,nykk)
c..
      write(isum,*)
      write(isum,*)  '      angular distribution x 1/5   '
      write(isum,*)  ' zeta    total      B + B      pi + N    pi+Delta
     1   rho+N    rho+Delta    pi+1520 '
      do nzeta = 0, 4
       zeta = nzeta * 0.2 + 0.1
       write(isum, 942) zeta, (angular(k,nzeta), k=0,6)
  942  format(f6.1,7e11.3)
      enddo
c..
      write(isum,*)
      write(isum,*) '                transverse mass  distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
 952  format( ' rapidity_lab',(23f11.3))
      do nplab = 0, nqtkk
        pp  = nplab * dqtk
        write(isum,954) pp, (sigpt(nplab,nrap),nrap=0,nykk)
 954    format(f11.3,23e11.3)
      enddo
c..
c      write(isum,*)
c      write(isum,*) 'transverse mass distribution  x  det. effic.'
c      write(isum,952)  (wcm(i), i = 0,nykk)
c      write(isum,*)
c      do nplab = 0, nqtkk
c        pp  = nplab * dqtk
c        write(isum,954) pp, (sigpt_det(nplab,nrap),nrap=0,nykk)
c      enddo
c..
      write(isum,*)
      write(isum,*) '            transverse momentum distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
      do ntramom = 0, ntmk
        tramom  = ntramom * dtmk
        write(isum,954) tramom, (sigtramom(ntramom,nrap),nrap=0,nykk)
      enddo
      write(isum,*) '            transverse momentum errors '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
      do ntramom = 0, ntmk
        tramom  = ntramom * dtmk
        write(isum,954) tramom, (sigtramom2(ntramom,nrap),nrap=0,nykk)
      enddo
c..
      if (reac .ne. reac0)  then
      write(isum,*)
      write(isum,*) 'tra. mom. distribution  x  det. effic.(Helitron)'
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
      do ntramom = 0, ntmk
        tramom  = ntramom * dtmk
        write(isum,954) tramom,
     &     (sigtramom_det_Hel(ntramom,nrap),nrap=0,nykk)
      enddo
c..
      write(isum,*)
      write(isum,*) 'tra. mom. distribution  x  det. effic.(CDC)'
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
      do ntramom = 0, ntmk
        tramom  = ntramom * dtmk
        write(isum,954) tramom,
     &     (sigtramom_det_CDC(ntramom,nrap),nrap=0,nykk)
      enddo
      endif  !  reac   check
c..
      write(isum,*)
      write(isum,*) 'density dependence'
      write(isum,*) 'density/rho0  no. of phi-s'
      do idens = 1,numdens
        dens = (idens-0.5)*ddens
        write(*,*) 'idens, dens, x_phi_dens(idens)',
     &     idens, dens, x_phi_dens(idens)
        write(isum,'(f12.3,e12.4)') dens, x_phi_dens(idens)
      end do
c----------
      end

************************************************************************
*                                                                      *
      subroutine pihwout(wref,wmin,wmax,isu)
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 thetal(0:9), wcm(0:9), phideg(0:9), dzeta(0:9), dangla(0:9)
      real*8 siginvm(0:50,0:9), sigazim(0:20,0:9), pcm(0:50,0:9)
      real*8 siginvp(0:50,0:9)
      real*8 v0(0:9), v1(0:9), v2(0:9), v3(0:9)
      real*8 sigpt(0:50,0:9)
      real*8 px_mean(0:9)
      integer  i,j,k, ii, ik, it, ith, iphik, nq2, nq3, irun, inp
      integer  ithetal, nykk, npi, nphi, nrap, nplab
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0, ptra2, ptra
      real*8  xtrav, pzlab, thl, xmt, tlab, ph, rap, gew, sum, dangle
      real*8 dw, yta, ypr, dphi, plab, e0lab, thldeg
      real*8  prob, qq, dang1, dang2
      integer ntrav, npis, npisp, npis0, npism, neta, nomega, nrho0
c
      save  siginvp, siginvm
      save  sigazim, sigpt, v0, v1, v2, v3, px_mean, npis
      save  npisp, npis0, npism, neta, nomega, nrho0
c-------------------------------------
      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*30
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*50
      thetal(5) = pi/180.*60
      thetal(6) = pi
      ithetal   = 5
      do  i = 1, ithetal
	  dang1 = thetal(i) - thetal(i-1)
	  dang2 = thetal(i+1) - thetal(i)
	  dang1 = .5* min(dang1,dang2, dangle)
	  dangla(i) = dang1
          dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
	  if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
          dzeta(1) = 1.0 - cos(dangla(1))
      if (insys .eq. 0 )  then
       write(isum,*)' NO calculation in the lab system at the moment'
         return
      endif
         nykk = nyk
         if(nyk .gt. 6)    nykk = 6
         if(nyk .lt. 2)    nykk = 3
         dw  = (wmax - wmin) / (nykk-2)
         wcm(0) = - dw
         do  i = 1,nykk
         wcm(i) = wcm(i-1) + dw
         enddo
      iphik  = 4
      dphi  = 180./(iphik)
      phideg(0) = .0
      do  i = 1,iphik
      phideg(i) = phideg(i-1) + dphi
      enddo
      dqtk = qtmaxk / nqtk
c
      IF (isu .eq. 1)   then
c
      neta = 0
      nomega = 0
      nrho0 = 0
      npis = 0
      npism = 0
      npis0 = 0
      npisp = 0
      do  nplab = 0, nqtk
      do  it    = 0, ithetal
          siginvm(nplab,it) = .0
          siginvp(nplab,it) = .0
      enddo
      enddo
      do  nrap = 0, nykk
          v0     (nrap) = .0
          v1     (nrap) = .0
          v2     (nrap) = .0
          v3     (nrap) = .0
          px_mean(nrap) = .0
      do  nphi = 0, iphik
          sigazim(nphi, nrap) = .0
      enddo
      do  ntrav = 0, nqtk
          sigpt(ntrav, nrap) = .0
      enddo
      enddo
      ENDIF
c
c
      npi = 0
      tmass = pmass
      chlab = cosh(wref)
      shlab = sinh(wref)
*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        do 800 ii  = 1,maxp-1
          ik  = ii + inp
      if (ipi(1,ik) .eq. 0) goto  800
      if (ipi(1,ik) .eq. 5) nomega = nomega + 1
      if (ipi(1,ik) .eq. 3 .and. ipi(2,ik) .eq. 0
     1  .and. epi(ik) .gt. 0.6    ) nrho0 = nrho0 + 1
      if (ipi(1,ik) .eq. 2) neta = neta +1
      if (ipi(1,ik) .ne. 1) goto  800
      npi = npi + 1
      if (ipi(2,ik) .eq. 1) npisp = npisp + 1
      if (ipi(2,ik) .eq. 0) npis0 = npis0 + 1
      if (ipi(2,ik) .eq.-1) npism = npism + 1
      prob = 1. / real(num*isubs)
      px = ppi(1,ik)
      py = ppi(2,ik)
      pz = ppi(3,ik)
      ptra2 = px**2 + py**2
      ptra  = sqrt(ptra2)
      xtrav = sqrt(tmass**2+ptra2)
      p0    = sqrt(xtrav**2 + pz**2)
      pzlab =  p0 * shlab + pz * chlab
      e0lab  =  p0 * chlab + pz * shlab
      plab  =  sqrt(pzlab**2 + ptra2)
      tlab  =  sqrt(tmass**2 + pzlab**2 + ptra2) - tmass
      thl   = acos(pzlab / plab)
      thldeg = 180./pi * thl
      xmt   =  xtrav - tmass
      ph    =  px/ptra
      ph    =  min(1.,ph)
      ph    =  max(-1.,ph)
      ph    =  acos(ph)
      nphi =  nint(180./pi * ph/dphi)
      rap   =  pzlab / e0lab
      rap   =  0.5 * log((1.+rap)/(1.-rap))
      nrap  =   nint ((rap + 1.0*dw) / dw)
      nplab = nint(plab/dqtk)
      ntrav = nint(xmt/dqtk)
c----
c      if (npi .lt. 10)  then
c      write(isum, *) ' pion ', npi, pp,px,py,pz
c      write(isum, *) nphi, plab, ph, rap, nrap
c      write(isum, *) nplab, pzlab, thldeg, dangle
c      endif
      ith = 0
      do  it = 1, ithetal
        if (abs(thl - thetal(it)) .lt. dangla(it))  ith = it
      enddo
      if (ith .ne. 0 .and. nplab .le. nqtk) then
      if (ipi(2,ik) .eq. 1)
     1       siginvp(nplab,ith) = siginvp(nplab,ith) + prob
      if (ipi(2,ik) .eq. -1)
     1       siginvm(nplab,ith) = siginvm(nplab,ith) + prob
      endif
c--------
      if (nrap .le. nykk  .and. nrap .ge. 0)  then
          v0(nrap) = v0(nrap) + prob
          v1(nrap) = v1(nrap) + prob * cos(ph)
          v2(nrap) = v2(nrap) + prob * cos(2.*ph)
          v3(nrap) = v3(nrap) + prob * cos(3.*ph)
          px_mean(nrap) = px_mean(nrap) + px*prob
         if (nphi .le. iphik)
     1      sigazim(nphi, nrap) = sigazim(nphi, nrap) + prob
         if (ntrav .le. nqtk)
     1      sigpt(ntrav,nrap) = sigpt(ntrav,nrap) +  prob
      endif
  800 continue
 1000 continue
      npis = npis + npi
      if (isu .lt. isubs)   return
c
      write(isum,852)   wmax-wmin, wref, dqtk, dangleb
  852 format(/,' rapidity difference:',2f7.3, ' dqtk, dangleb', 2f7.3)
      write(isum,850) npisp, npis0, npism, neta, nomega, nrho0
  850 format(/,'    pion - production ','  number of pions +0- ',3i7
     1           , /,  ' etas',i7, ' omegas ',i7,'  rhos ',i7)
      do  nplab = 0, nqtk
          pp = (nplab+0.0) * dqtk
          p0 = sqrt(tmass**2+pp*pp)
          if (nplab .eq. 0) pp = dqtk / sqrt(24.)
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1                          (pp*sin(thetal(it)))   **2  )
          siginvp(nplab,it) = siginvp(nplab,it) * sqrt(tmass**2+pp**2) /
     1                       (dzeta(it)*pp**2 * 2.*pi * dqtk)
          siginvm(nplab,it) = siginvm(nplab,it) * sqrt(tmass**2+pp**2) /
     1                       (dzeta(it)*pp**2 * 2.*pi * dqtk)
        enddo
        do  nrap = 0,nykk
          sigpt(nplab,nrap) = sigpt(nplab,nrap) / (dw * qq * dqtk)
        enddo
      enddo
c
      do  nrap = 0, nykk
      do nphi = 0, iphik
      gew = 1.
      if (nphi .eq. 0)      gew = 2.
      if (nphi .eq. iphik)  gew = 2.
        sigazim(nphi, nrap) = sigazim(nphi, nrap) * gew
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
      do nphi = 0, iphik
        sum = sum + sigazim(nphi, nrap)
      enddo
      enddo
      do  nrap = 0, nykk
      do nphi = 0, iphik
        if (sum .gt. .0) then
           sigazim(nphi, nrap) = sigazim(nphi, nrap) / sum
        endif
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        sum = sum + v0(nrap)
        if (v0(nrap) .gt. .0) then
            px_mean(nrap) =     px_mean(nrap)/v0(nrap)
            v1     (nrap) = 2.* v1     (nrap)/v0(nrap)
            v2     (nrap) = 2.* v2     (nrap)/v0(nrap)
            v3     (nrap) = 2.* v3     (nrap)/v0(nrap)
        endif
      enddo
      do  nrap = 0, nykk
        if (sum .gt. .0) v0(nrap) = v0(nrap) / sum
      enddo
      write(isum,*) '  invariant pion spectra pi+ / pi-'
      write(isum,*)
      nq2 = ithetal/2
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900 format(30x, ' theta_lab ',/,(2x, 9f30.1))
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1 (pcm(nplab,it), siginvp(nplab,it),siginvm(nplab,it),it=1,nq2)
  910 format(f9.2, 10(f7.3,2e12.4))
      enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1 (pcm(nplab,it), siginvp(nplab,it),siginvm(nplab,it),
     1                                      it=nq3,ithetal)
      enddo
c....
      write(isum,*) '  azimuthal  pion distribution '
      write(isum,950)  (wcm(i), i = 0,nykk)
      write(isum,*)
  950 format(30x, ' rapidity_lab  ',/,(10x, 9f10.3))
      do i = 0, iphik
      write(isum,960) phideg(i), (sigazim(i,j),j=0,nykk)
  960 format(f10.1,10f10.4)
      enddo
c....
      write(isum,*)
      write(isum,972) sum
  972 format(15x,' rapidity  distribution --','pion number=',e12.4)
      write(isum,970)  (v0(j),j=0,nykk)
  970 format(10x,10f10.4)
      write(isum,*) '       p_transverse  distribution -- v1 - v2- v3 -'
      write(isum,970)  (px_mean(j),j=0,nykk)
      write(isum,970)  (v1     (j),j=0,nykk)
      write(isum,970)  (v2     (j),j=0,nykk)
      write(isum,970)  (v3     (j),j=0,nykk)
c..
      write(isum,*)
      write(isum,*) '                 transverse  m_t distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
  952 format( ' rapidity_lab',(11f12.3))
      do nplab = 0, nqtk
      pp  = nplab * dqtk
      write(isum,954) pp, (sigpt(nplab,nrap),nrap=0,nykk)
  954 format(f12.3,10e12.4)
      enddo
c----------
      end
************************************************************************
*                                                                      *
      subroutine prohwout(wref,wmin,wmax,isu)
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 thetal(0:9), wcm(0:9), phideg(0:9), dzeta(0:9), dangla(0:9)
      real*8 siginv(0:50,0:9), sigazim(0:20,0:9), pcm(0:50,0:9)
      real*8 v0(0:9), v1(0:9), v2(0:9), v3(0:9), v1pt(0:5), v0pt(0:5)
      real*8 sigpt(0:50,0:9)
      real*8 px_mean(0:9)
      integer  i,j,k, ii, ik, it, ith, iphik, nq2, nq3, irun, inp
      integer  ithetal, nykk, npr, nphi, nrap, nplab
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0, ptra2, ptra
      real*8  xtrav, pzlab, thl, xmt, tlab, ph, rap, gew, sum, dangle
      real*8 dw, yta, ypr, dphi, plab, e0lab, thldeg
      real*8  prob, qq, dang1, dang2
      integer ntrav, nprs, nptra
c
      save  siginv, sigazim, sigpt, v0, v1, v2, v3, px_mean, nprs
c-------------------------------------
      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*32
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*48
      thetal(5) = pi/180.*56
      thetal(6) = pi
      ithetal   = 5
      do  i = 1, ithetal
	  dang1 = thetal(i) - thetal(i-1)
	  dang2 = thetal(i+1) - thetal(i)
	  dang1 = .5* min(dang1,dang2, dangle)
	  dangla(i) = dang1
          dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
	  if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
          dzeta(1) = 1.0 - cos(dangla(1))
      if (insys .eq. 0 )  then
       write(isum,*)' NO calculation in the lab system at the moment'
         return
      endif
         nykk = nyk
         if(nyk .gt. 6)    nykk = 6
         if(nyk .lt. 2)    nykk = 3
         dw  = (wmax - wmin) / (nykk-2)
         wcm(0) = - dw
         do  i = 1,nykk
         wcm(i) = wcm(i-1) + dw
         enddo
      iphik  = 4
      dphi  = 180./(iphik)
      phideg(0) = .0
      do  i = 1,iphik
      phideg(i) = phideg(i-1) + dphi
      enddo
      dqtk = qtmaxk / nqtk
c
      IF (isu .eq. 1)   then
      nprs = 0
      do  nplab = 0, nqtk
      do  it    = 0, ithetal
          siginv(nplab,it) = .0
      enddo
      enddo
      do  nrap = 0, nykk
          v0     (nrap) = .0
          v1     (nrap) = .0
          v2     (nrap) = .0
          v3     (nrap) = .0
          px_mean(nrap) = .0
      do  nphi = 0, iphik
          sigazim(nphi, nrap) = .0
      enddo
      do  ntrav = 0, nqtk
          sigpt(ntrav, nrap) = .0
      enddo
      enddo
      do  nptra = 0, 5
        v0pt(nptra) = .0
	v1pt(nptra) = .0
      enddo
      ENDIF
c
c
      npr = 0
      tmass = rmass
      chlab = cosh(wref)
      shlab = sinh(wref)
*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxb
        do 800 ii  = 1,maxb
          ik  = ii + inp
      if (id(1,ik) .ne. 1) goto  800
      if (id(2,ik) .ne. 1) goto  800
      npr = npr + 1
      prob = 1. / real(num*isubs)
      px = p(1,ik)
      py = p(2,ik)
      pz = p(3,ik)
      ptra2 = px**2 + py**2
      ptra  = sqrt(ptra2)
      xtrav = sqrt(tmass**2+ptra2)
      p0    = sqrt(xtrav**2 + pz**2)
      pzlab =  p0 * shlab + pz * chlab
      e0lab  =  p0 * chlab + pz * shlab
      plab  =  sqrt(pzlab**2 + ptra2)
      tlab  =  sqrt(tmass**2 + pzlab**2 + ptra2) - tmass
      thl   = acos(pzlab / plab)
      thldeg = 180./pi * thl
      xmt   =  xtrav - tmass
      ph    =  px/ptra
      ph    =  min(1.,ph)
      ph    =  max(-1.,ph)
      ph    =  acos(ph)
      nphi =  nint(180./pi * ph/dphi)
      rap   =  pzlab / e0lab
      rap   =  0.5 * log((1.+rap)/(1.-rap))
      nrap  =   nint ((rap + 1.0*dw) / dw)
      nplab = nint(plab/dqtk)
      ntrav = nint(xmt/dqtk)
      nptra = nint(ptra*10.)
c----
c      if (npr .lt. 10)  then
c      write(isum, *) ' proton ', npr, prob,px,py,pz
c      write(isum, *) nphi, plab, ph, rap, nrap
c      write(isum, *) nplab, pzlab, thldeg, dangle
c      endif
      ith = 0
      do  it = 1, ithetal
        if (abs(thl - thetal(it)) .lt. dangla(it))  ith = it
      enddo
      if (ith .ne. 0 .and. nplab .le. nqtk)
     1      siginv(nplab,ith) = siginv(nplab,ith) + prob
c--------
      if (nrap .le. nykk  .and. nrap .ge. 0)  then
          v0(nrap) = v0(nrap) + prob
          v1(nrap) = v1(nrap) + prob * cos(ph)
          v2(nrap) = v2(nrap) + prob * cos(2.*ph)
          v3(nrap) = v3(nrap) + prob * cos(3.*ph)
          px_mean(nrap) = px_mean(nrap) + px*prob
         if (nphi .le. iphik)
     1      sigazim(nphi, nrap) = sigazim(nphi, nrap) + prob
         if (ntrav .le. nqtk)
     1      sigpt(ntrav,nrap) = sigpt(ntrav,nrap) + prob
      endif
      if (nptra .lt. 6 .and. nrap .lt. nykk/2 .and. nrap .ge. 0)
     1   then
           v0pt(nptra) = v0pt(nptra) + prob
           v1pt(nptra) = v1pt(nptra) + prob * cos(ph)
       endif
			
  800 continue
 1000 continue
      nprs = nprs + npr
      if (isu .lt. isubs)   return
c
      write(isum,852)   wmax-wmin, wref, dqtk, dangleb
  852 format(/,' rapidity difference:',2f7.3, ' dqtk, dangleb', 2f7.3)
      write(isum,850) nprs
  850 format(/,'    proton - spectra ','  number of protons ',i7)
      do  nplab = 0, nqtk
          pp = (nplab+0.0) * dqtk
          p0 = sqrt(tmass**2+pp*pp)
          if (nplab .eq. 0) pp = dqtk / sqrt(24.)
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1                          (pp*sin(thetal(it)))   **2  )
          siginv(nplab,it) = siginv(nplab,it) * sqrt(tmass**2+pp**2) /
     1                       (dzeta(it)*pp**2 * 2.*pi * dqtk)
        enddo
        do  nrap = 0,nykk
          sigpt(nplab,nrap) = sigpt(nplab,nrap) / (dw * qq * dqtk)
        enddo
      enddo
c
      do  nrap = 0, nykk
      do nphi = 0, iphik
      gew = 1.
      if (nphi .eq. 0)      gew = 2.
      if (nphi .eq. iphik)  gew = 2.
        sigazim(nphi, nrap) = sigazim(nphi, nrap) * gew
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
      do nphi = 0, iphik
        sum = sum + sigazim(nphi, nrap)
      enddo
      enddo
      do  nrap = 0, nykk
      do nphi = 0, iphik
        if (sum .gt. .0) then
           sigazim(nphi, nrap) = sigazim(nphi, nrap) / sum
        endif
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        sum = sum + v0(nrap)
        if (v0(nrap) .gt. .0) then
            px_mean(nrap) =     px_mean(nrap)/v0(nrap)
            v1     (nrap) = 2.* v1     (nrap)/v0(nrap)
            v2     (nrap) = 2.* v2     (nrap)/v0(nrap)
            v3     (nrap) = 2.* v3     (nrap)/v0(nrap)
        endif
      enddo
      do  nrap = 0, nykk
        if (sum .gt. .0) v0(nrap) = v0(nrap) / sum
      enddo
      do   nptra= 1,5
        if (v0pt(nptra) .gt. .0) v1pt(nptra) =v1pt(nptra)/v0pt(nptra)
      enddo
      write(isum,*) '  invariant proton spectra '
      write(isum,*)
      nq2 = ithetal/2
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900 format(30x, ' theta_lab ',/,(10x, 9f22.1))
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=1,nq2)
  910 format(f10.3, 10(f10.3,e12.4))
      enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=nq3,ithetal)
      enddo
c....
      write(isum,*) '  azimuthal  proton  distribution '
      write(isum,950)  (wcm(i), i = 0,nykk)
      write(isum,*)
  950 format(30x, ' rapidity_lab  ',/,(10x, 9f10.3))
      do i = 0, iphik
      write(isum,960) phideg(i), (sigazim(i,j),j=0,nykk)
  960 format(f10.1,10f10.4)
      enddo
c....
      write(isum,*)
      write(isum,972) sum
  972 format(15x,' rapidity  distribution --','proton number=',e12.4)
      write(isum,970)  (v0(j),j=0,nykk)
  970 format(10x,10f10.4)
      write(isum,*) '       p_transverse  distribution -- v1 - v2- v3 -'
      write(isum,970)  (px_mean(j),j=0,nykk)
      write(isum,970)  (v1     (j),j=0,nykk)
      write(isum,970)  (v2     (j),j=0,nykk)
      write(isum,970)  (v3     (j),j=0,nykk)
c..
      write(isum,*)  ' v1  as a function of p_t '
      write(isum,956)  (0.1*nptra, v1pt(nptra),nptra = 1,5)
  956 format(f12.3, f10.4)
      write(isum,*)
      write(isum,*) '                      transverse  distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
  952 format( ' rapidity_lab',(11f12.3))
      do nplab = 0, nqtk
      pp  = nplab * dqtk
      write(isum,954) pp, (sigpt(nplab,nrap),nrap=0,nykk)
  954 format(f12.3,10e12.4)
      enddo
c----------
      end
************************************************************************
*                                                                      *
      subroutine pro_kaon(wref,wmin,wmax,isu)
*
*     calculates the correlation of a forward nucleon from
*     a deuteron  with  kaons
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      common /correlkaon/ ikairun(0:maxkaon)
      integer ikairun
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 sigpt(0:1,0:50,0:9), kaoreal(0:1)
      real*8 siginv(0:1,0:50,0:9), thetal(0:9), dangla(0:9)
      real*8  dzeta(0:9), pcm(0:50,0:9)
      real*8 dangle, dang1, dang2, thl
      integer ith, nq2, nq3, ithetal, nplab
      integer  i, ii, ik, it, irun, np, nz0, inp
      integer  nykk, npr(0:1), nrap, ntrav, kanums(0:1)
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0
      real*8  ptra2, ptra, e0lab, pzlab, xmt, plab, rap
      real*8 dw, wcm(0:30), yta, ypr
      real*8  prob, qq, xtrav
c
      save  sigpt,  npr, kanums, kaoreal
c-------------------------------------
      dqtk = qtmaxk / nqtk
      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*32
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*48
      thetal(5) = pi/180.*56
      thetal(6) = pi
      ithetal   = 5
      do  i = 1, ithetal
          dang1 = thetal(i) - thetal(i-1)
          dang2 = thetal(i+1) - thetal(i)
          dang1 = .5* min(dang1,dang2, dangle)
          dangla(i) = dang1
          dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
      if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
      dzeta(1) = 1.0 - cos(dangla(1))
         nykk = nyk
         if(nyk .gt. 6)    nykk = 6
         if(nyk .lt. 2)    nykk = 3
         dw  = (wmax - wmin) / (nykk-2)
         wcm(0) = - dw
         do  i = 1,nykk
         wcm(i) = wcm(i-1) + dw
         enddo
c
      IF (isu .eq. 1)   then
      npr(0) = 0
      npr(1) = 0
      kanums(0) = 0
      kanums(1) = 0
      kaoreal(0) = .0
      kaoreal(1) = .0
      do  nz0 = 0, 1
      do  nrap = 0, nykk
      do  ntrav = 0, nqtk
          sigpt(nz0, ntrav, nrap) = .0
      enddo
      enddo
      do  nplab = 0, nqtk
      do  it    = 0, ithetal
          siginv(nz0,nplab,it) = .0
      enddo
      enddo
      enddo
      ENDIF
c
c
      tmass = rmass
      chlab = cosh(wref)
      shlab = sinh(wref)
*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxb
        do 800 ii  = 1,maxb
          ik  = ii + inp
      if (id(1,ik) .ne. 1) goto  800
      nz0 = id(2,ik)
c     write(isum, *) ' proton at ',irun, ii, ik,
c    1                maxp, npr
      prob = 1. / real(num*isubs)
      px = p(1,ik)
      py = p(2,ik)
      pz = p(3,ik)
      ptra2 = px**2 + py**2
      ptra  = sqrt(ptra2)
      xtrav = sqrt(tmass**2+ptra2)
      if (ptra .gt. .20 ) goto 800
      p0    = sqrt(xtrav**2 + pz**2)
      pzlab =  p0 * shlab + pz * chlab
      e0lab  =  p0 * chlab + pz * shlab
      if (e0lab .lt. 1.0 * elab) goto 800
c----
      npr(nz0) = npr(nz0) + 1
      tmass = xkmas
      do  600  ik = 1, maxkaon
      if (ika(1,ik) .ne. 1)     goto 600
      if (ikairun(ik) .ne. irun) goto 600
      prob = pkao(4,ik) / real(num*isubs)
      kanums(nz0)  = kanums(nz0)  + 1
      kaoreal(nz0) = kaoreal(nz0) + prob
c
c      write(*,*)  '  prot_kaon check  ',isu, ik, ika(2,ik), prob
c
      px = pkao(1,ik)
      py = pkao(2,ik)
      pz = pkao(3,ik)
      ptra2 = px**2 + py**2
      ptra  = sqrt(ptra2)
      xtrav = sqrt(tmass**2+ptra2)
      p0    = sqrt(xtrav**2 + pz**2)
      pzlab =  p0 * shlab + pz * chlab
      e0lab  =  p0 * chlab + pz * shlab
      plab  =  sqrt(pzlab**2 + ptra2)
      thl   = acos(pzlab / plab)
      xmt   =  xtrav - tmass
      rap   =  pzlab / e0lab
      rap   =  0.5 * log((1.+rap)/(1.-rap))
      nplab = nint(plab/dqtk)
      nrap  =   nint ((rap + 1.0*dw) / dw)
      ntrav = nint(xmt/dqtk)
c----
      if (nrap .le. nykk  .and. nrap .ge. 0)  then
         if (ntrav .le. nqtk)
     1      sigpt(nz0,ntrav,nrap) = sigpt(nz0,ntrav,nrap) + prob
      endif
      ith = .0
      do  it = 1, ithetal
        if (abs(thl - thetal(it)) .lt. dangla(it))  ith = it
      enddo
c     write(isum,*) '  in correla  ',ik, nplab,ith,prob
      if (ith .ne. 0 .and. nplab .le. nqtk)
     1      siginv(nz0,nplab,ith) = siginv(nz0,nplab,ith) + prob

c---------------------
  600 continue
			
  800 continue
 1000 continue
      if (isu .lt. isubs)   return
c
      write(isum,852)   wmax-wmin, wref, dqtk, dangleb
  852 format(/,' rapidity difference:',2f7.3, ' dqtk, dangleb', 2f7.3)
      write(isum,850) (npr(i), kanums(i),i=0,1)
  850 format(/,' correlation of kaons with forward nucleons ',/,
     1         '   neutrons:',i7,' fict. kaons:',i7,/,
     2         '    protons:',i7,' fict. kaons:',i7 )
      do  nz0 = 0,1
      do  nplab = 0, nqtk
          pp = (nplab+0.0) * dqtk
          p0 = sqrt(tmass**2+pp*pp)
          if (nplab .eq. 0) pp = dqtk / sqrt(24.)
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1                          (pp*sin(thetal(it)))   **2  )
          siginv(nz0, nplab,it) = siginv(nz0, nplab,it) *
     1           sqrt(tmass**2+pp**2) / (dzeta(it)*pp**2 * 2.*pi * dqtk)
        enddo
        do  nrap = 0,nykk
          sigpt(nz0,nplab,nrap) = sigpt(nz0,nplab,nrap)
     1                               / (dw * qq * dqtk)
        enddo
      enddo
      enddo
c
c..
      do  nz0  = 0,1
        write(isum,*)
        write(isum,950) nz0, kaoreal(nz0)
  950 format( '    transverse  kaon distribution for',
     1        ' charge of forward nucleon:', i2,
     2        ' multiplicity:',e12.4)
        write(isum,952)  (wcm(i), i = 0,nykk)
        write(isum,*)
  952 format( ' rapidity_lab',(11f12.3))
        do np = 0, nqtk
          pp  = np * dqtk
          write(isum,954) pp, (sigpt(nz0,np,nrap),nrap=0,nykk)
  954     format(f12.3,10e12.4)
        enddo
c----------------------    invariant cross section  -------
      write(isum,*) '  invariant kaon spectra '
      write(isum,*)
      nq2 = ithetal/2
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900 format(30x, ' theta_lab ',/,(10x, 9f22.1))
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nz0,nplab,it),it=1,nq2)
  910 format(f10.3, 10(f10.3,e12.4))
      enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nz0,nplab,it),it=nq3,ithetal)
      enddo
      enddo
c....

      end
c=====================================
************************************************************************
      real*8 function detprob(ptra,rapi,system,reac)
*     detection prob. of a phi                                         *
*----------------------------------------------------------------------*
      implicit none
      include 'common'

      real*8 ptra,rapi
      character*3 system
      character*4 reac

      character*80 fileHel, fileCDCnini, fileCDCruru, fileCDC
      parameter(fileHel = 'buuinput/phiacHel.dat')
      parameter(fileCDCnini = 'buuinput/phiacCDCnini.dat')
      parameter(fileCDCruru = 'buuinput/phiacCDCruru.dat')

      integer iHel, iCDC, Hel, CDC, sys
      parameter(iHel=61, iCDC=62)
      parameter(Hel=0,CDC=1)

      character*8 sys_in, reac_in
      real*8 rapimin(Hel:CDC), rapimax(Hel:CDC)
      real*8 ptmin(Hel:CDC), ptmax(Hel:CDC)
      real*8 dpt, drap, xpt, ptdif, xrap, rapidif
      integer nrap(Hel:CDC), npt(Hel:CDC), irap, ipt
      integer ipt1, ipt2, irap1, irap2

      real*8 accept(Hel:CDC,0:30,0:30)
      integer opening

      data opening /0/
      save opening, accept, rapimin, rapimax, ptmin, ptmax,
     &   nrap, npt

c      write(*,*) 'In detprob; ptra,rapi,system,reac: ',
c     &   ptra,rapi,system,reac

      if (opening.eq.0) then
        opening = 1
        if(reac .eq. 'NiNi') then
          fileCDC = fileCDCnini
        else if(reac .eq. 'RuRu') then
          fileCDC = fileCDCruru
        else
          fileCDC = fileCDCnini ! the code shouldn't stop
          write(isum,*) 'No phi acceptance for this reaction'
        end if
c
        if(system.ne.'Hel' .and. system.ne.'CDC')
     &     write(isum,*) 'invalid name for detector system'
c
        open(iHel, file=fileHel, status='old')
        read(iHel,*) sys_in
        read(iHel,*) reac_in
        if(sys_in .ne. 'Helitron') then
          write(isum,*) 'phi accept. file for Helitron is invalid!'
        end if
        read(iHel,*) rapimin(Hel), rapimax(Hel), nrap(Hel)
        read(iHel,*) ptmin(Hel), ptmax(Hel), npt(Hel)
        do ipt = 0,npt(Hel)
          read(iHel,*) (accept(Hel,ipt,irap), irap = 0,nrap(Hel))
        end do
        close(iHel)
c
        open(iCDC, file=fileCDC, status='old')
        read(iCDC,*) sys_in
        read(iCDC,*) reac_in
        if(sys_in.ne.'CDC' .or. reac_in.ne.reac) then
          write(isum,*) 'phi accept. file for CDC is invalid!'
        end if
        read(iCDC,*) rapimin(CDC), rapimax(CDC), nrap(CDC)
        read(iCDC,*) ptmin(CDC), ptmax(CDC), npt(CDC)
        do ipt = 0,npt(CDC)
          read(iCDC,*) (accept(CDC,ipt,irap), irap = 0,nrap(CDC))
        end do
        close(iCDC)
c
        write(*,*) 'After readin in detprob'
      endif

      if(system.eq.'Hel') then
        sys = Hel
      else if(system.eq.'CDC') then
        sys = CDC
      else
        sys = 0                 ! code shouldn't stop
      end if

      dpt = (ptmax(sys)-ptmin(sys))/npt(sys)
      drap = (rapimax(sys)-rapimin(sys))/nrap(sys)

      if (ptra.lt.ptmin(sys) .or. ptra.gt.ptmax(sys) .or.
     &   rapi.lt.rapimin(sys) .or. rapi.gt.rapimax(sys)) then
        detprob = 0.
c        write(*,*) 'detprob; outside table!!'
      else
        xpt = (ptra-ptmin(sys))/dpt
        xrap = (rapi-rapimin(sys))/drap
        ipt = nint(xpt+0.5)
        irap = nint(xrap+0.5)
c        write(*,*) 'detprob; ipt, irap, accept: ',
c     &     ipt, irap, accept(sys,ipt,irap)
        if (ipt.lt.1 .or. ipt.gt.npt(sys) .or. irap.lt.1
     &     .or. irap.gt.nrap(sys)) write(*,*) 'Hiba in detprob !!'
        ptdif = xpt - ipt + 1.
        rapidif = xrap - irap + 1.
        ipt1 = max(ipt-1,0)
        ipt2 = min(ipt,npt(sys))
        irap1 = max(irap-1,1)
        irap2 = min(irap,nrap(sys))
c        write(*,*) 'ipt1, ipt2, irap1, irap2:', ipt1, ipt2, irap1, irap2
        detprob = accept(sys,ipt1,irap1)*(1.-ptdif)*(1.-rapidif) +
     &     accept(sys,ipt2,irap1)*ptdif*(1.-rapidif) +
     &     accept(sys,ipt1,irap2)*(1.-ptdif)*rapidif +
     &     accept(sys,ipt2,irap2)*ptdif*rapidif
      endif

      return
      end
************************************************************************
*                                                                      *
      subroutine ksi_out(wref,wmin,wmax,isu)
*                                                                      *
*         qq(i,j) /i=1-4/   -  i. coord. of the j. point in mom.space  *
*         qy(1,j)           -  the transverse mom belonging to j.      *
*         qy(2,j)           -  the rapidity       belonging to j.      *
*         qy(3,j)           -  angle in trans.dir.belonging to j.      *
*         iqq(0,j)          -  1 -> allowed by filter, 0-> not allowed;*
*         iqq(1,j)          -  number of n+n  with    pauli            *
*         iqq(2,j)          -  number of n+n  without pauli            *
*         iqq(3,j)          -  number of n+pi with    pauli            *
*         iqq(4,j)          -  number of n+pi without pauli            *
*         iqq(5,j)          -  number of resonance decay               *
*       variables:                                                     *
*         yref    - -1 * target rapidity in the frame                  *
*         nyk     - number of rapidities               (integer,input) *
*         nqtk    - number of transverse momenta       (integer,input) *
*         nfk     - number of angles                   (integer,input) *
*         isu     - loop number of runs                                *
*                                                                      *
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_kminu'
*----------------------------------------------------------------------*
      real*8  wref,wmin,wmax
      integer  isu
      real*8 thetal(0:9), wcm(0:9), phideg(0:9), dzeta(0:9), dangla(0:9)
      real*8 siginv(0:50,0:9), sigazim(0:20,0:9), pcm(0:50,0:9)
      real*8 v0(0:9), v1(0:9), v2(0:9), v3(0:9)
      real*8 sigpt(0:50,0:9)
      real*8 px_mean(0:9), multip(0:20)
      integer  i,j,ik, it, ith, iphik, nq2, nq3
      integer  ithetal, nykk, nksi, nphi, nrap, nplab
      real*8  dqtk, chlab, shlab, tmass, pp, px, py, pz, p0, ptra2, ptra
      real*8  xtrav, pzlab, thl, xmt, tlab, ph, rap, gew, sum, dangle
      real*8 dw, yta, ypr, dphi, plab, e0lab, thldeg
      real*8  prob, qq, dang1, dang2, ksimas
      integer ntrav, kanums, icollks,  nptra
c
      save  siginv, sigazim, sigpt, v0, v1, v2, v3, px_mean
      save  kanums, icollks, multip
c-------------------------------------
      ksimas = 1.321
      yta = wmin
      ypr = wmax
      dangle = pi/180.*dangleb
      thetal(0) = -pi
      thetal(1) = .0
      thetal(2) = pi/180.*32
      thetal(3) = pi/180.*40
      thetal(4) = pi/180.*48
      thetal(5) = pi/180.*56
      thetal(6) = pi
      ithetal   = 5
      do  i = 1, ithetal
	  dang1 = thetal(i) - thetal(i-1)
	  dang2 = thetal(i+1) - thetal(i)
	  dang1 = .5* min(dang1,dang2, dangle)
	  dangla(i) = dang1
          dzeta(i) = cos(thetal(i)-dang1)- cos(thetal(i)+dang1)
      enddo
	  if (thetal(1) .eq. 0 .and. thetal(2).gt. dangle)
     1          dangla(1) = dangle
          dzeta(1) = 1.0 - cos(dangla(1))
      if (insys .eq. 0 )  then
       write(isum,*)' NO calculation in the lab system at the moment'
         return
      endif
         nykk = nyk
         if(nyk .gt. 6)    nykk = 6
         if(nyk .lt. 2)    nykk = 3
         dw  = (wmax - wmin) / (nykk-2)
         wcm(0) = - dw
         do  i = 1,nykk
         wcm(i) = wcm(i-1) + dw
         enddo
      iphik  = 4
      dphi  = 180./(iphik)
      phideg(0) = .0
      do  i = 1,iphik
      phideg(i) = phideg(i-1) + dphi
      enddo
      dqtk = qtmaxk / nqtk
c
      IF (isu .eq. 1)   then
      kanums  = 0
      icollks = 0
      do  i = 0,20
          multip(i) = .0
      enddo
      do  nplab = 0, nqtk
      do  it    = 0, ithetal
          siginv(nplab,it) = .0
      enddo
      enddo

      do  nrap = 0, nykk
          v0     (nrap) = .0
          v1     (nrap) = .0
          v2     (nrap) = .0
          v3     (nrap) = .0
          px_mean(nrap) = .0
      do  nphi = 0, iphik
          sigazim(nphi, nrap) = .0
      enddo
      do  ntrav = 0, nqtk
          sigpt(ntrav, nrap) = .0
      enddo
      enddo
      ENDIF
c
      nksi = 0
      tmass = ksimas
      chlab = cosh(wref)
      shlab = sinh(wref)
      do  2000  ik = 1, max_ksi
      if (nx_ksi(0,ik) .ne. 1) goto  2000
      kanums  = kanums  + 1
      nksi = nksi + 1
      prob = p_ksi(4,ik) / real(num*isubs)
      i = nx_ksi(1,ik)
      if (i .eq. 1)                 multip(1) = multip(1) + prob
      if (i .eq. 10)                multip(2) = multip(2) + prob
      if (i .eq. 11)                multip(3) = multip(3) + prob
      if (i .eq. 5 .or. i.eq.6)     multip(4) = multip(4) + prob

c
c      write(*,*)  '  ksi_out check  ',isu, ik, ika(2,ik), prob
c
      px = p_ksi(1,ik)
      py = p_ksi(2,ik)
      pz = p_ksi(3,ik)
      ptra2 = px**2 + py**2
      ptra  = sqrt(ptra2)
      xtrav = sqrt(tmass**2+ptra2)
      p0    = sqrt(xtrav**2 + pz**2)
      pzlab =  p0 * shlab + pz * chlab
      e0lab  =  p0 * chlab + pz * shlab
      plab  =  sqrt(pzlab**2 + ptra2)
      tlab  =  sqrt(tmass**2 + pzlab**2 + ptra2) - tmass
      thl   = acos(pzlab / plab)
      thldeg = 180./pi * thl
      xmt   =  xtrav - tmass
      ph    =  px/ptra
      ph    =  min(1.,ph)
      ph    =  max(-1.,ph)
      ph    =  acos(ph)
      nphi =  nint(180./pi * ph/dphi)
      rap   =  pzlab / e0lab
      rap   =  0.5 * log((1.+rap)/(1.-rap))
      nrap  =   nint ((rap + 1.0*dw) / dw)
      nplab = nint(plab/dqtk)
      ntrav = nint(xmt/dqtk)
      nptra = nint(ptra*10.)
c     write(isum,*) ' ksi no',kanums, thl, ptra, pzlab, rap,
c    1                nplab, nrap, ntrav, nx_ksi(1,ik), prob
c----
      ith = .0
      do  it = 1, ithetal
        if (abs(thl - thetal(it)) .lt. dangla(it))  ith = it
      enddo
      if (ith .ne. 0 .and. nplab .le. nqtk)
     1      siginv(nplab,ith) = siginv(nplab,ith) + prob
c--------
      if (nrap .le. nykk  .and. nrap .ge. 0)  then
          v0(nrap) = v0(nrap) + prob
          v1(nrap) = v1(nrap) + prob * cos(ph)
          v2(nrap) = v2(nrap) + prob * cos(2.*ph)
          v3(nrap) = v3(nrap) + prob * cos(3.*ph)
          px_mean(nrap) = px_mean(nrap) + px*prob
         if (nphi .le. iphik)
     1      sigazim(nphi, nrap) = sigazim(nphi, nrap) + prob
         if (ntrav .le. nqtk)
     1      sigpt(ntrav,nrap) = sigpt(ntrav,nrap) + prob
      endif
 2000 continue
      if (isu .lt. isubs)   return
c
      write(isum,852)   wmax-wmin, wref, dqtk, dangleb
  852 format(/,' rapidity difference:',2f7.3,' dqtk, dangleb',2f7.3)
      write(isum,850)  kanums
  850 format(/,'    ksi - production ',/, ' fictive ksi_s:', i9)
      write(isum, 854) (multip(i), i=1,4)
  854 format(' multiplicities from  K-N', e12.4,
     1       ' piL',e12.4,'  piS', e12.4, ' K-Y',e12.4)

c
      do  nplab = 0, nqtk
          pp = (nplab+0.0) * dqtk
          p0 = sqrt(tmass**2+pp*pp)
          if (nplab .eq. 0) pp = dqtk / sqrt(24.)
          qq = pp + tmass
          if (nplab .eq. 0) qq = qq / 2.0
        do  it = 1,ithetal
          pz = pp * cos(thetal(it))
          pcm(nplab,it) = sqrt( (-p0*shlab + pz*chlab) **2 +
     1                          (pp*sin(thetal(it)))   **2  )
          siginv(nplab,it) = siginv(nplab,it) * sqrt(tmass**2+pp**2) /
     1                       (dzeta(it)*pp**2 * 2.*pi * dqtk)
        enddo
        do  nrap = 0,nykk
          sigpt(nplab,nrap) = sigpt(nplab,nrap) / (dw * qq * dqtk)
        enddo
      enddo
c
      do  nrap = 0, nykk
      do nphi = 0, iphik
      gew = 1.
      if (nphi .eq. 0)      gew = 2.
      if (nphi .eq. iphik)  gew = 2.
        sigazim(nphi, nrap) = sigazim(nphi, nrap) * gew
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
      do nphi = 0, iphik
        sum = sum + sigazim(nphi, nrap)
      enddo
      enddo
      do  nrap = 0, nykk
      do nphi = 0, iphik
        if (sum .gt. .0) then
           sigazim(nphi, nrap) = sigazim(nphi, nrap) / sum
        endif
      enddo
      enddo
c....
      sum = .0
      do  nrap = 0, nykk
        sum = sum + v0(nrap)
        if (v0(nrap) .gt. .0) then
            px_mean(nrap) =     px_mean(nrap)/v0(nrap)
            v1     (nrap) = 2.* v1     (nrap)/v0(nrap)
            v2     (nrap) = 2.* v2     (nrap)/v0(nrap)
            v3     (nrap) = 2.* v3     (nrap)/v0(nrap)
        endif
      enddo
      do  nrap = 0, nykk
        if (sum .gt. .0) v0(nrap) = v0(nrap) / sum
      enddo
      write(isum,*) '  invariant ksi_ spectra '
      write(isum,*)
      nq2 = ithetal/2
      nq3 = nq2+1
      write(isum,900)  (180./pi*thetal(i), i = 1,nq2)
  900 format(30x, ' theta_lab ',/,(10x, 9f22.1))
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=1,nq2)
  910 format(f10.3, 10(f10.3,e12.4))
      enddo
      write(isum,900)  (180./pi*thetal(i), i =nq3,ithetal )
      do  nplab = 0, nqtk
      write(isum,910) nplab*dqtk,
     1       (pcm(nplab,it), siginv(nplab,it),it=nq3,ithetal)
      enddo
c....
      write(isum,*)
      write(isum,*) '  azimuthal  ksi_ distribution '
      write(isum,950)  (wcm(i), i = 0,nykk)
      write(isum,*)
  950 format(30x, ' rapidity_lab  ',/,(10x, 9f10.3))
      do i = 0, iphik
      write(isum,960) phideg(i), (sigazim(i,j),j=0,nykk)
  960 format(f10.1,10f10.4)
      enddo
c....
      write(isum,*)
      write(isum,972) sum
  972 format(15x,' rapidity  distribution --','ksi_ number=',e12.4)
      write(isum,970)  (v0(j),j=0,nykk)
  970 format(10x,10f10.4)
      write(isum,*) '       p_transverse  distribution -- v1 - v2- v3 -'
      write(isum,970)  (px_mean(j),j=0,nykk)
      write(isum,970)  (v1     (j),j=0,nykk)
      write(isum,970)  (v2     (j),j=0,nykk)
      write(isum,970)  (v3     (j),j=0,nykk)
c..
c  956 format(f12.3, f10.4)
      write(isum,*)
      write(isum,*) '                  transverse mass distribution '
      write(isum,952)  (wcm(i), i = 0,nykk)
      write(isum,*)
  952 format( ' rapidity_lab',(11f12.3))
      do nplab = 0, nqtk
      pp  = nplab * dqtk
      write(isum,954) pp, (sigpt(nplab,nrap),nrap=0,nykk)
  954 format(f12.3,10e12.4)
      enddo
c----------
      end
************************************************************************
      function bin_num(mdilep,ninterval,minmas,maxmas,dinv)
************************************************************************
      implicit none
*----------------------------------------------------------------------*
      real*8  mdilep,minmas,maxmas,dinv
      integer  bin_num,ninterval
      dinv = (maxmas-minmas) / ninterval
      bin_num = int((mdilep-minmas)/dinv)+1
      return
      end
************************************************************************
      function dilmass(n_bin,ninterval,minmas,maxmas,dinv)
************************************************************************
      implicit none
*----------------------------------------------------------------------*
      real*8  dilmass,minmas,maxmas,dinv
      integer  n_bin,ninterval
      dinv = (maxmas-minmas) / ninterval
      dilmass = minmas + (n_bin-0.5)*dinv
      return
      end

