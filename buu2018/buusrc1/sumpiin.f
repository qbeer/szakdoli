
************************************************************************
*                                                                      *
      subroutine sumpiin(ipispe,bik,bit,bim,jfram,mang,anpi,
     &                   cang,wref,ydel)
*                                                                      *
*         variables:                                                   *
*            ipispe   - 0-> isospin average for spectra; 1-> only pi-; *
*            bink     - kinetic energy bin  (gev)                      *
*            binr     - rapidity bin                                   *
*            bint     - bin in trans.mom/pmass                         *
*            binm     - bin in momentum  (gev)                         *
*            iframe   - 0 -> spectra at angles of cms; 1 -> of lab.;   *
*            nangle   - number of angles for spectra (max.:5);         *
*            angle    - angles for spectra in degrees (max.:5);        *
*            dangle   - bin for angles;                                *
*            yref     - rapidity of cms in the lab.;                   *
*                                                                      *
************************************************************************
      implicit none
      real*8 bik,bit,bim,cang,wref,ydel,bink,bint,binm,dangle,yref
      real*8 binr,rat,px,py,pz,etot,ekin,pio,rap,ptr,ylab,teta,tim
      real*8 scalpro,ymax,pitrmi,pitrma,fect,val,e1,elab,facta,b0,e0
      real*8 tmin,tmax,b1
      integer ipispe,jfram,mang,maxbin,maxpr,maxpt,isosp,iframe,nangle
      integer j,i,nr,ii,jj,npipl,npize,npimi,npi,nbink,nbinr,nbint,nbinm
      integer num,isubs,ikin,ifit,imin,imax,ja,k,nt
      parameter     (maxbin =   50)
      parameter     (maxpr  =    5,  maxpt =    20)
      include"common"
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      save pitr,rappi,rapptpi,specpi
      save isosp,bink,binr,bint,binm,iframe,nangle,dangle,yref
      save angle
      save dsdp,ptsp
*----------------------------------------------------------------------*
      real*8         pitr(11,maxbin),epik(maxbin)
      real*8         rapptpi(10,-maxpr:maxpr,0:maxpt)
      real*8         rappi(10,-maxpr:maxpr)
      real*8         specpi(10,5,maxbin)
      real*8         dsdp(11,maxbin)
      real*8         ptsp(11,maxbin),pipt(maxbin)
      real*8         angle(5),anpi(5)
      real*8         timpisp(10)
*----------------------------------------------------------------------*
      isosp  = ipispe
      bink   = bik
      bint   = bit
      binm   = bim
      iframe = jfram
      nangle = mang
      dangle = pi * cang / 180.0
      yref   = wref
      binr   = ydel/8.0
      do 100  j = 1,nangle
        angle(j) = pi * anpi(j) / 180.0
 100  continue
      do 600  j = 1,10
        timpisp(j) = 0.0
        do 300  i = 1,maxbin
          pitr(j,i) = 0.0
          ptsp(j,i) = 0.0
          dsdp(j,i) = 0.0
          do 200  k = 1,5
            specpi(j,k,i) = 0.0
 200      continue
 300    continue
        do 500  nr = - maxpr , maxpr
          rappi(j,nr) = 0.0
          do 400  nt = 0 , maxpt
            rapptpi(j,nr,nt) = 0
 400      continue
 500    continue
 600  continue
      do 800  ii = 1,50
        etotpi(ii) = 0.0
        do 700  jj = 1 , 3
          ptotpi(jj,ii) = 0.0
 700    continue
 800  continue
      return
************************************************************************
*                                                                      *
      entry sumpion(tim ,npipl,npize,npimi)
*                                                                      *
************************************************************************
      timpisp(itip) = tim
      npi  = npipl + npize + npimi
      if(npi.le.0) return
      rat  = float(npimi) / float(npi)
      if(isosp .eq. 1) rat = 1.0
      do 1200 i = 1,maxppar
        if(ipi(1,i) .ne. 1)                                    goto 1200
        px = ppi(1,i)
        py = ppi(2,i)
        pz = ppi(3,i)
        etot = sqrt(pmass**2+px**2+py**2+pz**2)
        if((ipi(2,i).ne.-1) .and. (isosp.eq.1))                goto 1200
        ekin = etot - pmass
        pio  = sqrt(px**2+py**2+pz**2)
        rap  = 0.5 * log( (etot+pz)/(etot-pz) )
        ptr  = sqrt(px**2+py**2)
        nbink= nint( ekin / bink ) + 1
        nbinr= nint( rap/binr )
        nbint= nint(  ptr/pmass/bint )
        ylab = rap + yref
c       if(ptr.le.1.e-8) write(isum,'(''c:ptr:'',i5,e12.3)') i,ptr
c       if((nbint.lt.maxbin) .and. (ylab.ge.0.5) .and. (ylab.le.0.9))
        if(nbint.lt.maxbin)
     &      ptsp(itip,nbint+1) =ptsp(itip,nbint+1) + rat/ptr
        if(nbink.le.maxbin)
     &      dsdp(itip,nbink)=dsdp(itip,nbink)+rat/pio/ekin
        if((nbink.le.maxbin) .and. (ptr/pio.ge.cos(dangle)))
     &    pitr(itip,nbink)=pitr(itip,nbink) + rat/pio
        if(abs(nbinr).le.maxpr) rappi(itip,nbinr)=rappi(itip,nbinr)+rat
        if((nbint.le.maxpt) .and. (abs(nbinr).le.maxpr) )
     &    rapptpi(itip,nbinr,nbint)=rapptpi(itip,nbinr,nbint)+rat/ptr
        if(nangle .eq. 0)                                      goto 1200
        if(iframe .gt. 0) then
          pz = sqrt(pmass**2+px**2+py**2) * sinh(ylab)
          etot = sqrt(pmass**2+px**2+py**2+pz**2)
          pio  = sqrt(px**2+py**2+pz**2)
        end if
        nbinm= nint(pio/binm) + 1
        if(nbinm .gt. maxbin)                                  goto 1200
        teta = asin(ptr/pio)
        do 1100 j =1,nangle
           if(abs(teta-angle(j)) .le. dangle)
     &       specpi(itip,j,nbinm)=specpi(itip,j,nbinm) + rat*etot/pio**2
 1100   continue
 1200 continue
      pitr(itip,1)=pitr(itip,1) * 2.0
      dsdp(itip,1)=dsdp(itip,1) * 2.0
      ptsp(itip,1)=ptsp(itip,1) * 2.0
      do 1300 nr = - maxpr , maxpr
        rapptpi(itip,nr,0) = rapptpi(itip,nr,0) * 2.0
 1300 continue
      do 1400 j = 1, nangle
        specpi(itip,j,1)=specpi(itip,j,1) * 2.0
 1400 continue
      return
************************************************************************
*                                                                      *
      entry sumpioq(tim )
*                                                                      *
************************************************************************
      do 3200 i = 1,maxppar
        if(ipi(1,i) .ne. 1)                                    goto 3200
        px = ppi(1,i)
        py = ppi(2,i)
        pz = ppi(3,i)
        etot = sqrt(pmass**2+px**2+py**2+pz**2)
        ptotpi(1,itet) = ptotpi(1,itet) + px
        ptotpi(2,itet) = ptotpi(2,itet) + py
        ptotpi(3,itet) = ptotpi(3,itet) + pz
        etotpi(itet)   = etotpi(itet) + etot
 3200 continue
      return
************************************************************************
*                                                                      *
      entry sumpiout(scalpro,num,isubs,ikin,tmin,tmax,ifit)
*                                                                      *
*         variables:                                                   *
*            isum     - output unit                                    *
*            scalpro  - 2*b*db /number of test-part.                   *
*            num      - number of test-part.                           *
*            ikin     - number of points in the spectra to print out;  *
*            tmin     - minimal kin. energy for exp. fit               *
*            tmax     - maximal kin. energy for exp. fit               *
*                                                                      *
************************************************************************
      do 2100 i=1,ikin
        epik(i)     = float(i-1) * bink
 2100 continue
      imin    = nint( tmin / bink) + 1
      imax    = nint( tmax / bink) + 1
      if(ifit.gt.0) call fitexplms(pitr,epik,ikin,itip,e0,b0,imin,imax)
      if(ifit.gt.0) then
      do 2101 i=1,maxpt
        pipt(i)     = sqrt(float((i-1)**2)*bint**2 + 1.0) * pmass
        ptsp(itip,i)= ptsp(itip,i)/sqrt(pipt(i))
 2101 continue
      call fitexplms(ptsp,pipt,maxpt,itip,e1,b1,3,14)
      do 2102 i=1,maxpt
        ptsp(itip,i)= ptsp(itip,i)*sqrt(pipt(i))
        ptsp(11  ,i)= ptsp(11  ,i)*sqrt(pipt(i))
 2102 continue
      end if
      facta = 4.0 * pi * dangle * bink
      val  = 180.0 / pi
      write(isum,'(/''c:kinetic energy spectrum of transv. pions'')')
      if(isosp.eq.0) write(isum,'(''c: all pions'')')
      if(isosp.eq.1) write(isum,'(''c: pion-'')')
      write(isum,'(''c: transverse pions in cms, dfi(degree)='',f5.2)')
     &                           val * dangle
      write(isum,'(''n: kinetic energy of pions (gev)'')')
      write(isum,'(''n: e * dszigma/dp**3'')')
      write(isum,'(''c: fitted slope (mev) '',f8.2)') e0
      write(isum,'(''c: time(fm/c)'',f7.2,9x,9(3x,f7.2))')
     &    (timpisp(itip+1-i),i=1,itip)
      write(isum,'(
     &''n:   x       n,lst0   n(fit),dt0   n         n         n    '',
     &''     n         n         n         n         n         n'')')
      do 2200 i=1,ikin
        write(isum,'(g11.2,11e10.3)') bink*float(i-1),
     &  scalpro*pitr(itip,i)/facta,scalpro*pitr(11,i)/facta,
     &  (scalpro*pitr(itip-j,i)/facta,j=1,itip-1)
 2200 continue
*
      write(isum,'(/''c:plot of rapidity distribution of pions'')')
      write(isum,'(''c: rapidity of target in cms:'',f6.3)') yref
      write(isum,'(''n: rapidity'')')
      write(isum,'(''n: yield'')')
      write(isum,'(''c: time(fm/c)'',10(f7.2,3x))')
     &    (timpisp(itip+1-i),i=1,itip)
      write(isum,'(
     &''n:   x      n,lst0       n         n         n         n    '',
     &''     n         n         n         n         n'')')
      do 2300 i=-maxpr,maxpr
      write(isum,'(g12.2,10e10.3)')
     &  binr*float(i),(rappi(itip+1-j,i)/float(isubs*num)/binr,j=1,itip)
2300  continue
*
      facta = 2.0 * pi * binr * bint
      ymax  = float(maxpr) * binr
      pitrmi= 0.0
      pitrma= float(maxpt) * bint
      write(isum,2050) yref,pitrmi,pitrma,bint,-ymax,ymax,binr
2050  format(/'n:',1h','momentum distribution of pions',1h'/
     &'c:rapidity of target in cms',f6.3/
     &'n:rapidity '/
     &'n:transverse momenta/mass'/
     &'n2: x = ',f5.2,' to ',f5.2,' by ',f5.2,
     &  '; y = ',f5.2,' to ',f5.2,' by ',f5.2,';')
      write(isum,2150)
     &  ((scalpro*rapptpi(itip,i,j)/facta,i=-maxpr,maxpr),j=0,maxpt)
2150  format(21(11(e10.3,',')/))
*
      fect = 4.0 * pi * bink
      write(isum,'(/''c:kinetic energy spectrum of pions'')')
      if(isosp.eq.0) write(isum,'(''c: all pions'')')
      if(isosp.eq.1) write(isum,'(''c: pion-'')')
      write(isum,'(''c: d(sigma)/dp**3'')')
      write(isum,'(''n: kinetic energy of pions (gev)'')')
      write(isum,'(''n: dszigma/dp**3'')')
      write(isum,'(''c: fitted slope (mev) '',f8.2)') e0
      write(isum,'(''c: time(fm/c)'',f7.2,9x,9(3x,f7.2))')
     &    (timpisp(itip+1-i),i=1,itip)
      write(isum,'(
     &''n:   x       n,lst0   n(fit),dt0   n         n         n    '',
     &''     n         n         n         n         n         n'')')
      do 2201 i=1,ikin
        write(isum,'(g11.2,11e10.3)') bink*float(i-1),
     &  scalpro*dsdp(itip,i)/fect,scalpro*dsdp(11,i)/fect,
     &  (scalpro*dsdp(itip-j,i)/fect,j=1,itip-1)
 2201 continue
*
      facta = bint * pmass
      write(isum,'(/''c:transverse mom. spectrum of pions'')')
      if(isosp.eq.0) write(isum,'(''c: all pions'')')
      if(isosp.eq.1) write(isum,'(''c: pion-'')')
      write(isum,'(''n: transverse mom of pions (gev/c)'')')
      write(isum,'(''n: 1./pt * d(sigma)/dpt'')')
      write(isum,'(''c: fitted slope (mev) '',f8.2)') e1
      write(isum,'(''c: time(fm/c)'',f7.2,9x,9(3x,f7.2))')
     &    (timpisp(itip+1-i),i=1,itip)
      write(isum,'(
     &''n:   x       n,lst0   n(fit),dt0   n         n         n    '',
     &''     n         n         n         n         n         n'')')
      do 2202 i=1,maxpt
        write(isum,'(g11.2,11e10.3)') bint*float(i-1)*pmass,
     &  scalpro*ptsp(itip,i)/facta,scalpro*ptsp(11,i)/facta,
     &  (scalpro*ptsp(itip-j,i)/facta,j=1,itip-1)
 2202 continue
*
      if(nangle .gt. 0) then
        write(isum,'(/''c:inv. cross-sect. of pions at eb= '',f5.1,
     &              '' gev'')') elab
        if(isosp.eq.0) write(isum,'(''c: all pions'')')
        if(isosp.eq.1) write(isum,'(''c: pion-'')')
      if(iframe.eq.0) write(isum,'(''c: in cms, dfi='',f4.1)')dangle*val
      if(iframe.eq.1) write(isum,'(''c: in lab, dfi='',f4.1)')dangle*val
        write(isum,'(''n: momentum of pions (gev)'')')
        write(isum,'(''n: e * dszigma/dp**3'')')
        do 2500 ja=1,nangle
          facta = 4.0 * pi * sin(angle(ja)) * dangle * binm
          write(isum,'(/''c: angle(degree)='',f6.2)') val * angle(ja)
          write(isum,'(''c: time(fm/c)'',10(f7.2,3x))')
     &          (timpisp(itip+1-i),i=1,itip)
          write(isum,'(
     &''n:   x       n,lst0       n         n         n         n   '',
     &''      n         n         n         n         n'')')
          do 2400 i=1,ikin
          write(isum,'(g12.2,10e10.3)')
     & binm*float(i-1),(scalpro*specpi(itip+1-j,ja,i)/facta,j=1,itip)
 2400     continue
 2500   continue
      end if
*
      return
      end
