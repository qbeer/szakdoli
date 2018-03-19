************************************************************************
*                                                                      *
      subroutine mesonin(wmin,wmax)
*                                                                      *
*       purpose:    calculating vector meson production from:          *
*                                n+n  collision      (entry mesonbb)   *
*                                n+pi collision      (entry mesonpi)   *
*                   -------------------------------------------------  *
*                   initialization by subroutine mesonin               *
*                   -------------------------------------------------  *
*                   final cross section by entry mesonout              *
*                   -------------------------------------------------  *
*                                                                      *
*       method :    calculating and summing up the probabilities       *
*                      to produce a meson at a given momentum          *
*                      and than integrating over the meson angles      *
*                                                                      *
*     let j the index of the integration point in the dilepton momentum*
*          space, then qq and qy give the appopriate momenta in the    *
*                                                                      *
*         qq(i,j) /i=1-4/   -  i. coord. of the j. point in mom.space  *
*         qy(1,j)           -  the transverse mom belonging to j.      *
*         qy(2,j)           -  the rapidity       belonging to j.      *
*         qy(3,j)           -  angle in trans.dir.belonging to j.      *
*         sig(1,j)          -  for pi+pi                               *
*         sig(2,j)          -  for rho                                 *
*         sig(3,j)          -  for omega                               *
*         dsdm(1,j)         -  for pi+pi                               *
*         dsdm(2,j)         -  for n+n              ->n+n rho          *
*         dsdm(3,j)         -  for n*delta          ->n+n rho          *
*         dsdm(4,j)         -  for n+r, delta+delta ->n+n rho          *
*         dsdm(5,j)         -  for pi+n             ->n+n rho          *
*         dsdm(6,j)         -  for pi+r             ->n+n rho          *
*         dsdm(7,j)         -  for n+n              ->n+n omega        *
*         dsdm(8,j)         -  for n*delta          ->n+n omega        *
*         dsdm(9,j)         -  for n+r, delta+delta ->n+n omega        *
*         dsdm(10,j)        -  for pi+n             ->n+n omega        *
*         dsdm(11,j)        -  for pi+r             ->n+n omega        *
*         dsdm(12,j)        -  for total rho                           *
*         dsdm(13,j)        -  for total omega                         *
*       variables:                                                     *
*         yref    - -1 * target rapidity in the frame                  *
*         my      - number of rapidities               (integer,input) *
*         mqt     - number of transverse momenta       (integer,input) *
*         mf      - number of angles                   (integer,input) *
*                                                                      *
************************************************************************
      implicit none
*-----------------------------------------------------------------------

      include"cominput"
      include"common"
      real*8 wmax,wmin
      integer maxq,maxy,maxde,maxdm,isdm,nmas,ntotq,ntoty
      integer npiann,ii,jj,jma,ima,iqt,iy,nde,ij,iph,nq,nr,np,ire
      integer i1,i2,ireac,ixx,iyy,izz,jjm,id2,iz1,iz2

      real*8 p3,beta,xmasme,dsdm,sigma,sigma1,sigma2,sigma3,ymin,ymax
      real*8 factpi,sigpipi,dy,dqt,dphi,xmmas
      real*8 dmas,qtra,phi,tramass,xmesmas,pt,rapy,fi,densi,denst,gamma
      real*8 srt,etotal,xxx,yyy,zzz,sig0,rn,facreac,em12,em22,emm2,srtom
      real*8 srtro,sigom,sigro,ommasmax,romasmax,omnorm,ronorm,szigmom,s
      real*8 szigmro,traf,xmmas2,pmax2,pmax,ommasd,romasd,qbeta,eprim
      real*8 ppri2,ppri,sigm,ppr,fac,fac2,xx,yy,zz,r2,rr,qqx,qqy,qqz,qqe
      real*8 facta,e3,pbeta,transf,szigma,scal,yref,spr
      real*8 vol1,volm,event,yrap,phase,qqabs
      real*8 ede,pdx,pdy,pdz,dm2,gamm,trafa,betax,betay,betaz,q0,qq2
      real*8 rapydel,qbtra,trams,xyz,rap1,rap2,rap,dif1,dif2,qzpr
      real*8 q0pr,factu,phasv,sigmaom,sigmaro,densste
      integer ny,nqt,nf
      parameter     (maxq =60000)
      parameter     (maxy =3000)
      parameter     (maxde=20)
      parameter     (maxdm=35)
*----------------------------------------------------------------------*
      real*8     qq(4,maxq),qy(0:3,maxq),sig(3,maxy,maxde)
      integer iqq(3,maxy,maxde)
*----------------------------------------------------------------------*
      dimension p3(3),beta(3)
      dimension xmasme(maxdm)
      dimension dsdm(3,maxdm)
      dimension isdm(3,maxdm)
      dimension sigma(13,maxdm), sigma1(13),sigma2(13),sigma3(3)
*----------------------------------------------------------------------*
*          save for the whole                                          *
      save qy,qq,sig,iqq,nmas,ntotq,dy,dqt,dphi,ntoty,densste
      save xmasme,ny,nqt,nf
      save sigpipi
*----------------------------------------------------------------------*
      save ymin,ymax
      save sigma, sigma1,npiann
*----------------------------------------------------------------------*
      densste=0.2
      ymin   = wmin
      ymax   = wmax
      factpi = 10.0 * dispi**2 * pi
      sigpipi = sigpiom/factpi
      nmas   = max0(1,nmesma  )
      ny     = max0(1,nyk    )
      nqt    = max0(1,nqtk   )
      nf     = max0(1,nfk    )
      ntotq  = nmas*ny*nqt*nf
      ntoty  = nmas*ny*nqt
      write(*,*) 'vectmesin:', ntotq, dispi, pi
      npiann = 0
      if(ntoty .gt. maxy .or. ntotq.gt.maxq) then
        write(isum,'(''c: too many q-s are to be calculated for mes'')')
        write(*,'(''c: too many q-s are to be calculated for mes'')')
        stop
      end if
***                                                                  ***
***                                                                  ***
      dy     = (ymax-ymin) / float(ny)
      dqt    = qtmaxk / float(nqt)
      dphi   = 2. * pi / float(nf)
      ii = 0
      jj = 0
      jma= 0
      xmmas = xmesmi -dmesma
      do 1500 ima= 1,nmas
        dmas  = dmesma
        if(xmmas+dmesma.ge.xmesm2 .and. jma.lt.nmesm2) then
          dmas  = dmesm2
          jma   = jma+1
        endif
        xmmas = xmmas + dmas
        xmasme(ima) = xmmas
        do 1400 iqt= 1,nqt
          qtra = (float(iqt) - 0.5) * dqt
          do 1300 iy = 1,ny
            yrap = ymin + (float(iy) - 0.5) * dy
            jj = jj + 1
            do 1101 nde= 1,maxde
              do 1100 ij = 1,3
                iqq(ij,jj,nde) = 0
                sig(ij,jj,nde) = 0.0
 1100         continue
 1101       continue
            do 1200 iph = 0,nf-1
              ii = ii + 1
              phi = - pi + (float(iph)+0.5) * dphi
              qy(0,ii) = xmmas
              qy(1,ii) = qtra
              qy(2,ii) = yrap
              qy(3,ii) = phi
              tramass  = sqrt(xmmas**2+qtra**2)
              qq(1,ii) = qtra * cos(phi)
              qq(2,ii) = qtra * sin(phi)
              qq(3,ii) = tramass * sinh(yrap)
              qq(4,ii) = tramass * cosh(yrap)
 1200       continue
 1300     continue
 1400   continue
 1500 continue
      do ij=1,13
        sigma1(ij) = 0.0
        do ima=1,maxdm
          sigma(ij,ima) = 0.0
        enddo
      enddo
      do ij=1,3
        do ima=1,maxdm
          dsdm(ij,ima)   = 0.0
        enddo
      enddo
      return
*----------------------------------------------------------------------*
      entry pipiann(xmesmas,pt,rapy,fi,densi)
*       variables:                                                     *
*         xmesmas - invariant mass  of the pion pair      (real,input) *
*         pt      - transverse mom. of the pion pair      (real,input) *
*         rapy    - rapidity        of the pion pair      (real,input) *
*         fi      - azymuthal angle of the pion pair      (real,input) *
*----------------------------------------------------------------------*
c      write(*,*) 'pipiann',rho0,densste,dmas,dqt,dy,dphi
      npiann = npiann + 1
      if(xmesmas.lt.xmasme(1)-0.5*dmesma .or.
     &   xmesmas.gt.xmasme(nmas)+0.5*dmesma)                   return
      ima     = 0
 2100 ima     = ima+1
      if(xmesmas.gt..5*(xmasme(ima)+xmasme(ima+1)).and.ima.lt.nmas)
     &                                  goto 2100
      dmas   = dmesma
      if((ima.ne.1).and.(ima.ne.nmas))
     &       dmas=.5*(xmasme(ima+1)-xmasme(ima-1))
      nq     = nint((pt-qy(1,1))/dqt+1.0)
      nr     = nint((rapy-qy(2,1))/dy+1.0)
      np     = nint((fi-qy(3,1))/dphi+1.0)
      ii     = (ima-1)*nqt*ny*nf+(nq-1)*ny*nf+(nr-1)*nf+np
      jj     = (ima-1)*nqt*ny+(nq-1)*ny+nr
      if((ima.lt.1)  .or. (ima.gt.nmas))                   return
      if((nq.lt.1)  .or. (nq.gt.nqt ))                   return
      if((nr.lt.1)  .or. (nr.gt.ny  ))                   return
      if((np.lt.1)  .or. (np.gt.nf  ))                   return
      denst= densi/rho0
      nde  = nint(denst/densste+1.0)
      nde  = min(nde,maxde)
      sig(1,jj,nde) = sig(1,jj,nde)+sigpipi/(dmas*qy(1,ii)*dqt*dy*dphi)
      iqq(1,jj,nde) = iqq(1,jj,nde) + 1
      sigma(1,ima) = sigma(1,ima) + sigpipi/dmas
      sigma1(1) = sigma1(1) + sigpipi
      return
************************************************************************
      entry mesonbb(i1,i2,beta,gamma,srt,etotal,xxx,yyy,zzz,sig0)
*       variables:                                                     *
*         ireac   - 1->n+n; 2->n+d; 3->d+d 4->n+r reac.(integer,input) *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*         srt     - cms energy                                         *
*         etotal  - 2 particle energy in the frame                     *
*         xxx...zzz - coordinates of the event                         *
*         sig0    - cross section corresponds to calling this routine  *
*     cross-sections are taken from fit meson.fort(ppmeson)            *
*     for channels including delta we use a factor facreac multiplying *
*                  the n+n cross sections (taken from ko and randrup)  *
*              delta + n:     facreac=0.75                             *
*              delta + delta: facreac=0.5                              *
*----------------------------------------------------------------------*
c      write(*,*) 'mesonpro0 srt xmesmi rmass ncalmes',
c     &  srt,xmesmi,rmass,ncalmes
      if(srt.lt.xmesmi+2.*rmass)                                return
      if(rn(iseed) .gt. 1./float(ncalmes))                      return
***                                                                  ***
***         preparing    for the pauli blocking                      ***
***                                                                  ***
      call paulpro0(xxx,yyy,zzz)
***                                                                  ***
      ireac = 0
      ire   = 0
      if(i1+i2 .eq. 2) ireac = 1
      if(i1+i2 .eq. 3) ireac = 2
      if(i1+i2 .ge. 4) ireac = 4
      if(i1.eq.2 .and. i2.eq.2) ireac = 3
      if(i1+i2 .eq. 2) ire = 1
      if(i1+i2 .eq. 3) ire = 2
      if(i1+i2 .ge. 4) ire = 3
      if(ireac.le.0 .or. ireac.gt.4) write(6,'('' hiba,ireac,meson'')')
      if(ire.le.0 .or. ire.gt.3) write(6,'('' hiba,ire,meson'')')
      facreac = 1.0
      if(i1+i2 .eq. 3)          facreac = 0.75
      if(i1.eq.2 .and. i2.eq.2) facreac = 0.5
      em12  = rmass**2
      em22  = rmass**2
      emm2  = (rmass+rmass)**2
      s     = srt**2
      srtom = 2.0*rmass+omass-0.5*owidth
      srtro = 2.0*rmass+romas-0.5*rowidth
      sigom = 0.0
      sigro = 0.0
      if(srt.gt.srtom) sigom=0.36*(srt-srtom)**1.4/(1.25+(srt-srtom)**2)
      if(srt.gt.srtro) sigro=0.24*(srt-srtro)/(1.40+(srt-srtro)**2)
      ommasmax= amin1(srt-2.0*rmass,omass+0.5*owidth*gammesma)
      romasmax= amin1(srt-2.0*rmass,romas+0.5*rowidth*gammesma)
      omnorm= atan(2.0*(ommasmax-omass)/owidth) + atan(gammesma)
      ronorm= atan(2.0*(romasmax-romas)/rowidth) + atan(gammesma)
      if(ommasmax.lt.omass-0.5*owidth*gammesma) omnorm= 0.0
      if(romasmax.lt.romas-0.5*rowidth*gammesma) ronorm= 0.0
      szigmom = sigom*float(ncalmes)/sig0 * facreac
      szigmro = sigro*float(ncalmes)/sig0 * facreac
c      write(*,*) 'mesonpro:srt,xmesmi,szigmom,sigom',
c     & srt,xmesmi,szigmom,sigom,szigmro,sigro
      if(srt.ge.xmesmi) sigma1(ire+6) = sigma1(ire+6) + szigmom
      if(srt.ge.xmesmi) sigma1(13   ) = sigma1(13   ) + szigmom
      if(srt.ge.xmesmi) sigma1(ire+1) = sigma1(ire+1) + szigmro
      if(srt.ge.xmesmi) sigma1(12   ) = sigma1(12   ) + szigmro
c zm:
c      write(*,*) 'No. of rho mesons: ', sigma1(12)
      traf  = gamma / (gamma+1.0)
      ii    = 0
      jj    = 0
      do 4303 ima=1,nmas
        xmmas = qy(0,ii+1)
        xmmas2= xmmas**2
      pmax2 =.25*(s-(rmass+rmass+xmmas)**2)*(s-(rmass+rmass-xmmas)**2)/s
        if(pmax2 .le. 0.0)                                        return
c       pmax  = the maximal meson momentum
        pmax  = sqrt(pmax2)
        ommasd = 0.0
        romasd = 0.0
        if(abs(xmmas-omass) .lt. 0.5*gammesma*owidth)
     &    ommasd= 0.5*owidth/((xmmas-omass)**2+0.25*owidth**2)/omnorm
        if(abs(xmmas-romas) .lt. 0.5*gammesma*rowidth)
     &    romasd= 0.5*rowidth/((xmmas-romas)**2+0.25*rowidth**2)/ronorm
        sigma(ire+6,ima) = sigma(ire+6,ima) + szigmom * ommasd
        sigma(13   ,ima) = sigma(13   ,ima) + szigmom * ommasd
        sigma(ire+1,ima) = sigma(ire+1,ima) + szigmro * romasd
        sigma(12   ,ima) = sigma(12   ,ima) + szigmro * romasd
c-----------------------------------------------------------------------
c           density dependence
        ixx = nint(xxx)
        iyy = nint(yyy)
        izz = nint(zzz)
        denst = 0.0
      if(iabs(ixx).le.maxx.and.iabs(iyy).le.maxx .and.iabs(izz).le.maxz)
     &      denst = rhb(ixx,iyy,izz) / rho0
        nde  = nint(denst/densste+1.0)
        nde  = min(nde,maxde)
c***********************************************************************
        do 4302 iqt=1,nqt
        do 4301 iy =1,ny
          jj = jj + 1
        do 4300 iph=1,nf
          ii = ii + 1
*
*   lorentz-transformation in i1-i2-c.m. system
*
          qbeta = beta(1)*qq(1,ii) + beta(2)*qq(2,ii) + beta(3)*qq(3,ii)
*
*   eprim:    meson energy in i1-i2-c.m. system
*
          eprim = gamma * (qq(4,ii) - qbeta)
*   spr:      total energy**2 of the nucleon+nucleon system in their cms
          spr   = srt * (srt - 2. * eprim) + xmmas2
          if(spr .le. emm2)                                    goto 4300
          ppri2 = eprim**2 - xmmas2
          ppri  = sqrt(ppri2)
c       sigm is calculated by some phase space model (depend on pmax)
          sigm  = 3.*eprim/(pi*pmax**3)*(1.-ppri/pmax)
***                                                                  ***
***              looking for the pauli blocking                      ***
***                                                                  ***
*   ppr:      momentum of the nucleon in the nucleon+nucleon  cms
        ppr    = sqrt(amax1(0.,(0.25*(spr-em12-em22)**2-em12*em22)/spr))
        fac    = - 0.5 * (1. + (em12-em22)/spr)
        fac2   = ppr / (spr + (srt - eprim) * sqrt(spr))
 4100   continue
        xx     = rn(iseed) - 0.5
        yy     = rn(iseed) - 0.5
        zz     = rn(iseed) - 0.5
        r2     = xx**2 + yy**2 + zz**2
        if( (r2.lt.0.0001) .or. (r2.gt.0.25) )                 goto 4100
        rr     = sqrt(r2)
c
*   lorentz-transformation of meson momentum in i1-i2-c.m. system
        transf = gamma * (traf * qbeta - qq(4,ii) )
        qqx    = qq(1,ii) + beta(1) * transf
        qqy    = qq(2,ii) + beta(2) * transf
        qqz    = qq(3,ii) + beta(3) * transf
        qqe    = (qqx * xx + qqy * yy + qqz * zz) / rr
        facta  = fac + fac2 * qqe
*   p3:                    nucleon momentum in i1-i2-c.m. system
        p3(1)  = facta * qqx + xx/rr * ppr
        p3(2)  = facta * qqy + yy/rr * ppr
        p3(3)  = facta * qqz + zz/rr * ppr
        e3     = sqrt(em12+p3(1)**2+p3(2)**2+p3(3)**2)
*   p3:                    nucleon momentum in observable system
        pbeta  = beta(1)*p3(1) + beta(2)*p3(2) + beta(3)*p3(3)
        transf = gamma * (traf * pbeta + e3)
        p3(1)  = p3(1) + beta(1) * transf
        p3(2)  = p3(2) + beta(2) * transf
        p3(3)  = p3(3) + beta(3) * transf
        call paulpro1(i1,phase,p3(1),p3(2),p3(3))
        szigma = sigm * (1.0 - phase)
        p3(1)  = beta(1) * etotal - p3(1) - qq(1,ii)
        p3(2)  = beta(2) * etotal - p3(2) - qq(2,ii)
        p3(3)  = beta(3) * etotal - p3(3) - qq(3,ii)
        call paulpro1(i2,phase,p3(1),p3(2),p3(3))
        szigma = szigma * (1.0 - phase)
        sig(3,jj,nde)= sig(3,jj,nde) + szigma*szigmom*ommasd
        sig(2,jj,nde)= sig(2,jj,nde) + szigma*szigmro*romasd
        iqq(2,jj,nde)= iqq(2,jj,nde) + 1
 4300 continue
 4301 continue
 4302 continue
 4303 continue
c      write(*,*) 'end of mesonbb'
      return

*----------------------------------------------------------------------*
*                                                                      *
      entry vecmespi(ede,pdx,pdy,pdz,srt,xxx,yyy,zzz,sig0,id2,iz1,iz2)
*       variables:                                                     *
*         ede,pdx,- energy, momentum of the pi+baryon sys.(real,input) *
*         srt     - cms energy (mass of the pi+baryon sys.(real,input) *
*         xxx...zzz - coordinates of the event                         *
*         iseed   - seed for random number generator   (integer,input) *
*         beta    - velocity of the 2 part. cms in the frame           *
*         gamma   - gamma corresponds to beta                          *
*                                                                      *
*----------------------------------------------------------------------*
      if(iz1+iz2.ne.1 .and. iz1+iz2.ne.0)                    return
c      if(rn(iseed) .ge. dalinum)                            return
c      write(*,*)'in vecmespi id2,id1=',id1,id2
***                                                                  ***
***         preparing    for the pauli blocking                      ***
***                                                                  ***
      call paulpro0(xxx,yyy,zzz)
***                                                                  ***
      if(id2 .lt. 1) write(*,*) 'hiba vectmespi id2 =',id2
      if(id2 .eq. 1) ire = 1
      if(id2 .gt. 1) ire = 2
      dm2   = srt**2
      gamm  = ede / srt
      trafa = gamm / (gamm + 1.0)
      betax = pdx / ede
      betay = pdy / ede
      betaz = pdz / ede
c     write(6,*) jd, xmass,kk, id(1,kk),ipert
c      write(*,*) 'in wecmespi2:'
      srtom = rmass+omass-0.5*owidth
      srtro = rmass+romas-0.5*rowidth
      sigom = 0.0
      sigro = 0.0
      if(srt.gt.srtom)
     &      sigom=1.38*(srt-srtom)**1.6/(0.0011+(srt-srtom)**1.7)
      if(srt.gt.srtro)
     &      sigro=1.5 *(srt-srtro)**2.2/(0.018 +(srt-srtro)**3.5)
      if(id2.ne.2 .and. iz1.eq.0) sigom = 0.5 * sigom
      if(id2.eq.2 .and. (iz2.eq.2.or.iz2.eq.-1)) sigom=1.5*sigom*pidelka
      if(id2.eq.2 .and. (iz2.eq.1.or.iz2.eq.0).and.iz1.ne.0)
     &      sigom=0.5*sigom*pidelka
      if(id2.ne.2 .and. iz1.eq.0) sigro = 0.0
      if(id2.eq.2 .and. (iz2.eq.2.or.iz2.eq.-1)) sigro=1.5*sigro*pidelka
      if(id2.eq.2 .and. (iz2.eq.1.or.iz2.eq.0).and.iz1.ne.0)
     &      sigro=0.5*sigro*pidelka
      ommasmax= amin1(srt-rmass,omass+0.5*owidth*gammesma)
      romasmax= amin1(srt-rmass,romas+0.5*rowidth*gammesma)
      omnorm= atan(2.0*(ommasmax-omass)/owidth) + atan(gammesma)
      ronorm= atan(2.0*(romasmax-romas)/rowidth) + atan(gammesma)
      if(ommasmax.lt.omass-0.5*owidth*gammesma) omnorm= 0.0
      if(romasmax.lt.romas-0.5*rowidth*gammesma) ronorm= 0.0
      szigmom = sigom/sig0
      szigmro = sigro/sig0
c      write(*,*) 'mesonpro:srt,xmesmi,szigmom,sigom',
c     & srt,xmesmi,szigmom,sigom,szigmro,sigro
      if(srt.ge.xmesmi) sigma1(ire+9) = sigma1(ire+9) + szigmom
      if(srt.ge.xmesmi) sigma1(13   ) = sigma1(13   ) + szigmom
      if(srt.ge.xmesmi) sigma1(ire+4) = sigma1(ire+4) + szigmro
      if(srt.ge.xmesmi) sigma1(12   ) = sigma1(12   ) + szigmro
      ii    = 0
      jj    = 0
c      write(*,*) 'in vectmes', nqt,ny,nf
      do 3400 ima=1,nmas
        ii  = (ima-1) * nqt * ny * nf
        jj  = (ima-1) * nqt * ny
        xmmas = qy(0,ii+1)
        xmmas2= xmmas**2
        q0      = 0.5 * (dm2 + xmmas2 - rmass**2) / srt
        if(q0 .le. xmmas)                                    goto 3400
        ommasd = 0.0
        romasd = 0.0
        if(abs(xmmas-omass) .lt. 0.5*gammesma*owidth)
     &    ommasd= 0.5*owidth/((xmmas-omass)**2+0.25*owidth**2)/omnorm
        if(abs(xmmas-romas) .lt. 0.5*gammesma*rowidth)
     &    romasd= 0.5*rowidth/((xmmas-romas)**2+0.25*rowidth**2)/ronorm
        sigma(ire+9,ima) = sigma(ire+9,ima) + szigmom * ommasd
        sigma(13   ,ima) = sigma(13   ,ima) + szigmom * ommasd
        sigma(ire+4,ima) = sigma(ire+4,ima) + szigmro * romasd
        sigma(12   ,ima) = sigma(12   ,ima) + szigmro * romasd
        qq2     = q0**2 - xmmas2
        qqabs   = sqrt(qq2)

c***********************************************************************
c           density dependence
          ixx = nint(xxx)
          iyy = nint(yyy)
          izz = nint(zzz)
          denst = 0.0
          if(iabs(ixx).le.maxx .and. iabs(iyy).le.maxx .and.
     &       iabs(izz).le. maxz)
     &      denst = rhb(ixx,iyy,izz)/rho0
          nde  = nint(denst/densste+1.0)
          nde  = min(nde,maxde)
c***********************************************************************
        rapydel = 0.5 * log((1.+betaz)/(1.-betaz))
        do 3300 iqt=1,nqt
          do 3200 iy=1,ny
            jj = jj + 1
            do 3100 iph=1,nf
              ii = ii + 1
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
              qbtra = qq(1,ii) * betax + qq(2,ii) * betay
              trams = sqrt(qy(0,ii)**2 + qy(1,ii)**2)
              xyz   = cosh(rapydel) * (q0/gamm+qbtra)/trams
              if(abs(xyz) .le. 1.0)                           goto 3100
              rap1= rapydel - log(xyz + sqrt(xyz**2 - 1.0))
              rap2= rapydel + log(xyz + sqrt(xyz**2 - 1.0))
              dif1= abs(rap1-qy(2,ii))
              dif2= abs(rap2-qy(2,ii))
              if((dif1.gt.dy/2.).and.(dif2.gt.dy/2.))         goto 3100
              if(dif1.le.dy/2.) rap = rap1
              if(dif2.le.dy/2.) rap = rap2
              p3(1)  = pdx - qq(1,ii)
              p3(2)  = pdy - qq(2,ii)
              p3(3)  = pdz - qq(3,ii)
              call paulpro1(id2,phase,p3(1),p3(2),p3(3))
              qzpr = trams * sinh(rap)
              q0pr = trams * cosh(rap)
              factu = 4.*pi * qqabs * gamm  * dy * abs(qzpr- betaz*q0pr)
              phasv= dy*dqt*dphi*qy(1,ii)
              sigmaom= szigmom * ommasd / factu
              sigmaro= szigmro * romasd / factu
              sig(2,jj,nde) = sig(2,jj,nde)+sigmaro*(1.-phase)
              sig(3,jj,nde) = sig(3,jj,nde)+sigmaom*(1.-phase)
              iqq(3,jj,nde) = iqq(3,jj,nde) + 1
 3100       continue
 3200     continue
 3300   continue
 3400 continue
      return
*----------------------------------------------------------------------*
      entry mesonout(scal,yref)
*----------------------------------------------------------------------*
*  in sigma(i) the total x-section is calculated by sum of collisions
*  in sigma1(i) the total x-section is calculated by integrated over
*               dqt, dy, dphi and density
*  in sig(i,j,k) gives dsigma/dqt*qt /dy /dphi at density, at j. q
*                                         and  at i. channel
*  dsdm(i) gives the dsigma for channel i; i=6 is sum of all
*  angular filter is commented out
*  density dependence is not studied, although possible
*----------------------------------------------------------------------*
      do 5100 ima = 1,maxdm
        do 5100 jjm = 1,3
          dsdm(jjm,ima) = 0.0
          isdm(jjm,ima) = 0
          sigma3(jjm) = 0.0
 5100 continue
      do ij=1,13
        sigma2(ij) = 0.0
      enddo
      ii = 0
      jj = 0
      do 5900 ima= 1,nmas
        xmmas = qy(0,ii+1)
        dmas  = dmesma
      if(ima.ne.1.and.ima.ne.nmas) dmas=.5*(xmasme(ima+1)-xmasme(ima-1))
        do ij=1,13
          sigma2(ij) = sigma2(ij)+sigma(ij,ima)*dmas
        enddo
        do 5700 iqt= 1,nqt
          do 5500 iy = 1,ny
            jj = jj + 1
              ii = ii + nf
              vol1 = dy * dqt * dphi * qy(1,ii)
              volm = dmas * dy * dqt * dphi * qy(1,ii)
              do 5301 nde=1,maxde
              do 5300 ij =1,3
                dsdm(ij,ima)= dsdm(ij,ima) + sig(ij,jj,nde)*vol1
                sigma3(ij)  = sigma3(ij)   + sig(ij,jj,nde)*volm
 5300         continue
 5301         continue
 5500     continue
 5700   continue
 5900 continue
*
 100  format(//79('*')//'# vector meson production'//
     &      'c:bombarding energy:',f8.2,' GeV'/
     &      'c:impact parameter: ',f8.2,' fm')
 101  format(/'c: ny:',i4,5x,'nqt:',i4,5x,'nphi:',i4)
 102  format('c:rapidity interval:',f9.4,2x,'-',f8.4,3x,'lab. rapidity:'
     &      ,f6.3,3x,'qtmaxk:',f6.3)
 103  format('#  n+n',7x,'n+delta   n+res     n+pi      r+pi      all')
 104  format(6e10.3)
 105  format('# mass     n+n',9x,'n+delta     n+res       n+pi',
     &          '        r+pi')
 106  format(f6.3,5e12.3)
 107  format('n: mass pi+pi-> omega  lat.',6x,'rho',9x,'lat.      omega'
     &      ,6x,'lat.')
 108  format(50x,e10.3)
      write(isum,100) elab,b
      write(isum,101) ny,nqt,nf
      write(isum,102) ymin,ymax,yref,qtmaxk
      write(isum,'(/''c: pi+pi cross section:'',3e11.3,i6)')
     &      sigma1(1)*scal,sigma2(1)*scal,sigma3(1)*scal,npiann
*
c zm:
      write(isum,'(/''# The following data are cross section '//
     &     '# obtained from event numbers by multiplication by'//
     &     '# scal = '',e11.3)') scal
c end zm
        write(isum,'(/''# total cross sections [microb]'')')
        write(isum,'(''#  rho'')')
        write(isum,103)
        write(isum,104) (sigma1(ij)*scal,ij=2,6),sigma1(12)*scal
        write(isum,104) (sigma2(ij)*scal,ij=2,6),sigma2(12)*scal
        write(isum,108)  sigma3(2)*scal
        write(isum,'(''#  omega'')')
        write(isum,103)
        write(isum,104) (sigma1(ij)*scal,ij=7,11),sigma1(13)*scal
        write(isum,104) (sigma2(ij)*scal,ij=7,11),sigma2(13)*scal
        write(isum,108) sigma3(3)*scal
*
        write(isum,'(/'' rho mass spectra'')')
        write(isum,105)
        do ima = 1,nmas
          write(isum,106) xmasme(ima),(sigma(jj,ima)*scal,jj=2,6)
        enddo
*
        write(isum,'(/'' omega mass spectra'')')
        write(isum,105)
        do ima = 1,nmas
          write(isum,106) xmasme(ima),(sigma(jj,ima)*scal,jj=7,11)
        enddo
*
        write(isum,'(/''# vector meson mass spectra'')')
        write(isum,'(''# mass (GeV)'')')
        write(isum,'(''# dsigma (mikrob/GeV)'')')
        write(isum,107)
        do ij=1,nmas
        write(isum,'(f7.3,6e11.3)')
     &    xmasme(ij),
     &    sigma(1,ij)*scal,dsdm(1,ij)*scal,
     &    sigma(12,ij)*scal,dsdm(2,ij)*scal,
     &    sigma(13,ij)*scal,dsdm(3,ij)*scal
        enddo
      return
*----------------------------------------------------------------------*
      entry finalmes(event)
*----------------------------------------------------------------------*
      write(mmespri,'(4i8)') nmas,nqt,ny,maxde
 109  format('c: mass    qt     y     rho      pi+pi     rho     ',
     &      '  omega')
      write(mmespri,109)
      do 7000 jj = 1,ntoty
        ii = nf * jj
      do 7000 nde= 1,maxde
        denst= float(nde-1)*densste
      write(mmespri,'(4f7.3,2x,5e10.3)')qy(0,ii),qy(1,ii),qy(2,ii),
     & denst,(sig(ij,jj,nde)*dphi/event,ij=1,3)
 7000 continue
      return
      end
