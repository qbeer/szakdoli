
************************************************************************
*                                                                      *
      subroutine pauli(i,ntag,iseed,phase,xxx,yyy,zzz,px,py,pz)
*                                                                      *
*       variables:                                                     *
*          i        - number of particle (integer, input)              *
*          ntag     - flag which tells if phase-space is pauli-blocked *
*                     ntag =  0 => phase space open                    *
*                     ntag = -1 => phase space blocked                 *
*          iseed    - seed for random number generator (integer,input) *
*          phase    - phase space factor                (real, output) *
*          xxx,yyy  - middle point of two colliding nucleons           *
*             ,zzz                                                     *
*                                                                      *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      integer ihw,nthw
      common /nthwhw/  nthw
      common /counthw/ ihw
      integer i,ntag,iseed,maxind
      real*8 phase,xxx,yyy,zzz,px,py,pz
      parameter (maxind = 10000)
      real*8 paulf(maxind),count,sek,v00,pp,pf,pf2,rnx,ppx,ppy,ppz
      real*8 rad,zahl,pad,pad2,rn
      integer index(maxind), num,numpaul, ind,ix0,iy0,iz0
      integer jx,kx,jy,ky,jz,kz,lx,ly,lz,ic,ie,in,j,nruns,npaul,m,ntot
      save  count, num, pad2, numpaul, ind, index, paulf
*-----------------------------------------------------------------------
      real*8  v(0:3), pmatt(0:3), ppart(0:3)
      ntag=0
      sek=0.0
      phase = .0
      ix0 = nint(xxx)
      iy0 = nint(yyy)
      iz0 = nint(zzz)
!	return
c -----------------
      goto  8000    !   old version
c--------------------------
c       new version
c-------------------------
*
c     write(*,*)  '  begin   pauli ', ix0,iy0,iz0,rhob_4(0,ix0,iy0,iz0)
      if( abs(ix0).lt.maxx .and. abs(iy0).lt.maxx .and.
     &    abs(iz0).lt.maxz)                             then  ! grid
         v(0)  =  rhob_4(0,ix0,iy0,iz0)
         if ( v(0) .lt. .01)         goto 7000
         v(1)  = -rhob_4(1,ix0,iy0,iz0)
         v(2)  = -rhob_4(2,ix0,iy0,iz0)
         v(3)  = -rhob_4(3,ix0,iy0,iz0)
         v00  =  rhob_4(0,ix0,iy0,iz0) / rhob_4(6,ix0,iy0,iz0)
c           ppart(0) = p0 corrected in next line Gyuri
           ppart(0) = sqrt(e(i)**2+px**2+py**2+pz**2)
           ppart(1) = px
           ppart(2) = py
           ppart(3) = pz
           call lorentz_hw(v, ppart, pmatt)
           pp = pmatt(1)**2 + pmatt(2)**2 + pmatt(3)**2
!            pp = ppart(1)**2 + ppart(2)**2 + ppart(3)**2		!!!????????????  neu pmatt = ppart
           pf = (3.0  *v00*3./(16.*pi))**.333 * (2.*pi*.197)
c   the factor 3.0  estimates that the fermi sphere is 33%  populated
           pf2 = pf * pf
           if (pp .gt. pf2) goto 7000
           rnx  =  rn(iseed)
           phase = 0.666 * rnx
           if ( rnx .gt. .666)  ntag = -1
           goto    7000
       else
           goto    7000
       endif
ccc  ---------------------
c        now old version    - NOT USED
ccc  --------------------
 8000 continue
c     write(*,*)  '  in  pauli ', numpaul,ip,ir,pad2
      do 100 j=1,maxpar,numpaul
*
        if(id(1,j).ne.1) goto 100
*
         ppx=p(1,j)
         ppy=p(2,j)
         ppz=p(3,j)
*
        if((px-ppx)**2+(py-ppy)**2+(pz-ppz)**2.gt.pad2)
     &                                           goto 100
*
        jx=nint(r(1,j))
        kx=ix0-jx
           if(abs(kx).gt.ir) goto 100
        jy=nint(r(2,j))
        ky=iy0-jy
           if(abs(ky).gt.ir) goto 100
        jz=nint(r(3,j))
        kz=iz0-jz
           if(abs(kz).gt.ir) goto 100
*
        lx=nint(float(2*ip+1)*(r(1,j)-float(jx)))
        if(abs(lx) .eq. ip+1) lx = lx/abs(lx) * ip
        ly=nint(float(2*ip+1)*(r(2,j)-float(jy)))
        if(abs(ly) .eq. ip+1) ly = ly/abs(ly) * ip
        lz=nint(float(2*ip+1)*(r(3,j)-float(jz)))
        if(abs(lz) .eq. ip+1) lz = lz/abs(lz) * ip
*
        ic=1+(lz+ip)+(ly+ip)*(2*ip+1)+(lx+ip)*(2*ip+1)**2
        ie=1+(kz+ir)+(ky+ir)*(2*ir+1)+(kx+ir)*(2*ir+1)**2
*
        if(dm(ic,ie).le.1.0e-12) dm(ic,ie) = 0.0
        sek=sek+dm(ic,ie)*float(num)
*
  100 continue
*
c       write(*,*)  ' sek  in pauli ',num,sek, count
      phase = amin1(1.0,sek/count*float(numpaul))
      if (phase .gt. rn(iseed)) ntag = -1
*
      return
 7000 continue
!       write(20,*)  '  end of pauli ',v00, pf, pp, phase, ntag
      return
*
*-----------------------------------------------------------------------
      entry paulpro0(xxx,yyy,zzz)
*                                                                      *
*       variables:                                                     *
*          i        - number of particle (integer, input)              *
*          xxx,yyy  - middle point of two colliding nucleons           *
*             ,zzz                                                     *
*                                                                      *
*                                                                      *
*-----------------------------------------------------------------------
      ix0 = nint(xxx)
      iy0 = nint(yyy)
      iz0 = nint(zzz)
      ind = 0
*
      do 200 j=1,maxpar,numpaul
*
        if(id(1,j).ne.1) goto 200
*
        jx=nint(r(1,j))
        kx=ix0-jx
           if(abs(kx).gt.ir) goto 200
        jy=nint(r(2,j))
        ky=iy0-jy
           if(abs(ky).gt.ir) goto 200
        jz=nint(r(3,j))
        kz=iz0-jz
           if(abs(kz).gt.ir) goto 200
*
        lx=nint(float(2*ip+1)*(r(1,j)-float(jx)))
        if(abs(lx) .eq. ip+1) lx = lx/abs(lx) * ip
        ly=nint(float(2*ip+1)*(r(2,j)-float(jy)))
        if(abs(ly) .eq. ip+1) ly = ly/abs(ly) * ip
        lz=nint(float(2*ip+1)*(r(3,j)-float(jz)))
        if(abs(lz) .eq. ip+1) lz = lz/abs(lz) * ip
*
        ic=1+(lz+ip)+(ly+ip)*(2*ip+1)+(lx+ip)*(2*ip+1)**2
        ie=1+(kz+ir)+(ky+ir)*(2*ir+1)+(kx+ir)*(2*ir+1)**2
*
        if(dm(ic,ie).le.1.0e-12) dm(ic,ie) = 0.0
*
        index(ind+1) = j
        paulf(ind+1) = dm(ic,ie)*float(num)
        ind          = ind+1
        if(ind.eq.maxind) write(6,'(''too small indexblock in pauli'')')
        if(ind.eq.maxind) then
        write(50,*) "stop HS 49"
        stop
      endif
*
  200 continue
*
      return
*
*-----------------------------------------------------------------------
      entry paulpro1(i,phase,px,py,pz)
*                                                                      *
*       variables:                                                     *
*          i        - number of particle (integer, input)              *
*          ntag     - flag which tells if phase-space is pauli-blocked *
*                     ntag =  0 => phase space open                    *
*                     ntag = -1 => phase space blocked                 *
*          iseed    - seed for random number generator (integer,input) *
*          phase    - phase space factor                (real, output) *
*          xxx,yyy  - middle point of two colliding nucleons           *
*             ,zzz                                                     *
*                                                                      *
*                                                                      *
*-----------------------------------------------------------------------
      sek=0.0
      phase=0.0
      if(ind.eq.0)                                                return
*
      do 300 in=1,ind
*
         j  =index(in)
*
         ppx=p(1,j)
         ppy=p(2,j)
         ppz=p(3,j)
*
        if((px-ppx)**2+(py-ppy)**2+(pz-ppz)**2.gt.pad2)
     &                                           goto 300
*
        sek=sek+paulf(in)
*
  300 continue
*
      phase = amin1(1.0,sek/count*float(numpaul))
*
      return
*
*-----------------------------------------------------------------------
      entry paulin(nruns,npaul,m)
*
*       this entry computes the length of a phase-space-cell in
*       coordinate- and momentum-space for the pauli-blocking
*
*     variables:
*                nruns - number of testparticles per nucleon
*                m     - output unit
*                ntot  - total number of test particles
*                npaul - ntot/npaul testparticle is taking into account
*
*     this entry is used to initialize the value for pad, and
*     count for the use in pauli
*-----------------------------------------------------------------------
*
      num    = nruns
      numpaul= npaul
*      rad    = 3.0
      rad = float(ir)
      zahl   = 1./6.
      pad    = (3.0*zahl/16.0/pi)**(1./3.)*(2.*pi*.197)/rad
      pad2   = pad**2
      count  = zahl*float(num)
*
      if (m .ne. 0)
     &write(m,'(/''c:'',
     &43x,''==== initialization for subroutine pauli =============''/
     &''c:'',
     &47x,'' num of part incl. in pauli / total   = '',i5/''c:'',
     &47x,'' num of part per phase-space-cell     = '',f10.5/''c:'',
     &47x,'' length of cell in coord-space (fm)   = '',f10.5/''c:'',
     &47x,'' radius of cell in momen-space (1/fm) = '',f10.5/
     &       )')  numpaul,count, rad, pad/.197
*
      return
      end
