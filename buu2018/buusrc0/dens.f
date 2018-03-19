************************************************************************
*                                                                      *
      subroutine dens(minnum,maxnum,num,nesc)
*                                                                      *
*       purpose:     calculation of nuclear density from spatial       *
*                    distribution of testparticles                     *
*                    and average momentum in spacial cell              *
*                    125*125 points gaussian smearing                  *
*                    width of gaussian is 1.0 fm, cutoff=sqrt(5.0) fm  *
*                                                                      *
*       variables (all input, all integer)                             *
*         minnum  -  first testparticle treated in one run for density *
*         maxnum  -  last testparticle treated in one run for density  *
*         num     -  number of testparticles per nucleon               *
*         nesc    -  number of escaped particles      (integer,output) *
*         isp     -  0-> no, 1-> calculate <p> for suumry              *
*                                                                      *
************************************************************************
      implicit none
      include"common"

      integer num,i,ib,ic,ie,jy,kx,ky,kz,irun,ix,iy,iz,jx
      integer minnum,maxnum,nesc,jz
      real*8 gauss, dcut ,rd, xi,yi,zi,sek,xj,yj,zj,rsqr

      write(*,*)'begin of dens'
*
*-----------------------------------------------------------------------
      do 340 iz = -maxz,maxz
      do 240 iy = -maxx,maxx
      do 140 ix = -maxx,maxx
            rhb(ix,iy,iz) = 0.0
  140 continue
  240 continue
  340 continue
*
*-----------------------------------------------------------------------
      nesc  = 0
*
      do 600 irun=1,num
      do 400 i = minnum+(irun-1)*maxb,maxnum+(irun-1)*maxb
        if(id(1,i) .eq. 0) goto 400
*
        ix = nint( r(1,i) )
        iy = nint( r(2,i) )
        iz = nint( r(3,i) )

      if(abs(ix).gt.maxx.or.abs(iy).gt.maxx.or.abs(iz).gt.maxz) then
        nesc = nesc + 1
      else
*
        kx=nint(dble(2*ip+1)*(r(1,i)-dble(ix)))
        if(abs(kx) .eq. ip+1) kx = kx/abs(kx) * ip
        ky=nint(dble(2*ip+1)*(r(2,i)-dble(iy)))
        if(abs(ky) .eq. ip+1) ky = ky/abs(ky) * ip
        kz=nint(dble(2*ip+1)*(r(3,i)-dble(iz)))
        if(abs(kz) .eq. ip+1) kz = kz/abs(kz) * ip
        ic=1+(kz+ip)+(ky+ip)*(2*ip+1)+(kx+ip)*(2*ip+1)**2
        ib=0
        do 700 jx=ix-iq,ix+iq
        do 700 jy=iy-iq,iy+iq
        do 700 jz=iz-iq,iz+iq
        ib=ib+1
        if(cm(ic,ib).gt.0.0.and.
     &     abs(jx).le.maxx.and.abs(jy).le.maxx.and.
     &     abs(jz).le.maxz) then
*
        rhb(jx,jy,jz)=rhb(jx,jy,jz)+cm(ic,ib)
***      proton and neutron densities
         if( (id(1,i).eq.1) .and. (id(2,i).eq.1)) then
            rhob_4(4,jx,jy,jz)=rhob_4(4,jx,jy,jz)+cm(ic,ib)
         end if

         if( (id(1,i).eq.1) .and. (id(2,i).eq.0)) then
           rhob_4(5,jx,jy,jz)=rhob_4(5,jx,jy,jz)+cm(ic,ib)
         end if


!        radiale Dichteverteilung VON NUKLEONEN im Targetkern!
!       if (time.eq.0.0 .and. id(1,i).eq.1) then
!       write(86,*) r(1,i),r(2,i),r(3,i)	!write in Kpl_A*.dat

!       endif
!
!           ix = nint(xxx)
!           iy = nint(yyy)
!           iz = nint(zzz)
! 	
! 	if (zzz.le.0.0) then         !only front-hemisphere z<0
!
!       if(iabs(ix).le.maxx .and. iabs(iy).le.maxx .and. iabs(iz).le.maxz)

!
! 		if (time.eq.5.0) then
! 			do i = 0,10			!maxx
! 			do j = 0,10			!maxy = maxx
! 			do k= 0,10			!maxz
! 			deny(i,j,k) = rhb(i,j,k) / rho0
! 			enddo
! 			enddo
!                       enddo
!
!
!                    endif

*

        end if
  700   continue
*
      end if
  400 continue
  600 continue
*-----------------------------------------------------------------------
      write(*,*)'end of dens'
*
      return
*
************************************************************************
*                                                                      *
      entry denin(num)
*
*     ip ; inner mesh for dens
*     iq ; outer mesh for dens
*     ir ; outer mesh for pauli and average momentum
      ip=2
      iq=2
      ir=3
*                                                                      *
************************************************************************
*
      gauss=1.0
      dcut=5.0
*
      rd=1.0/dble(ip*2+1)
      do 101 ix=-ip,ip
      do 101 iy=-ip,ip
      do 101 iz=-ip,ip
        ic=1+(iz+ip)+(iy+ip)*(2*ip+1)+(ix+ip)*(2*ip+1)**2
        xi=dble(ix)*rd
        yi=dble(iy)*rd
        zi=dble(iz)*rd
        ib=0
        sek=0.0
        do 201 jx=-iq,iq
        do 201 jy=-iq,iq
        do 201 jz=-iq,iq
          ib=ib+1
          xj=dble(jx)
          yj=dble(jy)
          zj=dble(jz)
          rsqr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
          if(rsqr.gt.dcut) then
            cm(ic,ib)=0.0
          else
            cm(ic,ib)=exp(-rsqr/2./gauss**2)
          end if
          sek=sek+cm(ic,ib)
  201   continue
        ib=0
        do 301 jx=-iq,iq
        do 301 jy=-iq,iq
        do 301 jz=-iq,iq
          ib=ib+1
          cm(ic,ib)=cm(ic,ib)/sek/dble(num)
  301   continue
  101 continue
*
      do 121 ix=-ip,ip
      do 121 iy=-ip,ip
      do 121 iz=-ip,ip
        ic=1+(iz+ip)+(iy+ip)*(2*ip+1)+(ix+ip)*(2*ip+1)**2
        ie=0
        do 221 jx=-ir,ir
        do 221 jy=-ir,ir
        do 221 jz=-ir,ir
          ie=ie+1
          sek=0.0
          do 231 kx=jx-1,jx+1
            if(abs(kx).gt.iq) goto 231
            do 232 ky=jy-1,jy+1
              if(abs(ky).gt.iq) goto 232
              do 233 kz=jz-1,jz+1
                if(abs(kz).gt.iq) goto 233
                ib=1+(kz+iq)+(ky+iq)*(2*iq+1)+(kx+iq)*(2*iq+1)**2
                sek=sek+cm(ic,ib)
  233         continue
  232       continue
  231     continue
          dm(ic,ie)=sek
c     write(6,*) ic,ie,sek
  221   continue
  121 continue
      return
      end
