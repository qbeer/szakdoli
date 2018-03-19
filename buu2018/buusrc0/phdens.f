
************************************************************************
*                                                                      *
      subroutine phdens(minnum,maxnum,out,fac)
*                                                                      *
*       purpose:   determine phase-space-density in z-pz-plane         *
*       variables:                                                     *
*         minnum - number of first pseudoparticle      (integer,input) *
*         maxnum - number of last pseudoparticle       (integer,input) *
*         num    - number of pseudoparticles/nucleon   (integer,input) *
*         out    - number of pseudopart. out of range (integer,output) *
*         fac    - factor                                 (real,input) *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      integer out,i1,i2,ix,iy,iz,kx,ky,kz,ic,ib,jx,jy,jz
      integer minnum,maxnum,i
      real*8 fac
*
      do 200 i2 = -24,24
        do 100 i1 = -maxz,maxz
          phrho(i1,i2) = 0.0
  100   continue
  200 continue
*
      do 300 i = minnum,maxnum
        ix = nint( r(1,i) )
        iy = nint( r(2,i) )
        iz = nint( r(3,i) )
        i2 = nint( fac * p(3,i) )
        if( abs(iz).gt.maxz .or. abs(i2).gt.24) then
          out = out + 1
        else
*
        kx=nint(float(2*ip+1)*(r(1,i)-float(ix)))
        if(abs(kx) .eq. ip+1) kx = kx/abs(kx) * ip
        ky=nint(float(2*ip+1)*(r(2,i)-float(iy)))
        if(abs(ky) .eq. ip+1) ky = ky/abs(ky) * ip
        kz=nint(float(2*ip+1)*(r(3,i)-float(iz)))
        if(abs(kz) .eq. ip+1) kz = kz/abs(kz) * ip
        ic=1+(kz+ip)+(ky+ip)*(2*ip+1)+(kx+ip)*(2*ip+1)**2
        ib=0
        do 710 jx=-iq,iq
        do 710 jy=-iq,iq
        do 710 jz=iz-iq,iz+iq
        ib=ib+1
*
        if(cm(ic,ib).gt.0.0.and.abs(jz).le.maxz)
     &  phrho(jz,i2)=phrho(jz,i2)+cm(ic,ib)
*
  710   continue
*
        end if
  300 continue
*
      return
      end
