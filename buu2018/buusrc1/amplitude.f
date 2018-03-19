************************************************************************
      subroutine amplread
      implicit none
c  +-----------------------------------------------------------+
c  |  self = dens*(tot(3,i),tot(4,i)) + (tot(5,i),tot(6,i))    |
c  |  sqrt(s) = 0.141 + (i-1)*0.001                            |
c  +-----------------------------------------------------------+

      real ampnr(0:8,2000),ampel(0:8,2000)
      real*8 totom(7,2000),totro(7,2000)
      COMMON/self_stored/ totom,totro
      include 'common'
      real*8 pifa,pifa2,dmasss,width,pm2,dmr2,srt,srt2,srt0
      real*8 ppion0,ppion02
      real*8 realp,gamma,rej1,ppion2,ppion,rei1,reint,gaint,gaint0,dens
      integer i,j
*----------------------------------------------------------------------*
      open(amplit,file='buuinput/totro+n_ro+nf.0.out',status='old')
      do i=1,27
        read(amplit,*)
      end do
      do i=1,1000
        read(amplit,*) (ampnr(j,i),j=0,8)
      end do
      close(amplit)

c      open(amplit,file='buuinput/totme+n_me+nf.0.out',status='old')
c      do i=1,27
c        read(amplit,*)
c      end do
c      do i=1,999
c        read(amplit,*) (ampel(j,i),j=0,8)
c      end do
c      write(*,*) 'ampel read: ',(ampel(j,999),j=0,8)
c      close(amplit)

      dens = 1.
      pifa = 2.0*pi
      pifa2=pifa*pifa
      dmasss = 0.779
      width  = 0.141
      pm2=pmass*pmass
c      dmo2=omass**2
      dmr2=dmasss**2
      do i=1,999
        srt0 = ampnr(0,i)
c        if(abs(ampnr(0,i)-ampel(0,i)).gt.1.e-3) write(*,*) 'hiba srt',
c     &  ampnr(0,i), ampel(0,i)
        srt = srt0 - rmass
        srt2 = srt**2
c        totom(1,i)=
c     &      hbc/(4.0*pi)*rmass/srt0*(0.333*ampel(1,i)+0.6667*ampel(3,i))
c        totom(2,i)=
c     &      hbc/(4.0*pi)*rmass/srt0*(0.333*ampel(2,i)+0.6667*ampel(4,i))
        totro(1,i) = hbc/(4.0*pi)*rmass/srt0*(
     &               0.1111*ampnr(1,i)+0.2222*ampnr(3,i)+
     &               0.2222*ampnr(5,i)+0.4444*ampnr(7,i))
        totro(2,i) = hbc/(4.0*pi)*rmass/srt0*(
     &               0.1111*ampnr(2,i)+0.2222*ampnr(4,i)+
     &               0.2222*ampnr(6,i)+0.4444*ampnr(8,i))
c        totom(3,i) = 4.0*pi*(1.0+srt/rmass)*totom(1,i)*dens*hbc**2
c        totom(4,i) = 4.0*pi*(1.0+srt/rmass)*totom(2,i)*dens*hbc**2
        totro(3,i) = 4.0*pi*(1.0+srt/rmass)*totro(1,i)*dens*hbc**2
        totro(4,i) = 4.0*pi*(1.0+srt/rmass)*totro(2,i)*dens*hbc**2
c        totom(5,i) = 1.0e0*(srt2-dmo2+totom(3,i))/
c     &    ((srt2-dmo2+totom(3,i))**2+(omass*gamome+totom(4,i))**2)
c        totom(6,i) = 1.0e0*(omass*gamome+totom(4,i))/
c     &    ((srt2-dmo2+totom(3,i))**2+(omass*gamome+totom(4,i))**2)
c        totom(7,i) = 1.0e0*(omass*gamome)/
c     &    ((srt2-dmo2)**2+(omass*gamome)**2)
        ppion02=0.25e0*(dmr2-4.0*pm2)
        ppion0 = sqrt(ppion02)
        realp=srt2-dmr2
        gamma=0.0
        if(srt.gt.2.0*pmass) then
          rej1=ppion0/dmasss*log((2.0*pm2-dmr2+2.0*dmasss*ppion0)/
     &                           (2.0*pm2-dmr2-2.0*dmasss*ppion0))
          ppion2=0.25e0*(srt2-4.0e0*pm2)
          if(ppion2.lt.0.0e0) ppion2=0.0e0
          ppion = sqrt(ppion2)
          rei1=ppion/srt*log((2.0e0*pm2-srt2+2.e0*srt*ppion)/
     &                       (2.0e0*pm2-srt2-2.e0*srt*ppion))
          reint=ppion2*(rei1-rej1)/pifa2
          gaint =ppion2 *ppion /srt/pifa
          gaint0=ppion02*ppion0/dmasss/pifa
          gamma = width*gaint/gaint0
          realp =(srt2-dmr2)+dmasss*width*reint/gaint0
          totro(5,i) = dmasss*width*reint/gaint0
        else
          totro(5,i) = 0.
        end if
        totro(6,i) = dmasss*gamma

c        totro(5,i) = 1.0e0*(realp+totro(3,i))/
c     &    ((realp+totro(3,i))**2+(dmasss*gamma+totro(4,i))**2)
c        totro(6,i) = 1.0e0*(dmasss*gamma+totro(4,i))/
c     &    ((realp+totro(3,i))**2+(dmasss*gamma+totro(4,i))**2)
c        totro(7,i) = 1.0e0*(dmasss*gamma)/
c     &    (realp**2+(dmasss*gamma)**2)
c        spectv(1,i) = srt
c        spectv(2,i) = totom(6,i)/pi
c        spectv(3,i) = totro(6,i)/pi
c        write(42,102) (spectv(l,i),l=1,3)
c        write(52,102) srt,(totom(l,i),l=1,7)

c          write(63,*) i, srt, srt0, gaint,gaint0,totro(6,i)
      end do

c     do i=1,1,111
c       write(63,*) 'ampnr read: ',i,(ampnr(j,i),j=0,8)
c       write(63,*) 'ampnr calc: ',i,0.141+(i-1)*0.001,(totro(j,i),j=3,6)
c     end do

c 102  format(f7.3,8(1x,e10.3))
      return
      end
