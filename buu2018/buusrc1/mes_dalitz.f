************************************************************************
*
      subroutine mes_dalitz(final)
*
************************************************************************
      implicit none
      include 'common'
      include 'cominput'
      include 'com_cont_epair'
      real*8 xxx,yyy,zzz,srt,px,py,pz
      real*8 gamma,beta(3),dt0,energy,prob
      integer irun,inp,ii,i1,id1
      integer id2
      logical final

      write(*,*) 'In mes_dalitz'
      call f77flush()

*     loop over all parallel runs
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        do 800 ii  = 1,maxp-1
          i1  = ii + inp
          if(ipi(1,i1) .eq. 0)                                 goto 800

          id1 = ipi(1,i1)
          id2 = ipi(2,i1)
          px  = ppi(1,i1)
          py  = ppi(2,i1)
          pz  = ppi(3,i1)
          xxx = rpi(1,i1)
          yyy = rpi(2,i1)
          zzz = rpi(3,i1)
          srt = epi(i1)
          prob=1.0
ccc call direct decays
          if(ivmesdil.eq.1.and.(id2.eq.0 .and. id1.eq.3 .or. id1.eq.5))
     &      call vectmes_dilep(id1,srt,px,py,pz,xxx,yyy,zzz,prob,1)
c          if(id1.eq.5) write(*,*) 'mes_dalitz call omega',srt
ccc call dalitz dec.
          if( .not.( ipi(1,i1).eq.1 .or. ipi(1,i1).eq.2 .or.
     &       ipi(1,i1).eq.5))                                  goto 800

          energy  = sqrt(srt**2+ppi(1,i1)**2+ppi(2,i1)**2+ppi(3,i1)**2)
          gamma = energy/srt
          beta(1) = ppi(1,i1)/energy
          beta(2) = ppi(2,i1)/energy
          beta(3) = ppi(3,i1)/energy
c       write(*,*) ' in mes_dalitz ', i1,id1,id2, px1,srt
          dt0 = dt/gamma
          if(final) dt0 = 50.0
          if((id1.eq.1.and.id2.eq.0) .or. id1.eq.2 .or. id1.eq.5) then 
            call meson_dalitz(id1,srt,dt0,beta,gamma,xxx,yyy,zzz)
          end if

 800    continue
 1000 continue
c       write(*,*) 'end of mes_dalitz2'
       call f77flush()
      return
      end

