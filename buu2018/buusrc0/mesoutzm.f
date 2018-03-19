************************************************************************
*                                                                      *
      subroutine mesout_in
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"

      integer i, npart(6)
      real*8 nkaons(2), nphi_rhoN, nphi_rhoD, norma

c      write(*,*) 'in mesout_in', imestimpri

      write(mmestimpri,*) '# Time dependence of the no. of mesons'
      write(mmestimpri,*) '# phi,eta,rho,sigma: No. of part.'
      write(mmestimpri,*) '# kaon,phi: total prob. of exist.'
      write(mmestimpri,*) '# time is in fm/c'
      write(mmestimpri,*)
      write(mmestimpri,*) '#time          pi         eta',
     &     '         rho         sig       omega',
     &     '        kaon   kaon(pert)        phi',
     &     '   phi_rhoN    phi_rhoD'

*----------------------------------------------------------------------*
      entry mesout
*----------------------------------------------------------------------*

      do i=1, 6
         npart(i) = 0
      end do

      nkaons(1) = 0.
      nkaons(2) = 0.
      nphi_rhoN = 0.
      nphi_rhoD = 0.

      norma = 1./real(num*isubs)

c      write(*,*) 'after initialization'

      do i=1, maxppar
         if(ipi(1,i).gt.6) write(*,*) 'ipi(1,i)', ipi(1,i)
         if(ipi(1,i).gt.0) npart(ipi(1,i)) = npart(ipi(1,i)) + 1
      end do

      do i=1, maxkaon
         if(ika(1,i).gt.0) nkaons(ika(1,i)) =
     &        nkaons(ika(1,i)) + pkao(4,i)
         if(ika(1,i).eq.2) then
            if(ika(2,i).eq.8) nphi_rhoN = nphi_rhoN + pkao(4,i)
            if(ika(2,i).eq.9) nphi_rhoD = nphi_rhoD + pkao(4,i)
         end if
      end do

      write(*,*) '*************************'
      write(*,*) 'No. of mesons:'
      write(*,*) 'pi:    ', npart(1)
      write(*,*) 'eta:   ', npart(2)
      write(*,*) 'rho:   ', npart(3)
      write(*,*) 'sigma: ', npart(4)
      write(*,*) 'omega: ', npart(5)
      write(*,*) 'K+:    ', npart(6)

      write(*,*) 'Total prob. of existence of K+/-:'
      write(*,*) 'kaon:  ', nkaons(1)
      write(*,*) 'phi:   ', nkaons(2)

      write(*,*) 'Total prob. of phi from rho: ',nphi_rhoN+nphi_rhoD

      if(imestimpri.eq.1)
     &     write(mmestimpri,'(f6.1,10e12.3)') time,
     &     (npart(i)*norma,i=1,6),nkaons(1)*norma,nkaons(2)*norma,
     &     nphi_rhoN*norma,nphi_rhoD*norma

      return
      end
