
************************************************************************
*                                                                      *
      subroutine ptest(num,iseed,io,radta)
*                                                                      *
*       variables:                                                     *
*          num      - number of test particle/nucleon  (integer, input)*
*          iseed    - seed of random number            (integer, input)*
*          io       - output unit                      (integer, input)*
*          radta    - radius of target                    (real, input)*
*                                                                      *
************************************************************************
      implicit none
      integer num,iseed,io,ii,jj,irph,ntag
      real*8 radta,dd,rsq,phase,ftot,fins,fout
      include"common"
      real*8 fd(0:100,2)
*
      dd=0.04
      do 100 ii=0,100
      fd(ii,1)=0.0
      fd(ii,2)=0.0
  100 continue
*
      do 200 jj=1,num
         do 300 ii=1+(jj-1)*maxb,jj*maxb
         if(id(1,ii).eq.0)   goto 300
         rsq=sqrt(r(1,ii)**2+r(2,ii)**2+r(3,ii)**2)
         call pauli(ii,ntag,iseed,phase,r(1,ii),r(2,ii),r(3,ii),
     &                                 p(1,ii),p(2,ii),p(3,ii))
         irph=nint(phase/dd)
         if(irph.lt.100) then
         if(rsq.lt.radta-1.0) fd(irph,1)=fd(irph,1)+1.0
         if(rsq.ge.radta-1.0) fd(irph,2)=fd(irph,2)+1.0
         end if
  300 continue
  200 continue
      write(io,600)
  600 format(/'c:plot of phase factor distribution'/
     &       'n: f(x,p)'/
     &       'n: count'/
     &       'n: x   y(tot),lh0  y(r<radius-1),dbh0 y(r>radius-1),mgh0')
      do 400 ii=0,80
      ftot=fd(ii,1)+fd(ii,2)
      fins=fd(ii,1)
      fout=fd(ii,2)
      if(ii.eq.0) then
      write(io,'(7e13.5)') float(ii)*dd,ftot,fins,fout
      else
      write(io,'(7e13.5)') float(ii)*dd-dd/2.,ftot,fins,fout
      end if
  400 continue
      return
      end
