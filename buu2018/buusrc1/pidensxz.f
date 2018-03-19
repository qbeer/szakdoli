
************************************************************************
*                                                                      *
      subroutine pidensxz(tim )
*                                                                      *
************************************************************************
      implicit none
      include"common"
*----------------------------------------------------------------------*
      integer lplost,nx,nz,ii,jj
      real*8 tim
*----------------------------------------------------------------------*
      integer      idexz(-maxx:maxx,-maxz:maxz)
*----------------------------------------------------------------------*
      lplost = 0
      do 200 nx = -maxx,maxx
      do 100 nz = -maxz,maxz
        idexz(nx,nz) = 0
 100  continue
 200  continue
      do 300 ii = 1,maxppar
        if(ipi(1,ii) .ne. 1)                                    goto 300
        nx   = nint(rpi(1,ii))
        nz   = nint(rpi(3,ii))
        if((iabs(nx) .gt. maxx) .or. (iabs(nz) .gt. maxz)) then
          lplost = lplost + 1
        else
          idexz(nx,nz) = idexz(nx,nz) + 1
        end if
 300  continue
      write(isum,1000) tim , lplost
1000  format('n:',1h','xz density of pions',1h'/
     &'c: time = ',f65.2,' fm/c'/
     &'c: ',i5,'   pions are lost'/
     &'n:x (fm)'/
     &'n:z (fm)'/
     &'n2: y = -24.0 to 24.0 by 1.0;',
     &  '  x = -20.0 to 20.0 by 1.0;')
      write(isum,1100) ((idexz(ii,jj),ii=-maxx,maxx),jj=-maxz,maxz)
1100  format(49(41(i2,',')/))
      return
      end
