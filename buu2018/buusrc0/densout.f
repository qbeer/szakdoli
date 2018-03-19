
************************************************************************
*                                                                      *
      subroutine densir(tim,iflag)
*                                                                      *
************************************************************************
      implicit none
      integer maxtim,iflag,itid,nx,nz,i,j,l,ny
      real*8 tim
      parameter (maxtim = 50)
      include"common"
*----------------------------------------------------------------------*
      save densxy, densxz, densyz, itid, timeide
      real*8 densxy(0:maxtim,-maxx:maxx,-maxx:maxx)
      real*8 densxz(0:maxtim,-maxx:maxx,-maxz:maxz)
      real*8 densyz(0:maxtim,-maxx:maxx,-maxz:maxz)
      real*8 timeide(maxtim)
*----------------------------------------------------------------------*
      if(iflag .eq.0) itid = 0
      itid=min0(itid+1,maxtim)
      timeide(itid) = tim
        do nx = - maxx , maxx
        do nz = - maxz , maxz
          densxz(itid,nx,nz) = rhb(nx,0,nz)
          densyz(itid,nx,nz) = rhb(0,nx,nz)
        enddo
        do ny = - maxx , maxx
          densxy(itid,nx,ny) = rhb(nx,ny,0)
        enddo
        enddo
      return
************************************************************************
*                                                                      *
      entry densout
*                                                                      *
************************************************************************
*
      do i=1,itid
        write(mdenpri,'(''c:time: '',f6.2)') timeide(i)
        write(mdenpri,'(''n:z''/''n:x'')')
        write(mdenpri,2060) -maxx,maxx,-maxz,maxz
2060    format(/'n:',1h','density profile in xz plane',1h'/
     &     'n2: y = ',i4,' to ',i4,' by 1;',
     &         'x = ',i4,' to ',i4,' by 1;')
        do j=-maxx,maxx
          write(mdenpri,2160) (densxz(i,j,l),l=-maxz,maxz)
2160      format(7(7(e10.3,',')/))
        enddo
        write(mdenpri,'(''n:z''/''n:y'')')
        write(mdenpri,2061) -maxx,maxx,-maxz,maxz
2061    format(/'n:',1h','density profile in yz plane',1h'/
     &     'n2: y = ',i4,' to ',i4,' by 1;',
     &         'x = ',i4,' to ',i4,' by 1;')
        do j=-maxx,maxx
          write(mdenpri,2160) (densyz(i,j,l),l=-maxz,maxz)
        enddo
        write(mdenpri,'(''n:x''/''n:y'')')
        write(mdenpri,2062) -maxx,maxx,-maxx,maxx
2062    format(/'n:',1h','density profile in xy plane',1h'/
     &     'n2: y = ',i4,' to ',i4,' by 1;',
     &         'x = ',i4,' to ',i4,' by 1;')
        do j=-maxx,maxx
          write(mdenpri,2163) (densxy(i,l,j),l=-maxx,maxx)
2163      format(6(7(e10.3,',')/))
        enddo
      enddo
*
      return
      end


