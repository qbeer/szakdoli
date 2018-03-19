
************************************************************************
*                                                                      *
      subroutine resdyn(jresdalitz)
*                                                                      *
*       purpose:    calculating and writing out lifetime, mass,        *
*                   denstity, collision number, free path distributions*
*                   -------------------------------------------------  *
*                                                                      *
*       variables:                                                     *
*                                                                      *
************************************************************************
      implicit none
      integer jresdalitz,iresdalitz,jj,kl,kk,id61,id62,idres
      integer ireal,ixx,iyy,izz,iendel,icq,jd,ile,idj,kdi,kt0,km,iel,kml
      integer numdt,numdm,i,j
      real*8 dlbin,tim0,dendel,px,py,pz,distan,delife,ener,gamm0,delife0
      real*8 xnum,emd,gam,bwdist
      include"common"
      save    iresdalitz
************************************************************************
*      save    dlbin
      save    numco,spac,idista
      save    delwi, idelwi
*-----------------------------------------------------------------------
      parameter( dlbin = 2.0)
      real*8   spac(3,maxpar)
      real*8   delwi(40,20) , dummyf(9)
      integer idelwi(40,20)
      integer numco(1:nres,50),idista(1:nres,50)
      iresdalitz=jresdalitz
*-----------------------------------------------------------------------
      do 437 jj =1,80
        do 436 kl =1,nres
          idelt(kl,jj) = 0
          idelmas(kl,jj) = 0
          if(jj.le.50) idedel(kl,jj) = 0
          if(jj.le.50) idedel(kl+3,jj) = 0
          if(jj.le.50) numco(kl,jj) = 0
          if(jj.le.50) idista(kl,jj) = 0
  436   continue
  437 continue
      do jj =1,40
        do kl =1,20
          delwi(jj,kl) = 0.0
         idelwi(jj,kl) = 0
        enddo
      enddo

      return
*-----------------------------------------------------------------------
      entry resprod(kk,tim0,id61,id62,idres,ireal)
*-----------------------------------------------------------------------
*
*                   ireal = 0 : r -> r
*                   ireal = 1 : n -> r
*                   ireal = 2 : n + pi -> r
*                   ireal = 3 : r -> n + eta
*                   ireal = 4 : r -> n + rho
*                   ireal = 5 : r -> n + sigma
*
*     store the time, position and the density at the creation point
*         spac(3,j) is the position of the resonance j at creation
*         dtim(j)   is the time of the creation of resonance j
*         ndedel(j)   is the density(*10) at creation of resonance j
*         idedel(1,j)   is number of res. with density j/10 at creation
*         numco(i,j)  is number of res. i with coll. number of parents j
*-----------------------------------------------------------------------
*         time
                  dtim(kk)= time + tim0
*         position
                  spac(1,kk)= r(1,kk)
                  spac(2,kk)= r(2,kk)
                  spac(3,kk)= r(3,kk)
*         density
                  ixx = nint(r(1,kk))
                  iyy = nint(r(2,kk))
                  izz = nint(r(3,kk))
                  if(abs(ixx).le.maxx.and.abs(iyy).le.maxx
     +               .and. abs(izz).le.maxz) then
                    dendel = rhb(ixx,iyy,izz)/rho0
                  else
                    dendel= 0.0
                  end if
                  iendel = nint(10.0*dendel)
                  if(iendel .lt.  1) iendel = 1
                  if(iendel .gt. 50) iendel = 50
                  idedel(idres,iendel) = idedel(idres,iendel) + 1
                  ndedel(kk) = iendel
*         collision number
                  icq        = min(iabs(id61) + iabs(id62),30)
                  numco(idres,icq) = numco(idres,icq) + 1
                  if(icq.le.1) write(*,'(''hiba resprod'',5i6,f8.2)')
     &                    id61,id62,idres,ireal,kk,time
      return
*-----------------------------------------------------------------------
      entry resabs(kk,tim0,px,py,pz,ireal)
*-----------------------------------------------------------------------
*
*                   ireal = 0 : r -> r
*                   ireal = 1 : r -> n
*                   ireal = 2 : r -> n + pi
*                   ireal = 3 : r -> n + eta
*                   ireal = 4 : r -> n + rho
*                   ireal = 5 : r -> n + sigma
*
      jd = id(1,kk)-1
      if(jd.lt.1) write(*,'(''hiba: resabs: jd lt 1'',i8)') ireal
*-----------------------------------------------------------------------
*               meson collision chain
c***         the mesonic collision number is calculated    *************
              if(ireal.eq.1 .and. id(7,kk).ge.1) then
                ile = 2*id(7,kk)-1
                idj = min(iabs(id(4,kk)) , 50)
                mlife(ile,idj) = mlife(ile,idj) + 1
                id(4,kk) = 0
              endif
c***********************************************************************
c***********      delta free path distribution   (path bin = 0.5 fm)
                  distan = sqrt( (r(1,kk)-spac(1,kk))**2 +
     &               (r(2,kk)-spac(2,kk))**2 + (r(3,kk)-spac(3,kk))**2 )
                  kdi= min( nint(2.*(distan))+1, 50)
**************st 4.4
                  if(jd.gt.0) then
                   idista(jd,kdi) = idista(jd,kdi) + 1
                  end if
******************st 4.4
*
c***********      delta lifetime distribution   (lifetime bin = 0.5 fm)
                  delife = time + tim0 - dtim(kk)
                  ener = sqrt(e(kk)**2+p(1,kk)**2+p(2,kk)**2+p(3,kk)**2)
                  gamm0= ener / e(kk)
                  delife0= delife / gamm0
                  kt0= nint(dlbin*delife0)+1
                  if(kt0 .lt. 1) kt0 = 1
                  if(kt0 .ge.80) kt0 = 80
                  idelt(jd,kt0)  = idelt(jd,kt0) + 1
                  dtim(kk)  = 0.0
*
c**********       delta mass distribution          mass bin = 0.025 GeV
                  km = nint(40.0 * (e(kk)-1.050) )
                  if (km .le. 80) idelmas(jd,km) = idelmas(jd,km) + 1
c********** density
                  ixx = nint(r(1,kk))
                  iyy = nint(r(2,kk))
                  izz = nint(r(3,kk))
                  if(abs(ixx).le.maxx.and.abs(iyy).le.maxx
     +               .and. abs(izz).le.maxz) then
                    dendel = rhb(ixx,iyy,izz)/rho0
                  else
                    dendel= 0.0
                  end if
                  iendel = nint(10.0*dendel)
                  if(iendel .lt.  1) iendel = 1
                  if(iendel .gt. 50) iendel = 50
                  idedel(jd+3,iendel) = idedel(jd+3,iendel) + 1
c********** delta width (mass,density)
                if(jd.eq.1) then
                  iel = nint(2.0*dendel)
                  if(iel .lt.  1) iel = 1
                  if(iel .gt. 20) iel = 20
                  kml = km
                  if(kml .lt.  1) kml = 1
                  if(kml .gt. 40) kml = 40
                  delwi(kml,iel) = delwi(kml,iel) + delife0
                  idelwi(kml,iel) =idelwi(kml,iel) + 1
                endif
c********** dilepton
c                  if((ideldil.eq.2) .and. (kt0.gt.1) .and.
c     &                 ((id(2,kk)+2)/2.eq.1)) then
c                    ww = 1.0
c                    if(ideldil.eq.2) call deldil(iseed,delife,kk,
c     &           id(1,kk),e(kk),r(1,kk),r(2,kk),r(3,kk),px,py,pz,ww,0,
c     &                                            ndedel(kk),iendel)
c                  end if
*
      return
*-----------------------------------------------------------------------
      entry resdyout(xnum)
*
*=======================================================================
      numdt = 0
      numdm = 0
      do 438 jj =1,80
        numdt = numdt + idelt(1,jj)
        numdm = numdm + idelmas(1,jj)
  438 continue
ccc_hw
      goto  77777
ccc_hw
      write(isum,'(/''c:delta proper lifetime distribution'')')
      write(isum,'(''c: total number of deltas:'',f8.3)') numdt*xnum
      write(isum,'(''n: time (fm/c)'')')
      write(isum,'(''n: number of deltas:dn/dt(c/fm)'')')
      write(isum,'(''n:   x      y(delta),lt0     y(n*1),dt0   '',
     &                         ''y(n*2),mt0'')')
      write(isum,'(f8.2,3f14.3)') (float(i-1)/dlbin,
     &  (dlbin*float(idelt(j,i))*xnum,j=1,3), i=1,40)
*
      write(isum,'(/''c:delta pathlength distribution'')')
      write(isum,'(''c: total number of deltas:'',f8.3)') numdt*xnum
      write(isum,'(''n: path length (fm)'')')
      write(isum,'(''n: number of deltas:dn/ddist(1/fm)'')')
      write(isum,'(''n:   x      y(delta),lt0   y(n1440)   y(1535)'')')
      write(isum,'(f8.2,3f14.3)')
     &  (float(i-1)*0.5,(float(2*idista(j,i))*xnum,j=1,3),i=1,40)
*
      write(isum,'(/''c:delta mass distribution'')')
      write(isum,'(''c: total number of deltas:'',f8.3)') numdm*xnum
      write(isum,'(''n: mass (gev)'')')
      write(isum,'(''n: number of deltas:dn/dm(1/gev)'')')
      write(isum,'(''n:   x        y(delta),l0    y(n*1),d0  '',
     &             ''   y(n*2),m0'')')
      write(isum,'(f8.3,3f14.3)')
     & (1.05+float(i)/40.,(float(40*idelmas(j,i))*xnum,j=1,3), i=1,38)
      write(isum,'(/''c:res final mass distribution'')')
      write(isum,'(''n: mass (gev)'')')
      write(isum,'(''n: number of deltas:dn/dm(1/gev)'')')
      write(isum,'(''n:   x        y(delta),l0    y(n*1),d0  '',
     &             ''   y(n*2),m0'')')
      write(isum,'(f8.3,3f14.3)')
     & (1.05+float(i)/40.,(float(40*idelmas(j,i))*xnum,j=4,6), i=1,30)
*
*
      write(isum,'(/''c:density distribution of resonance creation'')')
      write(isum,'(''c: total number of deltas:'',f8.3)') numdm*xnum
      write(isum,'(''n: density (rho/rho0)'')')
      write(isum,'(''n: number of resonances:dn/drho)'')')
      write(isum,'(''n:   x        y(delta),l0    y(n*1),d0  '',
     &             ''   y(n*2),m0'')')
      write(isum,'(f8.3,3f14.3)')
     & (float(i)*0.1,(float(40*idedel(j,i))*xnum,j=1,3), i=1,35)
      write(isum,'(/''c:density distribution of final res. creation'')')
      write(isum,'(''n: density (rho/rho0)'')')
      write(isum,'(''n: number of resonances:dn/drho)'')')
      write(isum,'(''n:   x        y(delta),l0    y(n*1),d0  '',
     &             ''   y(n*2),m0'')')
      write(isum,'(f8.3,3f14.3)')
     & (float(i)*0.1,(float(40*idedel(j,i))*xnum,j=4,6), i=1,35)
*
*********
      write(isum,'(/''c:resonance creation vs coll. number'')')
      write(isum,'(''n: coll.number '')')
      write(isum,'(''n: number of resonances'')')
      write(isum,'(''n:   x  y(delta),l0   y(n*1),m0  y(n*2),d0'')')
      write(isum,'(i6,3i10)')
     & (i,(numco(j,i),j=1,3), i=1,30)
*
*********
      do i=1,40
        do j=1,20
          if(idelwi(i,j).ne.0) delwi(i,j)= delwi(i,j)/float(idelwi(i,j))
        enddo
      enddo
      write(isum,'(/''c:delta lifetime distribution'')')
      write(isum,'(''c: total number of deltas:'',f8.3)') numdm*xnum
      write(isum,'(''n: mass (gev)'')')
      write(isum,'(''n: average width(gev)'')')
      write(isum,'(''n:  x     y(0.5),l0  y(1.0),d0    '',
     &             ''y(1.5),m0    y(2.0),p0   y(2.5),r0  y(3.0),b0'')')
      write(isum,'(f8.3,6f12.3)')
     &              (1.05+float(i)/40.,(delwi(i,j),j=1,6), i=1,40)
      do i=1,40
        do j=1,20
          if(idelwi(i,j).ne.0) delwi(i,j)= hbc/delwi(i,j)
        enddo
      enddo
      write(isum,'(/''c:delta width distribution'')')
      write(isum,'(''c: total number of deltas:'',f8.3)') numdm*xnum
      write(isum,'(''n: mass (gev)'')')
      write(isum,'(''n: average lifetime(fm/c)'')')
      write(isum,'(''n:  x      y(0.0),lt0  y(0.5),l0  y(1.0),d0    '',
     &             ''y(1.5),m0    y(2.0),p0   y(2.5),r0  y(3.0),b0'')')
      do i=1,20
        emd = 1.050 + 0.025 * float(i)
        gam =bwdist(emd**2,1,0,3,0,0,1,0,dummyf)
        write(isum,'(f8.3,7f12.3)') emd,gam,(delwi(i,j),j=1,6)
      enddo
*
*=======================================================================
*
ccc_hw  ................i
77777 continue
      return
      end
