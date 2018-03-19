************************************************************************
*                                                                      *
      subroutine gradukao(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                   gradx,grady,gradz,gradpx, gradpy, gradpz)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))                          *
*       variables:                                                     *
*         nt                  - time step  <0 : etot = effective mass  *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy of meson      (real,output)     *
*         id2                 - charge of meson      (integer, input)  *
*         gradx, grady, gradz - gradient of u        (real,output)     *
*                                                                      *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"

      integer id2
      real*8    gradx, grady, gradz, cpot
      real*8    rx, ry, rz, px, py, pz, etot, gradpx, gradpy, gradpz
      real*8    emfox, emfoy, emfoz
      integer kx, ky, kz, nt
      real*8 v00, sikn, fpi, sigmakn, gamma
      real*8  v0, vx, vy, vz, dvxdy, dvxdz, dvydx, dvydz, dvzdx, dvzdy
      real*8  rotx, roty, rotz, rho00, rhoxp, rhoxm, rhoyp, rhoym,
     1      rhozp, rhozm, rhox, rhoy, rhoz, pj, wstar,
     2      dndx, dndy, dndz, djxdt, djydt, djzdt, pp2, cx, cy, cz

      data   fpi, sigmakn    /.093, .20/   ! effective chiral potential
c----------

      gradx = 0.0
      grady = 0.0
      gradz = 0.0
      etot  = sqrt(xkmas**2 + px**2 + py**2 + pz**2)
      if (nt .gt. 0) then
        gradpx= px / etot
        gradpy= py / etot
        gradpz= pz / etot
      endif
c
      v00   = 3./(8.*fpi**2) *  hbc**3
      sikn  = sigmakn / fpi**2 * hbc**3
*----------------------------------------------------------------------*
*       call routine emfoca to calculate the electro-magnetic force    *
*       emfox, emfoy, emfoz  ( link to coulomb-routine )               *

      if (ipou .eq. 1  .and.  id2 .eq. 1 .and. nt .ge. 0) then
           call emfoca(rx, ry, rz, px, py, pz, etot, id2, emfox, emfoy,
     &                 emfoz,ncont,cpot)
           gradx = -emfox
           grady = -emfoy
           gradz = -emfoz
           write(*,*)   '  kaon - coulomb not properly installed '
          		write(*,*) "stop HS 29"
           stop
      endif
      if ( ikaonpot .gt. 0)  then
         kx  = nint(rx)
         ky  = nint(ry)
         kz  = nint(rz)
      if( abs(kx).lt.maxx .and. abs(ky).lt.maxx .and.
     &    abs(kz).lt.maxz)                             then  ! grid
         pp2 = px*px + py*py + pz*pz
         v0  =  rhob_4(0,kx,ky,kz)
         gamma = rhob_4(6,kx  ,ky  ,kz  )
         rho00= v0 / gamma
         vx  =  rhob_4(1,kx,ky,kz)
         vy  =  rhob_4(2,kx,ky,kz)
         vz  =  rhob_4(3,kx,ky,kz)
c
         pj   =  px*vx + py*vy +pz*vz
         wstar = xkmas**2 + pp2 - sikn*rho00/gamma  -2.*v00*pj
     1                    + (v00*v0)**2
         if (wstar .lt. .02) wstar = .02
         wstar = sqrt (wstar)
         etot = wstar + v00 * v0
         if (nt .lt. 0) then
             etot = etot**2-pp2
             if (etot .lt. .1) etot = .25
             etot = sqrt(etot)
             return
         endif
c
         cx  = (px - v00*vx) / wstar
         cy  = (py - v00*vy) / wstar
         cz  = (pz - v00*vz) / wstar
         dvxdy = .5*(rhob_4(1,kx  ,ky+1,kz  )-rhob_4(1,kx  ,ky-1,kz  ))
         dvxdz = .5*(rhob_4(1,kx  ,ky  ,kz+1)-rhob_4(1,kx  ,ky  ,kz-1))
         dvydx = .5*(rhob_4(2,kx+1,ky  ,kz  )-rhob_4(2,kx-1,ky  ,kz  ))
         dvydz = .5*(rhob_4(2,kx  ,ky  ,kz+1)-rhob_4(2,kx  ,ky  ,kz-1))
         dvzdx = .5*(rhob_4(3,kx+1,ky  ,kz  )-rhob_4(3,kx-1,ky  ,kz  ))
         dvzdy = .5*(rhob_4(3,kx  ,ky+1,kz  )-rhob_4(3,kx  ,ky-1,kz  ))
         rotx  = v00 * (cy * (dvydx-dvxdy) + cz * (dvzdx-dvxdz))
         roty  = v00 * (cz * (dvzdy-dvydz) + cx * (dvxdy-dvydx))
         rotz  = v00 * (cx * (dvxdz-dvzdx) + cy * (dvydz-dvzdy))
c
cc         if (nt .gt. 16 .and. nt .lt. 19)  write(*,*)  '  gradu kao 1 '
         rhoxp= rhob_4(0,kx+1,ky  ,kz  )/rhob_4(6,kx+1,ky  ,kz  )
         rhoxm= rhob_4(0,kx-1,ky  ,kz  )/rhob_4(6,kx-1,ky  ,kz  )
         rhoyp= rhob_4(0,kx  ,ky+1,kz  )/rhob_4(6,kx  ,ky+1,kz  )
         rhoym= rhob_4(0,kx  ,ky-1,kz  )/rhob_4(6,kx  ,ky-1,kz  )
         rhozp= rhob_4(0,kx  ,ky  ,kz+1)/rhob_4(6,kx  ,ky  ,kz+1)
         rhozm= rhob_4(0,kx  ,ky  ,kz-1)/rhob_4(6,kx  ,ky  ,kz-1)
         rhox = .5 * (rhoxp - rhoxm)
         rhoy = .5 * (rhoyp - rhoym)
         rhoz = .5 * (rhozp - rhozm)
c
         dndx  = .5 * (rhob_4(0,kx+1,ky,kz) - rhob_4(0,kx-1,ky,kz))
         dndy  = .5 * (rhob_4(0,kx,ky+1,kz) - rhob_4(0,kx,ky-1,kz))
         dndz  = .5 * (rhob_4(0,kx,ky,kz+1) - rhob_4(0,kx,ky,kz-1))
c
         djxdt = (rhob_4(1,kx,ky,kz) - rhob_4(7,kx,ky,kz) )/ dt
         djydt = (rhob_4(2,kx,ky,kz) - rhob_4(8,kx,ky,kz) )/ dt
         djzdt = (rhob_4(3,kx,ky,kz) - rhob_4(9,kx,ky,kz) )/ dt
c
         gradx = - rotx + v00**2*rho00/wstar*rhox
     1           - 0.5*sikn/wstar/gamma * rhox + v00 * (dndx + djxdt)
         grady = - roty + v00**2*rho00/wstar*rhoy
     1           - 0.5*sikn/wstar/gamma * rhoy + v00 * (dndy + djydt)
         gradz = - rotz + v00**2*rho00/wstar*rhoz
     1           - 0.5*sikn/wstar/gamma * rhoz + v00 * (dndz + djzdt)
c
         gradpx  = cx
         gradpy  = cy
         gradpz  = cz
      endif               !              grid range
      endif               !              pot kaon
c
      return
      end
************************************************************************
*                                                                      *
      subroutine gradukaon(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                   gradx,grady,gradz,gradpx, gradpy, gradpz)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))    with  a simple
*                        U = n(a+becp(-cp)) ansatz                     *
*       variables:                                                     *
*         nt                  - time step  <0 : etot = effective mass  *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy o meson       (real, input)     *
*         id2                 - charge of meson      (integer, input)  *
*         gradx, grady, gradz - gradient of u        (real,output)     *
*                                                                      *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      common /nthwhw/  nthw
      integer nthw
      common /counthw/ ihw
      integer ihw

      integer id2, ncont0
      real*8    gradx, grady, gradz, cpot
      real*8    rx, ry, rz, px, py, pz, etot, gradpx, gradpy, gradpz
      real*8    emfox, emfoy, emfoz
      integer kx, ky, kz, nt
      real*8    dndx, dndy, dndz, rho_cr, aa
      real*8  pp, cx, cy, cz, dudp, fac, u0, p0, gamma, v00
      real*8  v(0:3), ppart(0:3), pmatt(0:3), pderiv(0:3)
      real*8  akpl, bkpl, ckpl, ak,bk,ck
      real*8  akpl1, bkpl1, ckpl1
      real*8  akpl2, bkpl2, ckpl2
      real*8  akpl3, bkpl3, ckpl3
      real*8  akpl4, bkpl4, ckpl4
      real*8  akpl5, bkpl5, ckpl5

      data   akpl1, bkpl1, ckpl1 /+.375,  .0, .0/
      data   akpl2, bkpl2, ckpl2 /+.147,  .0, 0./
      data   akpl4, bkpl4, ckpl4 /+.263,  .0, 0./
      data   akpl5, bkpl5, ckpl5 /+.0,  .0, 0./
      data   akpl3, bkpl3, ckpl3 /+.0,  .300, 4.0/
c----------
      etot  = xkmas
      if (id2 .gt. 0  .and.  ikaonpot    .eq. 0) goto 612
cc
      if (ikaonpot .eq. 1) then
         akpl = akpl1
         bkpl = bkpl1
         ckpl = ckpl1
      endif
      if (ikaonpot .eq. 5) then
         akpl = akpl5
         bkpl = bkpl5
         ckpl = ckpl5
      endif
      if (ikaonpot .eq. 2) then
         akpl = akpl2
         bkpl = bkpl2
         ckpl = ckpl2
      endif
      if (ikaonpot .eq. 3) then
         akpl = akpl3
         bkpl = bkpl3
         ckpl = ckpl3
      endif
      if (ikaonpot .eq. 4) then
         akpl = akpl4
         bkpl = bkpl4
         ckpl = ckpl4
      endif
c      rho_cr = 0.3
!       rho_cr = 0.17
      rho_cr = 0.16
      gradx = 0.0
      grady = 0.0
      gradz = 0.0
      p0  =   sqrt(xkmas**2 + px**2 + py**2 + pz**2)
      gradpx= px / p0
      gradpy= py / p0
      gradpz= pz / p0
c
      kx  = nint(rx)
      ky  = nint(ry)
      kz  = nint(rz)
      if( abs(kx).lt.maxx .and. abs(ky).lt.maxx .and.
     &    abs(kz).lt.maxz)                             then  ! grid
c        if (nthw  .ge. 15) write(*,*) ' check kao2 ',kx,ky,kz,rx,ry,rz
         v(0)  =  rhob_4(0,kx,ky,kz)
         if ( v(0) .lt. .01)           goto 612
         v(1)  = -rhob_4(1,kx,ky,kz)
         v(2)  = -rhob_4(2,kx,ky,kz)
         v(3)  = -rhob_4(3,kx,ky,kz)
         gamma =  rhob_4(6,kx,ky,kz)
         v00   = v(0) / gamma
c
        if (id2 .gt. 0)  then
           ak = akpl
           bk = bkpl
           ck = ckpl
        endif
c      if (ihw .eq. 3197)
c    1      write(*,*)  ' gradukao2 ',kx,ky,kz, v
         aa = .5 / gamma
         dndx  = aa * (rhob_4(0,kx+1,ky,kz) - rhob_4(0,kx-1,ky,kz))
         dndy  = aa * (rhob_4(0,kx,ky+1,kz) - rhob_4(0,kx,ky-1,kz))
         dndz  = aa * (rhob_4(0,kx,ky,kz+1) - rhob_4(0,kx,ky,kz-1))
c
         pp = 1.0
         if(abs(bk) .le. .01) then
           ppart(0) = p0
           ppart(1) = px
           ppart(2) = py
           ppart(3) = pz
           call lorentz_hw(v, ppart, pmatt)
c          if (nt .eq. -1)  write(*,*) ' kao2 ', v, ppart, pmatt
           pp = sqrt(pmatt(1)**2 + pmatt(2)**2 + pmatt(3)**2)
           pp = max(pp, .001)
           pmatt(0) = .0
           call lorentz_hw(v, pmatt, pderiv)
           dudp = -ck*bk*v00 * exp(-ck*pp) / pp
c          if (ihw .eq. 3197 ) write(*,*) v, pmatt, pderiv
           cx = dudp * (pderiv(1) + pderiv(0) * px/p0)
           cy = dudp * (pderiv(2) + pderiv(0) * py/p0)
           cz = dudp * (pderiv(3) + pderiv(0) * pz/p0)
           u0 = (ak + bk*exp(-ck*pp) )*v00
         else
           cx = .0
           cy = .0
           cz = .0
           u0 = ak * v00
         endif
         if (abs(u0) .gt. rho_cr)  then
            aa = rho_cr / abs(u0)
            u0 = aa * u0
            cx = aa * cx
            cy = aa * cy
            cz = aa * cz
         endif
c              	write(isum,*) 'henny4'
!          write(20,*) "pot kplus",u0
         etot = xkmas + u0
         if (nt .lt.  0)  return
c
c        p0     = sqrt(px**2+py**2+pz**2+etot**2)
c        fac    = etot / p0
         gradpx = px / p0 + cx
         gradpy = py / p0 + cy
         gradpz = pz / p0 + cz
c
         gradx  = u0 * dndx
         grady  = u0 * dndy
         gradz  = u0 * dndz
      endif               !              grid range
c
  612 continue
      if(ipou .eq. 1    .and.  abs(id2) .eq. 1) then
c       ncont0 = ncont
c       ncont = 777
        call emfoca(rx, ry, rz, px, py, pz, etot, id2,
     &             emfox, emfoy,emfoz,ncont,cpot)
c       ncont = ncont0
        gradx = gradx  - emfox
        grady = grady  - emfoy
        gradz = gradz  - emfoz
      end if
c       if (rx**2+ry**2 .lt. 2.)
c    1   write(*,*)  ' in gradukao2 - coul ', id2, rz, gradz, emfoz
c
      return
      end

************************************************************************
*                                                                      *
      subroutine gradukao2(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                   gradx,grady,gradz,gradpx, gradpy, gradpz)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))    with  a simple
*                        U = n(a+becp(-cp)) ansatz                     *
*       variables:                                                     *
*         nt                  - time step  <0 : etot = effective mass  *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy o meson       (real, input)     *
*         id2                 - charge of meson      (integer, input)  *
*         gradx, grady, gradz - gradient of u        (real,output)     *
*                                                                      *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      common /nthwhw/  nthw
      integer nthw
      common /counthw/ ihw
      integer ihw

      integer id2, ncont0
      real*8    gradx, grady, gradz, cpot
      real*8    rx, ry, rz, px, py, pz, etot, gradpx, gradpy, gradpz
      real*8    emfox, emfoy, emfoz
      integer kx, ky, kz, nt
      real*8    dndx, dndy, dndz, rho_cr, aa
      real*8  pp, cx, cy, cz, dudp, fac, u0, p0, gamma, v00
      real*8  v(0:3), ppart(0:3), pmatt(0:3), pderiv(0:3)
      real*8  akminu, bkminu, ckminu, ak,bk,ck
      real*8  akminu1, bkminu1, ckminu1
      real*8  akminu2, bkminu2, ckminu2
      real*8  akminu3, bkminu3, ckminu3
      real*8  akminu4, bkminu4, ckminu4
      real*8  akminu5, bkminu5, ckminu5
      real*8  akminu6, bkminu6, ckminu6
      real*8  akminu7, bkminu7, ckminu7

      data   akminu1, bkminu1, ckminu1 /-.59, .0, .0/
      data   akminu6, bkminu6, ckminu6 /-.61, .0, .0/
      data   akminu2, bkminu2, ckminu2 /-.47, .0, .0/
      data   akminu3, bkminu3, ckminu3 /-.341, -.823, 2.5/!Sibirtsev-Cassing
      data   akminu7, bkminu7, ckminu7 /-.60, -.823, 2.5/
      data   akminu4, bkminu4, ckminu4 /+.310, -1.27, 2.8/
      data   akminu5, bkminu5, ckminu5 /-.0, -.000, 0.0/
c----------
c
      write(*,*) 'in Gradukao2',nt, rx, ry, rz, px, py, pz, etot, id2,
     &     gradx,grady,gradz,gradpx, gradpy, gradpz, i_kminu_pot
      etot  = xkmas
      if (id2 .le. 0  .and.  i_kminu_pot .eq. 0) goto 600
cc
      if (i_kminu_pot .eq. 1) then
         akminu = akminu1
         bkminu = bkminu1
         ckminu = ckminu1
      endif
      if (i_kminu_pot .eq. 2) then
         akminu = akminu2
         bkminu = bkminu2
         ckminu = ckminu2
      endif
      if (i_kminu_pot .eq. 3) then
         akminu = akminu3
         bkminu = bkminu3
         ckminu = ckminu3
      endif
      if (i_kminu_pot .eq. 4) then
         akminu = akminu4
         bkminu = bkminu4
         ckminu = ckminu4
      endif
      if (i_kminu_pot .eq. 5) then
         akminu = akminu5
         bkminu = bkminu5
         ckminu = ckminu5
      endif
      if (i_kminu_pot .eq. 6) then
         akminu = akminu6
         bkminu = bkminu6
         ckminu = ckminu6
      endif
      if (i_kminu_pot .eq. 7) then
         akminu = akminu7
         bkminu = bkminu7
         ckminu = ckminu7
      endif
!       rho_cr = 0.3
      rho_cr = 0.16
      gradx = 0.0
      grady = 0.0
      gradz = 0.0
      p0  =   sqrt(xkmas**2 + px**2 + py**2 + pz**2)
      gradpx= px / p0
      gradpy= py / p0
      gradpz= pz / p0
c
      kx  = nint(rx)
      ky  = nint(ry)
      kz  = nint(rz)
      if( abs(kx).lt.maxx .and. abs(ky).lt.maxx .and.
     &    abs(kz).lt.maxz)                             then  ! grid
c        if (nthw  .ge. 15) write(*,*) ' check kao2 ',kx,ky,kz,rx,ry,rz
         v(0)  =  rhob_4(0,kx,ky,kz)
         if ( v(0) .lt. .01)           goto 600
         v(1)  = -rhob_4(1,kx,ky,kz)
         v(2)  = -rhob_4(2,kx,ky,kz)
         v(3)  = -rhob_4(3,kx,ky,kz)
         gamma =  rhob_4(6,kx,ky,kz)
         v00   = v(0) / gamma
c
        if (id2 .le. 0)  then
           ak = akminu
           bk = bkminu
           ck = ckminu
        endif
c      if (ihw .eq. 3197)
c    1      write(*,*)  ' gradukao2 ',kx,ky,kz, v
         aa = .5 / gamma
         dndx  = aa * (rhob_4(0,kx+1,ky,kz) - rhob_4(0,kx-1,ky,kz))
         dndy  = aa * (rhob_4(0,kx,ky+1,kz) - rhob_4(0,kx,ky-1,kz))
         dndz  = aa * (rhob_4(0,kx,ky,kz+1) - rhob_4(0,kx,ky,kz-1))
c
c		write(isum,*) 'henny3'
         pp = 1.0
         if(abs(bk) .le. .01) then
           ppart(0) = p0
           ppart(1) = px
           ppart(2) = py
           ppart(3) = pz
           call lorentz_hw(v, ppart, pmatt)
c          if (nt .eq. -1)  write(*,*) ' kao2 ', v, ppart, pmatt
           pp = sqrt(pmatt(1)**2 + pmatt(2)**2 + pmatt(3)**2)
           pp = max(pp, .001)
           pmatt(0) = .0
           call lorentz_hw(v, pmatt, pderiv)
           dudp = -ck*bk*v00 * exp(-ck*pp) / pp
c          if (ihw .eq. 3197 ) write(*,*) v, pmatt, pderiv
           cx = dudp * (pderiv(1) + pderiv(0) * px/p0)
           cy = dudp * (pderiv(2) + pderiv(0) * py/p0)
           cz = dudp * (pderiv(3) + pderiv(0) * pz/p0)
           u0 = (ak + bk*exp(-ck*pp) )*v00
         else
           cx = .0
           cy = .0
           cz = .0
           u0 = ak * v00
         endif
         if (abs(u0) .gt. rho_cr)  then
            aa = rho_cr / abs(u0)
            u0 = aa * u0
            cx = aa * cx
            cy = aa * cy
            cz = aa * cz
         endif
c              	write(isum,*) 'henny4'
         etot = xkmas + u0
!         write(20,*) "pot kminu",u0
         if (nt .lt.  0)  return
c
c        p0     = sqrt(px**2+py**2+pz**2+etot**2)
c        fac    = etot / p0
         gradpx = px / p0 + cx
         gradpy = py / p0 + cy
         gradpz = pz / p0 + cz
c
         gradx  = u0 * dndx
         grady  = u0 * dndy
         gradz  = u0 * dndz
      endif               !              grid range
c
  600 continue
      if(ipou .eq. 1    .and.  abs(id2) .eq. 1) then
c       ncont0 = ncont
c       ncont = 777
        call emfoca(rx, ry, rz, px, py, pz, etot, id2,
     &             emfox, emfoy,emfoz,ncont,cpot)
c       ncont = ncont0
        gradx = gradx  - emfox
        grady = grady  - emfoy
        gradz = gradz  - emfoz
      end if
c       if (rx**2+ry**2 .lt. 2.)
c    1   write(*,*)  ' in gradukao2 - coul ', id2, rz, gradz, emfoz
c
      return
      end
************************************************************************
*                                                                      *
      subroutine graduphi(nt, rx, ry, rz, px, py, pz, etot,
     &                   gradx,grady,gradz,gradpx, gradpy, gradpz)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))                          *
*              potantial: c * phimass * density/density_0              *
*       variables:                                                     *
*         nt                  - time step  <0 : etot = effective mass  *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy of meson      (real, output)    *
*         gradx, grady, gradz - gradient of u        (real, output)    *
*                                                                      *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"

      real*8  gradx, grady, gradz
      real*8  rx,ry,rz, p0,px,py,pz, etot, gradpx, gradpy, gradpz
      integer kx, ky, kz, nt
      real*8  v00, dndx, dndy, dndz, density, gamma, c_hatsuda
      real*8  c_2, c_3,c_4

      data   c_hatsuda   /.02175/   !  hatsuda kuwabara'95  more exact
c      data   c_hatsuda   /.0343/   !  hatsuda kuwabara'95  more exact  !!!! test_: -35MeV
      data   c_2   /.0098/   !  QCD sum rules min
      data   c_3   /.0330/   !  QCD sum rules max
c----------
c      write(*,*) 'graduphi',nt, rx, ry, rz, px, py, pz, etot,
c     &                   gradx,grady,gradz,gradpx, gradpy, gradpz
      call f77flush()
      gradx = 0.0
      grady = 0.0
      gradz = 0.0
      etot = xphimas
      p0  = sqrt(xphimas**2 + px**2 + py**2 + pz**2)
      if (nt .gt. 0) then
        gradpx= px / p0
        gradpy= py / p0
        gradpz= pz / p0
      endif
c
      if( iphi_pot .eq. 1)  then
        v00   = -c_hatsuda * xphimas / rho0
      elseif( iphi_pot .eq. 2)  then
        v00   = -c_2 * xphimas / rho0
      elseif ( iphi_pot .eq. 3)  then
        v00   = -c_3 * xphimas / rho0
      endif
c      write(*,*) 'graduphi 2',v00,etot,p0,gradpx,gradpy,gradpz
      call f77flush()
*----------------------------------------------------------------------*
      if ( iphi_pot .gt. 0)  then
        kx  = nint(rx)
        ky  = nint(ry)
        kz  = nint(rz)
        if( abs(kx).lt.maxx-1 .and. abs(ky).lt.maxx-1 .and.
     &    abs(kz).lt.maxz-1)                             then  ! grid
          density = rhob_4(0,kx,ky,kz)
          gamma =  rhob_4(6,kx,ky,kz)
          density = density / gamma
          etot = xphimas + v00 * density
          if (etot .lt. .5) etot = .5
          if (nt .lt. 0) then
            return
          endif
c
          dndx  = .5 * (rhob_4(0,kx+1,ky,kz) - rhob_4(0,kx-1,ky,kz))
          dndy  = .5 * (rhob_4(0,kx,ky+1,kz) - rhob_4(0,kx,ky-1,kz))
          dndz  = .5 * (rhob_4(0,kx,ky,kz+1) - rhob_4(0,kx,ky,kz-1))
c
          gradx =   v00 * dndx  / gamma
          grady =   v00 * dndy  / gamma
          gradz =   v00 * dndz  / gamma
c
          gradpx= px / p0
          gradpy= py / p0
          gradpz= pz / p0
        endif               !              grid range
      endif               !              pot kaon
c
      return
      end

************************************************************************
      subroutine graduxi(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                   gradx,grady,gradz,gradpx, gradpy, gradpz)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))                          *
*       variables:                                                     *
*         nt                  - time step  <0 : etot = effective mass  *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy o meson       (real, input)     *
*         id2                 - charge of meson      (integer, input)  *
*         gradx, grady, gradz - gradient of u        (real,output)     *
*                                                                      *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"

      integer id2
      real*8    gradx, grady, gradz
      real*8    rx,ry,rz, p0,px,py,pz, etot, gradpx, gradpy, gradpz
      integer kx, ky, kz, nt
      real*8 v00, dndx, dndy, dndz, density, gamma, c_hatsuda
      real*8 c_2, c_3,c_4,ksimas

      parameter (ksimas = 1.321)

      data   c_hatsuda   /.0160/   !  hatsuda kuwabara'95  more exact for hyperons    = -20.5 MeV
!       data   c_2   /.0098/   !  QCD sum rules min
!       data   c_3   /.0330/   !  QCD sum rules max
      data   c_2   /.0000/   !  = 0 MeV
      data   c_3   /.0454/   !  QCD sum rules max
      data   c_4   /.0757/   !  QCD sum rules max
c----------

      gradx = 0.0
      grady = 0.0
      gradz = 0.0
      etot = ksimas
      p0  = sqrt(ksimas**2 + px**2 + py**2 + pz**2)
      if (nt .gt. 0) then
        gradpx= px / p0
        gradpy= py / p0
        gradpz= pz / p0
      endif
c
      if( i_ksi_pot .eq. 1)  then
        v00   = -c_hatsuda * ksimas / 0.16      !     n_0   =  0.16
      elseif( i_ksi_pot .eq. 2)  then
        v00   = -c_2 * ksimas / 0.16      !     n_0   =  0.16
      elseif( i_ksi_pot .eq. 3)  then
        v00   = -c_3 * ksimas / 0.16      !     n_0   =  0.16
      elseif( i_ksi_pot .eq. 4)  then
        v00   = -c_4 * ksimas / 0.16      !     n_0   =  0.16
      endif
*----------------------------------------------------------------------*
      if ( i_ksi_pot .gt. 0)  then
         kx  = nint(rx)
         ky  = nint(ry)
         kz  = nint(rz)
      if( abs(kx).lt.maxx .and. abs(ky).lt.maxx .and.
     &    abs(kz).lt.maxz)                             then  ! grid
         density = rhob_4(0,kx,ky,kz)
         gamma =  rhob_4(6,kx,ky,kz)
         density = density / gamma
         etot = ksimas + v00 * density
         if (etot .lt. .5) etot = .5
         if (nt .lt. 0) then
             return
         endif
c
         dndx  = .5 * (rhob_4(0,kx+1,ky,kz) - rhob_4(0,kx-1,ky,kz))
         dndy  = .5 * (rhob_4(0,kx,ky+1,kz) - rhob_4(0,kx,ky-1,kz))
         dndz  = .5 * (rhob_4(0,kx,ky,kz+1) - rhob_4(0,kx,ky,kz-1))
c
         gradx =   v00 * dndx  / gamma
         grady =   v00 * dndy  / gamma
         gradz =   v00 * dndz  / gamma
c
         gradpx= px / p0
         gradpy= py / p0
         gradpz= pz / p0
      endif               !              grid range
      endif               !              pot xi
c
      return
      end

************************************************************************
*                                                                      *
      subroutine graduomega(nt, rx, ry, rz, px, py, pz, etot, id2,
     &                      gradx,grady,gradz,gradpx, gradpy, gradpz)
*                                                                      *
*       purpose:     determine grad(u(x,y,z))                          *
*       variables:                                                     *
*         nt                  - time step  <0 : etot = effective mass  *
*         rx, ry, rz          - coordinates of grid  (real, input)     *
*         px, py, pz          - momentum of meson    (real, input)     *
*         etot                - energy o meson       (real, input)     *
*         id2                 - charge of meson      (integer, input)  *
*         gradx, grady, gradz - gradient of u        (real,output)     *
*                                                                      *
************************************************************************
*
      implicit none
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"

      integer id2
      real*8    gradx, grady, gradz
      real*8    rx,ry,rz, p0,px,py,pz, etot, gradpx, gradpy, gradpz
      integer kx, ky, kz, nt
      real*8 v00, dndx, dndy, dndz, density, gamma, c_2
      data   c_2   /.160/   !  QCD sum rules min

c----------

      gradx = 0.0
      grady = 0.0
      gradz = 0.0
      etot = omass
      p0  = sqrt(omass**2 + px**2 + py**2 + pz**2)
      if (nt .gt. 0) then
        gradpx= px / p0
        gradpy= py / p0
        gradpz= pz / p0
      endif
c
      if( iphi_pot .eq. 1)  then
        v00   = -c_2 * omass / 0.16      !     n_0   =  0.16
      endif
*----------------------------------------------------------------------*
      if ( iphi_pot .gt. 0)  then
         kx  = nint(rx)
         ky  = nint(ry)
         kz  = nint(rz)
        if( abs(kx).lt.maxx .and. abs(ky).lt.maxx .and.
     &    abs(kz).lt.maxz)                             then  ! grid
          density = rhob_4(0,kx,ky,kz)
          gamma =  rhob_4(6,kx,ky,kz)
          density = density / gamma
          etot = omass + v00 * density
          if (etot .lt. .5) etot = .5
          if (nt .lt. 0) then
            return
          endif
c
          dndx  = .5 * (rhob_4(0,kx+1,ky,kz) - rhob_4(0,kx-1,ky,kz))
          dndy  = .5 * (rhob_4(0,kx,ky+1,kz) - rhob_4(0,kx,ky-1,kz))
          dndz  = .5 * (rhob_4(0,kx,ky,kz+1) - rhob_4(0,kx,ky,kz-1))
c
          gradx =   v00 * dndx  / gamma
          grady =   v00 * dndy  / gamma
          gradz =   v00 * dndz  / gamma
c
          gradpx= px / p0
          gradpy= py / p0
          gradpz= pz / p0
        endif               !              grid range
      endif               !              pot kaon
c
      return
      end
