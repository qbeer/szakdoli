c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=+1, XiN-potential
c
c    Xi^0 p, Sig^+ Sig0, Sig^+ Lam coupled channel potential
c
c      [3S1]         [3D1]
c    Xp, SS, SL,   Xp, SS, SL
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_XN_Ch_Plus_3SD(r,v)
       implicit   none
       real*8     r, v(6,6)

       real*8  vc1, vc2, vc3, wc1, wc2, wc3
       real*8  vt1, vt2, vt3, wt1, wt2, wt3
       real*8  wc11,wc12,wc13, wc21,wc22,wc23, wc31,wc32,wc33
       real*8  wt11,wt12,wt13, wt21,wt22,wt23, wt31,wt32,wt33
       real*8  A(3,3)

       real*8  temp(3,3,2,2)
       integer i,j,k,l, al,bt

       real*8  prmc1(7), prmc2(7), prmc3(7)
       common /potential_parameter_Vc/ prmc1, prmc2, prmc3

       real*8  prmt1(6), prmt2(6), prmt3(6)
       common /potential_parameter_Vt/ prmt1, prmt2, prmt3

       if(r.eq.0d0) then
        write(*,*)"pot_XN: r should not be zero."
        stop
       end if

c-----------------------------------------
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c   v1 = 10*,  v2 = 8a,  v3 = 10

       vc1 = prmc1(1)*DExp(-prmc1(2)*r**2)
     &     + prmc1(3)*DExp(-prmc1(4)*r**2)
     &     - prmc1(5)*( 1d0-DExp(-prmc1(6)*r**2) )**2
     &         *(DExp(-prmc1(7)*r)/r)**2

       vc2 = prmc2(1)*DExp(-prmc2(2)*r**2)
     &     + prmc2(3)*DExp(-prmc2(4)*r**2)
     &     - prmc2(5)*( 1d0-DExp(-prmc2(6)*r**2) )**2
     &         *(DExp(-prmc2(7)*r)/r)**2

       vc3 = prmc3(1)*DExp(-prmc3(2)*r**2)
     &     + prmc3(3)*DExp(-prmc3(4)*r**2)
     &     - prmc3(5)*( 1d0-DExp(-prmc3(6)*r**2) )**2
     &         *(DExp(-prmc3(7)*r)/r)**2

c-----------------------------------------
c   Vt = (1-Gauss)*Tensor + (1-Gauss)*Tensor
c   v1 = 10*,  v2 = 8a,  v3 = 10

       vt1 =  prmt1(1)*(1d0-DExp(-prmt1(2)*r**2))**2
     &               *(1d0+3d0/(prmt1(3)*r)+3d0/(prmt1(3)*r)**2)
     &               *DExp(-prmt1(3)*r)/r
     &      + prmt1(4)*(1d0-DExp(-prmt1(5)*r**2))**2
     &               *(1d0+3d0/(prmt1(6)*r)+3d0/(prmt1(6)*r)**2)
     &               *DExp(-prmt1(6)*r)/r

       vt2 =  prmt2(1)*(1d0-DExp(-prmt2(2)*r**2))**2
     &               *(1d0+3d0/(prmt2(3)*r)+3d0/(prmt2(3)*r)**2)
     &               *DExp(-prmt2(3)*r)/r
     &      + prmt2(4)*(1d0-DExp(-prmt2(5)*r**2))**2
     &               *(1d0+3d0/(prmt2(6)*r)+3d0/(prmt2(6)*r)**2)
     &               *DExp(-prmt2(6)*r)/r

       vt3 =  prmt3(1)*(1d0-DExp(-prmt3(2)*r**2))**2
     &               *(1d0+3d0/(prmt3(3)*r)+3d0/(prmt3(3)*r)**2)
     &               *DExp(-prmt3(3)*r)/r
     &      + prmt3(4)*(1d0-DExp(-prmt3(5)*r**2))**2
     &               *(1d0+3d0/(prmt3(6)*r)+3d0/(prmt3(6)*r)**2)
     &               *DExp(-prmt3(6)*r)/r

c-----------------------------------------
c   w1 = 10*,  w2 = 10,  w3 = 8a

      wc1 = vc1
      wc2 = vc3
      wc3 = vc2

      wt1 = vt1
      wt2 = vt3
      wt3 = vt2

c-----------------------------------------
c   Here 1st=XiN and 2nd=SigSig and 3rd=SigLam
c         (Xi first)                (Sig first)
c   and   1 = 10*,  2 = 10,  3 = 8a

       A(1,1) = +Sqrt( 1d0/3d0)
       A(1,2) = +Sqrt( 1d0/3d0)
       A(1,3) = -Sqrt( 1d0/3d0)

       A(2,1) =  Sqrt( 1d0/6d0)
       A(2,2) =  Sqrt( 1d0/6d0)
       A(2,3) =  Sqrt( 4d0/6d0)

       A(3,1) = -Sqrt( 1d0/2d0)
       A(3,2) =  Sqrt( 1d0/2d0)
       A(3,3) =  0d0

       wc11=A(1,1)*A(1,1)*wc1 + A(1,2)*A(1,2)*wc2 + A(1,3)*A(1,3)*wc3
       wc12=A(1,1)*A(2,1)*wc1 + A(1,2)*A(2,2)*wc2 + A(1,3)*A(2,3)*wc3
       wc13=A(1,1)*A(3,1)*wc1 + A(1,2)*A(3,2)*wc2 + A(1,3)*A(3,3)*wc3

       wc21=wc12
       wc22=A(2,1)*A(2,1)*wc1 + A(2,2)*A(2,2)*wc2 + A(2,3)*A(2,3)*wc3
       wc23=A(2,1)*A(3,1)*wc1 + A(2,2)*A(3,2)*wc2 + A(2,3)*A(3,3)*wc3

       wc31=wc13
       wc32=wc23
       wc33=A(3,1)*A(3,1)*wc1 + A(3,2)*A(3,2)*wc2 + A(3,3)*A(3,3)*wc3

       wt11=A(1,1)*A(1,1)*wt1 + A(1,2)*A(1,2)*wt2 + A(1,3)*A(1,3)*wt3
       wt12=A(1,1)*A(2,1)*wt1 + A(1,2)*A(2,2)*wt2 + A(1,3)*A(2,3)*wt3
       wt13=A(1,1)*A(3,1)*wt1 + A(1,2)*A(3,2)*wt2 + A(1,3)*A(3,3)*wt3

       wt21=wt12
       wt22=A(2,1)*A(2,1)*wt1 + A(2,2)*A(2,2)*wt2 + A(2,3)*A(2,3)*wt3
       wt23=A(2,1)*A(3,1)*wt1 + A(2,2)*A(3,2)*wt2 + A(2,3)*A(3,3)*wt3

       wt31=wt13
       wt32=wt23
       wt33=A(3,1)*A(3,1)*wt1 + A(3,2)*A(3,2)*wt2 + A(3,3)*A(3,3)*wt3

c-----------------------------------------

      temp(1,1,1,1) = wc11
      temp(1,1,1,2) = 2d0*Sqrt(2d0)*wt11
      temp(1,1,2,1) = 2d0*Sqrt(2d0)*wt11
      temp(1,1,2,2) = wc11 - 2d0*wt11

      temp(1,2,1,1) = wc12
      temp(1,2,1,2) = 2d0*Sqrt(2d0)*wt12
      temp(1,2,2,1) = 2d0*Sqrt(2d0)*wt12
      temp(1,2,2,2) = wc12 - 2d0*wt12

      temp(1,3,1,1) = wc13
      temp(1,3,1,2) = 2d0*Sqrt(2d0)*wt13
      temp(1,3,2,1) = 2d0*Sqrt(2d0)*wt13
      temp(1,3,2,2) = wc13 - 2d0*wt13

      temp(2,1,1,1) = wc21
      temp(2,1,1,2) = 2d0*Sqrt(2d0)*wt21
      temp(2,1,2,1) = 2d0*Sqrt(2d0)*wt21
      temp(2,1,2,2) = wc21 - 2d0*wt21

      temp(2,2,1,1) = wc22
      temp(2,2,1,2) = 2d0*Sqrt(2d0)*wt22
      temp(2,2,2,1) = 2d0*Sqrt(2d0)*wt22
      temp(2,2,2,2) = wc22 - 2d0*wt22

      temp(2,3,1,1) = wc23
      temp(2,3,1,2) = 2d0*Sqrt(2d0)*wt23
      temp(2,3,2,1) = 2d0*Sqrt(2d0)*wt23
      temp(2,3,2,2) = wc23 - 2d0*wt23

      temp(3,1,1,1) = wc31
      temp(3,1,1,2) = 2d0*Sqrt(2d0)*wt31
      temp(3,1,2,1) = 2d0*Sqrt(2d0)*wt31
      temp(3,1,2,2) = wc31 - 2d0*wt31

      temp(3,2,1,1) = wc32
      temp(3,2,1,2) = 2d0*Sqrt(2d0)*wt32
      temp(3,2,2,1) = 2d0*Sqrt(2d0)*wt32
      temp(3,2,2,2) = wc32 - 2d0*wt32

      temp(3,3,1,1) = wc33
      temp(3,3,1,2) = 2d0*Sqrt(2d0)*wt33
      temp(3,3,2,1) = 2d0*Sqrt(2d0)*wt33
      temp(3,3,2,2) = wc33 - 2d0*wt33

c-----------------------------------------
c  Because matrix for isospin to chage is the unit matrix.

       k = 0
       do i=1,2
        do al=1,3
         k = k + 1
         l = 0
         do j=1,2
          do bt=1,3
          l = l + 1
           v(k,l) = temp(al,bt,i,j)
          end do
         end do
        end do
       end do

       return
       end
c--------1---------2---------3---------4---------5---------6--------7--

