c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=0, XiN-potential
c
c    Xi^0 n, Xi^- p, Sig^+ Sig^-, Sig^0 Lam coupled channel potential
c
c      [3S1]          [3D1]
c    Xn,Xp,SS,SL,   Xn,Xp,SS,SL
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_XN_Ch_Zero_3SD(r,v)
       implicit   none
       real*8     r, v(8,8)

       real*8  vc1, vc2, vc3, wc1, wc2, wc3
       real*8  vt1, vt2, vt3, wt1, wt2, wt3
       real*8  A(3,3)
       real*8  wc11,wc12,wc13, wc21,wc22,wc23, wc31,wc32,wc33
       real*8  wt11,wt12,wt13, wt21,wt22,wt23, wt31,wt32,wt33
       real*8  tempc(4,4), tempt(4,4)
       real*8  C(4,4), vcij, vtij, vc(4,4), vt(4,4)
       real*8  temp(4,4,2,2)

       integer i,j,k,l, al,bt

       real*8  prmc1(7), prmc2(7), prmc3(7)
       common /potential_parameter_Vc/ prmc1, prmc2, prmc3

       real*8  prmt1(6), prmt2(6), prmt3(6)
       common /potential_parameter_Vt/ prmt1, prmt2, prmt3

       if(r.eq.0d0) then
        write(*,*)"pot_XN: r should not be zero."
        stop
       end if

c--------1---------2---------3---------4---------5---------6--------7--
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


c--------1---------2---------3---------4---------5---------6--------7--
c   For Isospin = 0

c   1st=XiN
c   vc = vc2,  vt = vt2

c--------1---------2---------3---------4---------5---------6--------7--
c   For Isospin = 1

c   w1 = 10*,  w2 = 10,  w3 = 8a

       wc1 = vc1
       wc2 = vc3
       wc3 = vc2

       wt1 = vt1
       wt2 = vt3
       wt3 = vt2

c   1st=XiN and 2nd=SigSig and 3rd=SigLam
c    and  1 = 10*,  2 = 10,  3 = 8a

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

c--------1---------2---------3---------4---------5---------6--------7--
c   1st=XiN(I=1), 2nd=XN(I=0), 3rd=SigSig(I=1), 4th=SigLam

       do i=1,4
        do j=1,4
         tempc(i,j) = 0d0
         tempt(i,j) = 0d0
        end do
       end do
       tempc(1,1) = wc11
       tempc(1,3) = wc12
       tempc(1,4) = wc13

       tempc(3,1) = wc21
       tempc(3,3) = wc22
       tempc(3,4) = wc23

       tempc(4,1) = wc31
       tempc(4,3) = wc32
       tempc(4,4) = wc33

       tempt(1,1) = wt11
       tempt(1,3) = wt12
       tempt(1,4) = wt13

       tempt(3,1) = wt21
       tempt(3,3) = wt22
       tempt(3,4) = wt23

       tempt(4,1) = wt31
       tempt(4,3) = wt32
       tempt(4,4) = wt33

       tempc(2,2) = vc2
       tempt(2,2) = vt2

c--------1---------2---------3---------4---------5---------6--------7--
c    1st=Xi0_n, 2nd=Xi-_p, 3rd=Sig+_Sig-, 4th=Sig0_Lam

       do i=1,4
        do j=1,4
         C(i,j) = 0d0
        end do
       end do
       C(1,1) =  1d0/Sqrt(2d0)
       C(1,2) =  1d0/Sqrt(2d0)
       C(2,1) =  1d0/Sqrt(2d0)
       C(2,2) = -1d0/Sqrt(2d0)
       C(3,3) =  1d0
       C(4,4) =  1d0

       do i=1,4
        do j=1,4
         vcij = 0d0
         vtij = 0d0
         do k=1,4
          do l=1,4
           vcij = vcij + C(k,i)*tempc(k,l)*C(l,j)
           vtij = vtij + C(k,i)*tempt(k,l)*C(l,j)
          end do
         end do
         vc(i,j) = vcij
         vt(i,j) = vtij
        end do
       end do

c-----------------------------------------
c    1st=3S1, 2nd=3D1

       do i=1,4
        do j=1,4
         temp(i,j,1,1) = vc(i,j)
         temp(i,j,1,2) = 2d0*Sqrt(2d0)*vt(i,j)
         temp(i,j,2,1) = 2d0*Sqrt(2d0)*vt(i,j)
         temp(i,j,2,2) = vc(i,j) - 2d0*vt(i,j)
        end do
       end do

       k = 0
       do i=1,2
        do al=1,4
         k = k + 1
         l = 0
         do j=1,2
          do bt=1,4
           l = l + 1
           v(k,l) = temp(al,bt,i,j)
          end do
         end do
        end do
       end do

       return
       end
c--------1---------2---------3---------4---------5---------6--------7--

