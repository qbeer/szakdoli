c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=0, XiN-potential
c
c    Xi^0 n, Xi^- p, Sig^+ Sig^-, Sig^0 Sig^0, Sig^0 Lam, Lam Lam
c                                           coupled channel potential
c
c    [1S0]
c    Xn, Xp, S+S-, S0S0, SL, LL
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_XN_Ch_Zero_1S0(r,v)
       implicit   none
       real*8     r, v(6,6)

       real*8  v1, v2, v3,  w1, w2, w3
       real*8  w11,w12,w13, w21,w22,w23, w31,w32,w33
       real*8  x11,x12,x13, x21,x22,x23, x31,x32,x33
       real*8  y11,y12, y21,y22, vij
       real*8  A(3,3), B(2,2), C(6,6), temp(6,6)
       integer i,j,k,l

       real*8  prm1(7), prm2(7), prm3(7)
       common /potential_parameter/ prm1, prm2, prm3

       if(r.eq.0d0) then
        write(*,*)"pot_XN: r should not be zero."
        stop
       end if

c--------1---------2---------3---------4---------5---------6--------7--
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c   v1 = 27,  v2 = 8s,  v3 = 1

       v1  = prm1(1)*DExp(-prm1(2)*r**2)
     &     + prm1(3)*DExp(-prm1(4)*r**2)
     &     - prm1(5)*( 1d0-DExp(-prm1(6)*r**2) )**2
     &         *(DExp(-prm1(7)*r)/r)**2

       v2  = prm2(1)*DExp(-prm2(2)*r**2)
     &     + prm2(3)*DExp(-prm2(4)*r**2)
     &     - prm2(5)*( 1d0-DExp(-prm2(6)*r**2) )**2
     &         *(DExp(-prm2(7)*r)/r)**2

       v3  = prm3(1)*DExp(-prm3(2)*r**2)
     &     + prm3(3)*DExp(-prm3(4)*r**2)
     &     - prm3(5)*( 1d0-DExp(-prm3(6)*r**2) )**2
     &         *(DExp(-prm3(7)*r)/r)**2


c--------1---------2---------3---------4---------5---------6--------7--
c   For Isospin = 0

c   w1 = 1,  w2 = 8s,  w3 = 27

       w1 = v3
       w2 = v2
       w3 = v1

c  1st=LamLam and 2nd=SigSig, 3rd=NXi
c   and 1 = 1, 2 = 8s, 3 = 27

       A(1,1) = -Sqrt( 5d0/40d0)
       A(1,2) = -Sqrt( 8d0/40d0)
       A(1,3) =  Sqrt(27d0/40d0)

       A(2,1) =  Sqrt(15d0/40d0)
       A(2,2) = -Sqrt(24d0/40d0)
       A(2,3) = -Sqrt( 1d0/40d0)

       A(3,1) =  Sqrt(20d0/40d0)
       A(3,2) =  Sqrt( 8d0/40d0)
       A(3,3) =  Sqrt(12d0/40d0)

       w11=A(1,1)*A(1,1)*w1 + A(1,2)*A(1,2)*w2 + A(1,3)*A(1,3)*w3
       w12=A(1,1)*A(2,1)*w1 + A(1,2)*A(2,2)*w2 + A(1,3)*A(2,3)*w3
       w13=A(1,1)*A(3,1)*w1 + A(1,2)*A(3,2)*w2 + A(1,3)*A(3,3)*w3

       w21=w12
       w22=A(2,1)*A(2,1)*w1 + A(2,2)*A(2,2)*w2 + A(2,3)*A(2,3)*w3
       w23=A(2,1)*A(3,1)*w1 + A(2,2)*A(3,2)*w2 + A(2,3)*A(3,3)*w3

       w31=w13
       w32=w23
       w33=A(3,1)*A(3,1)*w1 + A(3,2)*A(3,2)*w2 + A(3,3)*A(3,3)*w3

c  So that 1st=XiN, 2nd=SigSig, 3rd=LamLam

       x11 = w33
       x12 = w32
       x13 = w31

       x21 = w23
       x22 = w22
       x23 = w21

       x31 = w13
       x32 = w12
       x33 = w11

c--------1---------2---------3---------4---------5---------6--------7--
c   For Isospin = 1

c   Here 1st=XiN and 2nd=SigLam
c    and 1 = 27,  2 = 8s

       B(1,1) =  Sqrt( 2d0/5d0)
       B(1,2) = -Sqrt( 3d0/5d0)

       B(2,1) =  Sqrt( 3d0/5d0)
       B(2,2) =  Sqrt( 2d0/5d0)

       y11=B(1,1)*B(1,1)*v1 + B(1,2)*B(1,2)*v2
       y12=B(1,1)*B(2,1)*v1 + B(1,2)*B(2,2)*v2

       y21=y12
       y22=B(2,1)*B(2,1)*v1 + B(2,2)*B(2,2)*v2

c--------1---------2---------3---------4---------5---------6--------7--
c   1st=XiN(I=1), 2nd=XN(I=0), 3rd=SigSig(I=2), 4th=SigSig(I=0),
c   5th=SigLam,   6th=LamLam

       do i=1,6
        do j=1,6
         temp(i,j) = 0d0
        end do
       end do
       temp(1,1) = y11
       temp(1,5) = y12
       temp(5,1) = y21
       temp(5,5) = y22

       temp(2,2) = x11
       temp(2,4) = x12
       temp(2,6) = x13

       temp(4,2) = x21
       temp(4,4) = x22
       temp(4,6) = x23

       temp(6,2) = x31
       temp(6,4) = x32
       temp(6,6) = x33

       temp(3,3) = v1

c--------1---------2---------3---------4---------5---------6--------7--
c    1st=Xi0_n, 2nd=Xi-_p, 3rd=Sig+_Sig-, 4th=Sig0_Sig0, 5th=Sig0_Lam,
c    6th=LamLam

       do i=1,6
        do j=1,6
         C(i,j) = 0d0
        end do
       end do
       C(1,1) =  1d0/Sqrt(2d0)
       C(1,2) =  1d0/Sqrt(2d0)
       C(2,1) =  1d0/Sqrt(2d0)
       C(2,2) = -1d0/Sqrt(2d0)

       C(3,3) =  1d0/Sqrt(3d0)
       C(3,4) =  Sqrt(2d0/3d0)
       C(4,3) = -Sqrt(2d0/3d0)
       C(4,4) =  1d0/Sqrt(3d0)

       C(5,5) =  1d0
       C(6,6) =  1d0

       do i=1,6
        do j=1,6
         vij = 0d0
         do k=1,6
          do l=1,6
           vij = vij + C(k,i)*temp(k,l)*C(l,j)
          end do
         end do
         v(i,j) = vij
        end do
       end do

       return
       end
c--------1---------2---------3---------4---------5---------6--------7--

