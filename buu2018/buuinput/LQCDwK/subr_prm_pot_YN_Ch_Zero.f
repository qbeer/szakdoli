c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=0, YN-potential
c
c    Lam n, Sig^0 n, Sig^- p  coupled channel potential
c
c    [1S0],     [3S1],     [3D1]
c   Ln,Sn,Sp   Ln,Sn,Sp   Ln,Sn,Sp
c
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_YN_Ch_Zero(r,v)
       implicit   none
       real*8     r, v(9,9)

       real*8  v1, v2, vc1, vc2, vc3, vt1, vt2, vt3
       real*8  w11, w12, w21, w22
       real*8  vs(3,3), vc(3,3), vt(3,3)
       real*8  AS(2,2), AT(2,2), B(3,3), C(3,3)

       integer i,j,k,l

       real*8  prm1(7), prm2(7), prm3(7)
       common /potential_parameter/ prm1, prm2, prm3

       real*8  prmC1(7), prmC2(7), prmC3(7)
       common /potential_parameter_Vc/ prmC1, prmC2, prmC3

       real*8  prmT1(6), prmT2(6), prmT3(6)
       common /potential_parameter_Vt/ prmT1, prmT2, prmT3

c-----------------------------------------
c for Irr-rep to Spin-Singlet

       AS(1,1) =  Sqrt( 1d0/10d0)
       AS(1,2) =  Sqrt( 9d0/10d0)

       AS(2,1) =  Sqrt( 9d0/10d0)
       AS(2,2) = -Sqrt( 1d0/10d0)

c-----------------------------------------
c for Irr-rep to Spin-Triplet

       AT(1,1) =  Sqrt( 1d0/2d0)
       AT(1,2) =  Sqrt( 1d0/2d0)

       AT(2,1) =  Sqrt( 1d0/2d0)
       AT(2,2) = -Sqrt( 1d0/2d0)

c-----------------------------------------
c for Iso-spin to Chage

       B(1,1) = 1d0
       B(1,2) = 0d0
       B(1,3) = 0d0

       B(2,1) =  0d0
       B(2,2) =  Sqrt(1d0/3d0)
       B(2,3) = -Sqrt(2d0/3d0)

       B(3,1) =  0d0
       B(3,2) =  Sqrt(2d0/3d0)
       B(3,3) =  Sqrt(1d0/3d0)

c-----------------------------------------

       v1 =  prm1(1)*DExp(-prm1(2)*r**2)
     &     + prm1(3)*DExp(-prm1(4)*r**2)
     &     - prm1(5)*( 1d0-DExp(-prm1(6)*r**2) )**2
     &         *(DExp(-prm1(7)*r)/r)**2

       v2 =  prm2(1)*DExp(-prm2(2)*r**2)
     &     + prm2(3)*DExp(-prm2(4)*r**2)
     &     - prm2(5)*( 1d0-DExp(-prm2(6)*r**2) )**2
     &         *(DExp(-prm2(7)*r)/r)**2

c-----------------------------------------

       vc1 = prmC1(1)*DExp(-prmC1(2)*r**2)
     &     + prmC1(3)*DExp(-prmC1(4)*r**2)
     &     - prmC1(5)*( 1d0-DExp(-prmC1(6)*r**2) )**2
     &         *(DExp(-prmC1(7)*r)/r)**2

       vc2 = prmC2(1)*DExp(-prmC2(2)*r**2)
     &     + prmC2(3)*DExp(-prmC2(4)*r**2)
     &     - prmC2(5)*( 1d0-DExp(-prmC2(6)*r**2) )**2
     &         *(DExp(-prmC2(7)*r)/r)**2

       vc3 = prmC3(1)*DExp(-prmC3(2)*r**2)
     &     + prmC3(3)*DExp(-prmC3(4)*r**2)
     &     - prmC3(5)*( 1d0-DExp(-prmC3(6)*r**2) )**2
     &         *(DExp(-prmC3(7)*r)/r)**2

c-----------------------------------------

       vt1  = prmT1(1)*(1d0-DExp(-prmT1(2)*r**2))**2
     &               *(1d0+3d0/(prmT1(3)*r)+3d0/(prmT1(3)*r)**2)
     &               *DExp(-prmT1(3)*r)/r
     &      + prmT1(4)*(1d0-DExp(-prmT1(5)*r**2))**2
     &               *(1d0+3d0/(prmT1(6)*r)+3d0/(prmT1(6)*r)**2)
     &               *DExp(-prmT1(6)*r)/r


       vt2  = prmT2(1)*(1d0-DExp(-prmT2(2)*r**2))**2
     &               *(1d0+3d0/(prmT2(3)*r)+3d0/(prmT2(3)*r)**2)
     &               *DExp(-prmT2(3)*r)/r
     &      + prmT2(4)*(1d0-DExp(-prmT2(5)*r**2))**2
     &               *(1d0+3d0/(prmT2(6)*r)+3d0/(prmT2(6)*r)**2)
     &               *DExp(-prmT2(6)*r)/r

       vt3  = prmT3(1)*(1d0-DExp(-prmT3(2)*r**2))**2
     &               *(1d0+3d0/(prmT3(3)*r)+3d0/(prmT3(3)*r)**2)
     &               *DExp(-prmT3(3)*r)/r
     &      + prmT3(4)*(1d0-DExp(-prmT3(5)*r**2))**2
     &               *(1d0+3d0/(prmT3(6)*r)+3d0/(prmT3(6)*r)**2)
     &               *DExp(-prmT3(6)*r)/r

c-----------------------------------------

       w11 = AS(1,1)*AS(1,1)*v1 + AS(1,2)*AS(1,2)*v2
       w12 = AS(1,1)*AS(2,1)*v1 + AS(1,2)*AS(2,2)*v2

       w21 = w12
       w22 = AS(2,1)*AS(2,1)*v1 + AS(2,2)*AS(2,2)*v2

       C(1,1) = w22
       C(1,2) = w21
       C(1,3) = 0d0
       C(2,1) = w12
       C(2,2) = w11
       C(2,3) = 0d0
       C(3,1) = 0d0
       C(3,2) = 0d0
       C(3,3) = v1

       do i=1,3
        do j=1,3
         vs(i,j) = 0d0
         do k=1,3
          do l=1,3
           vs(i,j) = vs(i,j) + B(k,i)*C(k,l)*B(l,j)
          end do
         end do
        end do
       end do

c-----------------------------------------

       w11 = AT(1,1)*AT(1,1)*vc1 + AT(1,2)*AT(1,2)*vc2
       w12 = AT(1,1)*AT(2,1)*vc1 + AT(1,2)*AT(2,2)*vc2

       w21 = w12
       w22 = AT(2,1)*AT(2,1)*vc1 + AT(2,2)*AT(2,2)*vc2

       C(1,1) = w22
       C(1,2) = w21
       C(1,3) = 0d0
       C(2,1) = w12
       C(2,2) = w11
       C(2,3) = 0d0
       C(3,1) = 0d0
       C(3,2) = 0d0
       C(3,3) = vc3

       do i=1,3
        do j=1,3
         vc(i,j) = 0d0
         do k=1,3
          do l=1,3
           vc(i,j) = vc(i,j) + B(k,i)*C(k,l)*B(l,j)
          end do
         end do
        end do
       end do

c-----------------------------------------

       w11 = AT(1,1)*AT(1,1)*vt1 + AT(1,2)*AT(1,2)*vt2
       w12 = AT(1,1)*AT(2,1)*vt1 + AT(1,2)*AT(2,2)*vt2

       w21 = w12
       w22 = AT(2,1)*AT(2,1)*vt1 + AT(2,2)*AT(2,2)*vt2

       C(1,1) = w22
       C(1,2) = w21
       C(1,3) = 0d0
       C(2,1) = w12
       C(2,2) = w11
       C(2,3) = 0d0
       C(3,1) = 0d0
       C(3,2) = 0d0
       C(3,3) = vt3

       do i=1,3
        do j=1,3
         vt(i,j) = 0d0
         do k=1,3
          do l=1,3
           vt(i,j) = vt(i,j) + B(k,i)*C(k,l)*B(l,j)
          end do
         end do
        end do
       end do

c-----------------------------------------

       do i=1,9
        do j=1,9
         v(i,j) = 0d0
        end do
       end do

       do i=1,3
        do j=i,3
         v(i,j)     = vs(i,j)
         v(3+i,3+j) = vc(i,j)
         v(6+i,6+j) = vc(i,j)-2d0*vt(i,j)
        end do
       end do

       do i=1,3
        do j=1,3
         v(3+i,6+j) = 2d0*Sqrt(2d0)*vt(i,j)
        end do
       end do

       do i=2,9
        do j=1,i-1
         v(i,j) = v(j,i)
        end do
       end do

       return
       end
c--------1---------2---------3---------4---------5---------6--------7--
