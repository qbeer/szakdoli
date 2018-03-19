c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=-1, XiN-potential
c
c    Xi^- n, Sig^- Lam coupled channel potential
c
c    [1S0]
c    Xn, SL
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_XN_Ch_Mnus_1S0(r,v)
       implicit   none
       real*8     r, v(2,2)

       real*8  v1, v2, w11, w12, w21, w22
       real*8  A(2,2)

       real*8  prm1(7), prm2(7), prm3(7)
       common /potential_parameter/ prm1, prm2, prm3

       if(r.eq.0d0) then
        write(*,*)"pot_XN: r should not be zero."
        stop
       end if

c-----------------------------------------
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c   v1 = 27,  v2 = 8s

        v1 = prm1(1)*DExp(-prm1(2)*r**2)
     &     + prm1(3)*DExp(-prm1(4)*r**2)
     &     - prm1(5)*( 1d0-DExp(-prm1(6)*r**2) )**2
     &         *(DExp(-prm1(7)*r)/r)**2

        v2 = prm2(1)*DExp(-prm2(2)*r**2)
     &     + prm2(3)*DExp(-prm2(4)*r**2)
     &     - prm2(5)*( 1d0-DExp(-prm2(6)*r**2) )**2
     &         *(DExp(-prm2(7)*r)/r)**2

c-----------------------------------------
c   Here 1st=XiN and 2nd=SigLam
c   and   1 = 27,  2 = 8s

       A(1,1) =  Sqrt( 2d0/5d0)
       A(1,2) = -Sqrt( 3d0/5d0)

       A(2,1) =  Sqrt( 3d0/5d0)
       A(2,2) =  Sqrt( 2d0/5d0)

       w11=A(1,1)*A(1,1)*v1 + A(1,2)*A(1,2)*v2
       w12=A(1,1)*A(2,1)*v1 + A(1,2)*A(2,2)*v2

       w21=w12
       w22=A(2,1)*A(2,1)*v1 + A(2,2)*A(2,2)*v2

c-----------------------------------------
c  Because matrix for isospin to chage is the unit matrix.

       v(1,1) = w11
       v(1,2) = w12

       v(2,1) = w21
       v(2,2) = w22

       return
       end
c--------1---------2---------3---------4---------5---------6--------7--

