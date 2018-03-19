c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=-1, YN-potential
c
c    Sig^- n
c
c    [1S0],     [3S1],     [3D1]
c     Sn         Sn         Sn
c
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_YN_Ch_Minus(r,v)
       implicit   none
       real*8     r, v(3,3)

       real*8  v1, vc3, vt3

       real*8  prm1(7), prm2(7), prm3(7)
       common /potential_parameter/ prm1, prm2, prm3

       real*8  prmC1(7), prmC2(7), prmC3(7)
       common /potential_parameter_Vc/ prmC1, prmC2, prmC3

       real*8  prmT1(6), prmT2(6), prmT3(6)
       common /potential_parameter_Vt/ prmT1, prmT2, prmT3

c-----------------------------------------

       v1 =  prm1(1)*DExp(-prm1(2)*r**2)
     &     + prm1(3)*DExp(-prm1(4)*r**2)
     &     - prm1(5)*( 1d0-DExp(-prm1(6)*r**2) )**2
     &         *(DExp(-prm1(7)*r)/r)**2

c-----------------------------------------

       vc3 = prmC3(1)*DExp(-prmC3(2)*r**2)
     &     + prmC3(3)*DExp(-prmC3(4)*r**2)
     &     - prmC3(5)*( 1d0-DExp(-prmC3(6)*r**2) )**2
     &         *(DExp(-prmC3(7)*r)/r)**2

c-----------------------------------------

       vt3  = prmT3(1)*(1d0-DExp(-prmT3(2)*r**2))**2
     &               *(1d0+3d0/(prmT3(3)*r)+3d0/(prmT3(3)*r)**2)
     &               *DExp(-prmT3(3)*r)/r
     &      + prmT3(4)*(1d0-DExp(-prmT3(5)*r**2))**2
     &               *(1d0+3d0/(prmT3(6)*r)+3d0/(prmT3(6)*r)**2)
     &               *DExp(-prmT3(6)*r)/r

c-----------------------------------------

       v(1,1) = v1
       v(1,2) = 0d0
       v(1,3) = 0d0

       v(2,1) = V(1,2)
       v(2,2) = vc3
       v(2,3) = 2d0*Sqrt(2d0)*vt3

       v(3,1) = V(1,3)
       v(3,2) = V(2,3)
       v(3,3) = vc3 - 2d0*vt3

       return
       end
c--------1---------2---------3---------4---------5---------6--------7--
