c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=-1, XiN-potential : XmN1S0(2,2)
c    Xi^- n, Sig^- Lam coupled channel potential: XmN1S0
c    [1S0]
c
c   Ch=-1, XiN-potential :  XmN3SD(6,6)
c      Xi^- n, Sig^- Sig0, Sig^- Lam coupled channel potential
c      [3S1]         [3D1]
c     Xn, SS, SL,   Xn, SS, SL
c
c   Ch=+1, XiN-potential : XpN1S0(2,2)
c    Xi^0 p, Sig^+ Lam coupled channel potential
c    [1S0]
c    Xp, SL
c
c   Ch=+1, XiN-potential : XpN3SD(6,6)
c    Xi^0 p, Sig^+ Sig0, Sig^+ Lam coupled channel potential
c      [3S1]         [3D1]
c    Xp, SS, SL,   Xp, SS, SL
c
c   Ch=0, XiN-potential : XzN1S0
c    Xi^0 n, Xi^- p, Sig^+ Sig^-, Sig^0 Sig^0, Sig^0 Lam, Lam Lam
c    [1S0]
c    Xn, Xp, S+S-, S0S0, SL, LL
c
c   Ch=0, XiN-potential  : XzN3SD
c    Xi^0 n, Xi^- p, Sig^+ Sig^-, Sig^0 Lam coupled channel potential
c      [3S1]          [3D1]
c    Xn,Xp,SS,SL,   Xn,Xp,SS,SL
c
c   Ch=2, YN-potential  : Spp
c    Sig^+ p
c    [1S0],     [3S1],     [3D1]
c     Sp         Sp         Sp
c
c   Ch=-1, YN-potential  : YmN
c    Sig^- n
c    [1S0],     [3S1],     [3D1]
c     Sn         Sn         Sn
c
c   Ch=1, YN-potential  : YpN
c    Lam p, Sig^0 p, Sig^+ n  coupled channel potential
c    [1S0],     [3S1],     [3D1]
c   Lp,Sp,Sn   Lp,Sp,Sn   Lp,Sp,Sn
c
c   Ch=0, YN-potential  : YzN
c    Lam n, Sig^0 n, Sig^- p  coupled channel potential
c    [1S0],     [3S1],     [3D1]
c   Ln,Sn,Sp   Ln,Sn,Sp   Ln,Sn,Sp
c
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

      implicit   none

      real*8  rr,dr,hatmin,hatmax
      real*8  XmN1S0(2,2),XmN3SD(6,6),XpN1S0(2,2),XpN3SD(6,6),Spp(3,3)
      real*8  XzN1S0(6,6),XzN3SD(8,8),YmN(3,3),YpN(9,9),YzN(9,9)
      real*8 XN_Minus,XmNPot,XN_Plus,XpNPot,XN_Zeron,XznPot
      real*8 XN_Zerop,XzpPot,YN_2Plus,SppPot,YN_Minus,SmnPot
      real*8 YN_Plusp,SzpPot,YN_Zerop,SmpPot,YN_Plusn,SpnPot
      real*8 YN_Zeron,SznPot,YN_PlusL,SpLPot,YN_ZeroL,SzLPot
      real*8 pi4
      integer ii
      external XN_Minus,XN_Plus,XN_Zeron,XN_Zerop,YN_2Plus,YN_Minus
      external YN_Plusp,YN_Zerop,YN_Plusn,YN_Zeron,YN_PlusL,YN_ZeroL

      open(50,file='Xi_pot.dat')
      write(50,'("# r    XmN1S0   XmN3S1   XzN1S0   XzN3S1")')
      open(51,file='Sigma_pot.dat')
      write(51,'("# r     Spp0     Spp1     Szp0     Szp1     ",
     &   "Szn0     Szn1     LN0      LN1")')
      open(52,file='XiSigmaL_Dwave_pot.dat')
      write(52,'("# r     XzpD     XmpD     SppD     SzpD     ",
     &   "LND")')

      pi4 = 12.5664
      dr=0.1
      do ii=1,40
        rr=real(ii)*dr
        call prm_pot_XN_Ch_Mnus_1S0(rr,XmN1S0) 
        call prm_pot_XN_Ch_Mnus_3SD(rr,XmN3SD) 
        call prm_pot_XN_Ch_Plus_1S0(rr,XpN1S0) 
        call prm_pot_XN_Ch_Plus_3SD(rr,XpN3SD) 
        call prm_pot_XN_Ch_Zero_1S0(rr,XzN1S0) 
        call prm_pot_XN_Ch_Zero_3Sd(rr,XzN3SD) 
        call prm_pot_YN_Ch_2Plus(rr,Spp) 
        call prm_pot_YN_Ch_Minus(rr,YmN) 
        call prm_pot_YN_Ch_Plus(rr,YpN) 
        call prm_pot_YN_Ch_Zero(rr,YzN) 
        write(50,101)rr,XmN1S0(1,1),XmN3SD(1,1),XzN1S0(2,2),XzN3Sd(1,1)
        write(51,101)rr,Spp(1,1),Spp(2,2),YpN(2,2),YpN(5,5),
     &                  YpN(3,3),YpN(6,6),YpN(1,1),YpN(4,4)
        write(52,101)rr,XmN3SD(4,4),XzN3Sd(5,5),Spp(3,3),YpN(8,8),
     &                  YpN(7,7)
      end do
      hatmin=1.d-2
      hatmax=3.0
      call dqromb(XN_Minus,hatmin,hatmax,XmNPot)
      call dqromb(XN_Plus,hatmin,hatmax,XpNPot)
      call dqromb(XN_Zeron,hatmin,hatmax,XznPot)
      call dqromb(XN_Zerop,hatmin,hatmax,XzpPot)
      call dqromb(YN_2Plus,hatmin,hatmax,SppPot)
      call dqromb(YN_Minus,hatmin,hatmax,SmnPot)
      call dqromb(YN_Plusp,hatmin,hatmax,SzpPot)
      call dqromb(YN_Zerop,hatmin,hatmax,SmpPot)
      call dqromb(YN_Plusn,hatmin,hatmax,SpnPot)
      call dqromb(YN_Zeron,hatmin,hatmax,SznPot)
      call dqromb(YN_PlusL,hatmin,hatmax,SpLPot)
      call dqromb(YN_ZeroL,hatmin,hatmax,SzLPot)
      write(*,'("Xmn",f9.4,"  Xmp",f9.4,"  Xzp",f9.4,"  Xzn",f9.4,
     &   "  SpL",f9.4,"  SzL",f9.4)')pi4*XmNPot,pi4*XzpPot,pi4*XpNPot,
     &    pi4*XznPot,pi4*SpLPot,pi4*SzLPot 
      write(*,'("Spp",f9.4,"  Spn",f9.4,"  Smn",f9.4,"  Smp",f9.4,
     &"  Szp",f9.4,"  Szn",f9.4)')pi4*SppPot,pi4*SpnPot,pi4*SmnPot,
     &    pi4*SmpPot,pi4*SzpPot,pi4*SznPot
 101  format(f5.2,12f9.3)
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function XN_Minus(xx)
c    Xi^- n      
      implicit none
      real*8 xx,XN_Minus,XmN1S0(2,2),XmN3SD(6,6)
      call prm_pot_XN_Ch_Mnus_1S0(xx,XmN1S0) 
      call prm_pot_XN_Ch_Mnus_3SD(xx,XmN3SD) 
      XN_Minus=xx**2*0.25*(XmN1S0(1,1)+3.0*XmN3SD(1,1))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function XN_Plus(xx)
c    Xi^0 p      
      implicit none
      real*8 xx,XN_Plus,XpN1S0(2,2),XpN3SD(6,6)
      call prm_pot_XN_Ch_Plus_1S0(xx,XpN1S0) 
      call prm_pot_XN_Ch_Plus_3SD(xx,XpN3SD) 
      XN_Plus=xx**2*0.25*(XpN1S0(1,1)+3.0*XpN3SD(1,1))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function XN_Zeron(xx)
c    Xi^0 n
      implicit none
      real*8 xx,XN_Zeron,XzN1S0(6,6),XzN3SD(8,8)
      call prm_pot_XN_Ch_Zero_1S0(xx,XzN1S0) 
      call prm_pot_XN_Ch_Zero_3SD(xx,XzN3SD) 
      XN_Zeron=xx**2*0.25*(XzN1S0(1,1)+3.0*XzN3SD(1,1))
      return
      end
c--------1---------2---------3---------4---------5---------6--------7--

      function XN_Zerop(xx)
c    Xi^- p
      implicit none
      real*8 xx,XN_Zerop,XzN1S0(6,6),XzN3SD(8,8)
      call prm_pot_XN_Ch_Zero_1S0(xx,XzN1S0) 
      call prm_pot_XN_Ch_Zero_3SD(xx,XzN3SD) 
      XN_Zerop=xx**2*0.25*(XzN1S0(2,2)+3.0*XzN3SD(2,2))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_2Plus(xx)
c    Si^+ p
      implicit none
      real*8 xx,YN_2Plus,Spp(3,3)
      call prm_pot_YN_Ch_2Plus(xx,Spp)
      YN_2Plus=xx**2*0.25*(Spp(1,1)+3.0*Spp(2,2))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_Minus(xx)
c    Si^- n
      implicit none
      real*8 xx,YN_Minus,YmN(3,3)
      call prm_pot_YN_Ch_Minus(xx,YmN)
      YN_Minus=xx**2*0.25*(YmN(1,1)+3.0*YmN(2,2))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_Plusp(xx)
c    Si^0 p
      implicit none
      real*8 xx,YN_Plusp,YpN(9,9)
      call prm_pot_YN_Ch_Plus(xx,YpN)
      YN_Plusp=xx**2*0.25*(YpN(2,2)+3.0*YpN(5,5))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_Plusn(xx)
c    Si^+ n
      implicit none
      real*8 xx,YN_Plusn,YpN(9,9)
      call prm_pot_YN_Ch_Plus(xx,YpN)
      YN_Plusn=xx**2*0.25*(YpN(3,3)+3.0*YpN(6,6))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_PlusL(xx)
c    La p
      implicit none
      real*8 xx,YN_PlusL,YpN(9,9)
      call prm_pot_YN_Ch_Plus(xx,YpN)
      YN_PlusL=xx**2*0.25*(YpN(1,1)+3.0*YpN(4,4))
      return
      end
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_Zeron(xx)
c    Si^0 n
      implicit none
      real*8 xx,YN_Zeron,YzN(9,9)
      call prm_pot_YN_Ch_Zero(xx,YzN)
      YN_Zeron=xx**2*0.25*(YzN(2,2)+3.0*YzN(5,5))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_Zerop(xx)
c    Si^- p
      implicit none
      real*8 xx,YN_Zerop,YzN(9,9)
      call prm_pot_YN_Ch_Zero(xx,YzN)
      YN_Zerop=xx**2*0.25*(YzN(3,3)+3.0*YzN(6,6))
      return
      end
      
c--------1---------2---------3---------4---------5---------6--------7--

      function YN_ZeroL(xx)
c    La n
      implicit none
      real*8 xx,YN_ZeroL,YzN(9,9)
      call prm_pot_YN_Ch_Zero(xx,YzN)
      YN_ZeroL=xx**2*0.25*(YzN(1,1)+3.0*YzN(4,4))
      return
      end
      
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

c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=-1, XiN-potential
c
c    Xi^- n, Sig^- Sig0, Sig^- Lam coupled channel potential
c
c      [3S1]         [3D1]
c    Xn, SS, SL,   Xn, SS, SL
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_XN_Ch_Mnus_3SD(r,v)
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

c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=+1, XiN-potential
c
c    Xi^0 p, Sig^+ Lam coupled channel potential
c
c    [1S0]
c    Xp, SL
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_XN_Ch_Plus_1S0(r,v)
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

c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=2, YN-potential
c
c    Sig^+ p
c
c    [1S0],     [3S1],     [3D1]
c     Sp         Sp         Sp
c
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_YN_Ch_2Plus(r,v)
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
c--------1---------2---------3---------4---------5---------6--------7--
c   Ch=1, YN-potential
c
c    Lam p, Sig^0 p, Sig^+ n  coupled channel potential
c
c    [1S0],     [3S1],     [3D1]
c   Lp,Sp,Sn   Lp,Sp,Sn   Lp,Sp,Sn
c
c   Vc = Gauss + Gauss + ((1-Gauss)*Yukawa)^2
c
c      real*8 r : distance  [fm]
c      real*8 v : potential [MeV]
c
c--------1---------2---------3---------4---------5---------6--------7--

       subroutine prm_pot_YN_Ch_Plus(r,v)
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
       B(2,2) = -Sqrt(1d0/3d0)
       B(2,3) =  Sqrt(2d0/3d0)

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
c--------1---------2---------3---------4---------5---------6--------7--
c  Block data for LQCD YN and YY potentials measured with the Kei-conf
c
c                                        by Takashi Inoue, Aug 2016
c--------1---------2---------3---------4---------5---------6--------7--
       block data LQCD_YN_POT_Kei
       implicit   none

       real*8  prm1(7), prm2(7), prm3(7)
       common /potential_parameter/ prm1, prm2, prm3

       real*8  prmc1(7), prmc2(7), prmc3(7)
       common /potential_parameter_Vc/ prmc1, prmc2, prmc3

       real*8  prmt1(6), prmt2(6), prmt3(6)
       common /potential_parameter_Vt/ prmt1, prmt2, prmt3

       integer    i

c----------------------
c FL_27
       data (prm1(i),i=1,7) /
     &   1336.6549d0  ,  2.3742678d0  ,  1433.5221d0  ,  53.659831d0  ,
     &   141770.91d0  , 0.52964043d0  ,  2.5293523d0    /

c FL_8S
       data (prm2(i),i=1,7) /
     &   5638.7522d0  ,  12.951368d0  ,  411.74169d0  ,  1.2813190d0  ,
     &   1636.5727d0  ,  15.959756d0  ,  3.4365892d0  /

c FL_1
       data (prm3(i),i=1,7) /
     &  -966.67715d0  ,  3.1380113d0  , -84.682264d0  , 0.69326043d0  ,
     & -0.92626920d10 , 0.25353068d-2 ,  3.6762249d0  /

c----------------------
c FL_10*
       data (prmc1(i),i=1,7) /
     &   819.16615d0  ,  1.9129865d0  ,  1416.7240d0  ,  45.450617d0  ,
     &   51607.561d0  , 0.60079497d0  ,  2.1455409d0  /

c FL_8A
       data (prmc2(i),i=1,7) /
     &   223.30741d0  ,  39.544058d0  ,  193.89482d0  ,  1.6640125d0  ,
     &   242551.91d0  , 0.15265522d0  ,  2.2494186d0  /

c FL_10
       data (prmc3(i),i=1,7) /
     &   2180.6380d0  ,  7.9816661d0  ,  579.61863d0  ,  3.7219714d0  ,
     &   32656.853d0  ,  2.1772221d0  ,  4.1936891d0  /

c----------------------
c FL_10*
       data (prmt1(i),i=1,6) /
     &   11301.227d0  ,  9.4780722d0  ,  1.8860888d0  , -11202.034d0  ,
     &   9.4839043d0  ,  1.8782561d0  /

c FL_8A
       data (prmt2(i),i=1,6) /
     &  -15724.464d0  ,  3.2670958d0  ,  14.519991d0  ,  6.6975009d0  ,
     &   24.784160d0  ,  2.2239908d0  /

c FL_10
       data (prmt3(i),i=1,6) /
     &  -2109.3066d0  ,  16.472853d0  ,  1.2352148d0  ,  2092.6904d0  ,
     &   16.475335d0  ,  1.2303626d0  /
c----------------------
       end

c--------1---------2---------3---------4---------5---------6--------7--
c                END                           
c--------1---------2---------3---------4---------5---------6--------7--
c***********************************************************************
      subroutine dqromb(func,a,b,ss)
*
*      return ss as the integral of the function func from a to b.
*      integration is performed by romberg's method of order 2k,
*      where k=2 is simpson's rule
***********************************************************************
      implicit none
      real*8 func,a,b,ss,s,h,eps,dss
      external func
      integer jmax,jmaxp,k,km,j,l
      parameter(eps=1.d-6,jmax=90,jmaxp=jmax+1,k=10,km=k-1)
*      here eps is the fractional accuracy desired, as determined by
*      the extrapolation error estimate;
*      jmax limits the total number of steps;
*      k is the number of points used in the extrapolation
      dimension s(jmaxp),h(jmaxp)
      h(1)=1.
      do 11 j=1,jmax
        call trapzd(func,a,b,s(j),j)
        if (j.ge.k) then
	  l=j-km
          call polint(h(l),s(l),k,0.d0,ss,dss)
          if (dabs(dss).lt.eps*dabs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      write(*,*)'qrom:too many steps hiba:',abs(dss),eps*abs(ss)
      return
      end
************************************************************************
      subroutine trapzd(func,a,b,s,n)
      implicit none
      real*8 func,a,b,s,del,tnm,x,sum
      integer N,it,j
      save it
C     external func line may be wrong, but without it the compiler cries
c      external func
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
        it=1
      else
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
        it=2*it
      endif
      return
      end
***********************************************************************
      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      real*8 xa,ya,x,y,dy,c,d,dif,dift,ho,hp,w,den
      integer N,nmax,ns,I,m
      parameter (nmax=20)
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0) write(*,*) 'pause in qromb.f'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end

