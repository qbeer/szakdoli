************************************************************************
      subroutine lorentz(betax, betay, betaz, px, py, pz, en)
*                                                                      *
*     transforms the fourvector (en,px,py,pz) to a frame               *
*       moving with velocity (betax, betay, betaz)                     *
*       relative to the original frame                                 *
************************************************************************
c        be careful  it goes backwards, hw
      implicit none

      real*8   betax, betay, betaz, gamma
      real*8   px, py, pz, en
      real*8   trafo, bep, trans, pxz, pyz, pzz, enz

*      write(*,*)'begin of boost ', en**2 - px**2 -py**2 -pz**2

*---------------------------------------------------------------------*
*             evaluate gamma                                          *

      gamma = betax**2 + betay**2 + betaz**2
      if(gamma.ge.0.9999999) then
        write(*,*)'beta in lorentz 1 ', gamma,
     &     betax, betay, betaz, px, py, pz, en
        betax = betax/sqrt(gamma+0.0000001)
        betay = betay/sqrt(gamma+0.0000001)
        betaz = betaz/sqrt(gamma+0.0000001)
        gamma = 0.9999999
c        stop
      end if
      gamma = 1.0 - gamma

      if(gamma .gt. 0.0) then
        gamma = 1.0/sqrt(gamma)
      else
        write(*,*)'gamma in lorentz lt 0.0 ', gamma,
     &     betax, betay, betaz, px, py, pz, en
        stop
      end if
*                                                                     *
*---------------------------------------------------------------------*


            bep      = betax * px + betay * py + betaz * pz
            trafo    = gamma / (gamma + 1.0)
            trans    = trafo * bep - en
            pxz      = px + gamma * betax * trans
            pyz      = py + gamma * betay * trans
            pzz      = pz + gamma * betaz * trans
            enz      = gamma*(en - bep)

            px = pxz
            py = pyz
            pz = pzz
            en = enz
*      write(*,*)'end of boost ', en**2 - px**2 -py**2 -pz**2
      return
      end


      subroutine lorentz_hw(v,xi,xf)
ccccccccccccccccccccccccccccccccc     implicit  real*8 (a-h,p-z)
      implicit none
      real*8 v(0:3),xi(0:3),xf(0:3),eins,h1,h2,h3,vv,gam,xv,xtrans
c..................................................
c        transforms  x_final = L(v) * x_initial
c..................................................
      data  eins /1.e0/
      h1 = v(1)/v(0)
      h2 = v(2)/v(0)
      h3 = v(3)/v(0)
      vv   = h1**2 + h2**2 +h3**2
      gam = eins/sqrt(eins-vv)
c
      xv = xi(1)*h1 + xi(2)*h2 + xi(3)*h3
      xtrans = gam*(xi(0) + gam/(gam+eins)*xv)
c
      xf(1) = xi(1) + h1*xtrans
      xf(2) = xi(2) + h2*xtrans
      xf(3) = xi(3) + h3*xtrans
      xf(0) = gam * (xi(0) + xv)
c
      RETURN
      END


