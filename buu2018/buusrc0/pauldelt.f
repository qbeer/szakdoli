
************************************************************************
*                                                                      *
        subroutine pauldelt(pa,px,py,pz,pkf,iseed,phase,ntag,aa,bb)
*                                                                      *
*       purpose  : calulates pauliblocking for delta decay in          *
*                  local density approximation                         *
*                                                                      *
*       variables:                                                     *
*          pa       - momentum of outgoing pnucleon in delta rest frame*
*          ntag     - flag which tells if phase-space is pauli-blocked *
*                     ntag =  0 => phase space open                    *
*                     ntag = -1 => phase space blocked                 *
*          iseed    - seed for random number generator (integer,input) *
*          phase    - phase space factor                (real, output) *
*          px,py,pz - momenta of decaying delta                        *
*          pkf      - fermi momentum
*             differential c.s. for delta decay ds/dcos= aa + bb|cos|
*          aa       - diff. c.s parameter
*          bb       - diff. c.s parameter
*                                                                      *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      integer ntag,iseed
      real*8 phase,pdelt,px,py,pz,pa,pkf,aa,bb,thetmax,thet,rn
*-----------------------------------------------------------------------
      ntag=-1
      phase=1.0
      pdelt=sqrt(px**2+py**2+pz**2)
      if (pdelt + pa.lt.pkf)                                return
      if (pa.gt.pdelt+pkf.or.pdelt.gt.pa+pkf ) then
      phase=0.0
      ntag=0
      return
      else
      thetmax=abs(acos((-pkf**2+pa**2+pdelt**2)/(2.0*pdelt*pa)))
      thet=thetmax-pi/2.
      if (thet.lt.0.0)
     * phase= (aa*thetmax+bb*sin(thetmax))/(aa*pi+2*bb)
      if (thet.gt.0.0)
     * phase= (aa*thet-bb*sin(thetmax)+bb+aa*pi/2.+bb)/(aa*pi+2*bb)
      ntag=0
      if (phase.gt.rn(iseed)) ntag=-1
      endif
      return
      end
