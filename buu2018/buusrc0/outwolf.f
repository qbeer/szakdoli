c**********************************************************************c
      subroutine outdwolf(mcoll,mcol,time2,ictime,isubsqt,
     &                    mmesc,mmesa, mesc,mesa,mumes,lumes,
     &  mcopt,mcond,mcodd,mdpir,mdpiq,mrpiq,mept,mend,medd,ndis,
     &  mprt,mntr,mde2,mde1,mde0,mdem,mrpl,mrze,mqpl,mqze,
     &  mcnu,mcde,mcrn,mcqn,rati,suru,ener,resdens)

      implicit none

      include"common"
      include"coucom"
      include"cominput"

      integer dimlcoll
      parameter(dimlcoll = 3*nres+4+nres**2)
      integer dimlmesc
      parameter(dimlmesc = nres+6)

      real*8     suru(0:200),ener(0:200),rati(0:200),time2(0:200)
      real*8     resdens(-2:nres,0:200,50)
      integer  ndis(0:20,0:200)
      integer  lumes(18),mumes(18,0:200)
      integer  mcoll(-6:dimlcoll)
      integer  mmesc(dimlmesc),mmesa(dimlmesc)
      integer  mcol(-1:dimlcoll,0:200), mesc(dimlmesc,0:200),
     &             mesa(dimlmesc,0:200)
      integer  mend(0:200),medd(0:200),mept(0:200)
      integer  mprt(0:200),mntr(0:200),mrpl(0:200),mrze(0:200),
     &         mde2(0:200),mde1(0:200),mde0(0:200),mdem(0:200),
     &         mqpl(0:200),mqze(0:200),
     &         mcnu(0:200),mcde(0:200),mcrn(0:200),mcqn(0:200)
      integer ictime,i,isubsqt,j,ij,mcond,mcodd,mcopt,mdpir,mdpiq,mrpiq
      real*8 denu,xnum

      if(icoll.ne.1) then
      denu  = dt * float(nfreq*num*isubsqt)
      xnum  = 1.0/float(num*isubsqt)
*
      write(isum,'(/''c:number of collision in general''/
     &''n: time (fm/c)''/''n: collisions'')')
      write(isum,'('' c:time(fm/c)'',
     &''   all      blocked     elastic     proj+targ'')')
      write(isum,'(''c:total'',8i12)')
     & mcoll(-1),mcoll(0),mcoll(1),mcopt
      write(isum,'(f8.2,4f12.4)')
     &  (time2(i),
     &float(mcol(-1,i))*xnum,float(mcol(0,i))*xnum,
     &float(mcol(1,i))*xnum, float(mept(i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of inelastic collisions in general''/
     &''n: time (fm/c)''/''n: collisions'')')
      write(isum,'(''c:time(fm/c)'',
     &'' nn->nx   nx->nn  nx->nx1  nn->dd   dd->nd   dd->nx  xx->nx'',
     &'' x(n)x->dd'')')
      write(isum,'(''c:total'',8i9)')
     & mcoll(2)+mcoll(4)+mcoll(6),
     & mcoll(3)+mcoll(5)+mcoll(7),
     & mcoll(8)+mcoll(9)+mcoll(10)+mcoll(11)+mcoll(12)+mcoll(13),
     & mcoll(16),mcoll(17),mcoll(14)+mcoll(15),mcoll(19),
     & mcoll(20)+mcoll(21)
      write(isum,'(f8.2,8f9.4)')
     &  (time2(i),
     &float(mcol(2,i)+mcol(4,i)+mcol(6,i))*xnum,
     &float(mcol(3,i)+mcol(5,i)+mcol(7,i))*xnum,
     &float(mcol(8,i)+mcol(9,i)+mcol(10,i)+mcol(11,i)+mcol(12,i)
     &     +mcol(13,i))*xnum,
     &float(mcol(16,i))*xnum,float(mcol(17,i))*xnum,
     &float(mcol(14,i)+mcol(15,i))*xnum,float(mcol(19,i))*xnum,
     &float(mcol(20,i)+mcol(21,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of pion and eta collisions in general''/
     &''n: time (fm/c)''/''n: collisions'')')
      write(isum,'(
     &''c:time(fm/c) x->n+pi n+pi->x x+pi->x  x->x+pi'',
     &'' q->n+eta  n+eta->q '')')
      write(isum,'(''c:total'',10i9)')
     & mmesc(1)+mmesc(2)+mmesc(3),mmesa(1)+mmesa(2)+mmesa(3),
     & mmesc(4)+mmesc(5)+mmesc(6),mmesa(4)+mmesa(5)+mmesa(6),
     & mmesc(7),mmesa(7)
      write(isum,'(f8.2,6f9.4)')
     &  (time2(i),
     &float(mesc(1,i)+mesc(2,i)+mesc(3,i))*xnum,
     &float(mesa(1,i)+mesa(2,i)+mesa(3,i))*xnum,
     &float(mesc(4,i)+mesc(5,i)+mesc(6,i))*xnum,
     &float(mesa(4,i)+mesa(5,i)+mesa(6,i))*xnum,
     &float(mesc(7,i))*xnum,float(mesa(7,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of rho, sigma, omega collisions.''/
     &''n: time (fm/c)''/''n: collisions'')')
      write(isum,'(
     &''c: t    x>n+ro  n+ro>x  ro>2pi  2pi>rho'',
     &''  x>n+si  n+si>x  si>2pi  2pi>si o>pi+ro pi+ro>o'')')
      write(isum,'(''c:total'',i6,9i8)')
     & mmesc(8)+mmesc(9),mmesa(8)+mmesa(9),mmesc(12),mmesa(12),
     & mmesc(10)+mmesc(11),mmesa(10)+mmesa(11),mmesc(13),mmesa(13),
     & mmesc(14),mmesa(14)
      write(isum,'(f6.2,10f8.4)')
     &  (time2(i),
     &float(mesc(8,i)+mesc(9,i))*xnum,float(mesa(8,i)+mesa(9,i))*xnum,
     &float(mesc(12,i))*xnum,float(mesa(12,i))*xnum,
     &float(mesc(10,i)+mesc(11,i))*xnum,
     &float(mesa(10,i)+mesa(11,i))*xnum,
     &float(mesc(13,i))*xnum,float(mesa(13,i))*xnum,
     &float(mesc(14,i))*xnum,float(mesa(14,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of collision for deltas''/
     &''n: time (fm/c)''/''n: collisions'')')
      write(isum,'(
     &''c:time(fm/c) nn->nd  nd->nn   nn->dd   dd->nd   nd->nx  '',
     &''nx->nd  dd->nx x(n)x(d)->dd'')')
      write(isum,'(''c: total'',8i9)') mcoll(2),mcoll(3),mcoll(16),
     &  mcoll(17),mcoll(8)+mcoll(10),mcoll(9)+mcoll(11),
     &  mcoll(14)+mcoll(15),mcoll(20)+mcoll(21)
      write(isum,'(f8.2,8f9.4)')
     &  (time2(i),
     &float(mcol(2,i))*xnum,float(mcol(3,i))*xnum,
     &float(mcol(16,i))*xnum,float(mcol(17,i))*xnum,
     &float(mcol(8,i)+mcol(10,i))*xnum,float(mcol(9,i)+mcol(11,i))*xnum,
     &float(mcol(14,i)+mcol(15,i))*xnum,
     &float(mcol(20,i)+mcol(21,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of collision for n*1'')')
      write(isum,'(
     &''c:time(fm/c) nn->nr nr->nn   nr->nd   nd->nr   nr->nq'',
     &''   nq->nr   dd->nr'')')
      write(isum,'(''c: total'',7i9)') mcoll(4),mcoll(5),mcoll(9),
     &                mcoll(8),mcoll(12),mcoll(13),mcoll(14)
      write(isum,'(f8.2,7f9.4)')
     &  (time2(i),
     &float(mcol(4,i))*xnum,float(mcol(5,i))*xnum,
     &float(mcol(9,i))*xnum,float(mcol(8,i))*xnum,
     &float(mcol(12,i))*xnum,float(mcol(13,i))*xnum,
     &float(mcol(14,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of collision for n*2'')')
      write(isum,'(
     &''c:time(fm/c) nn->nq nq->nn   nq->nd   nd->nq   nq->nr'',
     &''   nr->nq   dd->nq'')')
      write(isum,'(''c: total'',7i9)') mcoll(6),mcoll(7),mcoll(11),
     &                mcoll(10),mcoll(13),mcoll(12),mcoll(15)
      write(isum,'(f8.2,7f9.4)')
     &  (time2(i),
     &float(mcol(6,i))*xnum,float(mcol(7,i))*xnum,
     &float(mcol(11,i))*xnum,float(mcol(10,i))*xnum,
     &float(mcol(13,i))*xnum,float(mcol(12,i))*xnum,
     &float(mcol(15,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:sources for n*2'')')
      write(isum,'(
     &''c:time(fm/c) nn       bb        pi     meson'')')
      write(isum,'(''c: total'',4i9)') mcoll(6)-mcoll(7),
     &          mcoll(10)+mcoll(12)+mcoll(15)-mcoll(11)-mcoll(13),
     & mmesa(3)+mmesa(5)+mmesa(6)-mmesc(3)-mmesc(5)-mmesc(6),
     & mmesa(9)+mmesa(11)-mmesc(9)-mmesc(11)
      write(isum,'(f8.2,4f9.4)')
     &  (time2(i),
     &float(mcol(6,i)-mcol(7,i))*xnum,
     &float(mcol(10,i)+mcol(12,i)+mcol(15,i)-mcol(11,i)-mcol(13,i))*
     &xnum,
     &float(mesa(3,i)+mesa(5,i)+mesa(6,i)-mesc(3,i)-mesc(5,i)-mesc(6,i))
     &*xnum,
     &float(mesa(9,i)+mesa(11,i)-mesc(9,i)-mesc(11,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of collision for pions''/
     &''n: time (fm/c)''/''n: collisions'')')
      write(isum,'(
     &''c:time d->n+pi n+pi->d x->n+pi n+pi->x'',
     &           '' x->x+pi x+pi->x ro->2pi 2pi->ro'',
     &           '' si->2pi 2pi->si om>pi+ro pi+ro>om'')')
      write(isum,'(''c: tot.'',12i8)') mmesc(1),mmesa(1),
     &mmesc(2)+mmesc(3),mmesa(2)+mmesa(3),mmesc(4)+mmesc(5)+mmesc(6),
     &mmesa(4)+mmesa(5)+mmesa(6),mmesc(12),mmesa(12),mmesc(13),mmesa(13)
     &,mmesc(14),mmesa(14)
      write(isum,'(f8.2,12f8.4)')
     &  (time2(i),
     &float(mesc(1,i))*xnum,float(mesa(1,i))*xnum,
     &float(mesc(2,i)+mesc(3,i))*xnum,float(mesa(2,i)+mesa(3,i))*xnum,
     &float(mesc(4,i)+mesc(5,i)+mesc(6,i))*xnum,
     &float(mesa(4,i)+mesa(5,i)+mesa(6,i))*xnum,
     &float(mesc(12,i))*xnum,float(mesa(12,i))*xnum,
     &float(mesc(13,i))*xnum,float(mesa(13,i))*xnum,
     &float(mesc(14,i))*xnum,float(mesa(14,i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of baryons''/
     &''n: time (fm/c)''/
     &''n: particles'')')
      write(isum,'(
     &''c:time    proton neutron delta++  delta+  delta0  '',
     &''delta-  nstar+  nstar0 nstar2+ nstar20'')')
      write(isum,'(f7.2,1x,10f8.3)')
     &  (time2(i),
     &float(mprt(i))*xnum,float(mntr(i))*xnum,float(mde2(i))*xnum,
     &float(mde1(i))*xnum,float(mde0(i))*xnum,float(mdem(i))*xnum,
     &float(mrpl(i))*xnum,float(mrze(i))*xnum,
     &float(mqpl(i))*xnum,float(mqze(i))*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of mesons''/
     &''n: time (fm/c)''/
     &''n: particles'')')
      write(isum,'(
     &''c:time    pi+     pi0     pi-     eta     '',
     &''ro+     ro0     ro-     sigma   omega'')')
      write(isum,'(''c: total'',i5,8i8)') (mumes(ij,ictime),ij=1,9)
      write(isum,'(''c: final'',i5,8i8)') (lumes(ij),ij=1,9)
      write(isum,'(f6.2,9f8.4)')
     &  (time2(i),
     &(float(mumes(ij,i))*xnum,ij=1,9),
     & i=1,ictime)
*
      write(isum,'(/''c:number of particles in central cells''/
     &''n: time (fm/c)''/
     &''n: particles'')')
      write(isum,'(
     &''c:time   nucleon   delta  nstar1  nstar2    pion '',
     &''    eta     rho   sigma res/nuc density energyd'')')
      write(isum,'(f7.2,1x,11f8.4)')
     &  (time2(i),
     &float(mcnu(i))*xnum,float(mcde(i))*xnum,float(mcrn(i))*xnum,
     &float(mcqn(i))*xnum,(float(mumes(ij,i))*xnum,ij=11,14),
     &rati(i),suru(i),ener(i)*xnum,
     & i=1,ictime)
*
      write(isum,'(/''c:number of collision/time in general''/
     &''n: time (fm/c)''/''n: collisions/time'')')
      write(isum,'(''c:time(fm/c)'',
     &'' all blocked elastic proj+targ'',
     &'' nn->nx nx->nn nr->nr1 d+d->nr nx->nxel xx->xxel'')')
      write(isum,'(''c:total'',10i8)')
     & mcoll(-1),mcoll(0),mcoll(1),mcopt,
     & mcoll(2)+mcoll(4)+mcoll(6),
     & mcoll(3)+mcoll(5)+mcoll(7),
     & mcoll(8)+mcoll(9)+mcoll(10)+mcoll(11)+mcoll(12)+mcoll(13),
     & mcoll(14)+mcoll(15),mcond,mcodd
      write(isum,'(f7.2,1x,10f8.4)')
     &  (time2(i),
     &float(mcol(-1,i)-mcol(-1,i-1))/denu,
     &float(mcol(0,i)-mcol(0,i-1))/denu,
     &float(mcol(1,i)-mcol(1,i-1))/denu,
     &float(mept(i)-mept(i-1))/denu,
     &float(mcol(2,i)+mcol(4,i)+mcol(6,i)
     &     -mcol(2,i-1)-mcol(4,i-1)-mcol(6,i-1))/denu,
     &float(mcol(3,i)+mcol(5,i)+mcol(7,i)
     &     -mcol(3,i-1)-mcol(5,i-1)-mcol(7,i-1))/denu,
     &float(mcol(8,i)+mcol(9,i)+mcol(10,i)+mcol(11,i)+mcol(12,i)
     &  -mcol(8,i-1)-mcol(9,i-1)-mcol(10,i-1)-mcol(11,i-1)-mcol(12,i-1)
     &     +mcol(13,i)-mcol(13,i-1))/denu,
     &float(mcol(14,i)+mcol(15,i)-mcol(14,i-1)-mcol(15,i-1))/denu,
     &float(mend(i)-mend(i-1))/denu,float(medd(i)-medd(i-1))/denu,
     & i=1,ictime)
*
      write(isum,'(/''c:number of collision/time for n*2'')')
      write(isum,'(
     &''c:time(fm/c) nx->nq  nq->nx  xpi->q   q->xpi'',
     &''   xme->q   q->xme   net->q   q->net'')')
      write(isum,'(''c: total'',8i9)') mcoll(6)+mcoll(10)+mcoll(15)+
     &                mcoll(12),mcoll(7)+mcoll(11)+mcoll(13),
     & mmesa(3)+mmesa(5)+mmesa(6),mmesc(3)+mmesc(5)+mmesc(6),
     & mmesa(9)+mmesa(11),mmesc(9)+mmesc(11),mmesa(7),mmesc(7)
      write(isum,'(f8.4,8f9.3)')
     &  (time2(i),
     &float(mcol(6,i)+mcol(10,i)+mcol(15,i)+mcol(12,i)
     &-mcol(6,i-1)-mcol(10,i-1)-mcol(12,i-1)-mcol(15,i-1))/denu,
     &float(mcol(7,i)+mcol(11,i)+mcol(13,i)
     &-mcol(7,i-1)-mcol(11,i-1)-mcol(13,i-1))/denu,
     &float(mesa(3,i)+mesa(5,i)+mesa(6,i)
     &-mesa(3,i-1)-mesa(5,i-1)-mesa(6,i-1))/denu,
     &float(mesc(3,i)+mesc(5,i)+mesc(6,i)
     &-mesc(3,i-1)-mesc(5,i-1)-mesc(6,i-1))/denu,
     &float(mesa(9,i)+mesa(11,i)-mesa(9,i-1)-mesa(11,i-1))/denu,
     &float(mesc(9,i)+mesc(11,i)-mesc(9,i-1)-mesc(11,i-1))/denu,
     &float(mesa(7,i)-mesa(7,i-1))/denu,
     &float(mesc(7,i)-mesc(7,i-1))/denu,
     & i=1,ictime)
*
      write(isum,'(/''c:sources/time for n*2'')')
      write(isum,'(
     &''c:time(fm/c) nn coll bb coll  pions   mesons'')')
      write(isum,'(''c: total'',4i9)') mcoll(6)-mcoll(7),
     &mcoll(10)+mcoll(15)+mcoll(12)-mcoll(11)-mcoll(13),
     & mmesa(3)+mmesa(5)+mmesa(6)-mmesc(3)-mmesc(5)-mmesc(6),
     & mmesa(9)+mmesa(11)-mmesc(9)-mmesc(11)
      write(isum,'(f8.2,4f9.3)')
     &  (time2(i),
     &float(mcol(6,i)-mcol(7,i)-mcol(6,i-1)+mcol(7,i-1))/denu,
     &float(mcol(10,i)+mcol(15,i)+mcol(12,i)-mcol(11,i)-mcol(13,i)
     &-mcol(10,i-1)-mcol(12,i-1)-mcol(15,i-1)+mcol(11,i-1)
     &+mcol(13,i-1))/denu,
     &float(mesa(3,i)+mesa(5,i)+mesa(6,i)
     &-mesa(3,i-1)-mesa(5,i-1)-mesa(6,i-1)-mesc(3,i)-mesc(5,i)-mesc(6,i)
     &+mesc(3,i-1)+mesc(5,i-1)+mesc(6,i-1))/denu,
     &float(mesa(9,i)+mesa(11,i)-mesa(9,i-1)-mesa(11,i-1)
     &-mesc(9,i)-mesc(11,i)+mesc(9,i-1)+mesc(11,i-1))/denu,
     & i=1,ictime)
*
      write(isum,'(/''c:plot of distribution of collision number''/
     &''n: time (fm/c)''/''n: distr. coll. num.''/
     &''c: time c(0)  c(1)  c(2)  c(3)  c(4)  c(5)  c(6)  c(7)  '',
     &''c(8)  c(9)  c(10) c(11) c(12) c(13) c(14) c(15) c(16) c(17) '',
     &''c(18) c(19)'')')
      write(isum,'(
     &''n: x     y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  y,l0 '',
     &'' y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  y,l0  '',
     &''y,l0  y,l0'')')
      write(isum,'(f7.2,20f6.2)')
     & (time2(i),(float(ndis(j,i))*xnum,j=0,19),i=0,ictime)
*
      end if
*
*=======================================================================
      call resdyout(xnum)
*
      write(isum,'(/''c:meson coll. number'')')
      write(isum,'(''c: # pion -> meson:'',i8)') pitomes
      write(isum,'(''c: # meson -> pion:'',i8)') mestopi
      write(isum,'(''c: # meson -> meson:'',i8)') mestomes
      write(isum,'(''c: total number of mesons:'',i8)')
     & mcoll(2)+mcoll(4)+mcoll(6)-mcoll(3)-mcoll(5)-mcoll(7)
     &-mcoll(14)-mcoll(15)-mdpir-mdpiq-mrpiq
      write(isum,'(''n: birth number'')')
      write(isum,'(''n: number of mesons'')')
      write(isum,'(''n:   x       y(pion),l0      y(eta),m0   '',
     &                       ''y(pion f.),d0    y(eta f.),p0'')')
      write(isum,'(i6,4i14)')
     & (i,(mlife(j,i),j=1,4), i=0,25)
*
      write(isum,'(/''c:meson free path'')')
      write(isum,'(''c: total number of mesons:'',i8)')
     & mcoll(2)+mcoll(4)+mcoll(6)-mcoll(3)-mcoll(5)-mcoll(7)
     &-mcoll(14)-mcoll(15)-mdpir-mdpiq-mrpiq
      write(isum,'(''n: free path(fm)'')')
      write(isum,'(''n: number of mesons'')')
      write(isum,'(''n:   x       y(pion),l0      y(eta),m0'')')
      write(isum,'(f6.2,2i14)')
     & (0.5*float(i),(mpath(j,i),j=1,2), i=0,30)
*
*
      write(isum,'(/''c:density dep. of meson creation'')')
      write(isum,'(''c: total number of mesons:'',i8)')
     & mcoll(2)+mcoll(4)+mcoll(6)-mcoll(3)-mcoll(5)-mcoll(7)
     &-mcoll(14)-mcoll(15)-mdpir-mdpiq-mrpiq
      write(isum,'(''n: density(rho0)'')')
      write(isum,'(''n: number of mesons'')')
      write(isum,'(''n:   x       y(pion),l0      y(eta),m0   '',
     &                       ''y(pion f.),d0    y(eta f.),p0'')')
      write(isum,'(f6.2,4i14)')
     & (float(i)*0.2,(mdens(j,i),j=1,4), i=0,20)
*
      write(isum,'(/''c:time dist of final meson creation'')')
      write(isum,'(''c: total number of mesons:'',i8)')
     & mcoll(2)+mcoll(4)+mcoll(6)-mcoll(3)-mcoll(5)-mcoll(7)
     &-mcoll(14)-mcoll(15)-mdpir-mdpiq-mrpiq
      write(isum,'(''n: time(fm/c)'')')
      write(isum,'(''n: number of mesons'')')
      write(isum,'(''n:   x       y(pion),l0      y(eta),m0'')')
      write(isum,'(i6,2i14)')
     & (i,(mbirt(j,i),j=1,2), i=1,50)
*
      write(isum,'(/''c:position of final meson creation'')')
      write(isum,'(''c: max. radius of the proj and targ.'',f8.3)')
     & radius
      write(isum,'(''n: distance(fm)'')')
      write(isum,'(''n: number of mesons'')')
      write(isum,'(''n:   x       y(pion),l0      y(eta),m0'')')
      write(isum,'(f6.2,2i14)')
     & (0.5*float(i),(mposi(j,i),j=1,2), i=1,20)
*
      write(isum,'(/''c:pauli blocking of delta decay '')')
      write(isum,'(''n: mass (gev)'')')
      write(isum,'(''n: blocked phase space distribution'')')
      write(isum,'(''n:   intervall  testparticle    localdensity'')')
      write(isum,'(f13.3,i10,i15)')
     & (float(i)/10.+0.05,iphapt(i),iphapc(i), i=0,9)
*
      write(isum,'(/''c:density dep. resonance life'')')
      write(isum,'(''n: density(rho0)'')')
      write(isum,'(''n: number of resonances'')')
      write(isum,'(''n:   x       y(nucl),l0      y(del),m0   '',
     &                       ''y(n1440),d0      y(n1535),p0'')')
      write(isum,'(f6.2,4f14.2)')
     & (float(i-1)*0.2,(resdens(j,0,i),j=1,4), i=1,18)


      return

      end


