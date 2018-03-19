
************************************************************************
*                                                                      *
      subroutine inecho(io    , massta, masspr, elab  , isubs,
     &            iseed , dt    , ntmax , icoll , num   , insys ,
     &            iavoid, iwidth,ipou  , mstapr,mstanu,msprpr,msprnu,
     &            idipi,iabso,iempi,ipiab, ipico,irescol, ipauli,dlife)
*                                                                      *
*     purpose:   echo of input                                         *
*                                                                      *
************************************************************************
*
      implicit none
      real*8 dlife, elab, dt
      integer io,massta,masspr,isubs,iseed,ntmax,icoll,num,insys,iavoid
      integer iwidth,ipou,mstapr,mstanu,msprpr,msprnu
      integer idipi,iabso,iempi,ipiab,irescol
      integer ipico,idelpot,ipauli
      include"common"
      write(io,666) massta,mstapr,mstanu,masspr,msprpr,msprnu
 666  format(/'c:',17x, 'input-echo for buu-program:' /
     &                    'c:',17x, '===========================' /
     &                    'c:',17x, 'mass of target:     ' ,i7/
     &                    'c:',17x, '        proton:     ' ,i7/
     &                    'c:',17x, '        neutron:    ' ,i7/
     &                    'c:',17x, 'mass of projectile: ' ,i7/
     &                    'c:',17x, '        proton:     ' ,i7/
     &                    'c:',17x, '        neutron:    ' ,i7/)
      write(io,667) dt,ntmax
 667  format(
     &                    'c:',17x, 'time step    "fm/c":' ,f7.2/
     &                    'c:',17x, 'number of t-steps:  ' ,i7)
      write(io,668) isubs,num,iseed
 668  format(
     &                    'c:',17x, 'subsequent runs.: ' ,i9/
     &                    'c:',17x, 'testparticle/nucl.: ' ,i7/
     &                    'c:',17x, 'iseed:       ' ,i14/)
      write(io,669)
 669  format(
     &                    'c:',17x, 'power of den.dep.int.=',f8.5/
     &                    'c:',17x, 't0  =',f8.1,'    t3  =',f8.1/
     &                    'c:',17x, 'smu =',f8.3,'    yv  =',f8.1/
     &                    'c:',17x, 'aaa =',f8.1,'    bbb =',f8.1/
     &                    'c:',17x, 'ccc =',f8.1/
     &                    'c:',17x, 'k   =',f8.1/)
      write(io,'(''c:'',15x,''================================'',
     &                      ''======='')')
      if (icoll .eq. 1) then
        write(io,'(''c:'',17x,''mean field only!'')')
      else if (icoll .eq. -1) then
        write(io,'(''c:'',17x,''cascade only!'')')
      else
        write(io,'(''c:'',17x,''full buu-theory used!'')')
      end if
c      if (idelpot.eq.0)
c     &  write(io,'(''c:'',17x,''same pot. for delta as for nucl.'')')
c      if (idelpot.eq.1)
c     &  write(io,'(''c:'',17x,''prop. pot. for delta to nucl.'')')
c      if (idelpot.eq.2)
c     &  write(io,'(''c:'',17x,''linear pot. for delta'')')
      if (iavoid .eq. 1) then
        write(io,'(''c:'',
     &             17x,''unimportant n-n collisions avoided.'')')
      else
        write(io,'(''c:'',
     &             17x,''all n-n collisions are allowed.'')')
      end if
      if(irescol.eq.0)
     &  write(io,'(''c:'',17x,''no nx->nx and dd->nx reactions'')')
      if(irescol.eq.1)
     &  write(io,'(''c:'',17x,''nx->nx and dd->nx reactions'')')
      if(iwidth.eq.0)
     &  write(io,'(''c:'',17x,''breit-wigner mass dist. for delta'')')
      if(iwidth.eq.1)
     &  write(io,'(''c:'',17x,''moniz width for for delta mass'')')
      if(iwidth.eq.2)
     &  write(io,'(''c:'',17x,''kitazoe width for for delta mass'')')
      if(idipi.eq.0)
     &  write(io,'(''c:'',17x,''no pions, only barions'')')
      if(idipi.eq.1)
     &  write(io,'(''c:'',17x,''pion emission only at the end'')')
      if(idipi.eq.2)
     &  write(io,'(''c:'',17x,''direct pions'')')
      if(iabso.eq.0)
     &  write(io,'(''c:'',17x,''no pion absorption'')')
      if(iabso.eq.1)
     &  write(io,'(''c:'',17x,''pions are absorbed (normal run)'')')
      if(iempi.eq.0)
     &  write(io,'(''c:'',17x,''no pion emission'')')
      if(iempi.eq.1)
     &  write(io,'(''c:'',17x,''pions are emitted (normal run)'')')
      if(ipiab.eq.0)
     &  write(io,'(''c:'',17x,''b-w pion absorption'')')
      if(ipiab.eq.1)
     &  write(io,'(''c:'',17x,''1/k**2 pion absorption'')')
      if(ipico.eq.0)
     &  write(io,'(''c:'',17x,''no pion-pion collision'')')
      if(ipico.eq.1)
     &  write(io,'(''c:'',17x,''pion-pion collision'')')
c      if(idnfac.eq.0)
c     &  write(io,'(''c:'',17x,''normal xsection for dn->nn'')')
c      if(idnfac.eq.1)
c     &  write(io,'(''c:'',17x,''sigma(dn->nn)=sigma0(dn->nn)*'',f5.2)')
c     &                                                         dnfac
c      if(idnfac.eq.2)
c     &  write(io,'(''c:'',17x,''sigma(dn->nn)=sigma0(dn->nn)*'',
c     &             ''(1+'',f5.2,''*rho/rho0)'')')         dnfac
      if(abs(dlife-1.0) .gt. 1.0e-4)
     &  write(io,'(''c:'',17x,''deltalife = deltalife *'',f6.3)')  dlife
      if (ipauli.eq.1) then
        write(io,'(''c:'',17x,''pauli testparticle in delta decay '')')
      end if
c      if(ggam .lt. 0.01)
c     &  write(io,'(''c:'',17x,''pauli correction in resonance width'')')
c      if(abs(ggam-1.0) .lt. 0.01)
c     &  write(io,'(''c:'',17x,''no pauli corr. in resonance width'')')
c      if( (abs(ggam-1.0) .gt. 0.01) .and. (abs(ggam) .gt. 0.01) )
c     &  write(io,'(''c:'',17x,''resonance width corr. by '',f6.2)') ggam
      write(io,'(''c:'',15x,''================================'',
     &                      ''======='')')
      return
      end
