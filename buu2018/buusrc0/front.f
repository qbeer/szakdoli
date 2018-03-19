
************************************************************************
*                                                                      *
      subroutine front(io,massta,masspr,elab,mstapr,mstanu,msprpr,
     &                 msprnu,mpion,b)
*                                                                      *
*       purpose:    writing front page of output                       *
*       variables:  io     - output-unit               (integer,input) *
*                   elab   - beam energy "gev/a"          (real,input) *
*                                                                      *
************************************************************************
*
      implicit none
      integer io,massta,masspr,mstapr,mstanu,msprpr,msprnu,mpion
      real*8 elab,b
      character *4 pich
      include"common"
      write(io,'(/''c:'',17x,''bbbbbbbb     uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbbbbbbbb    uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbb    bbb   uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbb    bbb   uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbbbbbbb     uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbbbbbbb     uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbb    bbb   uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbb    bbb   uuu    uuu   uuu    uuu'')')
      write(io,'(''c:'',17x,''bbbbbbbbb     uuuuuuuu     uuuuuuuu'')')
      write(io,'(''c:'',17x,''bbbbbbbb       uuuuuu       uuuuuu'')')
*
      write(io,'(//''c:'',67(''*'')/''c:*'',65x,''*''/''c:'',
     &             ''*'',10x,''reaction:'',46x,''*''/''c:'',
     &             ''*'',65x,''*'')')
      if(mpion .eq. 0) then
      write(io,'(''c:*'',14x,''mass'',i4,''('',i3,'','',i3,
     &             '') ==> mass'',i4,''('',i3,'','',i3,'')'',12x,''*''/
     &             ''c:*'',65x,''*''/''c:*'',15x,
     &             ''beam energy = '',f9.3,'' gev/nucleon'',15x,''*'')')
     &      masspr,msprpr,msprnu,massta,mstapr,mstanu, elab
      else if(mpion .ne. 0) then
      if(mpion .eq. 1) pich = 'pi+ '
      if(mpion .eq.-1) pich = 'pi- '
      if(mpion .eq. 2) pich = 'pi0 '
      write(io,'(''c:*'',20x,a4,
     &             '' ==> mass'',i4,''('',i3,'','',i3,'')'',19x,''*''/
     &             ''c:*'',65x,''*''/''c:*'',19x,
     &             ''beam energy = '',f9.3,'' gev'',19x,''*'')')
     &      pich,massta,mstapr,mstanu, elab
      end if
      write(io,'(''c:*'',15x,''impact parameter ='',f6.2,'' fm'',20x,
     &    ''*''/''c:*'',65x,''*''/''c:'',67(''*'')/)') b
      return
      end
