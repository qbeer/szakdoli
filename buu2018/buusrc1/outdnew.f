c**********************************************************************c
      subroutine outdat(isu,nt, lcoll, lmesc, lmesa, lmesa2,cres )

      implicit none

      include"common"
      include"coucom"
      include"cominput"

      common /nameresc/nameres(1:nres+3)
             character *8  nameres

      real*8     rboxdelta
      parameter(rboxdelta=1.0)
      integer  nmaxdim
      parameter(nmaxdim = 200 )
      integer  isu, nt, ionline
      parameter(ionline = 1)
      integer  i, repro, j, jw, jj
      real*8     resmass
      character *80 outfile

      character *2  nst1, de1, fi, lam, sigm
      character *1  endst

      integer dimlmesc
      parameter(dimlmesc = nres+6)
      integer lmesc(1:dimlmesc)
      integer lmesa(1:dimlmesc)
      integer lmescst(1:nmaxdim,1:dimlmesc)

      integer lcodim
      parameter(lcodim = 3*nres+6+nres**2)
      integer lcoll(-6:lcodim)
      integer lcollst(1:nmaxdim,-6:lcodim)

      integer lmesast(1:nmaxdim,1:dimlmesc)
      integer nuba(1:nres+3), nubast(1:nmaxdim,1:nres+3)
      integer numes(1:18), numesst(1:nmaxdim,1:18)
      integer nin, nout, ii
      integer elast(1:4), selast(1:nmaxdim,1:4)
      integer nrnn(1:5), snrnn(1:nmaxdim,1:5)
      integer nnnr(1:5), snnnr(1:nmaxdim,1:5)
      integer nndd, snndd(1:nmaxdim), ddnn, sddnn(1:nmaxdim)
      integer nrnr(1:7), snrnr(1:nmaxdim,1:7)
      integer ssta, sssta(1:nmaxdim)
      real*8    ncdelta(1:nmaxdim)
      real*8    w(1:5), test
      integer deltamass(0:38), ind , lmesa2, lmesa2st(1:nmaxdim)
      integer cres(1:nres,1:9), cresst(1:nmaxdim,1:nres,1:9)
      real*8    cenden(1:nmaxdim,1:2), j0, j1, j2, j3
      real*8    betlrfx, betlrfy, betlrfz
      real*8    norm

************** save for output **********************************
      save    lcollst, nubast, numesst, selast, snrnn, snnnr
      save    sddnn, snndd, sssta, snrnr, norm,cresst,lmescst
      save    cenden, ncdelta
*****************************************************************



******************************************************************
********* normalization factor **********************************
      norm = 1.0/(float(isubs)*float(num))

      write(*,*)'in outdat '
      write(*,*)'nt = ',nt
      write(*,*)'mod = ',mod(nt,1)
      write(isum,*)'************************************'
      write(isum,*)'************************************'
      write(isum,*)' TIME STEP = ', nt

      if(nt.gt. nmaxdim) then
        write(*,*)'stop in oud '
        write(*,*)'ntmax too high !!!!!!!!!!!!!!'
        stop
      end if

      if(nt.eq.1) then
        write(isum,*)'*****************************************'
        call croswwrite
        write(isum,*)'*****************************************'
      end if


c      if(isu.le.isubs .and. nt.le.nmaxdim) then

*******************************************************************
************** output of density-profiles************************
c          write(*,*) denout
      if(idenspr.eq.1)   then

        if(nt.eq.1 .or. mod(nt,10).eq.0) then

c 10     format(i4,1x,5(e8.4,1x))


        if(nt.lt.10) then
          write(outfile,11) denout,'.00',nt
 11       format(A,A,I1)
        else if(nt.lt.100) then
          write(outfile,12) denout,'.0',nt
 12       format(A,A,I2)
        else
          write(outfile,13) denout,'.',nt
 13       format(A,A,I3)
        end if

        write(*,*)'nach 1'

        open(102,file=outfile,status='unknown')



        do i = -maxz, maxz
          do j = -maxx, maxx


            write(102,*)i,j, rhob_4(0,j,0,i)

          end do
        end do


        end if
        end if
***************************************************************
*      build strings for particle names


 100  format(A,F5.0,A)
      endst = ')'
      nst1 = 'N('
      de1  = 'D('
      lam  = 'L('
      sigm = 'S('
      do i = 1,nres+1
        if(i.eq.1) then
          repro = 1
          resmass = rmass*1000.
          fi = nst1
        else if(i.eq.nres+2) then
          resmass = xlmas*1000.
          fi = lam
        else if(i.eq.nres+3) then
          resmass = xsmas*1000.
          fi = sigm
        else
          repro = resprop2(i-1,1)
          resmass = resprop1(i-1,1)*1000.
          if(repro.eq.1) then
            fi = nst1
          else
            fi = de1
          end if
        end if

      write(nameres(i),100)fi,resmass,endst



      end do

******************************************************************
********** baryon - coll bookkeeping *****************************


**    set storeing fields to zero

      if(isu.eq.1 .and. nt. eq.1) then

        do i = 1, nmaxdim
          ncdelta(i) = 0.0
          cenden(i,1) = 0.0
          cenden(i,2) = 0.0
          lmesa2st(i) = 0
          do j = 1, dimlmesc
            lmescst(i,j) = 0
            lmesast(i,j) = 0
          end do
          do j = -6,lcodim
            lcollst(i,j) =0
          end do

          do j = 1,18
            numesst(i,j) =0
          end do

          do j = 1, nres+3
            nubast(i,j) =0
          end do

          sssta(i) = 0
          snndd(i) = 0
          sddnn(i) = 0
          do j = 1, 5
            snrnn(i,j) = 0
            snnnr(i,j) = 0
          end do
          do j = 1,4
            selast(i,j) = 0
          end do
          do j = 1,7
            snrnr(i,j) = 0
          end do
          do j = 1, nres
            do jj = 1,6
              cresst(i,j,jj) = 0
            end do
          end do
        end do



      end if


******** add up the numbers for baryon-baryon coll.
      do i = -6,lcodim
        lcollst(nt,i) = lcollst(nt,i) + lcoll(i)
      end do

      do j = 1, nres
        do i = 1,6
          cresst(nt,j,i) = cresst(nt,j,i)+ cres(j,i)
        end do
      end do

***** determine the number of mesons and baryons and store them *

      do i = 1, 18
        numes(i) = 0
      end do


      do  i = 1,maxppar
        if(ipi(1,i) .ne. 0)then
*
          if(ipi(1,i) .eq. 1) then                       ! pion
            if(ipi(2,i).eq. 1) numes(1) = numes(1) + 1
            if(ipi(2,i).eq. 0) numes(2) = numes(2) + 1
            if(ipi(2,i).eq.-1) numes(3) = numes(3) + 1
          end if
          if(ipi(1,i) .eq. 2)  numes(4) = numes(4) + 1   ! eta
          if(ipi(1,i) .eq. 3) then                       ! rho
            if(ipi(2,i).eq. 1) numes(5) = numes(5) + 1
            if(ipi(2,i).eq. 0) numes(6) = numes(6) + 1
            if(ipi(2,i).eq.-1) numes(7) = numes(7) + 1
          end if
          if(ipi(1,i) .eq. 4)  numes(8) = numes(8) + 1   ! sigma
          if(ipi(1,i) .eq. 5)  numes(9) = numes(9) + 1   ! omega
          if(ipi(1,i) .eq. 6)  numes(10) = numes(10) + 1 ! kaon
          if(rpi(1,i)**2+rpi(2,i)**2+rpi(3,i)**2 .lt. 4.0) then
            if(ipi(1,i).eq. 1) numes(11)  = numes(11) + 1
            if(ipi(1,i).eq. 2) numes(12)  = numes(12) + 1
            if(ipi(1,i).eq. 3) numes(13)  = numes(13) + 1
            if(ipi(1,i).eq. 4) numes(14)  = numes(14) + 1
            if(ipi(1,i).eq. 5) numes(15)  = numes(15) + 1
            if(ipi(1,i).eq. 6) numes(16)  = numes(16) + 1
          end if

        end if



        if(ipi(1,i).gt.6 .or. ipi(1,i).lt.0 .or. abs(ipi(2,i)).gt.1 .or.
     &    ( (ipi(1,i).eq.2.or.ipi(1,i).eq.4.or.ipi(1,i).eq.5)
     &       .and.ipi(2,i).ne.0 ) )
     &        write(*,*) 'fehler in outdat mesno  ',
     &            ipi(1,i),ipi(2,i),ipi(5,i)


      end do



      do i = 1, nres+3
        nuba(i) = 0
      end do

      do i = 1, maxpar
        if(id(1,i).ne.0) then
          nuba(id(1,i)) = nuba(id(1,i))+1
        end if
      end do


************ count Deltas in central cell rboxdelta
      do i = 1, maxpar
        if(id(1,i).ne.0) then
          if(abs(r(1,i)).le.rboxdelta) then
            if(abs(r(2,i)).le.rboxdelta) then
              if(abs(r(3,i)).le.rboxdelta) then
                ncdelta(nt) = ncdelta(nt) + 1.0
              end if
            end if
          end if
        end if
      end do
***********************************************************



      do i = 1, 18
        numesst(nt,i) = numesst(nt,i)+ numes(i)
      end do


      do i = 1, nres+3
        nubast(nt,i) = nubast(nt,i)+ nuba(i)
      end do





****************** ONLINE - OUTPUT, if selected *****************
      if(ionline.eq.1) then

        write(isum,*)' ******* resonnance statistics *********'
        write(isum,*)' resonance         id         # '
        do i = 1,nres+1
          write(isum,*)nameres(i),i,nuba(i)
        end do

*
        write(isum,*)'****************************'
        write(isum,*)'baryon collision numbers '
        write(isum,*)'****************************'

        nin = 2
        nout = 1
        do i = 1, lcodim
         if(i.le. (nres+1)) then
      write(isum,*) lcoll(i) ,'   ',nameres(1), ' + ',nameres(i),
     &              ' --  '  ,nameres(1), ' + ',nameres(i),'el'

       else if(i.gt.(nres+1).and. i.le.(2*nres+1)) then
       write(isum,*)lcoll(i) ,'    ',nameres(1), ' + ',nameres(i-nres),
     &                   ' --  ',nameres(1),' + ',nameres(1)

       else if(i.gt.(2*nres+1).and. i.le.(3*nres+1)) then
       write(isum,*)lcoll(i),'     ',nameres(1),' + ',nameres(1),' -- '
     &              ,nameres(1),' + ',nameres(i-2*nres)

       else if(i.gt.(3*nres+1).and.i.le.(nres**2+3*nres+1)) then

         if(nout.gt.2.and. mod(nout,nres).eq.1) then
           nout = 1
           nin  = nin + 1
         end if
         nout = nout + 1
      write(isum,*)lcoll(i),'  ',nameres(1),' + ',nameres(nin),' --  ',
     &          nameres(1),' + ' ,nameres(nout)

       else if(i.eq.(nres**2+3*nres+2))then
        write(isum,*)lcoll(i),'  ',nameres(1),' + ',nameres(1),' -- ',
     &                          nameres(2),' + ',nameres(2)

       else if(i.eq.(nres**2+3*nres+3))then
        write(isum,*)lcoll(i),'  ',nameres(2),' + ',nameres(2),' -- ',
     &                          nameres(1),' + ',nameres(1)


       else if(i.eq.(nres**2+3*nres+4))then
                 write(isum,*)lcoll(i),'    s-state pion prod.'

         end if
       end do


        write(isum,*)' ******* meson-statistics *********'
        write(isum,*)'    id         # '
        write(isum,*)' pions : ', numes(1)+numes(2)+numes(3),numes(1)
     &                           , numes(2), numes(3)
        write(isum,*)'etas :   ',   numes(4)
        write(isum,*)'rhos :   ', numes(5)+numes(6)+numes(7),numes(5),
     &                         numes(6), numes(7)
        write(isum,*)'sigmas : ', numes(8)



        write(isum,*)' **********************************'
        write(isum,*)' *******      decays      *********'
        write(isum,*)
        write(isum,*)' res - mes + something '
        do i = 1, nres
          write(isum,*) nameres(i+1), lmesc(i)
          write(isum,*) (cres(i,j),j=1,6)

        end do



        write(isum,*)'************* pi absorption  '
        write(isum,*)' N pi       - D(1232) :', lmesa(1)
        write(isum,*)' N pi       - N(1440) :', lmesa(2)
        write(isum,*)' N pi       - N(1544) :', lmesa(3)
        write(isum,*)' D(1232) pi - N(1440) :', lmesa(4)
        write(isum,*)' D(1232) pi - N(1544) :', lmesa(5)
        write(isum,*)' N(1440) pi - N(1544) :', lmesa(6)
        write(isum,*)' N eta      - N(1544) :', lmesa(7)
        write(isum,*)' N rho      - N(1440) :', lmesa(8)
        write(isum,*)' N rho      - N(1544) :', lmesa(9)
        write(isum,*)' N sigma    - N(1440) :', lmesa(10)
        write(isum,*)' N sigma    - N(1544) :', lmesa(11)
        write(isum,*)' N meson    - Reson   :', lmesa(12)
        write(isum,*)' Reson meson- Reson   :', lmesa(13)
        write(isum,*)' sstate abs           :', lmesa2
*********************************************************
************ deltamsse

      do i = 0, 38
        deltamass(i) = 0
      end do


      do i = 1, maxpar
        if(id(1,i).eq.2) then
          ind = nint( (e(i)- 1.05)/0.025 )
          if(ind.ge.39) ind = 38
           deltamass(ind) = deltamass(ind) + 1
        end if
      end do


      write(isum,*)'********************************************'
      write(isum,*)'******* DELTA MASSES **********************'
      do i = 0, 38
           test = float(i)*0.025 + 1.05
        write(isum,*) test, deltamass(i)
      end do

      write(isum,*)'********************************************'

***    end of ionline if statement
       end if

*************************** end  of ONLINE OUTPUT *********************

*********** s e l e c t i v e  s u m a t i o n  o f  b a r y o n s *****

       do i = 1, 4
         elast(i) = 0
       end do
       do i = 1,5
         nrnn(i) = 0
         nnnr(i) = 0
       end do
         nndd = 0
         ddnn = 0
         ssta = 0
       do i = 1, 7
         nrnr(i) = 0
       end do


       nin = 2
       nout = 1
       do i = 1, lcodim
         if(i.le. (nres+1)) then
*        elast colls
           if(i.eq.1) elast(1) = lcoll(1)            ! NN
           if(i.eq.2) elast(2) = lcoll(2)            ! ND
           if(i.gt.2) elast(3) = elast(3) + lcoll(i) ! NR wo. ND
                      elast(4) = elast(4) + lcoll(i) ! total

         else if(i.gt.(nres+1).and. i.le.(2*nres+1)) then
*        N R - N N
           if(i.eq.nres+2) nrnn(1) = lcoll(i)        ! ND(1232)
           if(i.eq.nres+3) nrnn(2) = lcoll(i)        ! NN(1440)
           if(i.eq.nres+4) nrnn(3) = lcoll(i)        ! NN(1535)
           if(i.gt.nres+4) nrnn(4) = nrnn(4) + lcoll(i) !NR > 1535
                           nrnn(5) = nrnn(5) + lcoll(i) !NR totals

         else if(i.gt.(2*nres+1).and. i.le.(3*nres+1)) then
*        N N - N R
           if(i.eq.2*nres+2) nnnr(1) = lcoll(i)      ! ND(1232)
           if(i.eq.2*nres+3) nnnr(2) = lcoll(i)      ! NN(1440)
           if(i.eq.2*nres+4) nnnr(3) = lcoll(i)      ! NN(1535)
           if(i.gt.2*nres+4) nnnr(4) = nnnr(4)+ lcoll(i) !NR > 1535
                             nnnr(5) = nnnr(5)+ lcoll(i) !NR totals

         else if(i.gt.(3*nres+1).and.i.le.(nres**2+3*nres+1)) then
*        N R - N R'

           if(nout.gt.2.and. mod(nout,nres).eq.1) then
             nout = 1
             nin  = nin + 1
           end if
           nout = nout + 1
           if(nin.eq.2)  nrnr(1) = nrnr(1) + lcoll(i)
           if(nout.eq.2) nrnr(2) = nrnr(2) + lcoll(i)
           if(nin.eq.3)  nrnr(3) = nrnr(3) + lcoll(i)
           if(nout.eq.3) nrnr(4) = nrnr(4) + lcoll(i)
           if(nin.eq.4)  nrnr(5) = nrnr(5) + lcoll(i)
           if(nout.eq.4) nrnr(6) = nrnr(6) + lcoll(i)
                         nrnr(7) = nrnr(7) + lcoll(i)


         else if(i.eq.(nres**2+3*nres+2))then
*        N N - D D
           nndd = lcoll(nres**2+3*nres+2)

         else if(i.eq.(nres**2+3*nres+3))then
*        D D - N N
           ddnn = lcoll(nres**2+3*nres+3)

         else if(i.eq.(nres**2+3*nres+4))then
*         s-state pions
           ssta = lcoll(nres**2+3*nres+4)

         end if
       end do
*********************
*     density in central cell

       j0 = rhob_4(0,0,0,0)
       j1 = rhob_4(1,0,0,0)
       j2 = rhob_4(2,0,0,0)
       j3 = rhob_4(3,0,0,0)


       if(j0 .gt. 1.0e-6) then
         betlrfx = j1/j0
         betlrfy = j2/j0
         betlrfz = j3/j0
       else
c              write(*,*)'warning from potcalc j0 = ', j0
         betlrfx = 0.0
         betlrfy = 0.0
         betlrfz = 0.0
       end if

       cenden(nt,2) = j0
*         det. density in LRF
              if(j1**2+j2**2+j3**2.gt.j0**2) then
                 write(*,*) "hiba outdnew lorentz, mass<0",j0,j1,j2,j3
c                 stop
              end if
          call lorentz(betlrfx, betlrfy, betlrfz, j1, j2, j3, j0)

          cenden(nt,1) = j0

*********************************************************************
************    s u m a t i o n  o f  m e s o n s  ******************

       do i = 1, dimlmesc
         lmescst(nt,i) = lmesc(i)
       end do
       lmesa2st(nt) = lmesa2
       do i = 1, 13
         lmesast(nt,i) = lmesa(i)
       end do

***********************************************************************
***********************************************************************
       return

      entry outdout(lcoll,cres,lmesc)
******************* output of stuff ******************

      write(isum,*)'************************************************'
      write(isum,*)'********** NUMBER OF BARYONS *******************'
      write(isum,*)'************************************************'

      write(isum,*)'--------------- LEGEND--------------------------'
      write(isum,*)' resonance         id         # '
        do i = 1,nres+1
          write(isum,*)nameres(i),i,nuba(i)
        end do
      write(isum,*)
 201  format(6(1x,f8.3))
 202  format(5(1x,f8.3))
 203  format(3(1x,f8.3))
      write(isum,*)'time   numbers '
      write(isum,*)'time  1       2        3        4        5    '
      do i = 1, ntmax
        do j = 1,5
          jw = j
          w(jw)= float(nubast(i,j))*norm
        end do
      write(isum,201)float(i)*dt,w(1),w(2),w(3),w(4),w(5)
      end do

      write(isum,*)'time   numbers '
      write(isum,*)'time    6     7      8        9      10 '
      do i = 1, ntmax
        do j = 1,5
          jw = j + 5
          w(j)= float(nubast(i,jw))*norm
        end do
      write(isum,201)float(i)*dt,w(1),w(2),w(3),w(4),w(5)
      end do

      write(isum,*)'time   numbers '
      write(isum,*)'time    11     12       13     14  15 '
      do i = 1, ntmax
        do j = 1,5
          jw = j + 10
          w(j)= float(nubast(i,jw))*norm
        end do
      write(isum,201)float(i)*dt,w(1),w(2),w(3),w(4),w(5)
      end do

      write(isum,*)'************************************************'
      write(isum,*)'********** NUMBER OF MESONS  *******************'
      write(isum,*)'************************************************'

      write(isum,*)'time   tot# pi     pi+    pi0    pi- '
      do i = 1, ntmax
        do j = 1, 3
          jw = j + 1
          w(jw)= float(numesst(i,j))*norm
        end do
        w(1) = w(2)+w(3)+w(4)
        write(isum,202)float(i)*dt,w(1),w(2),w(3),w(4)
      end do

      write(isum,*)'time   tot# rho     rho+    rho0    rho- '
      do i = 1, ntmax
        do j = 5, 7
          jw  = j - 3
          w(jw)= float(numesst(i,j))*norm
        end do
        w(1) = w(2)+w(3)+w(4)
        write(isum,202)float(i)*dt,w(1),w(2),w(3),w(4)
      end do

      write(isum,*)'time   tot# eta    tot # sigma  '
      do i = 1, ntmax
        w(1) = float(numesst(i,4))*norm
        w(2) = float(numesst(i,8))*norm
        write(isum,203)float(i)*dt,w(1),w(2)
      end do

      write(isum,*) 'density in central cell: '
      do i = 1, ntmax
        write(isum,*)float(i)*dt,cenden(i,1),cenden(i,2), ncdelta(i)
      end do

      write(isum,*)
      write(isum,*)
      write(isum,*)
      write(isum,*)'**********************************************'
      write(isum,*)'********** C O L L I S I O N S ***************'
      write(isum,*)'**********************************************'

      write(isum,*)'************baryons***************************'
      write(isum,*)'                                              '
      write(isum,*)'****************************'
      write(isum,*)'baryon collision numbers '
      write(isum,*)'****************************'

      do i = 1, lcodim
        lcoll(i) = 0
      end do
      do j = 1, nmaxdim
        do i = 1, lcodim
          lcoll(i) = lcollst(j,i) + lcoll(i)
        end do
      end do


      nin = 2
      nout = 1
      do i = 1, lcodim
        if(i.le. (nres+1)) then
          write(isum,*) lcoll(i) ,'   ',nameres(1), ' + ',nameres(i),
     &              ' --  '  ,nameres(1), ' + ',nameres(i),'el'
        else if(i.gt.(nres+1).and. i.le.(2*nres+1)) then
        write(isum,*)lcoll(i) ,'    ',nameres(1), ' + ',nameres(i-nres),
     &                   ' --  ',nameres(1),' + ',nameres(1)
        else if(i.gt.(2*nres+1).and. i.le.(3*nres+1)) then
        write(isum,*)lcoll(i),'     ',nameres(1),' + ',nameres(1),' -- '
     &              ,nameres(1),' + ',nameres(i-2*nres)
        else if(i.gt.(3*nres+1).and.i.le.(nres**2+3*nres+1)) then
          if(nout.gt.2.and. mod(nout,nres).eq.1) then
            nout = 1
            nin  = nin + 1
          end if
          nout = nout + 1
       write(isum,*)lcoll(i),'  ',nameres(1),' + ',nameres(nin),' --  ',
     &          nameres(1),' + ' ,nameres(nout)
        else if(i.eq.(nres**2+3*nres+2))then
          write(isum,*)lcoll(i),'  ',nameres(1),' + ',nameres(1),' -- ',
     &                          nameres(2),' + ',nameres(2)

        else if(i.eq.(nres**2+3*nres+3))then
          write(isum,*)lcoll(i),'  ',nameres(2),' + ',nameres(2),' -- ',
     &                          nameres(1),' + ',nameres(1)
        else if(i.eq.(nres**2+3*nres+4))then
                 write(isum,*)lcoll(i),'    s-state pion prod.'
        end if
      end do

      write(*,*)'1'
      do i = 1, nres
        do ii = 1,6
          cres(i,ii) = 0
        end do
      end do
      do j = 1, nmaxdim
        do i = 1, nres
          do ii = 1,6
          cres(i,ii) = cresst(j,i,ii) + cres(i,ii)
          end do
        end do
      end do

      do i = 1, nres
        lmesc(i) = 0
      end do
      do j = 1, nmaxdim
        do i = 1, nres
          lmesc(i) = lmescst(j,i) + lmesc(i)
        end do
      end do

      write(isum,*)' **********************************'
      write(isum,*)' *******      decays      *********'
      write(isum,*)
      write(isum,*)' res - mes + something '
      do i = 1, nres
        write(isum,*) nameres(i+1), lmesc(i)
        write(isum,*) (cres(i,j),j=1,6)
      end do

      return
      end


