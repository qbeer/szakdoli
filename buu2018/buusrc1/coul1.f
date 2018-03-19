**** this fiel contains the routines needed for the calculation of ****
**** the electromagnetic potential ************************************
***  pois, fielca, rprimc, rhside, sort, iniit, bounds, cdens ********
***  the link of the subroutines to the BUU code is made via the ******
***  routine emfoca, which is stored in a separate file (emfoca.f)*****
************************************************************************
************************************************************************
      subroutine initcoul

      implicit none

      include"coucom"
      include"common"
      include"cominput"
      integer it, iix, iiy, iiz
      integer ncount, nloop, m, l, k

      ncount = 0
      do it = 1, maxpar
        if(id(1,it) .ne. 0) then
        iix = nint(r(1,it)/dgrid)
        iiy = nint(r(2,it)/dgrid)
        iiz = nint(r(3,it)/dgrid)

        if(abs(iix).gt.cmaxx-4 .or. abs(iiy).gt.cmaxx-4 .or.
     +     abs(iiz).gt.cmaxz-4) then
            ncount = ncount +1
        end if
        end if
      end do

      if(ncount .ge. 1) then
        write(*,*)'in initcoul ungeschaltet !!!!!!!!!!!'
        dgrid = 2.0*dgrid
        call cdenin
      end if

      return

      entry initcoulen

      dgrid = dgrids
      call cdenin

      return
c
      entry   initemrpot
c-------------------------  for reproducing consequent runs (isub)
c-------------------------  the elm-force should have the same value
c-------------------------  at beginning(=0)
          do nloop = 0, ncont
                 ncount = 1
          do m = -cmaxz,cmaxz
          do l = -cmaxx,cmaxx
          do k = -cmaxx,cmaxx
                 emrpot(nloop,k,l,m)  = 0.0
                 rprim (0,ncount) =  .0
                 ncount = ncount + 1
          end do
          end do
          end do
          end do
      return
      end

*****************************************************************

      subroutine pois(ncont)
      implicit none

      include"coucom"
*      integer kk, nrr
      integer  nit
*      parameter(kk   = 10)
*      parameter(nrr=2**(kk-1))
      parameter(nit= 200)

      integer i, k, l, m, ncount, nloop, nhelp, s
*      integer twopr

      real*8 w(inges), g(inges),  poth(0:inges)
      real*8 error, conpar, er, er1
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,tend,a(1:inges)
      real*8 c(1:inges), b(1:inges), rnenn
*      real*8 alpha, beta
*      real*8 sk(1:nrr,1:kk), rconpa(1:nrr), disc
*      real*8 alph(1:kk), ab, bet(1:kk)
      integer  nval, ncont, nhil, nhil2
      integer  kinput, j, nr, k1, ntest
      real*8 abw, zwi, rav, test
      integer tl, tk, tm


      real*8 t0, tin, tout
      t0 = 0.0
c      tin = secnds(t0)
      er1 = 0.0

*      write(*,*)'fester konvergenzpar. 2.1'
      conpar = 2.1


***  calculate potential values at boundaries***********************
            call bounds(ncont)

         pot(0) = 0.0

         ncount = 0
          do m = -cmaxz,cmaxz
            do l =-cmaxx,cmaxx
              do k =-cmaxx,cmaxx
                ncount = ncount +1
                i = ncount
                iz(i)  = m + cmaxz + 1
                iy(i)  = l + cmaxx + 1
                ix(i)  = k + cmaxx + 1
                pot(i) = 0.0
              end do
            end do
           end do

********************************************************************
**** calculation of convergence parameters *************************

*          if(ntest.eq.1) then
*
*            if(kinput .gt. kk-1) then
*             write(*,*)'problems with calculating the conv. parm.'
*             stop
*            end if
*            k1 = kinput+1
*            nr = 2**kinput
***   calculate alphas, betas
*              alph(1) = alpha
*              bet(1)  = beta
*            do j = 1,kinput
*              alph(j+1) = sqrt(alph(j)*bet(j))
*              bet(j+1)  = 0.5d0 * (alph(j)+bet(j))
*            end do

***   calculate the s's
*                sk(1,1) = sqrt(alph(k1)*bet(k1))
*              do j = 1,kinput
*                ab = alph(k1-j)*bet(k1-j)
*                twopr = 2**(j-1)
*                do m = 1, twopr
*                  disc        = dsqrt(sk(m,j)**2-ab)
*                  sk(2*m,j+1)  = sk(m,j) + disc
*                  sk(2*m-1,j+1)=ab/sk(2*m,j+1)
*                end do
*              end do
*              write(*,*)'break 2'
***   determine the r's
*              do k = 1, nr
*                rconpa(k) = sk(k,k1)
*                write(*,*)'conpar ',k,rconpa(k)
*              end do
*            end if


        call rprimc(ncont)

          nhil = inges/2+indx**2/4
          nhil2= (2*cmaxz+1) * (2*cmaxx+1)
c         write(*,*) ' pois2 ', inges,indx,ix(nhil),iy(nhil),
c    1    iz(nhil), pot(nhil), rprim(0,nhil2), rpri(nhil2)
**x*
********* loop over 0th dim or over 4-dim ***************************
***
         do nloop = 0,ncont
           if(loopte(nloop).eq.1) then
           error = 1000.0
           nhelp = 0

***  solve iteration equation (preparation)

        poth(0)= 0.0

*** store the correct potential and density in iteration variables ***
        call iniit(nloop)

        write(*,*)'eps(nloop) = ',eps(nloop)


        do while (error .gt. eps(nloop))
          nhelp = nhelp + 1
c          write(*,*)'iterationsschritt = ',nhelp, error
          if(nhelp .eq.nit) goto 10000


***   set value for convergence parameter
*        if(ntest.eq.1) then
*            ncount = mod(ntest,nrr)
*            conpar = rconpa(ncount)
*        end if

***   save potential for error-vector
          do k = 1, inges
            poth(k) = pot(k)
          end do

***   iteration-step starts
          do i = 1, 3

***     calculate right hand side
         call rhside(i,conpar)


***   solve one dim. equation

         if(i.eq.1) nval = indx
         if(i.eq.2) nval = indx
         if(i.eq.3) nval = indz

           do s = 1, inges

             if(mod(s,nval).eq.0) then
               c(s) = 0.0
             else
               c(s) = -1.0
             end if
             if(mod(s,nval).eq.1) then
               a(s) = 0.0
             else
               a(s) = -1.0
             end if
             b(s) = 2.0 + conpar
           end do
           w(1) = c(1) / b(1)
           g(1) = rhs(1)/b(1)
           do s = 2, inges
             rnenn = b(s) - a(s)*w(s-1)
             w(s) = c(s)/rnenn
             g(s)  = (rhs(s) - a(s)*g(s-1))/rnenn
           end do
***   solution of 1 dim. equation
         pot(inges) = g(inges)
         do s = inges-1, 1, -1
           pot(s) = g(s) - w(s)*pot(s+1)
         end do
c          write(*,*)'pois sort elott ',i

***   sort the vectors for next iteration-step
          call sort(i)

        end do


***   build error-vector  (maximum norm)
        er1 = 0.0
         do i = 1,inges
           er       = dabs(pot(i) - poth(i))
           error    = max(er1,er)
           er1      = error
         end do

       end do
*******   end of iteration *************************************


***   store the final potentials in the 'gloal vaector'
10000   continue
             rav = 0.0
             ncount = 1
          do m = -cmaxz,cmaxz
          do l = -cmaxx,cmaxx
          do k = -cmaxx,cmaxx
                 emrpot(nloop,k,l,m)  = sngl(pot(ncount))
                 rav = rav + emrpot(nloop,k,l,m)
                 ncount = ncount + 1
          end do
          end do
          end do
          rav = rav/ncount
          write(*,*)'average potential ',rav


***  test test
          abw = 0.0
          test= 0.0
           do m = -cmaxz,cmaxz
           do l = -cmaxx,cmaxx
           do k = -cmaxx,cmaxx
             t1 = -emrpot(nloop,k-1,l,m)
             t2 = 2.0 * emrpot(nloop,k,l,m)
             t3 = -emrpot(nloop,k+1,l,m)
             t4 = -emrpot(nloop,k,l-1,m)
             t5 = 2.0*emrpot(nloop,k,l,m)
             t6 = -emrpot(nloop,k,l+1,m)
             t7 = -emrpot(nloop,k,l,m-1)
             t8 = 2.0*emrpot(nloop,k,l,m)
             t9 = -emrpot(nloop,k,l,m+1)
             tend = t1+t2+t3+t4+t5+t6+t7+t8+t9
             zwi = abs(tend-12.56637061*chrho(nloop,k,l,m)/dgrid)
             if(zwi.gt.abw) then
               tl = l
               tm = m
               tk = k
               test =emrpot(nloop,tk,tl,tm)
             end if
             abw = max(abw,zwi)

           end do
           end do
           end do

           write(*,*)'max abw. ',abw, test
           write(*,*)tl, tk, tm


           do m = -cmaxzb,cmaxzb
           do l = -cmaxxb,cmaxxb
           do k = -cmaxxb,cmaxxb
             emrpot(nloop,k,l,m)=elmcon*emrpot(nloop,k,l,m)
           end do
           end do
           end do
           else if(loopte(nloop).eq.0) then
             write(*,*)'no current in the ',nloop,'-th component'
           else
             write(*,*)'problems with loopte '
           end if

        end do
****** end of dimension loop************************************


c         call fielca(ncont)


c         tout = secnds(t0)
c        tin = tout - tin
         write(*,*)' leaving pois = '
         return

      end


*******************************************************************
*******************************************************************

      subroutine rprimc(ncont)
***********************************************************************
***     this subroutine generates the vector rpri which is used as
***     an input vector for the iteration process.
***     the sorting is done in the way that is first needed by the
***     routine sort (x tridiagonal).
**********************************************************************

      implicit none
      include"coucom"

      integer i, k, l, m, ncount, ncont

      real*8 phix, phiy, phiz


      do i = 0, ncont
          ncount = 1
           do m = -cmaxz,cmaxz
            do l = -cmaxx,cmaxx
              do k = -cmaxx,cmaxx

                phix = 0.0
                phiy = 0.0
                phiz = 0.0

                if(m.eq.-cmaxz) phiz = dble(emrpot(i,k,l,-cmaxzb))
                if(m.eq. cmaxz) phiz = dble(emrpot(i,k,l, cmaxzb))
                if(l.eq.-cmaxx) phiy = dble(emrpot(i,k,-cmaxxb,m))
                if(l.eq. cmaxx) phiy = dble(emrpot(i,k, cmaxxb,m))
                if(k.eq.-cmaxx) phix = dble(emrpot(i,-cmaxxb,l,m))
                if(k.eq. cmaxx) phix = dble(emrpot(i, cmaxxb,l,m))

                rprim(i,ncount) = 12.56637061*dble(chrho(i,k,l,m)/
     +                            dgrid)+phix+phiy+phiz
                ncount = ncount + 1
              end do
            end do
          end do

        end do

        return
        end


************************************************************************
***********************************************************************


        subroutine iniit(nloop)
************************************************************************
***     this routine writes the potendials and rhs. to the fields
***     used by the iteration algorithm.
************************************************************************

        implicit none
        include"coucom"
        integer nloop, ncount, k, l, m
********************
**********************
        ncount = 1
          do m = -cmaxz,cmaxz
            do l = -cmaxx,cmaxx
              do k = -cmaxx,cmaxx

                ix(ncount)   = k + cmaxx + 1
                iy(ncount)   = l + cmaxx + 1
                iz(ncount)   = m + cmaxz + 1
                pot(ncount)  = dble(emrpot(nloop,k,l,m)/elmcon)
                rpri(ncount) = rprim(nloop,ncount)
                ncount       = ncount + 1
              end do
            end do
          end do

          return
          end

************************************************************************
************************************************************************

      subroutine rhside(nread,conpa)
      implicit none
      include"coucom"
      integer nread, id1,id2, id3
      integer ipotd
      integer ipot1y, ipot2y
      integer ipot1x, ipot2x
      integer ipot1z, ipot2z
      integer i
      real*8 conpa


********************************************************
***** I M P O R T A N T********************************
*******************************************************
         pot(0)  = 0.0


        if(nread .eq. 1) then

***    do derivative : x  tridiagonal

          do i = 1,inges

***    diagonal part
            id1 = ix(i)
            id2 = iy(i) - 1
            id3 = iz(i) - 1


            if (id1 .le. 0)  id1 = 0
            if (id2 .le. 0)  id2 = 0
            if (id3 .le. 0)  id3 = 0
            if (id1 .gt. indx)  id1 = 0
            if (id2 .gt. indx)  id2 = 0
            if (id3 .gt. indz)  id3 = 0

            ipotd = id1 + id2*indx + id3*insqxx


***    z-derivative
            ipot1z = 0
            if((iz(i)-1).gt.0) ipot1z=id1 + id2*indx + (id3-1)*insqxx

            ipot2z = 0
            if((iz(i)+1).le.indz) ipot2z=id1 + id2*indx +(id3+1)*insqxx

***    y-derivative
            ipot1y = 0
            if((iy(i)-1).gt.0) ipot1y=id1 + (id2-1)*indx + id3*insqxx

            ipot2y = 0
            if((iy(i)+1).le.indx) ipot2y=id1 + (id2+1)*indx + id3*insqxx

****    define rhs of iteration equn.
            rhs(i) = (conpa - 4.0)*pot(ipotd) + pot(ipot1z) +
     +               pot(ipot2z) + pot(ipot1y) + pot(ipot2y) +
     +               rpri(ipotd)

          end do

        else if(nread .eq. 2) then
***    do derivative : y tridiagonal

          do i = 1,inges
***    diagonal part
            id1 = ix(i) - 1
            id2 = iy(i)
            id3 = iz(i) - 1

            if (id1 .le. 0)  id1 = 0
            if (id2 .le. 0)  id2 = 0
            if (id3 .le. 0)  id3 = 0
            if (id1 .gt. indx)  id1 = 0
            if (id2 .gt. indx)  id2 = 0
            if (id3 .gt. indz)  id3 = 0


            ipotd = id2 + id3*indx + id1*insqxz

***    z-derivative
            ipot1z = 0
            if((iz(i)-1).gt.0) ipot1z=id2 + (id3-1)*indx + id1*insqxz

            ipot2z = 0
            if((iz(i)+1).le.indz) ipot2z=id2 + (id3+1)*indx + id1*insqxz

***    x-derivative
            ipot1x = 0
            if((ix(i)-1).gt.0) ipot1x=id2 + id3*indx + (id1-1)*insqxz

            ipot2x = 0
            if((ix(i)+1).le.indx) ipot2x=id2 + id3*indx + (id1+1)*insqxz

***    define rhs of iteration equn.
            rhs(i) = (conpa -4.0)*pot(ipotd) + pot(ipot1z) +
     +               pot(ipot2z) + pot(ipot1x) + pot(ipot2x) +
     +               rpri(ipotd)


          end do
        else if(nread .eq. 3) then
***    do derivative : z tridiagonal

          do i = 1,inges
***    diagonal part
            id1 = ix(i) - 1
            id2 = iy(i) - 1
            id3 = iz(i)

            if (id1 .le. 0)  id1 = 0
            if (id2 .le. 0)  id2 = 0
            if (id3 .le. 0)  id3 = 0
            if (id1 .gt. indx)  id1 = 0
            if (id2 .gt. indx)  id2 = 0
            if (id3 .gt. indz)  id3 = 0

            ipotd = id3 + id1*indz + id2*insqxz
*            write(3,*)'ipotd =',ipotd
***    y-derivative
            ipot1y = 0
            if((iy(i)-1).gt.0) ipot1y=id3 + id1*indz + (id2-1)*insqxz
*            write(3,*)'ipot1y =',ipot1y
            ipot2y = 0
            if((iy(i)+1).le.indx) ipot2y=id3 + id1*indz + (id2+1)*insqxz
*            write(3,*)'ipot2y =',ipot2y
***    x-derivative
            ipot1x = 0
            if((ix(i)-1).gt.0) ipot1x=id3 + (id1-1)*indz + id2*insqxz
*            write(3,*)'ipot1x =',ipot1x
            ipot2x = 0
            if((ix(i)+1).le.indx) ipot2x=id3 + (id1+1)*indz + id2*insqxz
*            write(3,*)'ipot2x =',ipot2x
***    define rhs of iteration equn.
            rhs(i) = (conpa - 4.0)*pot(ipotd) + pot(ipot1y) +
     +               pot(ipot2y) + pot(ipot1x) + pot(ipot2x) +
     +               rpri(ipotd)


          end do
        end if

        return
        end

************************************************************************
************************************************************************

      subroutine sort(nread)

      implicit none
      include"coucom"

      integer nread, i, k, ix1(1:inges), iy1(1:inges), iz1(1:inges)
      real*8 pot1(1:inges), rpri1(1:inges)

        if(nread .eq. 1) then
*     sort from x tridiagonal to y tridiagonal
          do i = 1,inges
            k = iy(i) + (iz(i)-1)*indx + (ix(i)-1)*insqxz
            pot1(k) = pot(i)
            ix1(k)  = ix(i)
            iy1(k)  = iy(i)
            iz1(k)  = iz(i)
            rpri1(k)= rpri(i)
          end do

        else if(nread .eq. 2) then
*     sort from y tridiagonal to z tridiagonal
          do i = 1,inges
            k = iz(i) + (ix(i)-1)*indz + (iy(i)-1)*insqxz
            pot1(k)= pot(i)
            ix1(k) = ix(i)
            iy1(k) = iy(i)
            iz1(k) = iz(i)
            rpri1(k)= rpri(i)
          end do

        else if(nread .eq. 3) then
*     sort from z tridiagonal to x tridiagonal
          do i = 1,inges
            k = ix(i) + (iy(i)-1)*indx + (iz(i)-1)*indx**2
            pot1(k)= pot(i)
            ix1(k) = ix(i)
            iy1(k) = iy(i)
            iz1(k) = iz(i)
            rpri1(k)= rpri(i)
          end do
        end if

       do i = 1,inges
         pot(i) = pot1(i)
         ix(i)  = ix1(i)
         iy(i)  = iy1(i)
         iz(i)  = iz1(i)
         rpri(i)= rpri1(i)

        end do

      return
      end

************************************************************************
************************************************************************



      subroutine bounds(ncont)
*****
**     berechung der bounds fuer die poisson-gleichung
**     es wird angenommen, dass die dichte- bzw. die stromverteilung auf
**     einem dreidimensionalen gitter gegeben ist.

      implicit none

      include"coucom"


      real*8 rcell
**
**    Felder fuer die verschiedenen momente (hat)

      real*8 mhat0(0:3), mhat1(0:3,1:3), mhat2(0:3,1:3,1:3)
      real*8 mhat3(0:3,1:3,1:3,1:3)
      real*8 dcms(0:3,1:3)
      real*8 m0(0:3), m1(0:3,1:3), m2(0:3,1:3,1:3)
      real*8 m3(0:3,1:3,1:3,1:3)
      real*8 xpri(1:3), rbig, rmin
      real*8 rabs, term1, term2, term3, term4, rmax

      integer i, j, k, l, m, n, o, p, ncont

**  felder auf null setzen

      do l = 0, ncont
        loopte(l) = 0
        mhat0(l) = 0.0
        m0(l)    = 0.0
        do i = 1, 3
          mhat1(l,i) = 0.0
          m1(l,i)    = 0.0
          dcms(l,i)   = 0.0
          do j = 1, 3
            mhat2(l,i,j) = 0.0
            m2(l,i,j)    = 0.0
            do k = 1, 3
              mhat3(l,i,j,k) = 0.0
              m3(l,i,j,k)    = 0.0
            end do
           end do
        end do
      end do

**  calculate the mhats
      do l = 0, ncont
        do i = -cmaxzb,cmaxzb
          do j = -cmaxxb, cmaxxb
            do k = -cmaxxb, cmaxxb

**    calculate xprime
              xpri(1) = float(k)*dgrid
              xpri(2) = float(j)*dgrid
              xpri(3) = float(i)*dgrid

***   evaluate mhats explicitely

              rcell = chrho(l,k,j,i)
              if(loopte(l).eq.0 .and. rcell.ne.0.0) loopte(l) = 1

              mhat0(l) = mhat0(l) + rcell

          do m = 1, 3
             mhat1(l,m)    = mhat1(l,m)    + rcell*xpri(m)
          do n = 1, 3
            mhat2(l,m,n)  = mhat2(l,m,n)  + rcell*xpri(m)*xpri(n)
            do o = 1, 3
            mhat3(l,m,n,o)= mhat3(l,m,n,o)+ rcell*xpri(m)*xpri(n)
     +  *xpri(o)
              end do
          end do
          end do
        end do
        end do
        end do
       end do
**** end of calculation of the mhats


******** set potential to zero and calculate pots at bounds **********

        do i = 0, ncont
          rbig = 0.0
          rmin = 1.0e+08
          rmax = 0.0

          do m = -cmaxzb, cmaxzb
            do l = -cmaxxb, cmaxxb
              do k = -cmaxxb, cmaxxb



               if( m.eq.-cmaxzb .or. m.eq.cmaxzb .or.
     +             l.eq.-cmaxxb .or. l.eq.cmaxxb .or.
     +             k.eq.-cmaxxb .or. k.eq.cmaxxb)  then


                  xpri(1) = float(k)*dgrid
                  xpri(2) = float(l)*dgrid
                  xpri(3) = float(m)*dgrid
********


**          absolute value of r
                  rabs = sqrt(xpri(1)**2 + xpri(2)**2 + xpri(3)**2)
                   if(rabs.lt.0.0) then
                     write(*,*) 'mistake while calc. of rabs'
                     stop
                   end if

**           calculate pots. at bounds.

                   term1 = mhat0(i)/rabs

                   term2 = 0.0
                   do n = 1,3
                     term2 = term2 + mhat1(i,n)*xpri(n)
                   enddo
                   term2 = term2/(rabs**3)

                   term3 = 0.0
                    do n = 1,3
                      do o = 1,3
                        if (n.ne.o) then
                          term3 = term3 + mhat2(i,n,o)*
     +                    (3.0*xpri(n)*xpri(o))
                        endif
                      end do
                    end do

                    term3 = term3/2.0/(rabs**5)

                    term4 = 0.0
                    do n = 1,3
                      do o = 1,3
                        do p = 1,3
                        if (n.ne.o) then
                          term4 = term4 + mhat3(i,n,o,p)*
     +                    5.0*xpri(n)*xpri(o)*xpri(p)
                        else
                          term4 = term4 + mhat3(i,n,o,p)*
     +                    (5.0*xpri(n)*xpri(o)*xpri(p)-
     +                     3.0*(rabs**2)*xpri(p))
                        endif
                        end do
                      end do
                    end do
                    term4 = term4/2.0/(rabs**7)
                    rmax = max(rmax,term4)
****************************************************************


                emrpot(i,k,l,m) = term1+term2+term3+term4
*                emrpot(i,k,l,m) = emrpot(i,k,l,m)*4.0*3.1415926
                rbig = max(rbig, abs(emrpot(i,k,l,m)))
                rmin = min(rmin, abs(emrpot(i,k,l,m)))
               end if

              end do
            end do
          end do
        eps(i) = 0.01*rmin
        end do

c 1000   continue
        return

        end


******************************************************************
***********************************************************************

************************************************************************
*
      subroutine cdens(nescc, nescpic,nt)
*
*       purpose:     calculation of electrical current j(mu) from
*                    the spatial distribution of the testparticles.
*                    two different smearing algorithms are implemneted:
*                      i) idec = 0 gaussian smearing (k. niita)
*                     ii) idec = 1 cub. spline (lattanzio)
*
*       variables (all input, all integer)
*         minnum   -  first testparticle treated in one run for density
*         maxnum   -  last testparticle treated in one run for density
*         num      -  number of testparticles per nucleon
*         nescc    -  number of escaped particles      (integer,output)
*
************************************************************************
      implicit none

      include"common"
      include"coucom"
      include"cominput"

      common /nthwhw/  nthw
      integer nthw

      integer ipc, iqc, idec
      real*8 hd, deltav, ccut

      parameter(idec= 0)  ! control for smearing
      parameter(ipc = 2)  !parameter for inner mesh
      parameter(iqc = 2)  !parameter for outer mesh
      parameter(hd  = 2.0)
      parameter(ccut= 5.0)
      real*8 ccm((2*ipc+1)**3, (2*iqc+1)**3)
      save ccm
      integer ixc, iyc, izc, inn
      integer nescc, nescpic
      integer kx, ky, kz, icc, ibc
      integer jx, jy, jz, i,nt
      real*8    zwi(0:3), psq, rdc, xi, yi, zi, sek
      real*8    xj, yj, zj, rsqr, quo
      real*8    wei, test
      real*8    x, y, z
      deltav = dgrid**3
      write(*,*)'in cdens'


*      set density and currents to zero
      do inn = 0, ncont
        do izc = -cmaxzb,cmaxzb
          do iyc = -cmaxxb,cmaxxb
            do ixc = -cmaxxb,cmaxxb
              chrho(inn,ixc,iyc,izc) = 0.0
            end do
          end do
        end do
      end do

*
*      counter for escaped particles (coulomb) to zero
      nescc   = 0
      nescpic = 0
      write(*,*)'vor baryon dens'

      do i = 1, maxpar
        if(id(1,i) .ne. 0) then
        ixc = nint( r(1,i)/dgrid )
        iyc = nint( r(2,i)/dgrid )
        izc = nint( r(3,i)/dgrid )

***     actual calculation of densities

        if(abs(ixc).gt.cmaxxb .or. abs(iyc).gt.cmaxxb .or.
     +     abs(izc).gt.cmaxzb) then
           nescc = nescc + 1
        else
*
           kx=nint(float(2*ipc+1)/dgrid*(r(1,i)-float(ixc)*dgrid))
           if(abs(kx) .eq. ipc+1) kx = kx/abs(kx) * ipc
           ky=nint(float(2*ipc+1)/dgrid*(r(2,i)-float(iyc)*dgrid))
           if(abs(ky) .eq. ipc+1) ky = ky/abs(ky) * ipc
           kz=nint(float(2*ipc+1)/dgrid*(r(3,i)-float(izc)*dgrid))
           if(abs(kz) .eq. ipc+1) kz = kz/abs(kz) * ipc


***          number in inner mesh
           icc=1+(kz+ipc)+(ky+ipc)*(2*ipc+1)+(kx+ipc)*(2*ipc+1)**2

***          counter in outer mesh
           ibc = 0

c--------------------------------------
c     if(nthw .lt. 3  .and. id(2,i) .ne. 0)
c    1    write(*,7777) i, nescc, ixc, iyc, izc, id(2,i)
c7777     format(' part in cdens ',2i7, 3i4, i2)
c--------------------------------------
***          outer mesh loop
           do jx=ixc-iqc,ixc+iqc
             do jy=iyc-iqc,iyc+iqc
               do jz=izc-iqc,izc+iqc
               ibc =ibc + 1

                if(ccm(icc,ibc).gt.0.0.and.
     +            abs(jx).le.cmaxxb.and.abs(jy).le.cmaxxb.and.
     +            abs(jz).le.cmaxzb) then
*

                    psq = p(1,i)**2+p(2,i)**2+p(3,i)**2
                    zwi(0) = sqrt((e(i)+upot(i))**2+psq)
                    zwi(1) = p(1,i)
                    zwi(2) = p(2,i)
                    zwi(3) = p(3,i)

                    do inn = 0, ncont

                     chrho(inn,jx,jy,jz)=chrho(inn,jx,jy,jz)+
     +               ccm(icc,ibc)*float(id(2,i))*zwi(inn)/zwi(0)
                    end do

                end if
*
               end do
             end do
           end do
*

      end if
      end if
***      end of loop over particles in one ensemble (baryons)
      end do

      write(*,*)'vor meson loop', chrho(0, 5, 5, 0)
***      start of meson loop
      do i = 1, maxppar
*
c        if(ipi(1,i) .ne. 1)    goto 401
        if(ipi(1,i) .eq. 0)    goto 401
        if(ipi(2,i) .eq. 0)    goto 401
*
        ixc = nint( rpi(1,i)/dgrid )
        iyc = nint( rpi(2,i)/dgrid )
        izc = nint( rpi(3,i)/dgrid )
*
        if(abs(ixc).gt.cmaxxb.or.abs(iyc).gt.cmaxxb.or.
     +                           abs(izc).gt.cmaxzb) then
        nescpic= nescpic + 1
        else
*
           kx=nint(float(2*ipc+1)/dgrid*(rpi(1,i)-float(ixc)*dgrid))
           if(abs(kx) .eq. ipc+1) kx = kx/abs(kx) * ipc
           ky=nint(float(2*ipc+1)/dgrid*(rpi(2,i)-float(iyc)*dgrid))
           if(abs(ky) .eq. ipc+1) ky = ky/abs(ky) * ipc
           kz=nint(float(2*ipc+1)/dgrid*(rpi(3,i)-float(izc)*dgrid))
           if(abs(kz) .eq. ipc+1) kz = kz/abs(kz) * ipc

***          number in inner mesh
           icc=1+(kz+ipc)+(ky+ipc)*(2*ipc+1)+(kx+ipc)*(2*ipc+1)**2

***          counter for the outer mesh
           ibc = 0
            do  jx=ixc-iqc,ixc+iqc
              do  jy=iyc-iqc,iyc+iqc
                do  jz=izc-iqc,izc+iqc
                 ibc =ibc + 1
                  if(ccm(icc,ibc).gt.0.0.and.
     +               abs(jx).le.cmaxxb.and.abs(jy).le.cmaxxb.and.
     +               abs(jz).le.cmaxzb) then


                    zwi(1) = ppi(1,i)
                    zwi(2) = ppi(2,i)
                    zwi(3) = ppi(3,i)
                 zwi(0) = sqrt((epi(i)+mpot(i))**2+zwi(1)**2+
     &                          zwi(2)**2+zwi(3)**2)

                    do inn = 0, ncont
                     chrho(inn,jx,jy,jz)=chrho(inn,jx,jy,jz)+
     +               ccm(icc,ibc)*float(ipi(2,i))*zwi(inn)/zwi(0)
                    end do
                  end if
               end do
             end do
           end do

*
      end if
 401   continue


       end do

***   total charge and center of charge in order to determine the     *
*     C-force acting on a prticle which has left the grid             *

      chatot = 0.0
      do i = 1,3
        rchar(i) = 0.0
      end do


      do izc = -cmaxzb,cmaxzb
        z = float(izc)*dgrid
        do iyc = -cmaxxb,cmaxxb
          y= float(iyc)*dgrid
          do ixc = -cmaxxb,cmaxxb
            x = float(ixc)*dgrid

            rchar(1) = rchar(1) + x*chrho(0,ixc,iyc,izc)
            rchar(2) = rchar(2) + y*chrho(0,ixc,iyc,izc)
            rchar(3) = rchar(3) + z*chrho(0,ixc,iyc,izc)
            chatot   = chatot   +   chrho(0,ixc,iyc,izc)

          end do
        end do
      end do

      if(chatot.lt. 1e-3) then
        write(*,*)'error in cdens total charge = ', chatot
        stop
      end if

      do i = 1, 3
        rchar(i) = rchar(i)/chatot
      end do

      write(*,*)'test test test cdens :' , chrho(0,5,5,0)
      write(*,*)'coul1 total charge = ',chatot, iseed
      write(*,*)'center of charge: ',rchar(1), rchar(2), rchar(3)




*

*-----------------------------------------------------------------------
*
      return
*      end
*
************************************************************************
*                                                                      *
      entry cdenin
*                                                                      *
************************************************************************
*

*
      write(*,*)'incdenin'
      write(*,*)
      write(*,*)'ipc =',ipc
      write(*,*)'iqc =',iqc
      write(*,*)'hd =',hd
      rdc = dgrid/float(ipc*2+1)

      do ixc=-ipc,ipc
        do iyc=-ipc,ipc
          do izc=-ipc,ipc
          icc=1+(izc+ipc)+(iyc+ipc)*(2*ipc+1)+(ixc+ipc)*(2*ipc+1)**2
          xi=float(ixc)*rdc
          yi=float(iyc)*rdc
          zi=float(izc)*rdc
          ibc = 0
          sek = 0.0

          do jx=-iqc,iqc
            do jy=-iqc,iqc
              do jz=-iqc,iqc
               ibc=ibc+1
               xj=float(jx)*dgrid
               yj=float(jy)*dgrid
               zj=float(jz)*dgrid
               rsqr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
               quo = rsqr/(hd**2)
               if(quo.gt.ccut) then
                ccm(icc,ibc)=0.0
               else
                ccm(icc,ibc)=wei(quo,idec)
               end if

               sek=sek+ccm(icc,ibc)
              end do
            end do
          end do
          ibc = 0
         do jx=-iqc,iqc
           do jy=-iqc,iqc
             do jz=-iqc,iqc
             ibc = ibc + 1
             ccm(icc,ibc)=ccm(icc,ibc)/sek/float(num)
             end do
           end do
         end do
          test = 0.0
         ibc = 0
         do jx=-iqc,iqc
           do jy=-iqc,iqc
             do jz=-iqc,iqc
             ibc = ibc + 1
              test = test +ccm(icc,ibc)
             end do
           end do
         end do
         test = test*float(num)
         if(abs(test-1.0).gt.0.1)then
           write(*,*)'normalization error'
          write(*,*)'test = ',test
         end if
          end do
        end do
      end do
      ibc = 0
          do jx=-iqc,iqc
            do jy=-iqc,iqc
              do jz=-iqc,iqc
                ibc = ibc + 1
               end do
               end do
               end do

*
      return
      end

      real*8 function wei(x,idec)

      implicit none

      integer idec
      real*8 x

      wei = 0.0
      if (idec.eq.0) then
        wei = exp(-x/2.0)
      else if(idec.eq.1) then
        x = sqrt(x)
        if(abs(x).le.1.0) then
          wei = 1.0 - 1.5*(x**2) + 0.75*(x**3)
        else if(abs(x).le.2.0 .and. abs(x).gt.1.0)then
          wei = 0.25*(2.0-x)**3
        else
          wei = 0.0
        end if
        x = x**2
      end if

      return
      end


************************************************************************
***********************************************************************
