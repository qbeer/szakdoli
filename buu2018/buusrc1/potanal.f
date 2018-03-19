c**********************************
      real*8 function  potanal(rho0, rho,plrf,id1,rhap)
c       for the threefluid model, uncomment the line in the if
c       loop ithree.eq.1
      implicit none
      include "cominput"
      include "resdata"

      real*8    nuclpot, plrf(1:3), rhap(1:3), strangeLpot
      integer  id1, idres,ii
      real*8    rho0, rho, p
      real*8    out
      logical deltares, nuclres, antipr,strange

      potanal = 0.0
      if(massta.le.1 .and. masspr.le.1) return
      out = 0.0

      if(id1.eq.-1) return
c      goto 10
      nuclres  = .false.
      deltares = .false.
      antipr   = .false.
      strange  = .false.
      idres = id1 - 1

      if(idres.eq.0) then
        nuclres = .true.
      else if(idres.eq.-2) then
        antipr = .true.        ! antiproton
      else if(idres.gt.nres .and. idres.le.nres+2) then
        strange = .true.        ! Lambda or Sigma !!!!
      else
        if(resprop2(idres,2).eq.1) then
          nuclres = .true.
        else if(resprop2(idres,2).eq.3) then
          deltares = .true.
        end if
      end if

********************* test **************
      p = sqrt(plrf(1)**2 + plrf(2)**2 + plrf(3)**2)
      if(p.gt.1.d20) write(*,*)'potanal p0',p,plrf(1),plrf(2),plrf(3)
      if(p.gt.1.d20) stop

cc      out1 = nuclpot(icomp,rho0, rho,p )
cc      out2 = nuclpot3(icomp,rho,plrf,rhap )

c      if(abs(out1-out2).gt.1.0e-06) then
c        write(*,*)'deviation in pot : ' , out1, out2
c        write(*,*)'plrf = ', plrf , p
c        write(*,*)'r = ', rhap
c        write(*,*)'rho = ', rho
c      end if

c      if((ideltapot.eq.0).or.(ideltapot.eq.1 .and. nuclres)) then
      if(nuclres .or. (deltares .and. ideltapot.eq.0) .or.
     &        (strange .and.istrangepot.eq.0)) then
c        if(ithree.eq.0) then
          p = sqrt(plrf(1)**2 + plrf(2)**2 + plrf(3)**2)
          if(p.gt.1.d20) write(*,*)'potanal p',p,plrf(1),plrf(2),plrf(3)
          out = nuclpot(icomp,rho0, rho,p )
c          write(*,*) 'potanal nuclpot',out,icomp,rho0,rho,p
c        else if(ithree.eq.1)then
c          out = nuclpot3(icomp,rho,plrf,rhap )
c        end if
      else if(deltares .and. ideltapot.eq.1) then
c        out = deltapot(rho,p)
        write(*,*)'ideltapot has to be zero'
        stop
      else if(strange .and. istrangepot.eq.1) then ! Lambda or Sigma
        out = 2./3. * nuclpot(icomp,rho0, rho,p )
      else if(strange .and. istrangepot.eq.2) then ! Lambda or Sigma
        out = strangeLpot(idres-nres,icomp,rho0,rho,p )
      end if
c 10    continue

c       write(600,*) out
c       write(601,*) out1, out2

      potanal = out
      return
      end

************************************************************************
      real*8 function  nuclpot(icomp,rho0, rho,pin )

*        this function provides the analytical momdep. potential

      implicit none
      include "potparam"

      real*8   rho0,  pi
      parameter(pi = 3.1415)

      integer icomp, isum
      parameter(isum   = 5)

      real*8      aval, bval, cval, rho, tauval
      real*8      pfermi0, pfermi, lambda, p
      real*8      pot, t1, t2, t3, t4, t5, t6, t7
      real*8      sky
      real*8      xtest, temp, zwi, pin
c      real*8      ebind, ener

      integer   isu

chw  test  hw
c     nuclpot = .0
c     return

      pot = 0.0
      if(icomp.gt.maxset) then
        write(*,*) 'problem in potanal '
        write(*,*)'icomp = ',icomp
        stop
      end if


        aval   = a(icomp)
        bval   = b(icomp)
        tauval = tau(icomp)
        cval   = c(icomp)
        lambda = la(icomp)


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*     convert p in units 1/fm                                        *
       p = pin/0.197
       pfermi0 = (3./2.*pi*pi*rho0)**(1./3.)
       pfermi  = (3./2.*pi*pi*rho )**(1./3.)

       if(p.lt.1.0e-10) then
*         write(*,*)'nuclpot p changed '
         p = 1.0e-10
       end if
*      lambda  = 1.58*pfermi0
*
*     further on everything is calculated in 1 / fm
*     the constnts a,b,c convert then the potential in gev


        sky = aval*(rho/rho0) + bval*(rho/rho0)**tauval

        if(icomp.eq.3 .or. icomp.eq.4 .or. icomp.eq.6) then
          pot = sky
        else if(icomp.le.2 .or. icomp.eq.5) then

*      determine weather the small momentum expansion has to be used
          xtest = 2.0*pfermi*p/(pfermi**2+lambda*2)

          if(xtest .gt. 1.0e-06)then

            t1 = pi*lambda**3
            t1 = t1*4.0/(2.0*pi)**3
            t2 = (pfermi**2+lambda**2-p**2)/(2.0*p*lambda)
            t3 = (p+pfermi)**2+lambda**2
            t4 = (p-pfermi)**2+lambda**2
            t5 = 2.0*pfermi/lambda
            t6 = (p+pfermi)/lambda
            t7 = (p-pfermi)/lambda
           pot = t2*log(t3/t4) + t5 - 2.0*(atan(t6)-atan(t7))
c           if(abs((atan(t6)-atan(t7))-
c     &  atan(2.0*pfermi*lambda/(lambda**2-pfermi**2+p**2))).gt.1.d-4) 
c     &  write(*,*)'nuclpot atan', atan(t6)-atan(t7),
c     &  atan(2.0*pfermi*lambda/(lambda**2+pfermi**2-p**2))
c           write(*,*) 'nuclpot 1',t6,t7,pot,p,pin
           pot = pot*t1
           pot = pot*2.0*cval/rho0
           pot = pot  + sky

          else if(xtest .le. 1.0e-06) then
           t1   = pi*lambda**3
           t1   = t1*4.0/(2.0*pi)**3
           t2   = (pfermi**2 + lambda**2 - p**2)/lambda
           t3   = 0.0
           temp = 2.0*pfermi/(pfermi**2+lambda**2)
           do isu = 0, isum
             zwi = p**float(2*isu) * temp**(float(2*isu+1))
             t3 = t3 + zwi/float(2*isu+1)
           end do
           t5 = 2.0*pfermi/lambda
           t6 = (p+pfermi)/lambda
           t7 = (p-pfermi)/lambda

           pot = t2*t3+t5-2.0*(atan(t6)-atan(t7))
c           write(*,*) 'nuclpot 2',t6,t7,pot,pfermi
           pot = pot*t1
           pot = pot*2.0*cval/rho0
           pot = pot + sky

        end if
      end if
      nuclpot = pot

      return
      end

************************************************************************
      function strangeLpot(id1,icomp,rho0,rho,pp )

*        this function provides the analytical momdep. potential

      implicit none

      real*8 strangeLpot,rho,pp,rho0,nuclpot
      integer id1,icomp
      

c      real*8 potL,potXzp,potXmp,potSpp,potSzp,potSmp
c      parameter(potL = -0.170318) !GeV fm^3
c      parameter(potXmn = 0.226251, potXmp = -0.003375) !GeV fm^3
c      parameter(potSpp=0.254552, potSpn=0.226251, potSzp=0.240402) !GeV fm^3

      strangeLpot = .0
      if(id1.eq.1) then
        strangeLpot = 0.47*nuclpot(icomp,rho0, rho,pp )
      else if(id1.eq.2)  then
        strangeLpot = -0.005*rho/rho0
      else if(id1.eq.3)  then
        strangeLpot = 0.01*rho/rho0
      end if

      return
      end

********************************************************************
      real*8 function  potmes( rho,p,id1)
*     do the meson potentials

      implicit none

      include"cominput"

      integer id1
      real*8     rho, out, p

      out  = 0.0

c      if(ipipot.eq.1 .and. id1.eq.1) then
c        out = pionpotmit(rho,p)
c      else if(ipipot.eq.2 .and. id1.eq.1) then
c        out = pionpotohn(rho,p)
c      end if


      potmes = out

      return
      end
******************************************************************
      real*8 function  mdpart(icomp,rho0, rho,pin )

*        this function provides the analytical momdep. potential

      implicit none

      include "potparam"

      real*8   rho0,  pi
      parameter(pi = 3.1415)

      integer icomp, isum
      parameter(isum   = 5)

      real*8      aval, bval, cval, rho, tauval
      real*8      pfermi0, pfermi, lambda, p
      real*8      pot, t1, t2, t3, t4, t5, t6, t7
      real*8      xtest, temp, zwi, pin
c      real*8      ebind, ener


      integer   isu

      pot = 0.0
      if(icomp.gt.maxset) then
        write(*,*) 'problem in potanal '
        write(*,*)'icomp = ',icomp
        stop
      end if

      aval   = a(icomp)
      bval   = b(icomp)
      tauval = tau(icomp)
      cval   = c(icomp)
      lambda = la(icomp)


*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*     convert p in units 1/fm                                        *
       p = pin/0.197
       pfermi0 = (3./2.*pi*pi*rho0)**(1./3.)
       pfermi  = (3./2.*pi*pi*rho )**(1./3.)

       if(p.lt.1.0e-10) then
*         write(*,*)'nuclpot p changed '
         p = 1.0e-10
       end if
*      lambda  = 1.58*pfermi0



*
*     further on everything is calculated in 1 / fm
*     the constnts a,b,c convert then the potential in gev


        if(icomp.eq.3 .or. icomp.eq.4 .or. icomp.eq.6) then
          write(*,*)'icomp in mdpart = ', icomp
          stop
        else if(icomp.le.2 .or. icomp.eq.5) then

*      determine weather the small momentum expansion has to be used
          xtest = 2.0*pfermi*p/(pfermi**2+lambda*2)

          if(xtest .gt. 1.0e-06)then

            t1 = pi*lambda**3
            t1 = t1*4.0/(2.0*pi)**3
            t2 = (pfermi**2+lambda**2-p**2)/(2.0*p*lambda)
            t3 = (p+pfermi)**2+lambda**2
            t4 = (p-pfermi)**2+lambda**2
            t5 = 2.0*pfermi/lambda
            t6 = (p+pfermi)/lambda
            t7 = (p-pfermi)/lambda
           pot = t2*log(t3/t4) + t5 - 2.0*(atan(t6)-atan(t7))
           pot = pot*t1
           pot = pot*2.0*cval/rho0

          else if(xtest .le. 1.0e-06) then
           t1   = pi*lambda**3
           t1   = t1*4.0/(2.0*pi)**3
           t2   = (pfermi**2 + lambda**2 - p**2)/lambda
           t3   = 0.0
           temp = 2.0*pfermi/(pfermi**2+lambda**2)
           do isu = 0, isum
             zwi = p**float(2*isu) * temp**(float(2*isu+1))
             t3 = t3 + zwi/float(2*isu+1)
           end do
           t5 = 2.0*pfermi/lambda
           t6 = (p+pfermi)/lambda
           t7 = (p-pfermi)/lambda

           pot = t2*t3+t5-2.0*(atan(t6)-atan(t7))
           pot = pot*t1
           pot = pot*2.0*cval/rho0

          end if

        end if

          mdpart = pot
c        write(*,*)'open up mdpart '
c        stop

      return
      end

*****************************************************************
*     the quick mdpart
******************************************************************
      real*8 function  mdpart1(icomp,rho0, rho,pin )

*        this function provides the analytical momdep. potential

      implicit none

      include "potparam"

      real*8   rho0,  pi, mdpartini
      parameter(pi = 3.1415)


      integer imaxp
      parameter(imaxp = 2000)
      integer imaxd
      parameter(imaxd =140)
      real*8     deltad
      parameter(deltad = 0.005)
      real*8     deltap
      parameter(deltap = 0.001)


      real*8 mdfield(0:imaxd, 0:imaxp)
      integer ilaufp, ilaufd
      real*8    laufp, laufd


      integer icomp, isum
      parameter(isum   = 5)

      real*8      aval, bval, cval, rho, tauval
      real*8      pfermi0, pfermi, lambda, p
      real*8      pot, t1, t2, t3, t4, t5, t6, t7
      real*8      xtest, temp, zwi, pin


      integer   isu

      save mdfield, a, b, c, la, tau

      pot = 0.0
      if(icomp.gt.maxset) then
        write(*,*) 'problem in potanal '
        write(*,*)'icomp = ',icomp
        stop
      end if


c      if(rho.le.1.0e-06) then
c        mdpart = 0.0
c        write(*,*)'open quick mdpart'
c        stop
c        return
c      end if

      ilaufd = nint(rho/deltad)
      ilaufp = nint(pin/deltap)

      if(ilaufd .gt. imaxd) ilaufd = imaxd
      if(ilaufp .gt. imaxp) ilaufp = imaxp

      mdpart1 = mdfield(ilaufd, ilaufp)
      write(*,*)'open quick mdpart'
      stop

      return


*************************************************
      entry  mdpartini(rho0, icomp)

      write(*,*)'vor mdpartini 1', icomp


        aval   = a(icomp)
        bval   = b(icomp)
        tauval = tau(icomp)
        cval   = c(icomp)
        lambda = la(icomp)




      do ilaufd = 0, imaxd
        do ilaufp = 0, imaxp
          mdfield(ilaufd, ilaufp) = 0.0
        end do
      end do
      write(*,*)'vor mdpartini 2'

      do ilaufd = 0, imaxd


       laufd = float(ilaufd)*deltad


        do ilaufp = 0, imaxp

          laufp = float(ilaufp)*deltap

****************** evaluate the mm-dep parts





*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
*     convert p in units 1/fm                                        *
       p = laufp/0.197
       pfermi0 = (3./2.*pi*pi*rho0)**(1./3.)
       pfermi  = (3./2.*pi*pi*laufd )**(1./3.)

       if(p.lt.1.0e-10) then
*         write(*,*)'nuclpot p changed '
         p = 1.0e-10
       end if
*      lambda  = 1.58*pfermi0



*
*     further on everything is calculated in 1 / fm
*     the constnts a,b,c convert then the potential in gev


        if(icomp.eq.3 .or. icomp.eq.4 .or. icomp.eq.6) then
          write(*,*)'icomp in mdpart = ', icomp
          stop
        else if(icomp.le.2 .or. icomp.eq.5) then

*      determine weather the small momentum expansion has to be used
          xtest = 2.0*pfermi*p/(pfermi**2+lambda*2)

          if(xtest .gt. 1.0e-06)then

            t1 = pi*lambda**3
            t1 = t1*4.0/(2.0*pi)**3
            t2 = (pfermi**2+lambda**2-p**2)/(2.0*p*lambda)
            t3 = (p+pfermi)**2+lambda**2
            t4 = (p-pfermi)**2+lambda**2
            t5 = 2.0*pfermi/lambda
            t6 = (p+pfermi)/lambda
            t7 = (p-pfermi)/lambda
           pot = t2*log(t3/t4) + t5 - 2.0*(atan(t6)-atan(t7))
           pot = pot*t1
           pot = pot*2.0*cval/rho0

          else if(xtest .le. 1.0e-06) then
           t1   = pi*lambda**3
           t1   = t1*4.0/(2.0*pi)**3
           t2   = (pfermi**2 + lambda**2 - p**2)/lambda
           t3   = 0.0
           temp = 2.0*pfermi/(pfermi**2+lambda**2)
           do isu = 0, isum
             zwi = p**float(2*isu) * temp**(float(2*isu+1))
             t3 = t3 + zwi/float(2*isu+1)
           end do
           t5 = 2.0*pfermi/lambda
           t6 = (p+pfermi)/lambda
           t7 = (p-pfermi)/lambda

           pot = t2*t3+t5-2.0*(atan(t6)-atan(t7))
           pot = pot*t1
           pot = pot*2.0*cval/rho0

          end if

        end if
          pot = pot

***************************
          mdfield(ilaufd, ilaufp) = pot
c       if (ilaufd/10*10 .eq. ilaufd
c    1      .and. ilaufp/10 * 10 .eq. ilaufp)
c    2      write(*,*) ' in mdpartini ', ilaufd, ilaufp, pot,p
        end do
      end do
          mdpartini = .0
          write(*,*)   '  end of mdpartini '
      return

      end



