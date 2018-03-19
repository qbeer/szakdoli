
************************************************************************
*                                                                      *
      subroutine inip(minnum,maxnum,z0,p3,p1,gamma,beta)
*                                                                      *
*       purpose:     providing initial conditions for momentum         *
*                    distribution of testparticles                     *
*       variables:   (all input)                                       *
*         minnum  - first testparticle treated in one run    (integer) *
*         maxnum  - last testparticle treated in one run     (integer) *
*         num     - number of testparticles per nucleon      (integer) *
*         z0      - displacement of center of nucleus in               *
*                   z-direction "fm"                            (real) *
*         p3,p1   - momentum-boost  "gev/c"                     (real) *
*         gamma   - relativistic gamma-factor                   (real) *
*         iseed   - seed for random-number generator         (integer) *
*         icoll   - =0 buu, =1 meanfield, =-1 cascade        (integer) *
*                                                                      *
************************************************************************
      implicit none

      real*8 rhow0
      parameter     (rhow0 = 0.168)
      include"common"
      include"cominput"
*
      real*8    ptot(3), rhows, rhowsn, rhowsp, pfermi, pfermin
      real*8    pfermip, px, py, pz, psqr, gamma, z0, p3, p1, rn
      real*8    pdens(3,35),pstep,pabs
      external rn
      integer irun, idir, ix, iy, iz, np, ij
      integer minnum, maxnum, npart, i
      real*8    beta
      integer init
      real*8 fmax, pdeut, aa, f, q, x, x23, y, seta, zeta
      logical  deuteron
*----------------------------------- ----------------------------------
c     radiu =1.124*(float(maxnum-minnum+1))**1./3.
*
*     loop over all parallel runs:
*
      write(*,*) ' inip',minnum,maxnum,z0,p3,p1,gamma,beta
      deuteron = maxnum .eq. minnum+1
      do i = 1,35
        do ij = 1,3
          pdens(ij,i) = 0.0
        enddo
      enddo

      do 1000 irun = 1,num
        init     = (irun-1)*maxb
*
        if(icoll.ne.-1 .and. minnum.lt.maxnum) then
          do 600 i = minnum+init,maxnum+init
            if(id(1,i).eq.0)      goto 600
            if (maxnum-minnum .eq. 1) then         !    deuteron
              if (init + minnum .eq. i) then         !    deuteron
                aa = 6.16                        !   Hulthen wave fct
                pdeut = 0.231 * .1973                     !  GeV/c
                fmax = .00162
   50           continue
                x = rn(iseed)
                y = (1.-x)**2
                x23 = x**(2./3.)
                f = (1.+2.*x)*y*y/((x23+y)*(x23+aa*aa*y))**2
                if (f  .lt. fmax*rn(iseed) ) goto 50
                q = pdeut * x**(1./3.)/(1-x)
                zeta = (2.*rn(iseed) - 1.)
                seta = sqrt(abs(1.-zeta**2))
                pz  = q* zeta
                x =  (6.283*rn(iseed))
                py =  q * seta * sin(x)
                px =  q * seta * cos(x)
                p(1,i) = px
                p(2,i) = py
                p(3,i) = pz
                p(1,i+1) = -px
                p(2,i+1) = -py
                p(3,i+1) = -pz
              endif
              goto 600
            endif
            ix=nint(r(1,i))
            iy=nint(r(2,i))
            iz=nint(r(3,i))
            if(ix.gt. maxx) ix=maxx
            if(ix.lt.-maxx) ix=-maxx
            if(iy.gt. maxx) iy=maxx
            if(iy.lt.-maxx) iy=-maxx
            if(iz.gt. maxz) iz=maxz
            if(iz.lt.-maxz) iz=-maxz
*
*     local fermi momentum by local thomas fermi
            rhows   = rhb(ix,iy,iz)
            rhowsn  = rhob_4(5,ix,iy,iz)
            rhowsp  = rhob_4(4,ix,iy,iz)
c            write(*,*) 'inip1 ',rhows,rhowsn,rhowsp
            pfermi  = 0.19732*(3.0/2.0*pi*pi * rhows )**(1./3.)
            pfermin = 0.19732*(3.0    *pi*pi * rhowsn)**(1./3.)
            pfermip = 0.19732*(3.0    *pi*pi * rhowsp)**(1./3.)
************************fuer ANKE ohne fermiimpuls***********
!		pfermi = 0.
!		pfermin = 0.
!		pfermip = 0.
*************************************************************

*
  500       continue
              px = 1.0 - 2.0 * rn(iseed)
              py = 1.0 - 2.0 * rn(iseed)
              pz = 1.0 - 2.0 * rn(iseed)

              psqr=px*px+py*py+pz*pz
              if (psqr .gt. 1.0) goto 500
*
              if(id(2,i).eq.0) then
                p(1,i) = pfermin * px
                p(2,i) = pfermin * py
                p(3,i) = pfermin * pz
              else if(id(2,i).eq.1)then
                p(1,i) = pfermip * px
                p(2,i) = pfermip * py
                p(3,i) = pfermip * pz
              end if
*
  600       continue
*
*         set total momentum to 0 in rest frame and boost
*
c            write(*,*) 'inip1 '
            do 700 idir = 1,3
              ptot(idir) = 0.0
  700       continue
            npart = 0
            do 900 i = minnum+init,maxnum+init
              npart = npart + 1
              do 800 idir = 1,3
                ptot(idir) = ptot(idir) + p(idir,i)
  800         continue
  900       continue
            do 950 i = minnum+init,maxnum+init
              do 925 idir = 1,3
                p(idir,i) = p(idir,i) - ptot(idir) / float(npart)
  925         continue

              pstep = 0.02
              pabs = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
              np = min(34, nint(pabs/pstep) ) + 1
             pdens(1,np)=pdens(1,np)+1./float(num)/(4.*pi*pabs**2*pstep)
              if(id(2,i).eq.1)
     &       pdens(2,np)=pdens(2,np)+1./float(num)/(4.*pi*pabs**2*pstep)
              if(id(2,i).eq.0)
     &       pdens(3,np)=pdens(3,np)+1./float(num)/(4.*pi*pabs**2*pstep)

*           booost
              p(3,i) = gamma * p(3,i) + p3
     &        * sqrt( 1.0 + (p(1,i)**2+p(2,i)**2+p(3,i)**2)/rmass**2 )
              if (deuteron) then        !       1.0 GeV/c     limit  !!
                if (p(3,i)-p3 .gt. 1.0)p(3,i) = p3 + 1.0
               if (p(3,i)-p3 .lt.-1.0)p(3,i) = p3 - 1.0
              endif
              p(1,i) = p(1,i) + p1

******** used to be a test, at least at 1 GeV/A no difference **********
c            ppart3 = p(3,i)
c            enpart = sqrt(rmass**2 + p(1,i)**2 + p(2,i)**2 + p(3,i)**2)
c            p(3,i) = gamma*(p(3,i) + beta*enpart)
**************************************************************************

              r(3,i) = (r(3,i)-z0) / gamma + z0

  950       continue
*
*     coherent boost for cascade
          else
            do 935 i = minnum+init,maxnum+init
              p(1,i) = p1
              p(2,i) = 0.0
              p(3,i) = p3

              r(3,i) = (r(3,i)-z0) / gamma + z0
  935       continue
          end if
 1000   continue

cc      if(minnum.eq.1) then
cc        write(isum,'(''#c: Momentum density for the target'')')
cc      else
cc        write(isum,'(''#c: Momentum density for the projectile'')')
cc      endif
cc      write(isum,'(''#x: momentum [GeV]'')')
cc      write(isum,'(''#y: rho(p)=dN/dp**3 [1/GeV**3]'')')
cc      write(isum,'(''#c:    p     rho(tot)     proton      neutron'')')
cc      do np=1,20
cc       write(isum,'(f8.3,3e12.4)') float(np-1)*pstep,(pdens(i,np),i=1,3)
cc      enddo
!              write(54,*) rhowsp
      return
      end

