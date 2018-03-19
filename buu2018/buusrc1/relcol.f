************************************************************************
*                                                                      *
      subroutine relcol(lcoll, nt)

*                                                                      *
*       purpose:    calculating the kinematics in a collision between  *
*                   two particles - relativistic formula used          *
*                   -------------------------------------------------  *
*                   in this sub production part is called !!!          *
*                   -------------------------------------------------  *
*                                                                      *
*       variables:                                                     *
*                                                                      *
*         lcoll   - number of collisions              (integer,output) *
*         deltar  - maximum spatial distance for which a collision     *
*                   still can occur                                    *
*                                                                      *
************************************************************************
      implicit none
      include"common"
      include"cominput"
      include 'com_kminu'

      common /nthwhw/  nthw
      integer nthw
      common /nameresc/nameres(1:nres+3)
      character *8  nameres

      real*8 sig0, bmaxd, sigd, deltar, bmax0, pirkaon
      real*8 pbarlab2,crosspbarann,rn
      integer lcodim
      parameter(lcodim = 3*nres+6+nres**2)

      parameter(bmax0=1.323142)
      parameter(sig0=55.0)
      parameter(bmaxd=1.6)
      parameter(sigd=80.4)
      parameter(deltar=2.0)
      parameter(pirkaon=0.8)

      integer idim
      parameter(idim = 4) !!!!!!!! has to be cons. with commsp

      real*8    pcm(3), beta(3), pcmi(3), pot1
      real*8    x1, y1, z1, px1, py1, pz1, em1, e1, e10, e20
      integer id11, id21, id61, ibar, ippbar
      real*8    x2, y2, z2, px2, py2, pz2, em2, e2
      real*8    psum2
      integer id12, id22, id62
      real*8    dx, dy, dz, rsqare, sigkaon
      real*8    s, srt, cutoff, srt_0, srt_p
      real*8    p12, p1dr, p2dr, a12, b12, c12, brel, bmax, sig,
     &        b21, t1, t2, xxx, yyy, zzz
      real*8    etotal, ptotx, ptoty, ptotz, gamma, traf
      real*8    p1beta, transf,  pcm2, prcm, e1cm, e2cm
      real*8    ecm, px, py, pz
      real*8    phase                             !needed in pauli
      real*8    j01, j02, help, meff1, meff2, help1, help2
      integer nt, i, test

      real*8    rpi3(1:3), ppi3(1:3),pnucl1(1:3),pnucl2(1:3)
      real*8    upi
      integer izpi
      logical testflag,sstate

      integer lcoll(-6:lcodim)
      integer irun, j1, j2 , i1, i2 , iblock, ntag
      integer kk, idres

c      integer ncputim
      integer icptimep
      common /cputime/ cptimepr, cputim0, cputim, clockmax
      real*8 cptimepr,cputim0, cputim, clockmax, second
      external second

*     stuff needed for the mom-dep. pots
      real*8 betlrfx(1:2), betlrfy(1:2), betlrfz(1:2)
      real*8 scapo1, scapo2, scapo3, scapo4,cutpot1, cutpot2
      real*8 pin, potanal, epicm,t0
      integer inp, ib

c      real*8 tin, tout
      real*8 x,y,z, deriv(0:4),j0me, j1me, j2me,j3me
      real*8 mbetlrfx, mbetlrfy,mbetlrfz
      integer qtest

      integer nuba(-1:nres+3)
      character *2  nst1, de1, fi, lam, sigm
      character *1  endst
      integer repro, rescheck1, rescheck2
      real*8 resmass
      real*8 crx, cry, crz, cpx, cpy, cpz, cpot
      real*8 cfox, cfoy, cfoz, etot
      integer cid2

      real*8 rhap1(1:3), rhap2(1:3) , plrf1(1:3), plrf2(1:3)


*-----------------------------------------------------------------------
*     initialization of counting variables
*
*      new meaning of countig vaiables !!!
*      built in : 05.02.94 (st)
*
*                i =  -6 ; number of bremsstrahlungs events. (lcbre)
*                i =  -5 ; number of colls target-projectile (lcopt)
*                i =  -4 ; d + d (elast.) (lcodd)
*                i =  -3 ; n + d (elast.) (lcond)
*                i =  -2 ; proton-neutron collision (elast.) (lcopn)
*                i =  -1 ; a collision has taken place
*                i =   0 ; nothing has happened
*                i =   1 ; elastic n-n collision
*                i =   2 ; n + n -> n + delta
*                i =   3 ; n + delta -> n + n
*                i =   4 ; n + n -> n + nstar
*                i =   5 ; n + nstar -> n + n
*                i =   6 ; n + n -> n + nstar2
*                i =   7 ; n + nstar2-> n + n
*                i =   8 ; n + d -> n + nstar1
*                i =   9 ; n + nstar1 -> n + d
*                i =  10 ; n + d -> n + nstar2
*                i =  11 ; n + nstar2 -> n + d
*                i =  12 ; n + nstar1 -> n + nstar2
*                i =  13 ; n + nstar2 -> n + nstar1
*                i =  14 ; d + d -> n + nstar1
*                i =  15 ; d + d -> n + nstar2
*
************************************************************************
c      write(*,*)'resabs auskommentiert !!!!!!!!!!!!!!!!!!! relcol'
      t0 = 0.0
c      tin = secnds(t0)
      do kk =-6,lcodim
        lcoll(kk) = 0
      end do
*      build strings for particle names

c            write(*,*)'relcol1 ', r(1,6424),r(2,6424),r(3,6424)
      write(*,*)'legende[5~ ++++++++++++++++++++++++'
 100  format(A,F5.0,A)
      endst = ')'
      nst1 = 'N('
      de1  = 'D('
      lam  = 'L('
      sigm = 'S('
      do i = 1,nres+3
        if(i.eq.1) then
          repro = 1
          resmass = rmass*1000.
          fi = nst1
        else if(i.eq.nres+2) then
          fi = lam
          resmass = xlmas*1000.
        else if(i.eq.nres+3) then
          fi = sigm
          resmass = xsmas*1000.
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

      do i = -1, nres+3
        nuba(i) = 0
      end do

      do i = 1, maxpar
        if(id(1,i).ne.0) then
          nuba(id(1,i)) = nuba(id(1,i)) + 1
        end if
      end do

      write(*,*)'******************************'
      write(*,*)'number of res after coll'
      do kk = 1, nres+3
        write(*,*)nameres(kk), nuba(kk)
      end do
      write(*,*) 'pbar    ', nuba(-1)
      write(*,*)'******************************'

      test = 0
************************************************************************
*    id -legend
*  id(1,i) = particle type (nucleon,delta, nstar...)
*  id(2,i) = particle charge
*  id(3,i) = last colliding partner of particle i
*  id(4,i) = for resonances: how many times the (abs.)meson was created
*  id(5,i) = number of pion absorption
*  id(6,i) = abs: number of barion collision; sign: target or projectile
*
************************************************************************
*
*
*     calculate the potentials used for the mom.dep stuff
*
*
      call dens_4

      write(*,*)   '  nach dens_4 vor  potcalc ', iseed
      write(*,*)   '  rhob ', rhob_4(0,0,0,0), rhob_4(0,5,5,0)
      call potcalc
      write(*,*)   ' ende  potcalc ', iseed

*
*
************************************************************************

*  loop over all parallel runs
      do 1000 irun = 1,num
        ippbar = 0
c        call otime(ncputim)
         cputim=second() - cputim0
         if(cputim.lt.cptimepr) icptimep = icptimep + 1
         cptimepr = cputim
c         cputim = cputim + float(icptimep) * clockmax
*
*  loop over all pseudoparticles 1 in the same run
        do 800 j1 = 2,maxb
          i1  = j1 + (irun - 1) * maxb
          if(id(1,i1).eq.0) goto 800
*         store positions, momenta, mass and ids of particle 1
          x1        = r(1,i1)
          y1        = r(2,i1)
          z1        = r(3,i1)
          px1       = p(1,i1)
          py1       = p(2,i1)
          pz1       = p(3,i1)
          em1       = e(i1)
          scapo1    = upot(i1)
          meff1     = em1 +  scapo1
          e1        = sqrt( meff1**2 + px1**2 + py1**2 + pz1**2 )
          e10       = sqrt( em1**2   + px1**2 + py1**2 + pz1**2 )
c          write(*,*) 'relcol 1',px1,py1,pz1,em1,e10 
c
          id61       = id(6,i1)
          id11       = id(1,i1)
          id21       = id(2,i1)
          betlrfx(1) = betlrfboo(i1,1)
          betlrfy(1) = betlrfboo(i1,2)
          betlrfz(1) = betlrfboo(i1,3)
          j01        = betlrfboo(i1,4)
          rhap1(1)   = x1
          rhap1(2)   = y1
          rhap1(3)   = z1

*
*  loop over all pseudoparticles 2 in the same run
          do 600 j2 = 1,j1-1
            i2  = j2 + (irun - 1) * maxb
            if(id(1,i2).eq.0) goto 600
            id62= id(6,i2)
*
c            write(*,*)'relcol 2', i1,i2,id61,id62,id(5,i1),id(5,i2)
*-----------------------------------------------------------------------
*   avoid first collisions within the same nucleus (aichelin):
            if((id61*id62.eq.iavoid) .and.
     &      ( (id(5,i1)+id(5,i2)) .eq. 0) )                    goto 400
c            write(*,*)' relc3a', i1,i2,id61,id62,id(5,i1),id(5,i2)
*
*   avoid second collisions for the same pairs (important):
            if ( (id(3,i1).eq.i2) .and. (id(3,i2).eq.i1) )     goto 400
*
*-----------------------------------------------------------------------
*   the following prescription (deltar) is not covariant, but useful
*   to save cpu time especially below 200 mev/u
*
*   we better perform a test on this feature , just to make sure !!!!!
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c            if(i2.ge.0) write(*,*)' relcol 3', i1,i2
            x2     = r(1,i2)
            dx     = x1 - x2
            if (abs(dx) .gt. deltar)                           goto 400
            y2     = r(2,i2)
            dy     = y1 - y2
            if (abs(dy) .gt. deltar)                           goto 400
            z2     = r(3,i2)
            dz     = z1 - z2
            if (abs(dz) .gt. deltar)                           goto 400
            rsqare = dx**2 + dy**2 + dz**2
            if (rsqare .gt. deltar**2)                         goto 400
*
c            if(i2.ge.0) write(*,*)' relcol 4', i1,i2
*   now particles are close enough to each other !
*-----------------------------------------------------------------------
            px2    = p(1,i2)
            py2    = p(2,i2)
            pz2    = p(3,i2)
            em2    = e(i2)
            scapo2 = upot(i2)
            meff2  = em2 + scapo2
            e2     = sqrt ( meff2**2 + px2**2 + py2**2 + pz2**2 )
            e20    = sqrt ( em2**2 + px2**2 + py2**2 + pz2**2 )
            e10    = sqrt( em1**2   + px1**2 + py1**2 + pz1**2 )
c            write(*,*) 'relcol 2',px2,py2,pz2,em2,e20 
c          write(*,*) 'relcol 1',px1,py1,pz1,em1,e10 
            id12   = id(1,i2)
            id22   = id(2,i2)
            betlrfx(2) = betlrfboo(i2,1)
            betlrfy(2) = betlrfboo(i2,2)
            betlrfz(2) = betlrfboo(i2,3)
            j02        = betlrfboo(i2,4)
            rhap2(1)   = x2
            rhap2(2)   = y2
            rhap2(3)   = z2

            if(nt.eq.-700 .and. i2.ge.0) write(*,*)' relcol  4'

c            if(id(1,i1).gt.4 .or. id(1,i1).lt.1 .or.
c     &         (id(1,i1).ne.2.and.(id(2,i1)+2)/2.ne.1) )
c     &         write(*,'(''hiba relcol1'',2i8,2f8.3)')
c     &                             id(1,i1),id(2,i1),e(i1),time
c            if(id(1,i2).gt.4 .or. id(1,i2).lt.1 .or.
c     &         (id(1,i2).ne.2.and.(id(2,i2)+2)/2.ne.1) )
c     &        write(*,'(''hiba relcol1'',2i8,2f8.3)')
c     &                             id(1,i2),id(2,i2),e(i2),time
*-----------------------------------------------------------------------
*   calculate s and invariant energy of particle pair
            psum2 = (px1+px2)**2 + (py1+py2)**2 + (pz1+pz2)**2
            s      = (e1+e2)**2 -  psum2
            srt    = sqrt(s)
            srt_p  = (e10+e20)**2 - psum2
            if(srt_p.le.(em1+em2)**2) then
              write(*,*) 'relcol hiba s < ',srt_p,(e10+e20)**2,psum2,
     &       px1,py1,pz1,em1,e10,px2,py2,pz2,em2,e20
              srt_p = (em1+em2)**2
              stop
            end if
            srt_p  = sqrt(srt_p) + upot(i1) + upot(i2)
c     if (irun .eq. 16)
c      write(*,*) ' in relcol sqrt ', srt, srt_p, (e10+e20)**2,psum2,
c     2      i1, i2, e1, e2, upot(i1), upot(i2)
c     &      ,e10,px1,py1,pz1,e20,px2,py2,pz2,em1,em2
*----------------------------------------------------------------------
*
*   low energy cutoff !
*           1.8966  =  rmas1 + rmas2 + 0.02 gev     (nucleon case)
*
            pin = 0.0
            plrf1(1) = 0.0
            plrf1(2) = 0.0
            plrf1(3) = 0.0

            plrf2(1) = 0.0
            plrf2(2) = 0.0
            plrf2(3) = 0.0

            cutpot1 = potanal(rho0,j01,plrf1,1,rhap1)
            cutpot2 = potanal(rho0,j02,plrf2,1,rhap2)
            cutoff  =  2.0 * rmass + cutpot1+cutpot2 +0.02
c         if(nt.eq.-700  .and. i2.ge.0) write(*,*)' relcol  6'
            if (srt .lt. cutoff)                               goto 400
c          if(nt.eq.-700  .and. i2.ge.0) write(*,*)' relcol  7'
*-----------------------------------------------------------------------
*

*   now there is enough energy available !
*
*-----------------------------------------------------------------------
*   are their impact parameter small enough?
*     (description according to Kodama)

              p12    = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
          if(nt.eq.-700   .and. i2.ge.0) write(*,*)' 8a ', p12
              p1dr   = px1 * dx + py1 * dy + pz1 * dz
             p2dr   = px2 * dx + py2 * dy + pz2 * dz
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)' 8b ',p1dr, p2dr
              a12    = 1.0 - ( meff1 * meff2 / p12 ) ** 2
              b12    = p1dr / meff1 - p2dr * meff1 / p12
              c12    = rsqare + ( p1dr / meff1 )**2
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'8c ', a12, b12,c12
             brel   = sqrt( abs(c12 - b12**2/a12) )
              if(id11+id12.ge.3) then
                bmax = bmaxd
                sig  = sigd
              else
                bmax = bmax0
                sig  = sig0
              end if
            if (brel .gt. bmax)                                goto 400
c          if(nt.eq.-700  .and. i2.ge.0) write(*,*)' relcol  9'
*-----------------------------------------------------------------------
*   average time-shift of the collision in the fixed frame
*   will particles get closest point in this time interval ?
              b21    = - p2dr / meff2 + p1dr * meff2 / p12
              t1     = ( p1dr / meff1 - b12 / a12 ) * e1 / meff1
              t2     = ( - p2dr / meff2 - b21 / a12 ) * e2 / meff2
            if( abs(t1+t2).gt.dt)                              goto 400
c          if(nt.eq.-700  .and. i2.ge.0) write(*,*)' relcol 10'

*needed for Gyuri's routines see below
              xxx=(x1+x2)/2.0
              yyy=(y1+y2)/2.0
              zzz=(z1+z2)/2.0

c              if(id11.eq.-1.or.id12.eq.-1) write(*,*)'pbar coll relcol'
c     &            ,i_JPsi
*-----------------------------------------------------------------------
*   lorentz-transformation in i1-i2-c.m. system
              etotal = e1 + e2
c          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'etotal = ', etotal
              ptotx  = px1+px2
              ptoty  = py1+py2
              ptotz  = pz1+pz2
              beta(1) = (px1+px2) / etotal
              beta(2) = (py1+py2) / etotal
              beta(3) = (pz1+pz2) / etotal
c          if(nt.eq.-700   .and. i2.ge.0) write(*,*)' relcol 11'
              gamma  = 1.0 / sqrt(1.0-beta(1)**2-beta(2)**2-beta(3)**2)
              traf   = gamma / (gamma + 1.)
c          if(nt.eq.-700 .and. i2.ge.0) write(*,*)' relcol 12'

*       transformation of momenta (p1c(i) = - p2c(i))
              p1beta = px1*beta(1) + py1*beta(2) + pz1*beta(3)
              transf = gamma * ( traf * p1beta - e1 )
              pcm(1) = beta(1) * transf + px1
              pcm(2) = beta(2) * transf + py1
              pcm(3) = beta(3) * transf + pz1
*              potndcm= gamma * potnd
              pcmi(1)= pcm(1)
              pcmi(2)= pcm(2)
              pcmi(3)= pcm(3)
              pcm2   = pcm(1)**2 + pcm(2)**2 + pcm(3)**2
              prcm   = sqrt (pcm2)
              if (prcm .le. 0.00001)                           goto 400
*
              e1cm   = sqrt ( meff1**2 + pcm2 )
              e2cm   = sqrt ( meff2**2 + pcm2 )
              if(abs(srt-(e1cm+e2cm)).gt. 1.0e-04) then
                write(*,*)' in relcol 772',srt, e1cm, e2cm, e1cm+e2cm
                stop
              end if
c          if(nt.eq.-700   .and. i2.ge.0) write(*,*)' relcol 13'
*-----------------------------------------------------------------------
*         call of routines which evaluate the cross sections

*-----------------------------------------------------------------
* JPsi production by Wolf, Balassa, Kovacs, Zetenyi (2017)
*-----------------------------------------------------------------
         if(i_JPsi .ge. 1) then
           srt_0 = srt_p - 2.* pot1
           write(*,*) 'in relcol JPsi cal', id11,id12,sig
           if(id11*id12.gt.0) call JPsi_pert_NN(i1,i2,id11,id12,beta,
     $          srt_0,xxx,yyy,zzz,sig,irun)
           if(id11*id12.lt.0) then
c             pbarlab2 = (0.5*srt_0**2-rmass**2)**2/rmass**2-rmass**2
c             if(pbarlab2.le.0)                                 goto 400
c     Landolt-Börnstein_1.12/b p.286
c             crosspbarann = 0.532 + 63.4*pbarlab2**(-0.71/2.0)
             crosspbarann = 25.
             if(rn(iseed) .le. crosspbarann/sig)  then
c              call JPsi_pert_pbarN(i1,i2,id11,id12,beta,
c     $          srt_0,xxx,yyy,zzz,sig,irun)
               call JPsi_pert_pbarN(i1,i2,id11,id12,beta,
     $          srt_0,xxx,yyy,zzz,crosspbarann,irun)
               if(id11*id12.eq.-1) ippbar = ippbar+1
             end if
           end if
           write(*,*) 'in relcol after JPsi cal', id11,id12,sig

           if(i_JPsi.eq.3.or.i_JPsi.eq.4.or.i_JPsi.eq.7.or.i_JPsi.eq.8)
     &                                                          goto 800
         endif
*-----------------------------------------------------------------------
         if(iDmes .ge. 1) then
           srt_0 = srt_p - 2.* pot1
c           write(*,*) 'in relcol Dmes cal', id11,id12,sig
           call Nbarn_DD_pert(i1,i2,id11,id12,beta,srt_0,xxx,yyy,zzz,
     $          sig,irun)
         endif
*-----------------------------------------------------------------------
*            kaon production (gyuri)
*      write(*,*)   '    begin kaondbb    inside  relcol ', ikaon
*-----------------------------------------------------------------------
           if(i_phi .eq. 1  .or. ikaon .eq. 1) then
             ibar = 1
             pot1 = 0.0
c         write(*,*)"relcol potcalci",id(1,i1),i1,ibar,beta(1),beta(2),
c     &           beta(3)
         
             call potcalc_i(i1, ibar, beta(1), beta(2), beta(3), pot1)
c         if (nt.eq.-700  ) write(*,*) ' call phi_dbb with ',
c     1     id11,id12,beta,gamma,srt,xxx,yyy,zzz,sig,id61,id62,
c     2     irun,i1,i2, etotal
             if(i_phi .eq. 1  ) then
             srt_0 = srt_p - 2.*pot1
c       write(*,*) ' before phi_dbb ', upot(i1),upot(i2) , pot1
             call phi_dbb(id11,id12,beta,gamma,srt_0,xxx,yyy,zzz,sig,
     1         id61,irun,i1,i2, etotal,id21,id22)
           endif
           if(i_kminu .eq. 1) then
             srt_0 = srt_p - 2.* pot1
c           write(*,*) ' before kminu_dbb ',upot(i1),srt ,srt_p,pot1
             call kminu_dbb(id11,id12,beta,gamma,srt_0,xxx,yyy,zzz,sig,
     1       irun,i1,i2)
           endif

*-----------------------------------------------------------------


c         if(ikaon .eq. 1  .and. brel .le. pirkaon
c     1     .and.   iphi_dec .ne. 1)                        then
c           sigkaon = 10.*pirkaon**2 *pi
c           srt_0 = srt_p - 1.66 * pot1
c           call kaondbb(id11,id12,beta,gamma,srt_0,xxx,yyy,zzz,
c     1       sigkaon,id61,id62, irun,i1,i2, etotal)
c         end if
       end if
*-----------------------------------------------------------------------
*            vector meson production (gyuri)
c      if(imeson .eq. 1) then
c        call mesonbb(id11,id12,beta,gamma,srt,etotal,
c     &              xxx,yyy,zzz,sig)
c      end if
*-----------------------------------------------------------------------
*            omega production (gyuri)
c      if(iomega .eq. 1) then
c        call omegabb(id11,id12,beta,gamma,srt,etotal,xxx,yyy,zzz,sig,
c     &              id61,id62)
c      end if
*-----------------------------------------------------------------------
cc        if(igamma.eq.1 .and. id21*id22.eq.0 .and. id21+id22.eq.1)
cc     &     call gamnp(i1,i2,beta,gamma,srt,
cc     &                          xxx,yyy,zzz,e(i1),e(i2),etotal,sig)
c*
*-----------------------------------------------------------------------
*            perturbative pion production (gyuri)
c      if(((ipertpi.eq.1).or.(idilper.eq.1)) .and. (id11*id12.eq.1)
c     &         .and. (srt.gt.2.016)) then
c                call pertupi(pcm,srt,iseed,sig,prcm,
c     &                       i1,i2,id21,id22,xxx,yyy,zzz,beta,gamma)
c      end if
c
*-----------------------------------------------------------------------
*            dilepton production from p-n-bremsstrahlung
      if(ibrems.ge.1 .and. id11+id12.eq.2 .and. id21*id22.eq.0 .and.
     &        id21+id22.eq.1) then
        call npbrems(srt,sig,beta,gamma,xxx,yyy,zzz)   ! new (zm)

        ecm=float(id(2,i1))*e1cm+float(id(2,i2))*e2cm
c        call bremslep(srt,beta, e(i1),e(i2),ecm,sig)
c        lcoll(-6) = lcoll(-6) + 1
      end if
*
*-----------------------------------------------------------------------
*            dilepton production from Drell-Yan
c      write(*,*) 'call drell Yan',iDrellYan,id11,id12,srt,sig,beta,
c     &    gamma,xxx,yyy,zzz
      if(iDrellYan.ge.1 .and. id11+id12.eq.0 .and. id11*id12.eq.-1) then
c        write(*,*) 'call Drell Yan',srt,sig,beta,gamma,xxx,yyy,zzz
        call Drell_Yan(srt,sig,beta,gamma,xxx,yyy,zzz) !
      end if
*     
*-----------------------------------------------------------------------
*           baryon dynamics
          if(nt.eq.-700   .and. i2.ge.0) write(*,*)' relcol 14'

       help=sqrt((e(i1)+upot(i1))**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)+
     +     sqrt((e(i2)+upot(i2))**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)

          help = help**2 -(p(1,i1)+p(1,i2))**2 -
     +          (p(2,i1)+p(2,i2))**2 - (p(3,i1)+p(3,i2))**2
          help2 = sqrt(help)
          if(abs(help2-srt).gt.0.0001) then
            write(*,*)'relcol prob '
            write(*,*)srt,em1,em2,id11,id12, help
            write(*,*)irun,i1,i2,j1,j2
            stop
          end if


          if(nt.eq.-700    .and. i2.ge.0)   write(*,*)'vor crosw'
        test = test + 1
c        write(*,*)'id11, id12', id11, id12
c        write(*,*)'i1, i2 ', i1, i2 , id(1,i1), id(1,i2)
c        write(*,*)'em 1, em2 ' , em1, em2
c        write(*,*)'ladungen ', id21, id22

        qtest = id21 + id22
        testflag = .false.

**************************************************************
c         write(*,*)' relcol 15'
c         if(nt.eq.-700   .and. i2.ge.0) write(*,*)' relcol 15'
            call crosw(pcm,srt,em1,em2,id11,id12,iblock,sig,
     &            cutoff, prcm,i1,i2,id21,id22,
     &            scapo1, scapo2, scapo3, scapo4,j01, j02,
     &            beta, betlrfx, betlrfy, betlrfz,testflag,nt,
     &            sstate,ppi3,pnucl1,pnucl2,upi,izpi,rpi3,
     &            rhap1, rhap2)
c         write(*,*)' relcol 15e'
******************************************************
          if(nt.eq.-700   .and. i2.ge.0)write(*,*)'nach crosw',
     &           iblock  , em1, em2
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'id11, id12',
     &          id11, id12


*---------------- check the charges ---------------------------------*
            if(.not.sstate) then
              if(qtest.ne.id21+id22) then
                write(*,*)'relcol prob mit ladung '
                write(*,*)qtest, iblock
                write(*,*)id11, id12
                write(*,*)id12, id22
                stop
              end if
            else if(sstate) then
              if(qtest.ne.id21+id22+izpi) then
                write(*,*)'relcol prob mit ladung '
                write(*,*)qtest, iblock
                write(*,*)id11, id12
                write(*,*)id12, id22
                write(*,*)'pion ', izpi
                stop
              end if
            end if

            if(id11.gt.1 .and. id11.lt.nres+2) then
              if(resprop2(id11-1,1).eq.1) then
                if(id21.lt.0 .or. id21.gt.1) then
                  write(*,*)'prob mit ladg in relcol 1',qtest
                  write(*,*)id11, id21, iblock, em1
                  write(*,*)id12, id22, iblock, em2
                  stop
                end if
              else if(resprop2(id11-1,1).eq.3) then
                if(id21.lt.-1 .or. id21.gt.2) then
                  write(*,*)'prob mit ladg in relcol 2',qtest
                  write(*,*)id11, id21, iblock, em1
                  write(*,*)id12, id22, iblock, em2
                  stop
                end if
              end if
            else if (id11.eq.1) then
              if(id21.lt.0 .or. id21.gt.1) then
                write(*,*)'prob mit ladg in relcol 1a', qtest
                write(*,*)id11, id21, iblock, em1
                write(*,*)id12, id22, iblock, em2
                stop
              end if
            else if (id11.eq.nres+2) then
              if(id21.ne.0) then
                write(*,*)'hiba mit ladg in relcol 1a lambda'
                write(*,*)id11, id21, iblock, em1
                write(*,*)id12, id22, iblock, em2
                stop
              end if
            else if (id11.eq.nres+3) then
              if(iabs(id21).gt.1) then
                write(*,*)'hiba mit ladg in relcol 1a sigma'
                write(*,*)id11, id21, iblock, em1
                write(*,*)id12, id22, iblock, em2
                stop
              end if
            end if

            if(id12.gt.1 .and. id12.lt.nres+2) then
              if(resprop2(id12-1,1).eq.1) then
                if(id22.lt.0 .or. id22.gt.1) then
                  write(*,*)'prob mit ladg in relcol 3',qtest
                  write(*,*)id11, id21, iblock, em1
                  write(*,*)id12, id22, iblock, em2
                  stop
                end if
              else if(resprop2(id12-1,1).eq.3) then
                if(id22.lt.-1 .or. id22.gt.2) then
                  write(*,*)'prob mit ladg in relcol 4',qtest
                  write(*,*)id11, id21, iblock, em1
                  write(*,*)id12, id22, iblock, em2
                  stop
                end if
              end if
            else if(id12.eq.1) then
              if(id22.lt.0 .or. id22.gt.1) then
                write(*,*)'prob mit ladg in relcol 3a',qtest
                write(*,*)id11, id21, iblock, em1
                write(*,*)id12, id22, iblock, em2
                stop
              end if
            else if (id12.eq.nres+2) then
              if(id22.ne.0) then
                write(*,*)'hiba mit ladg in relcol 3a lambda'
                write(*,*)id11, id21, iblock, em1
                write(*,*)id12, id22, iblock, em2
                stop
              end if
            else if (id12.eq.nres+3) then
              if(iabs(id22).gt.1) then
                write(*,*)'hiba mit ladg in relcol 3a sigma'
                write(*,*)id11, id21, iblock, em1
                write(*,*)id12, id22, iblock, em2
                stop
              end if
            end if

*--------------------------------------------------------------------*

***********************************************************************
*
*     nothing has happend
      if(iblock.eq.0)                                        goto 400
*
**********************************************************************
*
*   a collision has taken place !!
              lcoll(-1)=lcoll(-1)+1
*

*-----------------------------------------------------------------------
*        lorentz-transformation into lab frame (including the
*        possibility for s-state pion-production
                e1cm  = sqrt((em1+ scapo3)**2 + pnucl1(1)**2 +
     &                        pnucl1(2)**2 + pnucl1(3)**2)

                call lorentz(-beta(1),-beta(2),-beta(3),
     &               pnucl1(1),pnucl1(2),pnucl1(3),e1cm)

                p(1,i1) = pnucl1(1)
                p(2,i1) = pnucl1(2)
                p(3,i1) = pnucl1(3)


                e2cm  = sqrt((em2+ scapo4)**2 + pnucl2(1)**2 +
     &                        pnucl2(2)**2 + pnucl2(3)**2)

                call lorentz(-beta(1),-beta(2),-beta(3),
     &               pnucl2(1),pnucl2(2),pnucl2(3),e2cm)

                p(1,i2) = pnucl2(1)
                p(2,i2) = pnucl2(2)
                p(3,i2) = pnucl2(3)

          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x1'
*******   check the boost
                if(.not.sstate) then
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x2'
                  help = p(1,i1)+p(1,i2)
                  if(abs(help-ptotx).gt.1.0e-05) then
                    write(*,*)'boost problems in relcol1 '
                    write(*,*)i, p(1,i1),p(1,i2),ptotx
                    write(*,*)iblock
                    stop
                  end if

                  help = p(2,i1)+p(2,i2)
                  if(abs(help-ptoty).gt.1.0e-05) then
                    write(*,*)'boost problems in relcol2   '
                    write(*,*)i, p(2,i1),p(2,i2),ptoty
                    write(*,*)iblock
                    stop
                  end if


                  help = p(3,i1)+p(3,i2)
                  if(abs(help-ptotz).gt.1.0e-05) then
                    write(*,*)'boost problems in relcol3 '
                    write(*,*)i, p(3,i1),p(3,i2),ptotz
                    write(*,*)iblock
                    stop
                  end if
            if(nt.eq.-700 .and. i2.ge.0)
     &            write(*,*)'relcol x3', sstate

                help1 = (e1cm+e2cm)**2 -(pnucl1(1)+pnucl2(1))**2 -
     &                  (pnucl1(2)+pnucl2(2))**2 -
     &                  (pnucl1(3)+pnucl2(3))**2

                help1 = sqrt(help1)
                if(abs(help2-help1).gt.1.0e-04) then
                   write(*,*)'relcol prob 2. mal ', iblock
                   write(*,*)srt,em1,em2,id11,id12, help, help1
                   write(*,*)irun,i1,i2,j1,j2
                   stop
                 end if

                end if

          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x4'

*       check for pauli-blocking
*       was collision pauli-forbidden? if yes, ntag = -1
              ntag = 0


              if((massta+masspr.gt.2).and.(icoll.ne.-1)) then

                if(ipauli.eq.1.and.id(1,i1).eq.1) then
                   call pauli(i1,ntag,iseed,phase,xxx,yyy,zzz,
     &                         p(1,i1),p(2,i1),p(3,i1))
                end if
*
                if((ntag .ne. -1).and.
     &            ipauli.eq.1.and.id(1,i1).eq.1) then
                       call pauli(i2,ntag,iseed,phase,xxx,yyy,zzz,
     &                            p(1,i2),p(2,i2),p(3,i2))
                end if

              end if
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x5'
*-----------------------------------------------------------------------
*       collision was pauli-blocked
           if (ntag .eq. -1) then
*         set variable back to values of initial particles
              lcoll(0) = lcoll(0)+1
              p(1,i1) = px1
              p(2,i1) = py1
              p(3,i1) = pz1
              p(1,i2) = px2
              p(2,i2) = py2
              p(3,i2) = pz2
              em1     = e(i1)
              id11    = id(1,i1)
              id21    = id(2,i1)
           else

c          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x6'
*----------------------------------------------------------------------*
*      if coll. was s-state pion prod.                                 *

             if(sstate) then
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x7'
               epicm  = sqrt((pmass+ upi)**2 + ppi3(1)**2 +
     &                       ppi3(2)**2 + ppi3(3)**2)

               call lorentz(-beta(1),-beta(2),-beta(3),
     &              ppi3(1),ppi3(2),ppi3(3),epicm)

*           get a new pion id

               inp = (irun-1) * maxp
               ib =inp
  10           ib= ib+ 1
               if(ipi(1,ib) .ne. 0) goto 10
               if(ib .gt. (inp+maxp) ) then
                 write(isum,'(10x,''***  too many pions  ***'')')
                 write(*,'(''relcol  *** too many pions  ***'')')
                 stop
               end if

                rpi(1,ib) = rpi3(1)
                rpi(2,ib) = rpi3(2)
                rpi(3,ib) = rpi3(3)

                rpie(1,ib)= rpi3(1)
                rpie(2,ib)= rpi3(2)
                rpie(3,ib)= rpi3(3)
                rpie(4,ib)= j01 !der eingfachheit halber
c-hw               rpie(5,ib)= e(i1)
                rpie(6,ib)= time

                epi(ib)   = pmass
                mpot(ib)  = upi
                ppi(1,ib) = ppi3(1)
                ppi(2,ib) = ppi3(2)
                ppi(3,ib) = ppi3(3)
                ipi(1,ib) = 1
                ipi(2,ib) = izpi
                ipi(3,ib) = i1
                ipi(4,ib) = 1
                ipi(5,ib) = 1
                ipi(6,ib) = id(6,i1)
                ipi(7,ib) = id(5,i1)
                ipi(8,ib) = i2

************** store coulomb pot
                crx = rpi(1,ib)
                cry = rpi(2,ib)
                crz = rpi(3,ib)

                cpx = ppi(1,ib)
                cpy = ppi(2,ib)
                cpz = ppi(3,ib)

                cid2 = ipi(2,ib)
                etot = sqrt(pmass**2+cpx**2+cpy**2+cpz**2)
                call emfoca(crx,cry,crz,cid2,
     &                       cfox,cfoy,cfoz,ncont,cpot)

c-hw               rpie(7,ib) = cpot

*********************************************



                id(8,i1)  = ib
                id(8,i2)  = ib
*------------------ LRF-properties for the pion ----------------------*
                    x     = rpi(1,ib)
                    y     = rpi(2,ib)
                    z     = rpi(3,ib)

                    if(isplipi .eq. 1) then
                      call splinint1(x,y,z,deriv)
                    else if(isplipi.eq. 0) then
                      call linint1(x,y,z,deriv)
                    end if
                    j0me    = deriv(0)
                    j1me    = deriv(1)
                    j2me    = deriv(2)
                    j3me    = deriv(3)

                    if(j0me .gt. 1.0e-6) then
                      mbetlrfx = j1me/j0me
                      mbetlrfy = j2me/j0me
                      mbetlrfz = j3me/j0me
                    else
                      mbetlrfx = 0.0
                      mbetlrfy = 0.0
                      mbetlrfz = 0.0
                    end if
              if(j1me**2+j2me**2+j3me**2.gt.j0me**2) then
                 write(*,*) "hiba relcol4 lorentz, negative mass",
     &              j1me,j2me,j3me,j0me
c                 stop
              end if
           call lorentz(mbetlrfx,mbetlrfy,mbetlrfz,j1me,j2me,j3me,j0me)
                    betlrfbom(ib,1) = mbetlrfx
                    betlrfbom(ib,2) = mbetlrfy
                    betlrfbom(ib,3) = mbetlrfz
                    betlrfbom(ib,4) = j0me
*----------------------------------------------------------------------*

          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x8'
              else
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x9'
                id(8,i1) = 0
                id(8,i2) = 0
              end if
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x10'

*       collison was allowed by pauli-principle
                if(id(6,i1)*id(6,i2).lt. 0)       lcoll(-5)=lcoll(-5)+1
                if(iblock.eq.1) then
                  if(id(2,i1)+id(2,i2).eq. 1)     lcoll(-2)=lcoll(-2)+1
                  if((em1.gt..94.and.em2.lt..94).or.
     &               (em1.lt..94.and.em2.gt..94)) lcoll(-3)=lcoll(-3)+1
                  if(em1.gt..94.and.em2.gt..94)   lcoll(-4)=lcoll(-4)+1
                end if
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x11'

           if(iblock.gt.3*nres + nres**2 + 4)
     &           write(*,*)'hiba relcol ibloc>21',iblock
                if(iblock.ge.1) lcoll(iblock) = lcoll(iblock) +1
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x12'
*--------------------------------------------------------------------
*      evaluate different reaction channels
*          (still in the else branch of ntag-if-statement)
*--------------------------------------------------------------------
*                       n+n -> n+r (resonance production)
                if((iblock.gt.2*nres+1).and.(iblock.le.(1+3*nres))) then
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x13'
                  if(id11 .ne. id(1,i1)) then
                     kk=i1
                     idres = id11
                  end if
                  if(id12 .ne. id(1,i2)) then
                     kk=i2
                     idres = id12
                  end if
                  idres = idres-1
                  if(idres.lt.1) then
                   write(6,'(''warning relcol idres lt 1'',i8)')idres
                  end if
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'x13a'
                  call resprod(kk,0.5*(t1+t2),id61,id62,idres,1)
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x14'
                end if
*----------------------------------------------------------------------
*                      n+r -> n+n  (resonance absorption)
                if((iblock.gt.(nres+1)).and.(iblock.le.(1+2*nres))) then
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x15'
                  if(id(1,i1) .ne. 1) then
                    kk = i1
                    idres = id11 - 1
                  else if(id(1,i2) .ne. 1) then
                    kk = i2
                    idres = id12 - 1
                  end if

                    if(kk.eq.i1) then
                      px = px1
                      py = py1
                      pz = pz1
                    else
                      px = px2
                      py = py2
                      pz = pz2
                    end if

                  call resabs(kk,0.5*(t1+t2),px,py,pz,1)
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x16'
                end if

*--------------------------------------------------------------------------
*                    n+r -> n+r'
                if((iblock.gt.(1+3*nres)).and.
     &             (iblock.le.(1+3*nres+nres**2)))then
          if(nt.eq.-700 .and. i2.ge.0) then
              write(*,*)'relcol x17',id(1,i1),id(1,i2),i1,
     &              i2,iblock
          end if
                  if(id(1,i1) .ne. 1) then
                    kk = i1
c                    idres = id11 - 1
                    rescheck1 = id(1,i1)
                  else if(id(1,i2) .ne. 1) then
                    kk = i2
c                    idres = id12 - 1
                    rescheck1 = id(1,i2)
                  end if
c                  write(*,*)'kk = ', kk
                  if(kk.eq.i1) then
                    px = px1
                    py = py1
                    pz = pz1
                  else
                    px = px2
                    py = py2
                    pz = pz2
                  end if
*
                  call resabs(kk,0.5*(t1+t2),px,py,pz,0)

c                  if(idres.lt.1) then
c                    write(*,'(''warning  relcol2 idres lt 1'',5i8)')
c     &                    idres,id11,id(1,i1),id12,id(1,i2),i1,i2
c                  else
                     if(id11 .ne. 1) then
                        kk = i1
                        idres = id11 - 1
                        rescheck2 = id11
                     else if(id12 .ne. 1) then
                        kk = i2
                        idres = id12 - 1
                        rescheck2 = id12
                     end if
                     if(idres.lt.1) then
                       write(*,'(''warning  relcol2 idres lt 1'',5i8)')
     &                  idres,id11,id(1,i1),id12,id(1,i2),i1,i2
                     end if

                    call resprod(kk,0.5*(t1+t2),id61,id62,idres,0)

                  if(rescheck1.eq.rescheck2) then
                    write(*,*)'error in NR NRt '
                    write(*,*)id11, id12
                    write(*,*)id(1,i1),id(1,i2)
                    write(*,*)i1, i2
                    write(*,*) iblock
                  end if

          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x18'
                end if
*-----------------------------------------------------------------------
*                store number of colliding particles
                id(3,i1)= i2
                id(3,i2)= i1
c          write(*,*) 'relcol 1 stor elott',px1,py1,pz1,em1,e10 
c        if (i1 .eq. 156)
c    1   write(*,*) '  relcol - partners ',irun, i1,(id(i,i1),i=1,8),
c    1                  i2,(id(i,i2),i=1,8)
                px1     = p(1,i1)
                py1     = p(2,i1)
                pz1     = p(3,i1)
c          write(*,*) 'relcol 1 store utan',px1,py1,pz1,em1,e10 
*                 store new masses and ids
                if(iblock.ge.1) then
                  e(i1)   = em1
                  e(i2)   = em2
                  id(1,i1)= id11
                  id(1,i2)= id12
                  id(2,i1)= id21
                  id(2,i2)= id22
                  upot(i1)= scapo3
c                   scapo1 = scapo3
                  upot(i2)= scapo4
c                   scapo2 = scapo4
                  meff1   = em1 + scapo3
c                  if (id(1,i2).eq.2) then
c                    write(*,*) 'delta from NN - mass:',e(i2)
c                  end if
c          write(*,*) 'relcol upot',scapo3,scapo4,id11,id12,id21
                end if
                id(6,i1)= id61+ id61/abs(id61)
                id(6,i2)= id62+ id62/abs(id62)
*                end of storing
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x19'
*---------------------------------------------------------------------
*                set values of particle 1 (not changed in inner loop)
                e1      = sqrt( meff1**2 + px1**2 + py1**2 + pz1**2 )
                id61    = id(6,i1)
                id21    = id(2,i1)
                id11    = id(1,i1)
*---------------------------------------------------------------------
*               end of ntag if-statement
           end if
          if(nt.eq.-700 .and. i2.ge.0)  write(*,*)'relcol x20'
*-------------------------------------------------------------------
c            if(id(1,i1).gt.4 .or. id(1,i1).lt.1 .or.
c     &         (id(1,i1).ne.2.and.(id(2,i1)+2)/2.ne.1) ) then
c               write(6,'(''hiba relcol3'',3i8,2f8.3)')
c     &              id(1,i1),id(2,i1), iblock, e(i1),time
c            end if
c            if(id(1,i2).gt.4 .or. id(1,i2).lt.1 .or.
c     &        (id(1,i2).ne.2.and.(id(2,i2)+2)/2.ne.1) )  then
c              write(6,'(''hiba relcol3'',3i8,2f8.3)')
c     &              id(1,i2),id(2,i2),iblock, e(i2),time
c            end if
          if(nt.eq.-700 .and. i2.ge.0) write(*,*)'relcol x21'
*-----------------------------------------------------------------------
*              end of do-loops
*              nothing happend label
c     if (i1 .eq. 6280)
c    1          write(*,*)' end loop 2 relcol',irun, i1,i2,e1,upot(i1)
  400       continue
*              loop over particle 2
c          write(*,*) 'relcol 1 600 elott',px1,py1,pz1,em1,e10 
  600     continue
*              loop over particle 1
  800   continue
*              loop over parallel runs
        if(ippbar.gt.0 .and. ippbarcoll.eq.1) then
          do i1 = (irun - 1) * maxb + 1,irun*maxb
            if(id(1,i1).eq.-1) id(1,i1)=0
          end do
        end if
 1000 continue
***********************************************************************

      write(*,*)'******************************+'
      write(*,*)'dichte bei 0,0,0   ', rhob_4(0,0,0,0)
      write(*,*)'tot number of colls = ', lcoll(-1)
      write(*,*)'******************************+'

      do i = 1, nres+3
        nuba(i) = 0
      end do

      do i = 1, maxpar
        if(id(1,i).ne.0) then
          nuba(id(1,i)) = nuba(id(1,i)) + 1
        end if
      end do

      write(*,*)'******************************'
      write(*,*)'number of res after coll'
      do i = 1, nres+3
        write(*,*)nameres(i), nuba(i)
      end do
      write(*,*) 'pbar    ', nuba(-1)
      write(*,*)'******************************'
c************************************AlmasiG_begin************************************
c      open(70,FILE='barmesout.dat',access='append')
c      write(70,*)'******************************'
c      write(70,*)'number of res after coll'
c      do i = 1, nres+3
c        write(70,*)nameres(i), nuba(i)
c      end do
c      close(70)
c************************************AlmasiG_end**************************************

      write(*,*)
c       tout = secnds(t0)
c     tin = tout - tin
c     write(*,*)'time ellapsed in relcol = ',tin,'  sec.'

*
*              end of collision term
      return
      end
*
