************************************************************************
*                                                                      *
      subroutine mesdec(lropi,lsipi,lompi,flag,lmesa2)
*                                                                      *
************************************************************************
      implicit none
      common /nthwhw/  nthw
      integer nthw
c
      real*8 ede, ema2, pm2, gamma, dt0, w, gam, sgam, bwmes, rnxx
      real*8 xx,     epion
      real*8   dendel, path, rn
      logical gridte,flag
      integer lropi, lsipi, lompi,irun, inp, ii, jj, id1, id2, ib
      integer iendel, ixx, iyy, izz, npat, idj
      real*8 x1,y1,z1, x2,y2,z2,dist,rmesdec2,gtot,g(3),rho2,
     &     enuk1,enuk2,srts,pcmf,  beta(3),traf, p1beta,transf,
     &      fac,rho,phase, vrel, vrel0
      real*8 matrix(3,1000)
      integer j1,i3,i4,i1,i2,qpi,istore(2,0:1),inucl(0:1),ka,imat,qtot,
     &     index, ntag,i,lmesa2, iii, iix, iiy, iiz

      real*8    px1, px2, pxp, py1, py2, pyp, pz1, pz2, pzp, eps
      real*8    e1, e2, ep, meff1, meff2, meffp, trmass,ymass
      real*8    scapo1, scapo2, scapop, em1, em2, emp
      real*8    betlrfx(2), betlrfy(2), betlrfz(2), j01, j02
      real*8    u3, u4, pfincmx, pfincmy, pfincmz, srtshelp
      real*8    meff3, meff4, help1, help2, e3cm, e4cm
      real*8    etotal, ptotx, ptoty, ptotz , pbx, pby, pbz
      real*8    ptest(3), etotb, rhon, rhop, wolfc
      real*8    crx, cry, crz, cpx, cpy, cpz, cpot
      real*8    cfox, cfoy, cfoz, etot, rhap1(1:3), rhap2(1:3)
      real*8 density, rself, iself, jrho1, jrho2, jrho3, zmass
      real*8    ee0,eex, vvx(0:3), vvy(0:3), betx(0:3)
      integer cid2
      logical negp
      save  matrix
      include"common"
      include"cominput"
c      include"com_cont_epair"
*----------------------------------------------------------------------*
      lropi= 0
      lompi =0
      lsipi= 0
      lmesa2 = 0
      pm2  = pmass**2
      wolfc = 2.57672

*   loop over all parallel runs*
      do 1000 irun = 1,num
        inp = (irun-1) * maxp
        do 800 ii  = 1,maxp-1

          jj = inp + ii
          if((ipi(1,jj) .eq. 2).or.(ipi(1,jj).lt.1))          goto 800
          if ((ipi(1,jj) .eq. 3).and.
     1                          (epi(jj) .lt. 2.*pmass+.001)) goto 800
          if ((ipi(1,jj) .eq. 5).and.
     1                          (epi(jj) .lt. 3.*pmass+.001)) goto 800

          id1 = ipi(1,jj)
          id2 = ipi(2,jj)

          ema2     = epi(jj)**2
          ede      = sqrt(ema2+ppi(1,jj)**2+ppi(2,jj)**2+ppi(3,jj)**2)
c       write(*,*) ' in mesdec ' , ii, irun, jj, id1, ema2, ede
c       call f77flush()
          gamma    = ede / epi(jj)
          dt0      = dt / gamma

          gam = 0.0


          if(id1.eq.1 .and. isstatea.eq.1) then

*********             two-body absorption
***         2 Nukleonen suchen,deren Abstand (nichtrel.,da bedeutungslos)  zum
***         Pion kleiner als rmesdec(cominput,fuer rho0) ist

            x1=rpi(1,jj)
            y1=rpi(2,jj)
            z1=rpi(3,jj)

            iix = nint(x1)
            iiy = nint(y1)
            iiz = nint(z1)
            if(abs(iix).lt.maxx .and. abs(iiy).lt.maxx .and.
     &        abs(iiz).lt.maxz) then
              rho=rhob_4(0,iix,iiy,iiz)
              jrho1 = rhob_4(1,iix,iiy,iiz)
              jrho2 = rhob_4(2,iix,iiy,iiz)
              jrho3 = rhob_4(3,iix,iiy,iiz)
              rhop = rhob_4(4,iix,iiy,iiz)
              rhon = rhob_4(5,iix,iiy,iiz)
            else
              rho = 0.0
              jrho1 = 0.0
              jrho2 = 0.0
              jrho3 = 0.0
              rhop =  0.0
              rhon =  0.0
            end if
               if(rho.lt.1e-03) goto 800
            rmesdec2=rmesdec*(rho0/rho)**(1./3.)
            if(rmesdec2.gt.5.0) then
              rmesdec2 = 5.0
            end if

            inucl(0)=0    ! count protons found
            inucl(1)=0    ! cound neutrons found
            do i3=0,1
              do i4=1,2
                istore(i4,i3)=0 !array i.o. to store id of found nuleons
              end do
            end do
            j1 = 0

 17         j1 = j1+1
            if(j1.gt.maxb) goto 19

            i1  = j1 + (irun - 1) * maxb  !id of nucleon
*           allow only N-pi collisions
            if(id(1,i1).ne.1)                                goto 17

*           coll. forbidden, if the parent res. of the pion
*           is the coll. partner of the pion AND if the
*           last meson contact of the nucleon was the coll. pion
            if((ipi(3,jj).eq.i1) .and. (id(8,i1).eq.jj))     goto 17

*           coll forbidden, if nucleon is the parent res of coll pion
*           AND if the nucl is associated with the 2pi partner
c            if((ipi(3,jj).eq.i1) .and. (id(8,i1).eq.ipi(8,jj))
c     &         .and. (id(8,i1).ne.0))                        goto 17

*           coll forbidden, if nucleon was the last coll partner of
*           the parent of the pion and vice versa
c            if((ipi(7,jj).eq.i1) .and. (id(3,i1).eq.ipi(3,jj))
c     &         .and. (id(3,i1).ne.0))                        goto 17


            x2=r(1,i1)
            y2=r(2,i1)
            z2=r(3,i1)

            qpi=ipi(2,jj)
            dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

            if(dist.le.rmesdec2) then

*              write(217,*)id(2,i1),inucl(1),inucl(0)
              if(((id(2,i1).eq.0).and.(inucl(0).le.1)).or.
     &           ((id(2,i1).eq.1).and.(inucl(1).le.1))) then
                 inucl(id(2,i1))=inucl(id(2,i1))+1
                 istore(inucl(id(2,i1)),id(2,i1))=i1
              end if
              if((inucl(1).eq.2).and.(inucl(0).eq.2)) goto 19
            end if
            goto 17
 19         continue

            if(inucl(1)+inucl(0).lt.2) then
*              2 nucleons could not be found
c              write(217,*)'inucl <2',inucl(1),inucl(0)
              goto 800
            else if((inucl(1)+inucl(0)).eq.2) then
*              exactly 2 nucleons have been found
c              write(217,*)'inucl =2',inucl(1),inucl(0)
            else if((inucl(1)+inucl(0)).eq.3) then
*              three nucleons have been found
c              write(217,*)'inucl =3',inucl(1),inucl(0)
            else if((inucl(1)+inucl(0)).eq.4) then
c              write(217,*)'inucl =4',inucl(1),inucl(0)
            end if
*           at least two nucleons have been found
            epion=sqrt(epi(jj)**2+ppi(1,jj)**2+ppi(2,jj)**2+ppi(3,jj)
     &           **2)
            gtot=0.
            do ka=1,3  ! loop over all three possible channels
              imat=0
              g(ka)=0.
              if(ka.eq.1) then
*pp
                i1=istore(1,1)
                i2=istore(2,1)
*                rho2=ratp**2*rho**2
                rho2=rhop**2
                if(qpi.eq.0) imat=1
                if(qpi.eq.-1) imat=2
              else if(ka.eq.2) then
*pn
                i1=istore(1,1)
                i2=istore(1,0)
c                rho2=ratn*ratp*rho**2
                rho2=rhop*rhon
                if(abs(qpi).eq.1) imat=1
                if(qpi.eq.0) imat=3
              else if(ka.eq.3) then
*nn
                i1=istore(1,0)
                i2=istore(2,0)
*                rho2=ratn**2*rho**2
                rho2=rhon**2
                if(qpi.eq.0) imat=1
                if(qpi.eq.1) imat=2
              end if
              if((i1.ne.0).and.(i2.ne.0)) then

                qtot=id(2,i1)+id(2,i2)+qpi
                if((qtot.lt.3).and.(qtot.gt.-1)) then
                  enuk1=sqrt(e(i1)**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
                  enuk2=sqrt(e(i2)**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
                  srts=sqrt((epion+enuk1+enuk2)**2-
     &                 (ppi(1,jj)+p(1,i1)+p(1,i2))**2-
     &                 (ppi(2,jj)+p(2,i1)+p(2,i2))**2-
     &                 (ppi(3,jj)+p(3,i1)+p(3,i2))**2)
                  if(srts.lt.(2.*rmass+pmass)) then
                    write(*,*)'problems in mesdec srts',
     &                         srts-2.*rmass-pmass
c                    stop
                  end if
                  index=min(max(nint((srts-(2.*rmass+pmass))*1000.),1),
     &                 1000)
                  fac=1.
                  if(imat.eq.2) fac=2.
                  if((imat.le.0).or.(imat.ge.4)) then
                    write(*,*)'problems in mesdec imat',imat
                    stop
                  end if
                  pcmf=sqrt(srts**2/4.-rmass**2)
                  g(ka)=fac*matrix(imat,index)*pcmf/srts/4./pi*rho2/8./
     &                 enuk1/enuk2/epion* wolfc * hbc**6
c            write(*,*) ' g-1-1 ', imat, index,
c    1        fac,matrix(imat,index),
c    1        pcmf,srts,pi,rho2,enuk1,enuk2,epion
*                  write(217,*)srts,ka,g(ka)
                  gtot=g(ka)+gtot
                end if
              end if
            end do
c            write(*,*) ' g-1-2-3 ',g(1),g(2),g(3),gtot
*      Monte-Carlo decision if pion is going to be absorbed

            w=exp(-dt*gtot/hbc)
*dt anstatt dt0, da Zerfallsbreite im Laborsystem berechnet wurde
            rnxx=rn(iseed)
            if(rnxx.gt.w) then
*Pion wird eventuell absorbiert
*wuerfeln des Kanals
              rnxx=rn(iseed)
              if(rnxx.le.(g(1)/gtot)) then
                ka=1
              else if(rnxx.le.(g(1)+g(2))/gtot) then
                ka=2
              else
                ka=3
              end if
              if(ka.eq.1) then
*pp
                i1=istore(1,1)
                i2=istore(2,1)
              else if(ka.eq.2) then
*pn
                i1=istore(1,1)
                i2=istore(1,0)
              else if(ka.eq.3) then
*nn
                i1=istore(1,0)
                i2=istore(2,0)
              end if
              if((i1.eq.0).or.(i2.eq.0)) then
                write(*,*)'problems in mesdec nach Kanal wuerfeln',
     &               g(1),g(2),g(3),gtot,ka,rnxx
              end if
*---------------------------------------------------------------------*
*             determine the kinematics of outgoing nucleons           *
*             same as relcol.f
*     nucleon 1
          px1        = p(1,i1)
          py1        = p(2,i1)
          pz1        = p(3,i1)
          rhap1(1)   = r(1,i1)
          rhap1(2)   = r(2,i1)
          rhap1(3)   = r(3,i1)
          em1        = e(i1)
          scapo1     = upot(i1)
          meff1      = em1 + scapo1
          e1         = sqrt( meff1**2 + px1**2 + py1**2 + pz1**2 )
          betlrfx(1) = betlrfboo(i1,1)
          betlrfy(1) = betlrfboo(i1,2)
          betlrfz(1) = betlrfboo(i1,3)
          j01        = betlrfboo(i1,4)
*     nucleon 2
          px2        = p(1,i2)
          py2        = p(2,i2)
          pz2        = p(3,i2)
          rhap2(1)   = r(1,i2)
          rhap2(2)   = r(2,i2)
          rhap2(3)   = r(3,i2)
          em2        = e(i2)
          scapo2     = upot(i2)
          meff2      = em2 + scapo2
          e2         = sqrt( meff2**2 + px2**2 + py2**2 + pz2**2 )
          betlrfx(2) = betlrfboo(i2,1)
          betlrfy(2) = betlrfboo(i2,2)
          betlrfz(2) = betlrfboo(i2,3)
          j02        = betlrfboo(i2,4)
*     pion
          pxp        = ppi(1,jj)
          pyp        = ppi(2,jj)
          pzp        = ppi(3,jj)
          emp        = epi(jj)
          scapop     = mpot(jj)
          meffp      = epi(jj) + mpot(jj)
          ep         = sqrt( meffp**2 + pxp**2 + pyp**2 + pzp**2 )

          srts = (e1 + e2 + ep)**2 - (px1 + px2+ pxp)**2 -
     &           (py1 + py2 + pyp)**2 - (pz1 + pz2 + pzp)**2
          srts = sqrt(srts)
          srtshelp = srts

*    lorentz-transformation parameters for boost in system with
*     \vec{p1+p2+pp} = 0
              etotal = e1 + e2  + ep
              ptotx  = px1+ px2 + pxp
              ptoty  = py1+ py2 + pyp
              ptotz  = pz1+ pz2 + pzp

              beta(1) = ptotx/etotal
              beta(2) = ptoty/etotal
              beta(3) = ptotz/etotal
c          write(*,*) ' in mesdec 12 ',etotal, e1,e2,ep, i1,i2,jj,
c    1                  px1,px2,pxp

              gamma  = 1.0 / sqrt(1.0-beta(1)**2-beta(2)**2-beta(3)**2)
              traf   = gamma / (gamma + 1.)

***   now pretend like the calc is done in the p= 0 frame and call momiter
***   with isotrop parameters
              u3      = 0.0
              u4      = 0.0
              pfincmx = 0.0
              pfincmy = 0.0
              pfincmz = 0.0
              negp    = .false.
c       write(*,*) ' mesdec1  :  vor momiter ',irun, ii, jj, id1
c      write(*,*)  srts,1,rmass, rmass, 1, 1,
c    &                      rho0, j01, j02,
c    &                      betlrfx, betlrfy, betlrfz, beta,
c    &                      u3, u4, pfincmx, pfincmy, pfincmz, negp,
c    &                      1, 0, 0, 0.0, 0.0, 0.0, 0,iseed,0.0,
c    &                      rhap1,rhap2

              call momiter(srts,1,rmass, rmass, 1, 1,
     &                      rho0, j01, j02,
     &                      betlrfx, betlrfy, betlrfz, beta,
     &                      u3, u4, pfincmx, pfincmy, pfincmz, negp,
     &                      1, 0, 0, 0.0, 0.0, 0.0, 0,iseed,0.0,
     &                      rhap1,rhap2)

              if(abs(srtshelp-srts).gt.1.0e-04) then
                write(*,*)'in mesdec :prob with momiter :'
                write(*,*)srts, srtshelp
                stop
              end if

              if(negp) then
                write(*,*)'in mesdec mom-iteration is wrong 1 '
                write(*,*) srts,1,rmass, rmass, 1, 1,
     &                      rho0, j01, j02,
     &                      betlrfx, betlrfy, betlrfz, beta,
     &                      u3, u4, pfincmx, pfincmy, pfincmz, negp,
     &                      1, 0, 0, 0.0, 0.0, 0.0, 0,iseed,0.0,
     &                      rhap1,rhap2
               goto 800
                stop
              end if
              meff3 = em1 + u3
              meff4 = em2 + u4

*            pfinx, pfiny, pfinz contains the momenta in the p=0 frame
*            boost back

            help1 = sqrt(meff3**2+pfincmx**2+pfincmy**2+pfincmz**2)
            help2 = sqrt(meff4**2+pfincmx**2+pfincmy**2+pfincmz**2)
            help1 = help1 + help2
            if(abs(srtshelp-help1) .gt. 1.0e-04) then
                write(*,*)'in mesdec :prob nach momiter :'
                write(*,*)help1, srtshelp
                stop
            end if


*----------------------------------------------------------------------
*             lorentz-transformation into lab frame
            e3cm  = sqrt(meff3**2+pfincmx**2+pfincmy**2+pfincmz**2)
     &
            p1beta  = pfincmx*beta(1)+pfincmy*beta(2)+pfincmz*beta(3)
              transf  = gamma * ( traf * p1beta + e3cm )
              p(1,i1) = beta(1) * transf + pfincmx
              p(2,i1) = beta(2) * transf + pfincmy
              p(3,i1) = beta(3) * transf + pfincmz
*
              e4cm  = sqrt (meff4**2+pfincmx**2+pfincmy**2+pfincmz**2)
              transf  = gamma * (-traf * p1beta + e4cm)
              ptest(1)= beta(1) * transf - pfincmx
              ptest(2)= beta(2) * transf - pfincmy
              ptest(3)= beta(3) * transf - pfincmz
              p(1,i2) = ptotx - p(1,i1)
              p(2,i2) = ptoty - p(2,i1)
              p(3,i2) = ptotz - p(3,i1)
              upot(i1) = u3
              upot(i2) = u4
c              write(*,*) 'mesdec',u3,u4
*             test boost
              do i =1 ,3
                if(abs(ptest(i)-p(i,i2)).gt.1.0e-06) then
                  write(*,*)'prob with boost in mesdec '
                  write(*,*)i, ptest(i), p(i,i2)
c                  stop
                end if
              end do

              help1 = sqrt(meff3**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
              help2 = sqrt(meff4**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
              help1 = (help1+help2)**2 - (p(1,i1)+p(1,i2))**2 -
     &                (p(2,i1)+p(2,i2))**2 - (p(3,i1)+p(3,i2))**2
              help1 = sqrt(help1)
              if(abs(srtshelp-help1).gt.1.0e-06) then
                write(*,*)'mesdec prob with boost '
                write(*,*) help1, srtshelp
c                stop
              end if


*-------------------------------------------------------------------------
*       check for pauli-blocking
*       was collision pauli-forbidden? if yes, ntag = -1
              ntag = 0
              if(ipauli.eq.1)
     &         call pauli(i1,ntag,iseed,phase,r(1,i1),r(2,i1),r(3,i1),
     &                         p(1,i1),p(2,i1),p(3,i1))

*
              if(ntag .ne. -1 .and. ipauli.eq.1)then
                call pauli(i2,ntag,iseed,phase,r(1,i2),r(2,i2),r(3,i2),
     &                            p(1,i2),p(2,i2),p(3,i2))
              end if

*-----------------------------------------------------------------------
*       sstate decay  was pauli-blocked
              if (ntag .eq. -1) then
*         set variable back to values of initial particles
                p(1,i1) = px1
                p(2,i1) = py1
                p(3,i1) = pz1
                p(1,i2) = px2
                p(2,i2) = py2
                p(3,i2) = pz2
                upot(i1)= scapo1
                upot(i2)= scapo2
                goto 800
              end if

                qtot=id(2,i1)+id(2,i2)+qpi
*Verteilen der Ladungen
                if(qtot.eq.1) then
                  rnxx=rn(iseed)
                  if(rnxx.lt.0.5) then
                    id(2,i1)=1
                    id(2,i2)=0
                  else
                    id(2,i1)=0
                    id(2,i2)=1
                  end if
                else
c          write(*,*) i1, i2, id(2,i1), id(2,i2), inucl(0),inucl(1)
                  id(2,i1)=nint(qtot/2.)
                  id(2,i2)=nint(qtot/2.)
                end if
                ipi(1,jj)=0
                id(3,i1)=i2
                id(3,i2)=i1
                id(8,i1)= 0
                id(8,i2)= 0
                id(7,i1)= 1
                id(7,i2)= 1
                id(6,i1)=id(6,i1)+id(6,i1)/abs(id(6,i1))
                id(6,i2)=id(6,i2)+id(6,i2)/abs(id(6,i2))
                lmesa2=lmesa2+1
                if(qtot.ne.(id(2,i1)+id(2,i2))) then
                  write(*,*)'charge prob in mesdec'
                  write(*,*)qtot, id(2,i1), id(2,i2)
                  stop
                end if

            end if
            goto 800


***************** end of s state

          else if(id1.ge.3 .and. id1 .le.5) then

            x1=rpi(1,jj)
            y1=rpi(2,jj)
            z1=rpi(3,jj)

            iix = nint(x1)
            iiy = nint(y1)
            iiz = nint(z1)
            if(abs(iix).lt.maxx .and. abs(iiy).lt.maxx .and.
     &        abs(iiz).lt.maxz) then
              rho=rhob_4(0,iix,iiy,iiz)
              jrho1 = rhob_4(1,iix,iiy,iiz)
              jrho2 = rhob_4(2,iix,iiy,iiz)
              jrho3 = rhob_4(3,iix,iiy,iiz)
              rhop = rhob_4(4,iix,iiy,iiz)
              rhon = rhob_4(5,iix,iiy,iiz)
            else
              rho = .0
              jrho1 = 0.0
              jrho2 = 0.0
              jrho3 = 0.0
            endif
c
c            write(*,*) '  rhos ',id1,rho,jrho1,jrho2,jrho3
            if (abs(rho) .lt. 1.e-4)  then
              betx(0) = 1.0
              betx(1) = 0.0
              betx(2) = 0.0
              betx(3) = 0.0
              jrho1 = 0.0
              jrho2 = 0.0
              jrho3 = 0.0
            else
              betx(0) = rho
              betx(1) = -jrho1
              betx(2) = -jrho2
              betx(3) = -jrho3
            end if

c            if(id1.eq.3) gam=bwmes(ema2,1,idec2pi,iresmode,0,0,iwidth,0)
            if(id1.eq.3 .or. id1.eq.5) then
             ee0 = epi(jj)
             eex = sqrt(epi(jj)**2 + ppi(1,jj)**2 + ppi(2,jj)**2
     1            + ppi(3,jj)**2)
             vvx(0) = eex / ee0
             vvx(1) = ppi(1,jj) / ee0
             vvx(2) = ppi(2,jj) / ee0
             vvx(3) = ppi(3,jj) / ee0
             call lorentz_hw(betx, vvx, vvy)
c      write(*,*) ' in mesdec before bwmes ',jj,id1, epi(jj), ee0,vvx,
c    1              eex,ppi(1,jj),ppi(2,jj),ppi(3,jj), betx, vvy
             density = (rho**2 - jrho1**2 - jrho2**2 - jrho3**2)
             if (density .lt. .0) then
                     density = .0
             else
                     density = sqrt(density)
             endif
c             write(*,*)  ' in mesdec vvy ', density, vvy
             vrel = sqrt(vvy(1)**2+vvy(2)**2+vvy(3)**2) / vvy(0)
             if (id1.eq.3) then
               vrel0 = .0
               call self_rho(vrel0,density,ee0,rself,iself,sgam)
               ymass=romas+0.5*rself/romas
               if (icbro .eq. 2) trmass = .0
c----------------------------
               zmass = 1.5 * (em1 - ymass)/(em1-trmass)
               if (zmass .gt. .8)  zmass = .8
               if (icbro.eq.2 .and. em1.le.mrho_lim) zmass = .0
c--------------------------------------------------
               gam = sgam
             endif
             if(id1.eq.5)then
               vrel0 = .0
               call self_omega(vrel0,density,ee0,rself,iself,sgam)
               gam = sgam
               ymass=omass+0.5*rself/omass
               trmass=2.99*pmass
               zmass =  .0   !   1.5*(ee0-ymass)/(ee0-trmass)
c              if (zmass .gt. 0.8)  zmass = 0.8
             endif
             if (id1 .eq. 3  .or. id1.eq.5) then
c      at the moment omegas are not treated because  fixed sgam
ccc                  Danielewicz description
               if(ireslife.eq.1)
     1           gam = (sgam/2. + 2.0*(ee0-ymass)**2/sgam)
     2               /(1. - zmass)
c              if(ireslife.eq.2)
ccc                baz description
c    1          gam = (sgam/4. + (ee0-ymass)**2/sgam)
c    2               /(1. - zmass)
c         wolf prescription
c     1        gam = 0.5 * sgam + (sgam/2. + 2.0*(ee0-ymass)**2/sgam)
c     2               /(1. - zmass)
c             write(*,*) ' in mesdec after bwmes ',jj,id1
             end if
            end if
            if(id1.eq.4) gam=bwmes(ema2,2,idec2pi,iresmode,0,0,iwidth,0)
c           w = exp(-dt0 * gam / hbc / dlife)    !  life = input = 1.0
            eps = dt0 * gam / hbc
            if (abs(eps) .gt. 1.e-3) then
               w =  exp(-eps)
            else
               w = 1. - eps
            endif
            rnxx = rn(iseed)
C======================================================================
c---------------   no decay of  rho + omega mesons -------------
c            if (id1 .eq. 3 .or. id1 .eq. 5) then
c               w = 2.0
c            endif
c----------------------------------
c----------------------------------
c 732        format(i8,i3,i5,' mdv-G', 7f10.4)
c           if (id1 .eq. 5 ) then
c             write(53,732) jj, id1, nthw, epi(jj), density, vrel,
c             write(momepri,732) jj, id1, nthw, epi(jj), density, vrel,
c    1        rself,iself,gam,sgam     !   momepri=11
c           endif
c----------------------------------
c           if (id1.eq.3 .and. id2.eq.0) then
c             write(54,732) jj, id1, nthw, epi(jj), density, vrel,
c             write(mmespri,732) jj, id1, nthw, epi(jj), density, vrel,
c    1        rself,iself,gam,sgam    !  mmespri=12
c           endif
c----------------------------------
C=========================================================================
c
c=======================================================================
c               decay is treated   in mes_dilep
c-----------------------------------------------------------------------
c          if ((id1.eq.3 .and. id2.eq.0) .or. id1 .eq. 5 ) goto 800
c-----------------------------------------------------------------------
c          write(*,*)  '  d1,jj, gamma, w, rnxx  ',  id1,jj,gam,w,rnxx
            if ((id1 .eq. 3  .and.  ee0 .lt. 2.*pmass+.001)
     1         .or. (id1 .eq. 5  .and.  ee0 .lt. 3.*pmass+.001)) then
                 if (rnxx .gt. w) then
c                  meson decays
                   write(*,*) ' in mesdec - rho mass < 2mpi ',jj, id1
c                  ipi(1,jj) = 0
                 endif
                 goto 800
            endif

            if((rnxx .gt. w).or.flag) then
*                                                                      *
*            meson decays                                              *
*                                                                      *
             if(id1.eq.3  .or. id1.eq.5) then
c              write(48,730) jj, id1, nthw, epi(jj)
c  730          format(' meson decays  ', i8,i3,i5,f10.4)
               iii = (id1-1)/2
c               m_final(iii) = m_final(iii) + epi(jj)
c               t_final(iii) = t_final(iii) + dt*nthw
c               m2_final(iii) = m2_final(iii) + epi(jj)**2
c               n_final(iii) = n_final(iii) + 1
             endif
c             pa2= .25 * ema2 - pm2
c              pa = sqrt(pa2)
c***********************************************************************
              ib =inp
  10          ib= ib+ 1
              if(ipi(1,ib) .ne. 0) go to 10
              if(ib .gt. (inp+maxp) ) then
                write(isum,'(10x,''***  too many pions  ***'')')
                stop
                goto 1000
              else
***********************************************************************
**          determine the momenta of the outgoing pions          ******
***         always assumed that rho and sigma do not feel a potential *



*    lorentz-transformation parameters for boost in system with
*     \vec{p1} = 0 ; restframe of decaying particle
                etotal = epi(jj)
                ptotx  = ppi(1,jj)
                ptoty  = ppi(2,jj)
                ptotz  = ppi(3,jj)
                rhap1(1) = rpi(1,jj)
                rhap1(2) = rpi(2,jj)
                rhap1(3) = rpi(3,jj)
                rhap2(1) = rpi(1,jj)
                rhap2(2) = rpi(2,jj)
                rhap2(3) = rpi(3,jj)

                etotb  = sqrt(etotal**2+ptotx**2+ptoty**2+ptotz**2)

                beta(1) = ptotx/etotb
                beta(2) = ptoty/etotb
                beta(3) = ptotz/etotb

                gamma  =1.0/sqrt(1.0-beta(1)**2-beta(2)**2-beta(3)**2)
                traf   = gamma / (gamma + 1.)

***   now pretend like the calc is done in the p= 0 frame and call momiter
***   with isotrop parameters
                u3      = 0.0
                u4      = 0.0
                pfincmx = 0.0
                pfincmy = 0.0
                pfincmz = 0.0
                negp    = .false.
                srtshelp= etotal
                em1     = pmass
                em2     = pmass

                betlrfx(1) = betlrfbom(jj,1)
                betlrfy(1) = betlrfbom(jj,2)
                betlrfz(1) = betlrfbom(jj,3)
                betlrfx(2) = betlrfbom(jj,1)
                betlrfy(2) = betlrfbom(jj,2)
                betlrfz(2) = betlrfbom(jj,3)
                j01        = betlrfbom(jj,4)
                j02        = betlrfbom(jj,4)
c
c       write(*,*) ' mesdec2  :  vor momiter ',jj, etotal,id1
c
                call momiter(etotal,1,em1, em2, 1, 1,
     &                      rho0, j01, j02,
     &                      betlrfx, betlrfy, betlrfz, beta,
     &                      u3, u4, pfincmx, pfincmy, pfincmz, negp,
     &                      3, 0, 0, 0.0, 0.0, 0.0, 0,iseed,0.0,
     &                      rhap1,rhap2)

                if(abs(srtshelp-etotal).gt.1.0e-06) then
                  write(*,*)'in mesdec 2 :prob with momiter :'
                  write(*,*)etotal, srtshelp
c                  stop
                end if

                if(negp) then
                  write(*,*)'in mesdec mom-iteration is wrong 2'
                  write(*,*) etotal,1,em1, em2, 1, 1,
     &                      rho0, j01, j02,
     &                      betlrfx, betlrfy, betlrfz, beta,
     &                      u3, u4, pfincmx, pfincmy, pfincmz, negp,
     &                      3, 0, 0, 0.0, 0.0, 0.0, 0,iseed,0.0,
     &                      rhap1,rhap2
                  stop
                end if
                meff3 = em1 + u3
                meff4 = em2 + u4

*            pfinx, pfiny, pfinz contains the momenta in the p=0 frame
*            boost back

                help1 = sqrt(meff3**2+pfincmx**2+pfincmy**2+pfincmz**2)
                help2 = sqrt(meff4**2+pfincmx**2+pfincmy**2+pfincmz**2)
                help1 = help1 + help2
                if(abs(srtshelp-help1) .gt. 1.0e-05) then
                    write(*,*)'in mesdec 2 :prob nach momiter :'
                    write(*,*)help1, srtshelp
c                    stop
                end if

c                 write(*,*)'mesdec 2'
*----------------------------------------------------------------------
*             lorentz-transformation into lab frame
                e3cm  = sqrt(meff3**2+pfincmx**2+pfincmy**2+pfincmz**2)
                e4cm =  sqrt(meff4**2+pfincmx**2+pfincmy**2+pfincmz**2)

c                p1beta= pfincmx*beta(1)+pfincmy*beta(2)+pfincmz*beta(3)
c                transf  = gamma * ( traf * p1beta + e3cm )
c                ppi(1,ib) = beta(1) * transf + pfincmx
c                ppi(2,ib) = beta(2) * transf + pfincmy
c                ppi(3,ib) = beta(3) * transf + pfincmz
*
c                write(*,*)'mesdec 3', ib

c                transf  = gamma * (-traf * p1beta + e4cm)
c                ptest(1)= beta(1) * transf - pfincmx
c                ptest(2)= beta(2) * transf - pfincmy
c                ptest(3)= beta(3) * transf - pfincmz

                pbx = pfincmx
                pby = pfincmy
                pbz = pfincmz
              if(pbx**2+pby**2+pbz**2.gt.e3cm**2) then
                 write(*,*) "hiba mesdec1 lorentz, negative mass",
     &                e3cm,pbx,pby,pbz
c                 stop
              end if

            call lorentz(-beta(1),-beta(2),-beta(3),pbx,pby,pbz,e3cm)

                ppi(1,ib) = pbx
                ppi(2,ib) = pby
                ppi(3,ib) = pbz

                pbx = -pfincmx
                pby = -pfincmy
                pbz = -pfincmz

              if(pbx**2+pby**2+pbz**2.gt.e4cm**2) then
                 write(*,*) "hiba mesdec2 lorentz, negative mass",
     &                e4cm,pbx,pby,pbz
c                 stop
              end if
            call lorentz(-beta(1),-beta(2),-beta(3),pbx,pby,pbz,e4cm)

                ptest(1) = pbx
                ptest(2) = pby
                ptest(3) = pbz

                ppi(1,jj) = ptotx - ppi(1,ib)
                ppi(2,jj) = ptoty - ppi(2,ib)
                ppi(3,jj) = ptotz - ppi(3,ib)
                mpot(ib) = u3
                mpot(jj) = u4
                epi(jj)  = pmass
                epi(ib)  = pmass

*         check first boost
                pbx = ptotx
                pby = ptoty
                pbz = ptotz
                e4cm = sqrt(etotal**2+ptotx**2+ptoty**2+ptotz**2)
            call lorentz(beta(1),beta(2),beta(3),pbx,pby,pbz,e4cm)


*             test boost
                do i =1 ,3
                  if(abs(ptest(i)-ppi(i,jj)).gt.1.0e-06) then
                    write(*,*)'prob with boost in mesdec '
                    write(*,*)i, ptest(i), ppi(i,jj)
c                    stop
                  end if
                end do

c                write(*,*)'mesdec 5'
                help1 = sqrt(meff3**2+ppi(1,ib)**2+ppi(2,ib)**2
     &                    +ppi(3,ib)**2)
                help2 = sqrt(meff4**2+ppi(1,jj)**2+ppi(2,jj)**2
     &                    +ppi(3,jj)**2)
                help1 = (help1+help2)**2 - (ppi(1,ib)+ppi(1,jj))**2 -
     &             (ppi(2,ib)+ppi(2,jj))**2 - (ppi(3,ib)+ppi(3,jj))**2
                help1 = sqrt(help1)
                if(abs(srtshelp-help1).gt.1.0e-06) then
                  write(*,*)'pionco prob with boost '
                  write(*,*) help1, srtshelp
c                  stop
                end if


c                write(*,*)'in mesdec 2pi coll 1'

c                if(id1.eq.3) lropi = lropi + 1
c                if(id1.eq.4) lsipi = lsipi + 1
c  20            xx       = 1. - 2. * rn(iseed)
c                yy       = 1. - 2. * rn(iseed)
c                zz       = 1. - 2. * rn(iseed)
c                rr       = sqrt( xx**2 + yy**2 + zz**2 )
c                if((rr .lt. 0.001) .or. (rr .gt. 1.) ) go to 20
c                px       = pa * xx / rr
c                py       = pa * yy / rr
c                pz       = pa * zz / rr

***   the outgoing pion momentum : (px,py,pz) in the meson res. cms  ****
*-----------------------------------------------------------------------
*                 path length distribution
c        write(*,*) ' in mesdec - path ',jj, ipi(1,jj), epi(jj), ede,
c     1               rpie(1,jj),rpie(2,jj),
c     2               rpie(3,jj),rpi(1,jj),rpi(2,jj),rpi(3,jj)
                path=sqrt((rpie(1,jj)-rpi(1,jj))**2+
     &              (rpie(2,jj)-rpi(2,jj))**2+(rpie(3,jj)-rpi(3,jj))**2)
                npat=min(30,nint(2.0*path))
                mpath(id1,npat)=mpath(id1,npat) + 1
*-----------------------------------------------------------------------
*               lifetime distribution
                path=max( 0.0, time - rpie(6,jj))
                npat=min(50,nint(2.0*path))
                melet(id1,npat)=melet(id1,npat) + 1
*-----------------------------------------------------------------------
c           meson collision number

                idj = min(ipi(4,jj) , 50)
                mlife(id1,idj) = mlife(id1,idj) + 1
                mestopi  = mestopi + 1
************************************************************************
c****        store the position, density, mass, and time at creation *
                rpie(1,ib)= rpi(1,jj)
                rpie(2,ib)= rpi(2,jj)
                rpie(3,ib)= rpi(3,jj)
                rpie(1,jj)= rpi(1,jj)
                rpie(2,jj)= rpi(2,jj)
                rpie(3,jj)= rpi(3,jj)
                ixx = nint(rpi(1,jj))
                iyy = nint(rpi(2,jj))
                izz = nint(rpi(3,jj))
                if(abs(ixx).lt.maxx .and. abs(iyy).lt.maxx
     &                   .and.abs(izz).lt.maxz) then
                  gridte = .true.
                else
                  gridte = .false.
                endif
                if(gridte) then
                  rpie(4,ib)= rhb(ixx,iyy,izz)
                  rpie(4,jj)= rhb(ixx,iyy,izz)
                else
                  rpie(4,ib) = 0.0
                  rpie(4,jj) = 0.0
                end if
c-hw               rpie(5,ib)= epi(jj)
                rpie(6,ib)= time
c-hw                rpie(5,jj)= epi(jj)
                rpie(6,jj)= time

                rpi(1,ib)= rpi(1,jj)
                rpi(2,ib)= rpi(2,jj)
                rpi(3,ib)= rpi(3,jj)


c***********************************************************************
c                epion      = sqrt(pm2 + px**2 + py**2 + pz**2 )
c                bex      = - ppi(1,jj) / ede
c                bey      = - ppi(2,jj) / ede
c                bez      = - ppi(3,jj) / ede
c                bep      = bex * px + bey * py + bez * pz
c                gamma    = ede / epi(jj)
c                trafo    = gamma / (gamma + 1.0)
c                trans    = trafo * bep - epion

c                ppi(1,ib)  = px + gamma * bex * trans
c                ppi(2,ib)  = py + gamma * bey * trans
c                ppi(3,ib)  = pz + gamma * bez * trans
c                epi(ib)    = pmass

c                ppi(1,jj)  = ppi(1,jj) - ppi(1,ib)
c                ppi(2,jj)  = ppi(2,jj) - ppi(2,ib)
c                ppi(3,jj)  = ppi(3,jj) - ppi(3,ib)
c                epi(jj)    = pmass
*
c***********************************************************************
*           charge assignment
                  xx = rn(iseed)
                  ipi(2,ib) = 0
                  if(id1 .eq. 3 .or. id1.eq.5) then
                    if(ipi(2,jj) .eq. 0) then
                       if(xx.gt.0.5) then
                        ipi(2,ib) = -1
                        ipi(2,jj) = 1
                      else
                        ipi(2,ib) = 1
                        ipi(2,jj) = -1
                      endif
                    else
                      if(xx.gt.0.5) then
                        ipi(2,ib) = ipi(2,jj)
                        ipi(2,jj) = 0
                      endif
                    endif
                  else if(id1 .eq. 4) then
                    if(xx .lt. 0.33333) then
                      ipi(2,ib) = 0
                      ipi(2,jj) = 0
                    else if(xx.ge.0.33333 .and. xx.lt.66667) then
                      ipi(2,ib) = 1
                      ipi(2,jj) = -1
                    else
                      ipi(2,ib) = -1
                      ipi(2,jj) = 1
                    end if
                  end if
c              end if
*
                ipi(1,ib)= 1
                ipi(3,ib)= ipi(3,jj)
                ipi(4,ib)= 1
                ipi(5,ib)= -ipi(1,jj)
                ipi(6,ib)=jj

                ipi(1,jj)= 1
                ipi(4,jj)= 1
                ipi(5,jj)= ipi(5,ib)
                ipi(6,jj)=ib

******************* store the coulomb potential at creation time *****

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


                crx = rpi(1,jj)
                cry = rpi(2,jj)
                crz = rpi(3,jj)

                cpx = ppi(1,jj)
                cpy = ppi(2,jj)
                cpz = ppi(3,jj)

                cid2 = ipi(2,jj)
                etot = sqrt(pmass**2+cpx**2+cpy**2+cpz**2)
                call emfoca(crx,cry,crz,cid2,
     &                       cfox,cfoy,cfoz,ncont,cpot)
c-hw               rpie(7,jj) = cpot


*
************************************************************************
c              density of meson creation      *************************
                if(gridte) then
                  dendel = rhb(ixx,iyy,izz)/rho0
                else
                  dendel = 0.0
                end if
                iendel = nint(5.0*dendel)
                if(iendel .gt. 50) iendel = 50
                mdens(1,iendel) = mdens(1,iendel) + 2
************************************************************************
cc             end of enough particle loop                  ************
              end if
c              end of decay loop for res - pi pi
            end if
c              end of id test if statement
          end if
  800   continue
 1000 continue
        write(*,*) 'end mesdec '

      return

      entry mesdecini
*Einlesen der Matrixelemente von N N - N N pi (s-wave)
      open(msstread,file='buuinput/SSTATE.dat',status='old')
      do i=1,1000
        read(msstread,*) srts,matrix(1,i),matrix(2,i),matrix(3,i)
      end do
      close(msstread)
      return
      end
