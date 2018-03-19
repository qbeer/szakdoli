*----------------------------------------------------------------------
      subroutine mommass(srt,id1, id2, newmass1, newmass2,
     &                   massflag,srtinfree)
*-----------------------------------------------------------------------
*        evaluate mass distribution

      implicit none

      include"common"
      include"cominput"


      real*8    srt, newmass1, newmass2
      integer id1, id2, ireact
*     ireact = 0 : Delta, delta final channel
*     ireact= 1 : use Dimitriev to determin the masses
*           = other than 0,1: see indices of NNXSEC in file resdata
      logical massflag, loopflag
      real*8    bwdist,  ratio(1:9), normali, pfin, resmass
      integer ireso, indexs, indexslor, idsum, indexmass
      integer loopcount, indexf
      real*8    rmint, tmass, help, help1, x1,umass
      real*8    d1mass, d2mass, srt2d0, srt2de, x1d, x2d, rn
      integer nsave, nrandom
      real*8    sigtot, sumdsdm, dsdm,srtinfree

c      include"resdata"
c      include"resdata1"
c      newmass1 = 0.0
c      newmass2 = 0.0
      loopcount= 0
      massflag = .false.
      loopflag = .false.
      indexmass= -1
      indexs   = -1
      tmass    = 0.0
      nsave = 0

      if(id1.eq.1 .and. id2.eq.1) then
        write(*,*)' aus mommass 3: only nucleons '
        stop
      else if(id1.eq.1 .and. id2.ne.1) then
        ireso   = id2 - 1
        ireact = 1
      else if(id2.eq.1 .and. id1.ne.1) then
        ireso   = id1 - 1
        ireact = 1
      else if(id2.eq.2 .and. id1.eq.2) then
        ireso   = id1 - 1
        ireact = 0
      else
        write(*,*)'****************************************'
        write(*,*)'****************************************'
        write(*,*)'****************************************'
        write(*,*)'id problems in mommass ', id1, id2
        massflag =.true.
        return
      end if

      if(srt.lt.(2.0*rmass+pmass).and. ireact.gt.0) then
        write(*,*)'hiba aus mommass 1: warning srt =' , srt
        return
      else if(srt.lt.(2.0*rmass+2.0*pmass).and. ireact.eq.0) then
        write(*,*)'hiba aus mommass 2: warning srt =' , srt
        return
      end if

*     check for Dimitriev:
      idsum = id1 + id2
      if(idsum.eq.3) then
        ireact = 2
      end if


*     pick the normalzation for the Monte-Carlo decision


      indexs = nint( (srt-sigs0)/delsigs )
      if(indexs.lt.0) indexs = 0
      if(indexs.gt.nsmax-1) then
c        write(*,*)'mommass vigyazz s big', indexs,indexslor , srt
        indexs = nsmax-1
      end if

      indexslor = nint( (srt-sigs0)/delsiglor )
      if(indexslor.lt.0) indexslor = 0
      if(indexslor.gt.nsmaxlor-1) then
c        write(*,*)'mommass 2 vigyazz s big', indexslor , srt
        indexslor = nsmaxlor-1
      end if

      if(ireact.eq.1) then
        rmint   = srt - 2.0*rmass - pmass
        resmass  = rmass + pmass +0.5*rmint
*     do the MC-decision for reactions with NR in the final channel
*     (not DIMITIEV)
        normali = nrmassmax(2,ireso,indexslor)
        help    = nrmassmax(2,ireso,indexslor+1)
        if(help.gt.normali) normali = help

        if(normali.lt.1.0e-08) then
          write(*,*)'normali small', normali
        end if
        loopflag  =.true.
        loopcount = 0

        do while(loopflag)
          loopcount = loopcount + 1
          tmass     = rmass + pmass +rn(iseed)*rmint
          pfin      = ((srt**2-tmass**2+rmass**2)**2)/(4.0*srt**2)
     &                  -rmass**2
          if(pfin.lt.0.0) then
            write(*,*)'mommass fpin = ', pfin
            stop
          end if
          pfin      = sqrt(pfin)
          help  = bwdist(tmass**2,ireso,idec2pi,3,0,0,iwidth,2,ratio)
          help1 = pfin*help*tmass
          x1    = rn(iseed)

          if(normali.lt.1.0e-08) then
            write(*,*)'normali small', normali
          end if

          if(help1/normali.gt.1.05) then
c            write(*,*)'prob with norm  1', help1, normali, indexs, srt
c            stop
          end if

          if(x1.lt.help1/normali) then
            if(normali.lt.1.0e-08) then
              write(*,*)'after normali ', normali
            end if
            loopflag =.false.
            resmass  = tmass
            massflag = .true.
          end if

          if(loopcount.gt.100) then
            loopflag = .false.
            resmass  = rmass + pmass +0.5*rmint
            massflag = .true.
            write(*,*)'vigyazz: mommass1 loopcount.gt.100'
          end if
        end do
      else if(ireact.eq.0) then
*     do the MC-decision for reactions with DD in the final channel
        srt2d0    = 2.0*rmass + 2.0*pmass
        srt2de    = srt - srt2d0
        normali   = ddmassmax(2,indexs)
        help      = ddmassmax(2,indexs+1)
        if(help.gt.normali) normali = help
        loopflag  =.true.
        loopcount = 0
        d1mass   = rmass + pmass +0.3*srt2de
        d2mass   = rmass + pmass +0.3*srt2de

        do while(loopflag)
          loopcount = loopcount + 1

 111      x1d = rn(iseed)
          x2d = rn(iseed)
          if(x1d+x2d .ge.1) goto 111

          tmass = rmass + pmass + x1d*srt2de
          umass = rmass + pmass + x2d*srt2de
          pfin      = ((srt**2-tmass**2+umass**2)**2)/(4.0*srt**2)
     &                  -umass**2
          if(pfin.lt.0.0) then
            write(*,*)'mommass pfin 2 = ', pfin
           goto 111
c            stop
          end if
          pfin      = sqrt(pfin)
          help  = bwdist(tmass**2,ireso,idec2pi,3,0,0,iwidth,2,ratio)
          help1 = bwdist(umass**2,ireso,idec2pi,3,0,0,iwidth,2,ratio)
          help1 = pfin*help*help1*tmass*umass
          x1    = rn(iseed)

          if(normali.lt.1.0e-08) then
            write(*,*)'kleines normali1 ', normali
          end if

          if(help1/normali.gt.1) then
           write(*,*)'vigyazz prob with norm 2',help1,normali,indexs,srt
c            stop
          end if
          if(x1.lt.help1/normali) then
            if(normali.lt.1.0e-08) then
              write(*,*)'kleines normali2 ', normali
            end if
            loopflag =.false.
            d1mass   = tmass
            d2mass   = umass
            massflag = .true.
          end if

          if(loopcount.gt.100) then
            loopflag = .false.
            d1mass   = rmass + pmass +0.3*srt2de
            d2mass   = rmass + pmass +0.3*srt2de
            massflag = .true.
            write(*,*)'vigyazz: mommass2 loopcount.gt.100'
          end if
        end do

      else if(ireact.eq.2) then
*     do the MC-decision for reactions with ND DIMITRIEV

        indexf    = nint( (srt-sigs0)/delsigs )
 10     continue

        if(indexf.gt.nsmax) then
c          write(*,*)'**********************************'
c          write(*,*)'***********crosw *****************'
c          write(*,*)'maximaler wurz s index ', indexs , srt
c          write(*,*)'**********************************'
          indexf = nsmax
        end if

        sigtot    = dimixsec(indexf)

        if(sigtot.lt.1.0e-5) then
          write(*,*)'kleines sigtot ', sigtot
          nsave = nsave + 1
          if(nsave.gt.5) then
            write(*,*)'warning from mommass nsave ', nsave
            stop
          end if
          indexf    = nint( (srtinfree-sigs0)/delsigs )
          goto 10
        end if

        nrandom = 0
 20     continue
        nrandom = nrandom + 1
        sumdsdm   = 0.0
        x1        = rn(iseed)
        loopflag  =.true.
        loopcount = 0

        do while(loopflag)
          dsdm    = dmdimi(1,indexf,loopcount)
          sumdsdm = sumdsdm + dsdm*delmassdi

          help    = sumdsdm/sigtot

          if(x1.le. help) then
            indexmass= loopcount
            loopflag = .false.
            massflag = .true.
          end if
          if(loopcount .eq. nmassmax .and. nrandom.lt.5)goto 20
          if(loopcount .eq. nmassmax.and. nrandom.eq.5) then
            write(*,*)'hiba:problem hier stimmt was nicht in mommass :'
            write(*,*)' srt = ', srt
            write(*,*)'indexs,loopcount  = ', indexs, loopcount
            write(*,*)'sigtot = ', sigtot
            write(*,*)'dsdm, dsdmsum = ', dsdm, sumdsdm
            indexmass = 0
            loopflag = .false.
            massflag = .true.
*            stop
          end if
          loopcount = loopcount + 1
        end do
*     we found a mass, now distribute equally in the 2 MeV bin

        if(indexmass.eq.0) then
          if((srt-dimimass0-rmass).lt.0.0) then
            write(*,*)'prob in mommass 1'
            stop
          end if
          resmass = dimimass0 + rn(iseed)*(srt-dimimass0-rmass)
        else
*        resmass = dimimass0 + (float(indexmass-1)+rn(iseed))*delmassdi
         resmass = dimimass0 + (float(indexmass)-0.5)*delmassdi+
     &            rn(iseed)*delmassdi
         if(resmass+rmass .gt.srtinfree) then
           resmass = srtinfree - rmass - 1.0e-03
         end if
          if(resmass.lt.dimimass0) then
            resmass = dimimass0 + rn(iseed)*(srt-dimimass0-rmass)
          end if
        end if
      end if

      if(.not.massflag) then
        write(*,*)'at this place massflag better be true'
        stop
      end if

      if(massflag) then
        if(ireact.eq.0) then
          newmass1 = d1mass
          newmass2 = d2mass
        else
          if(id1.eq.1) then
            newmass1 = rmass
            newmass2 = resmass
          else if(id2.eq.1) then
            newmass1 = resmass
            newmass2 = rmass
          end if
        end if
      end if

      return
      end
