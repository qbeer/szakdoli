
************************************************************************
       subroutine resmasdi(srt, massmin, mass1, idd, dem)
*
*        evaluate mass distribution for decay to stabile + resonance
*
*        srt     : total energy of the system
*        mass1   : the stabile decay product
*        massmin : the minimal mass for the resonance
*        idd     : positive for baryon resonances, the number is the same
*                       as in bwdist
*                  negative for meson resonances, the number is the same
*                       as in bwmes
*        dem     : the mass of the resonance decay product
*-----------------------------------------------------------------------

       implicit none
       include"common"
       include"cominput"

       real*8    srt, s, mass1, massmin, dem, rn, bwd, bwdist, bwmes
       real*8    p30, p32, rmint, dem2, dummyf(9)
       integer idd, ncount
       logical lcount

c       write(*,*) srt, massmin, mass1, idd, dem

       lcount = .true.
       ncount = 0
       rmint = srt-massmin-mass1
       s = srt**2
       p30= 0.25*(s-massmin**2+mass1**2)**2/s-mass1**2
       do while(lcount)
         ncount = ncount+1
         dem = massmin + rn(iseed) * rmint
         dem2 = dem**2
        if(idd.ge.0) bwd=bwdist(dem2,idd,idec2pi,iresmode,0,0,iwidth,3,
     &     dummyf)
         if(idd.lt.0)
     &          bwd=bwmes(dem2,abs(idd),idec2pi,iresmode,0,0,iwidth,3)
         if(masdis.eq.1) then
           p32= 0.25*(s-dem2+mass1**2)**2/s-mass1**2
           bwd = sqrt(max(0.0,p32/p30))*bwd
         endif
         if(rn(iseed).le. bwd) then
           lcount =.false.
         end if
         if(ncount.ge.100) then
           dem = massmin + 0.5*rmint
           lcount =.false.
         end if
       end do

c       write(*,*) ' in resmasdi ', srt, massmin, mass1, idd, dem
*
        return
        end



