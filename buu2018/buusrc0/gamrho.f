************************************************************************
*                                                                      *
      subroutine gamrhoin(num,iseed,massta,egamma,delt,dengamro,
     &               masspec,drhom)
*                                                                      *
*       purpose:    calculating the dilepton production from:          *
*                                                                      *
************************************************************************
      implicit none
      real*8 romass, gammro,egamma,delt,drhom,rhomin,rhomax,rom2,pm2
      real*8 dilmass,rogam,roform,pmrho,betaz,gamma,ppio,prob,dlife
      real*8 dzed,xx,yy,zz,rr,ppx,ppy,ppz,epion,rn
      integer num,iseed,massta,maxmas,ii,jj,i1,imas,lstep
      parameter     (romass = 0.775, gammro = 0.149)
      parameter     (maxmas = 30)
      include 'common'
*----------------------------------------------------------------------*
      real*8     dengamro(10000)
      integer  masspec(3,maxmas)
*----------------------------------------------------------------------*
       rhomin= 2*pmass
       rhomax= egamma
       rom2  = romass**2
       pm2   = pmass**2
       do ii=1,maxmas
         masspec(1,ii) = 0
         masspec(2,ii) = 0
         masspec(3,ii) = 0
       enddo
       do 5000 ii=1,num
         jj=(ii-1)*maxp
         i1=(ii-1)*massta + max0(1,nint(float(massta)*rn(iseed)))
 1000    dilmass= rhomin + rn(iseed) * (rhomax-rhomin)
         rogam = ((dilmass**2-4.0*pm2)/(rom2-4.0*pm2))**1.5
         rogam = rogam * gammro
         roform = 1.0/(4.0*((dilmass-romass)/rogam)**2+1.0)
         if(rn(iseed) .gt. roform)       goto 1000
         imas = nint((dilmass-2.0*pmass)/drhom)+1
         imas = min0(imas,maxmas)
         masspec(1,imas) = masspec(1,imas) + 1
         pmrho = sqrt(egamma**2-dilmass**2)
         betaz = pmrho / dilmass
         gamma = egamma/dilmass
         ppio=sqrt(0.25 * dilmass**2 - pm2)
         prob = exp(-delt*rogam/hbc)
         lstep = -1
 2000    lstep = lstep + 1
         if(rn(iseed).lt.prob)       goto 2000
         dlife = float(lstep) * delt
         dzed = dlife* betaz
 3000    xx = rn(iseed)
         yy = rn(iseed)
         zz = rn(iseed)
         rr = sqrt(xx**2+yy**2+zz**2)
         if(rr .lt.0.01 .or. rr.gt.1.0) goto 3000
         ppx = xx/rr * ppio
         ppy = yy/rr * ppio
         ppz = zz/rr * ppio
***   the outgoing pion momentum : (ppx,ppy,ppz) in the rho cms    ***
         epion      = 0.5 * dilmass
*
*----------------------------------------------------------------------*
*
*  ipi(1,i) = particle type (pion, eta, ...)
*  ipi(2,i) = particle charge
*  ipi(3,i) = the serial number for the parent resonance of meson i
*  ipi(4,i) = how many times the meson was created + for pions - for eta
*  ipi(5,i) = type of the parent resonance (id(1,ipi(3,i))))
*
*----------------------------------------------------------------------*
          ppi(1,jj+1)  = ppx
          ppi(2,jj+1)  = ppy
          ppi(3,jj+1)  = gamma*(ppz + betaz * epion)
          rpi(1,jj+1)  = r(1,i1)
          rpi(2,jj+1)  = r(2,i1)
          rpi(3,jj+1)  = r(3,i1) + dzed
          ipi(1,jj+1)  = 1
          ipi(2,jj+1)  = 1
          ipi(3,jj+1)  = i1
          ipi(4,jj+1)  = 1
          ipi(5,jj+1)  = 1
*
          ppi(1,jj+2)  = -ppx
          ppi(2,jj+2)  = -ppy
          ppi(3,jj+2)  = gamma*(-ppz + betaz * epion)
          rpi(1,jj+2)  = r(1,i1)
          rpi(2,jj+2)  = r(2,i1)
          rpi(3,jj+2)  = r(3,i1) + dzed
          ipi(1,jj+2)  = 1
          ipi(2,jj+2)  = -1
          ipi(3,jj+2)  = i1
          ipi(4,jj+2)  = 1
          ipi(5,jj+2)  = 1
       dengamro(ii)=rhb(nint(r(1,i1)),nint(r(2,i1)),nint(r(3,i1)))/rho0
c       write(6,*) rho0
*
 5000 continue
      return
      end
c     debug unit(6),subchk,trace,init(numbpio)
c     debug unit(6),subchk,trace
c     end debug
************************************************************************
*                                                                      *
      subroutine gamrhout(num,drhom,densste,dengamro,masspec)
*                                                                      *
*       purpose:    calculating the dilepton production from:          *
*                                                                      *
************************************************************************
      implicit none
      real*8 romass,gammro,drhom,densste,pm2,rom2,rogam,roform,dilmass
      real*8 px1,py1,pz1,e1,px2,py2,pz2,e2,xx,yy,zz
      real*8 s,srt
      real*8 rn
      integer maxde,num,maxmas,ii,ni,nde,numpio,jj,jjn,imas
      parameter     (romass = 0.775, gammro = 0.149)
      parameter     (maxde = 20)
      parameter     (maxmas= 30)
*----------------------------------------------------------------------*
      real*8   dengamro(10000)
      integer  indpi(10), numbpio(0:4,maxde),idensdis(maxde)
      integer  masspec(3,maxmas),masspecd(maxmas,maxde)
      integer  masspe(maxmas)
*----------------------------------------------------------------------*
      include 'common'

       pm2   = pmass**2
       rom2   = romass**2
       do ii=0,4
       do ni=1,maxde
         numbpio(ii,ni) = 0
       enddo
       enddo
       do ii=1,maxmas
         do ni=1,maxde
           masspecd(ii,ni) = 0
         enddo
       enddo
       do 5000 ii=1,num
c-----------------------------------------------------------------------
c           density dependence
         nde  = nint(dengamro(ii)/densste+1.0)
         nde  = min(nde,maxde)
         numpio = 0
         do 4000 jj=1,maxp
          jjn=(ii-1)*maxp+jj
            if(ipi(1,jjn) .eq.1) then
              numpio = numpio + 1
              indpi(numpio) = jjn
            endif
 4000    continue
         numpio=min0(numpio,3)
         numbpio(numpio,nde) = numbpio(numpio,nde) + 1
         idensdis(nde) = idensdis(nde) + 1
         if(numpio.eq.2) then
           px1=ppi(1,indpi(1))
           py1=ppi(2,indpi(1))
           pz1=ppi(3,indpi(1))
           e1 =sqrt(pm2+px1**2+py1**2+pz1**2)
           px2=ppi(1,indpi(2))
           py2=ppi(2,indpi(2))
           pz2=ppi(3,indpi(2))
           e2 =sqrt(pm2+px2**2+py2**2+pz2**2)
           s  =(e1+e2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2
           srt = sqrt(s)
           imas = nint((srt-2.0*pmass)/drhom)+1
           imas = min0(imas,maxmas)
           masspec(2,imas) = masspec(2,imas) + 1
           masspecd(imas,nde) = masspecd(imas,nde) + 1
           if(ipi(3,indpi(1)).eq.ipi(3,indpi(2)))
     &       numbpio(4,nde) = numbpio(4,nde) + 1
           if(ipi(3,indpi(1)).eq.ipi(3,indpi(2)))
     &       masspec(3,imas) = masspec(3,imas) + 1
           if(ipi(4,indpi(1)).eq.1 .and. ipi(4,indpi(2)).eq.1)
     &       masspe(imas) = masspe(imas) + 1
         endif
 5000 continue
        write(isum,'(''c:rhomass spectra(eredetigamma)'')')
        do jj = 1,maxmas
          dilmass= 2.0 * pmass + float(jj) * drhom
          rogam = ((dilmass**2-4.0*pm2)/(rom2-4.0*pm2))**1.5
          rogam = rogam * gammro
          roform = 1.0/(4.0*((dilmass-romass)/rogam)**2+1.0)
          write(isum,'(2f15.4)') dilmass, roform
        enddo
        write(isum,'(''c:rhomass spectra(total)'')')
        do jj = 1,maxmas
          dilmass= 2.0 * pmass + float(jj-1) * drhom
          write(isum,'(f15.4,4i15)') dilmass, (masspec(ii,jj), ii=1,3),
     &              masspe(jj)
        enddo
        do nde=1,6
        write(isum,'(''c:rhoms spect(dens) '',f8.2)') densste*float(nde)
        do jj = 1,maxmas
          dilmass= 2.0 * pmass + float(jj-1) * drhom
          write(isum,'(f15.4,i15)') dilmass, masspecd(jj,nde)
        enddo
        enddo
        write(isum,'(''c:den dist for rho '')')
        do nde=1,maxde
         write(isum,'(f15.4,6i8)') densste*float(nde-1),idensdis(nde),
     &    (numbpio(ii,nde),ii=0,4)
        enddo
*
      return
      end
