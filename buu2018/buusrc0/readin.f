
      subroutine readin
*----------------------------------------------------------------------*
*                                                                      *
*    this routine reads the input from the main control file           *
*                                                                      *
*     a l l variables are stored in the common block :                 *
*           cominp                                                     *
*                                                                      *
*----------------------------------------------------------------------*
      implicit none
      character*80 totoutp, barout, pioout, cpiout, cpiout1, delout,
     #  dilout, mesout, mestimout, omeout, kaoout, ppiout, gamout, fsout
     # ,spfout,JPsiout
      include"common"
      include"cominput"
      include"com_kminu"
      include"com_pert"
      include"com_cont_epair"
      integer i,ijp2,ijp1
      integer icbro1,ivecmatt1
      real rhomshift1,rhowbroad1,omemshift1,omewbroad1
      NAMELIST /adatok/
     # totoutp, barout, pioout, cpiout, cpiout1, delout,
     # dilout, mesout, mestimout, omeout, kaoout,
     # ppiout, gamout, denout, fsout, spfout,JPsiout,
     # MASSTA, MSTAPR, MASSPR, MSPRPR, MPION, ELAB,
     # B, DT, NTMAX, ISUBS, NUM, NUMPAUL,
     # iseed, icoll, INSYS, IPAULI, icomp, ipou, istrangepot,ireslife,
     # IAVOID,IRESCOL,mescol,ippbarcoll,masdis,IWIDTH,iresmode,idec2pi,
     # dLIFE,IDIPI, IABSO, IEMPI, IPIAB, IPICO,
     # nbound, boxx, boxz, NTCELL, CELL,
     # NTHER,DELR, DELP, NPARTCEL, NFREQ,
     # TIMINSP, IFIT, NFRBA, IBASPE, IPARTI, TMINBA, TMAXBA,
     # BINKBA, BINTBA, BINMBA, DANGLEB, IFRAMEB, NANGLEB, ANBA,
     # NFRPI, IPISPE, TMINPI, TMAXPI, BINKPI, BINTPI, BINMPI,
     # DANGLEP, IFRAMEP, NANGLEP, ANPI, NFRPIRA,
     # idilepton, ivmesdil,  IBREMS, IpiNbrems,Iresdalitz, Imesdalitz, 
     # IFILT,IDILTRA,NDLMASlow,NDLMAS,DLMASL,DDLMASlow,DDLMAShigh,
     # NY,NQT,NF,QTMAXI, YMAXSCAL, YMINSCAL, NBREMS,npiNbrems,
     # IDELTAD, NDALITZ, IRAND, IDAL, NDILTIM, TIMINDI,
     # IDILDELM, IFORMRO, DENSSTE, DIMASMI, IPERTPI, NDMAS, ISZOG, MQT,
     # MY, MF, QTMAXP, NPERTPI, IDILPER, NDILPER, NDelMAS, ISZOGD,
     # ikaon, ikaonpot, ikaoncr, pidelka,
     # PKAONNUM, NQTK, NYK, nfk, YMAXSK,YMINSK,
     # qtmaxk, NPL, PLSTEP, ANGKAON, NPHIKAON, timkao,
     + i_phi, iphi_pot, max_phi, iphi_cr, pphi_num, 
     + iphi_col, iphi_dec, i_epair, i_kminu, i_kminu_pot, i_kminu_cr,
     + i_kminu_coll, i_ksi, i_ksi_pot,
     + i_JPsi,JPsi_scale_factor,iJPsimat,i_charm_matt_dec,
     & JPsi_massshift,JPsi_widthshift,
     + JP_dlmas,JP_dlmasl,JP_ddlmas,JP_y,JP_qt,JP_f,JP_qtmaxi,
     & JP_ymaxscal,JP_yminscal,iDrellYan,iDmes,IR_dlmasl,IR_ddlmas,
     & IR_qtmaxi,IR_yminscal,IR_ymaxscal,IR_dlmas,IR_y,IR_qt,IR_f,
     # IOMEGA, IMESON, NMESMA, DMESMA, XMESMI, XMESM2, NMESM2,
     # DMESM2, NCALMES, SIGPIOM, GAMMESMA,
     & ivecmatt,icbro, rhomshift,rhowbroad,omemshift,omewbroad,
     & igamma, igamsys, pgamma,
     # pgammad, nplgam, plstepga, plmingam, nphigam, nanglgam,
     # anglegam,nqtg, nyg, nfg, qtmaxg, yminsg, ymaxsg, ngamtim,timinga,
     # ifserat ,ipipot, ireapl, isplipi, iorder,ideltapot,
     # INUCLPRI, NNUCLPRI, IPIONPRI, icpipri, NPIONPRI, ideltpri,
     # IDILPRI,IMESPRI, IMESTIMPRI, IOMEPRI, ikaopri, igampri, IPpiPRI,
     # ITHERMO, IDENSPR, icelpri,ncont,ihades,TIMECPU
*=======================================================================

      common/com_matter/ icbro1,ivecmatt1,rhomshift1,rhowbroad1,
     &                 omemshift1,omewbroad1

      include"resdata1"

      rdist=2.9     ! intial distance of the centers in z direction: r1+r2+rdist
      delpi=3.0     ! for bigger distances the collision is thrown away
      dispi=2.5     ! maximal b_max in fm for pions
c      smu=2.175     ! range of yukawa potential (1/fm)
      rmesdec = 2.0
      isnndd  = 0
      isddnn  = 0
      isnrnrp = 0
      isstatep= 0
      isstatea = 0
      isdeltaonly = 0
      inotwopi = 0
      
******************** end of test variables ************************
c         mass     width     BR_p+p_   BR_e+e-  BR_JPsi+X            
      data((JPsi_prop(ijp1,ijp2),ijp2=1,5),ijp1=1,3)
     &  /3.0969, 0.0000929,  2.12D-3,  0.16 , 0.,
     &   3.6861, 0.000296,   2.88D-4,  0.00789, 0.61,
     &   3.7731, 0.027,      1.50D-4,  0.96D-5, 0.004/

c         mass     width     BR_e+e-
      data((Dmes_prop(ijp1,ijp2),ijp2=1,3),ijp1=1,2)
     &  /1.8696, 0.00019,  0.165,
     &   1.8648, 0.00048,  0.065/

      data(sigJPsib(ijp1),ijp1=1,3)
     &  /4.18, 7.6, 7.6/
      data(sigJPsim(ijp1),ijp1=1,3)
     &  /4.18, 7.6, 7.6/

      id_JPsi(1) = 2 ! J/Psi
      id_JPsi(2) = 3 ! Psi'(3686)
      id_JPsi(3) = 4 ! Psi(2s)(3773)
      id_Dmes(1) = 5 ! D+
      id_Dmes(2) = 6 ! D-
      id_Dmes(3) = 7 ! D0+ cdbar
      id_Dmes(4) = 8 ! D0- cbard

********************end of test variables ************************

      write(*,*) 'vor einlesen'
c      open(22,file='filename.dat',status='old',iostat=i)
c      read(22,fnames)
c      write(6,fnames)
c      close(22)
      read(5,adatok)
      ikaondi = ikaon                   !  hw
      JP_masl = JP_dlmasl
      JP_dmas = JP_ddlmas
      if (i_phi .eq. 0) iphi_dec = 0
      write(*,*) 'beolvasas utan', rdist, delpi, dispi
      call inputwri(0)
      open(isum,file=totoutp,status='new',iostat=i)
      if(i.ne.0) open(isum,file=totoutp,status='old',iostat=i)
      call inputwri(1)
      if(inuclpri.eq.1)
     &      open(mnucpri,file=barout,status='unknown')
      if(ipionpri.eq.1)
     &      open(mpiopri,file=pioout,status='unknown')
      if(idilpri.eq.1)
     &      open(mdilpri,file=dilout,status='unknown')
      if(ikaopri.eq.1)
     &      open(mkaopri,file=kaoout,status='unknown')
      if(iomepri.eq.1)
     &      open(momepri,file=omeout,status='unknown')
      if(imespri.eq.1)
     &      open(mmespri,file=mesout,status='unknown')
      if(imestimpri.eq.1)
     &      open(mmestimpri,file=mestimout,status='unknown')
*      if(ithermo.eq.1)
*     &      open(mterpri,file=terout,status='unknown')
      if(ippipri.eq.1)
     &      open(mppipri,file=ppiout,status='unknown')
      if(icpipri.eq.1 )
     &      open(mcpipri,file=cpiout,status='unknown')
      if(icpipri.eq.1 )
     &      open(mcpipri1,file=cpiout1,status='unknown')
      if(ifserat.ge.1 )
     &      open(mfserat,file=fsout,status='unknown')
      open(mspfpri,file=spfout,status='unknown')
      if(i_JPsi.ge.1 )
     &      open(m_JPsipri,file=JPsiout,status='unknown')
*--------------------------------------------------------------------

      return
      end
