c--------1---------2---------3---------4---------5---------6--------7--
c  Block data for LQCD YN and YY potentials measured with the Kei-conf
c
c                                        by Takashi Inoue, Aug 2016
c--------1---------2---------3---------4---------5---------6--------7--
       block data LQCD_YN_POT_Kei
       implicit   none

       real*8  prm1(7), prm2(7), prm3(7)
       common /potential_parameter/ prm1, prm2, prm3

       real*8  prmc1(7), prmc2(7), prmc3(7)
       common /potential_parameter_Vc/ prmc1, prmc2, prmc3

       real*8  prmt1(6), prmt2(6), prmt3(6)
       common /potential_parameter_Vt/ prmt1, prmt2, prmt3

       integer    i

c----------------------
c FL_27
       data (prm1(i),i=1,7) /
     &   1336.6549d0  ,  2.3742678d0  ,  1433.5221d0  ,  53.659831d0  ,
     &   141770.91d0  , 0.52964043d0  ,  2.5293523d0    /

c FL_8S
       data (prm2(i),i=1,7) /
     &   5638.7522d0  ,  12.951368d0  ,  411.74169d0  ,  1.2813190d0  ,
     &   1636.5727d0  ,  15.959756d0  ,  3.4365892d0  /

c FL_1
       data (prm3(i),i=1,7) /
     &  -966.67715d0  ,  3.1380113d0  , -84.682264d0  , 0.69326043d0  ,
     & -0.92626920d10 , 0.25353068d-2 ,  3.6762249d0  /

c----------------------
c FL_10*
       data (prmc1(i),i=1,7) /
     &   819.16615d0  ,  1.9129865d0  ,  1416.7240d0  ,  45.450617d0  ,
     &   51607.561d0  , 0.60079497d0  ,  2.1455409d0  /

c FL_8A
       data (prmc2(i),i=1,7) /
     &   223.30741d0  ,  39.544058d0  ,  193.89482d0  ,  1.6640125d0  ,
     &   242551.91d0  , 0.15265522d0  ,  2.2494186d0  /

c FL_10
       data (prmc3(i),i=1,7) /
     &   2180.6380d0  ,  7.9816661d0  ,  579.61863d0  ,  3.7219714d0  ,
     &   32656.853d0  ,  2.1772221d0  ,  4.1936891d0  /

c----------------------
c FL_10*
       data (prmt1(i),i=1,6) /
     &   11301.227d0  ,  9.4780722d0  ,  1.8860888d0  , -11202.034d0  ,
     &   9.4839043d0  ,  1.8782561d0  /

c FL_8A
       data (prmt2(i),i=1,6) /
     &  -15724.464d0  ,  3.2670958d0  ,  14.519991d0  ,  6.6975009d0  ,
     &   24.784160d0  ,  2.2239908d0  /

c FL_10
       data (prmt3(i),i=1,6) /
     &  -2109.3066d0  ,  16.472853d0  ,  1.2352148d0  ,  2092.6904d0  ,
     &   16.475335d0  ,  1.2303626d0  /
c----------------------
       end

c--------1---------2---------3---------4---------5---------6--------7--
c                END                           
c--------1---------2---------3---------4---------5---------6--------7--