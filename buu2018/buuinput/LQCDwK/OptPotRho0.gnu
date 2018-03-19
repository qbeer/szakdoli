set term postscript eps enhanced color 26
#set term pdf
set out "OptPotRho0.eps"
reset

#set encoding iso

#set linestyle 1 linetype 1 lw 4
#set linestyle 2 linetype 2 lw 4
#set linestyle 4 linetype 4 lw 4
#set linestyle 6 linetype 6 lw 4
#set linestyle 8 linetype 8 lw 4
#set linestyle 12 linetype 12 lw 4
#set linestyle 20 linetype 1 lw 6

#lambda=.125
#fact=mu2+LLambda/3.0+0.5*lambda*kkip
#ss11(x)=x*2*sqrt(fact)

set xrange [0.0:1000]
#set yrange [-30.0:200]
set xlabel "p [MeV/c]"
set ylabel "U({/Symbol r}_0,p) [MeV]"
#set ylabel "d{/Symbol s}^{dilepton}/dM"
#set xtics("" -900,"-800" -800,"" -700,"-600" -600,"" -500,"-400" -400,"" -300,"-200" -200,"" -100,"0" 0)
#set xtics("" 0.05,"0.1" 0.1,"" 0.15,"0.2" 0.2,"" 0.25,"0.3" 0.3,"" 0.35,"0.4" 0.4,"" 0.45,"0.5" 0.5)
#set xtics("10^7" 10000000,"" 5000000,"" 2000000,"10^6" 1000000,"" 500000,"" 200000,"10^5" 100000,"" 50000,"" 20000,"10^4" 10000,"" 5000,"" 2000,"10^3" 1000,"" 500,"" 200,"10^2" 100,"" 50,"" 20,"10" 10,"" 5,"" 2,"1" 1,"" 0.5,"" 0.2,"10^{-1}" 0.1,"" 0.05,"" 0.02,"10^{-2}" 0.01,"" 0.0000001,"10^{-8}" 0.00000001,"" 0.000000001)
#set ytics("10^7" 10000000,"" 5000000,"" 2000000,"10^6" 1000000,"" 500000,"" 200000,"10^5" 100000,"" 50000,"" 20000,"10^4" 10000,"" 5000,"" 2000,"10^3" 1000,"" 500,"" 200,"10^2" 100,"" 50,"" 20,"10" 10,"" 5,"" 2,"1" 1,"" 0.5,"" 0.2,"10^{-1}" 0.1,"" 0.05,"" 0.02,"10^{-2}" 0.01,"" 0.005,"" 0.002,"10^{-3}" 0.001,"" 0.0000001,"10^{-8}" 0.00000001,"" 0.000000001)
#set logscale y
#set logscale x
set size 0.8,1.0
set key at 900.,-15
#set title "{26 Symbol p}^+{26 Symbol p}^- annihilation cross section"
     
plot "spectrum_hyperon_SNM_Kei_939_AV18_rho1000.bin" i 0  us (197*column(1)):4 tit "{/Symbol L}" w l ls 1, \
     "spectrum_hyperon_SNM_Kei_939_AV18_rho1000.bin" i 0  us (197*column(1)):8 tit "{/Symbol S}" w l ls 2, \
     "spectrum_xi_SNM_Kei_939_AV18_rho1000.bin" i 0  us (197*column(1)):4 tit "{/Symbol X}" w l ls 3, \
     "EOSMMr10.dat" i 0  us (1000*column(1)):(470*column(2)) tit "0.47* N" w l ls 5
set out
