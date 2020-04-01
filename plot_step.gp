reset
set title 'Initial Wave Packet'
set xlabel 'x'
set xrange[-5:5]
set xtics 2.5
set ticslevel 0.000001

plot 'wp_t0.log' u 1:2 w l t 'Real part' lc rgb 'cyan' 
replot 'wp_t0.log' u 1:3 w l t 'Imaginary part' lc rgb 'green'
replot 'wp_t0.log' u 1:4 w l t 'Wave function squared' lc rgb 'black' 
replot 'wp_t0.log' u 1:5 w l t 'Potential' lc rgb 'red'
pause -1
reset
