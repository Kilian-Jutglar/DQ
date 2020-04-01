## Estructura per fer un gif a gnuplot
reset

set terminal gif size 1000,600 animate delay 1 loop 0 optimize
set output 'euler_pbc.gif'
datafile = 'wp_trajectory.log'
stats datafile
numblo=STATS_blocks

set title 'Initial Wave Packet'
set xrange[-5:5]
set yrange[-0.2:2.2]
set xlabel 'x'
set xtics 2.5
##set ticslevel 0.000001

do for [i=0:numblo-2] {
plot datafile index i u 1:2 w l t 'Wave function squared' lc rgb 'black' ,datafile index i u 1:3 w l t 'Potential' lc rgb 'red'
}

