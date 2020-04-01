reset
#set yrange [0:1]
plot 'RT.log' u 1:2 w l t 'Reflexion'
replot 'RT.log' u 1:3 w l t 'Transmission'
replot 'RT.log' u 1:($2+$3) w l t 'Norma total'
set xlabel 'time (a.u.)
pause -1
plot 'der_RT.log' u 1:2 w l t 'Der Reflexion'
replot 'der_RT.log' u 1:3 w l t 'Der Transmission'
pause -1
