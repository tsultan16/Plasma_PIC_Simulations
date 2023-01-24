set pm3d map
set palette color

nx = 32
ny = 32

set dgrid3d nx,ny

nt=600
tskip=1
N=1#4*nx*ny

set key font ",12"

set autoscale



set terminal gif animate delay 5 size 2560,1440
set size square
set output 'particles.gif'
do for [j=0:nt/tskip - 1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/nt)
  set title "Time Step = #".(i+1)
  set xlabel "x"
  set ylabel "y"
  set xrange [-0.1:6.3]
  set yrange [-0.1:6.3]
  #plot "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" ,\
  #"particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS"  
  plot "particles_i.txt" every ::0::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.2 lc rgb "blue" title "IONS" ,\
  "particles_e.txt" every ::0::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.2 lc rgb "red" title "ELECTRONS"  
  
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  unset xrange

  
}
unset output
reset

exit

set terminal gif animate delay 20 size 2560,1440
set output 'phase_space.gif'
do for [j=0:nt/tskip - 1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/nt)



  set xlabel "x"
  set ylabel "vx"
  #set xrange [-0.1:6.3]
  #set yrange [-0.2:0.2]
  plot "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:3 with point pointtype 7 pointsize 0.5 lc rgb "blue" title "IONS" ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:3 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange
  
  
}
unset output
reset


stats 'energy.txt' using 2
KE_min = STATS_min
KE_max = STATS_max 


  set terminal png size 2560,1440
  filename="energy_momentum.png"
  set output filename
  
  set multiplot layout 4,2 columnsfirst 
  
  set xlabel "t"
  set ylabel "KE"
  set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:2 with lines lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE"
  set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:3 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "t"
  set ylabel "KE+ESE"
  set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:4 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
    
  set xlabel "t"
  set ylabel "Px"
  set yrange [-0.005:0.005]
  plot "energy.txt" using 1:5 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "t"
  set ylabel "vd1"
  #set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:6 with lines lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "t"
  set ylabel "vspread1"
  #set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:7 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
 
  set xlabel "t"
  set ylabel "vd2"
  #set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:8 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
    
  set xlabel "t"
  set ylabel "vspread2"
  #set yrange [-0.005:0.005]
  plot "energy.txt" using 1:9 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  unset multiplot
  unset output




stats 'fields_horcut.txt' using 2
rho_min = STATS_min
rho_max = STATS_max 
stats 'fields_horcut.txt' using 3
phi_min = STATS_min
phi_max = STATS_max 
stats 'fields_horcut.txt' using 4
Ex_min = STATS_min
Ex_max = STATS_max 
stats 'fields_horcut.txt' using 5
Ey_min = STATS_min
Ey_max = STATS_max 

set terminal gif animate delay 30 size 2560,1440
set output 'fields_horcut.gif'
do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting fields. % Completed =".((100*i)/nt)

  
  set title "Time Step = #".(i+1)
  set multiplot layout 4,1 columnsfirst
  
  set xlabel "x"
  set ylabel "RHO"
  set yrange [rho_min:rho_max]
  plot "fields_horcut.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "PHI"
  set yrange[phi_min:phi_max]
  plot "fields_horcut.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "Ex"
  set yrange[Ex_min:Ex_max]
  plot "fields_horcut.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  
  set xlabel "x"
  set ylabel "Ey"
  set yrange[Ey_min:Ey_max]
  plot "fields_horcut.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  
  unset multiplot
}

unset output




set terminal gif animate delay 20 size 2560,1440
set output 'fields.gif'
do for [j=0:nt/tskip - 1] {

  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_t=".i.".txt"
  set title "Time Step = #".j
  set multiplot layout 1,4 rowsfirst  

  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3 title 'RHO'
  unset title
  unset xlabel
  unset ylabel

  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:4 title 'PHI'
  unset xlabel
  unset ylabel

  set label 1 'Ex' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:5 
  unset title
  unset xlabel
  unset ylabel

  set label 2 'Ey' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:6
  unset xlabel
  unset ylabel
  unset multiplot
  
}
unset output



reset