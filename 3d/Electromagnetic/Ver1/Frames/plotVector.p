set pm3d map
set palette color

nx = 50
ny = 50

set dgrid3d nx,ny

pi = 3.14159265359

nt=200
tskip=2
N=500

Lx = 1.0
Ly = 1.0
dx = Lx/nx
dy = Ly/ny
xmin = 0.5*dx
xmax = xmin+Lx
ymin = 0.5*dy
ymax = ymin+Ly

set key font ",12"

set autoscale


 
 
stats 'fields_horcut.txt' using 4
E_min = STATS_min
E_max = STATS_max 
stats 'fields_horcut.txt' using 3
ne_min = STATS_min
ne_max = STATS_max 
stats 'fields_vercut.txt' using 2
B_min = STATS_min
B_max = STATS_max 


set terminal gif animate delay 10 size 2560,1440
set output 'fields_horcut.gif'
do for [j=0:nt/tskip-1] {

  i = j*tskip
  print "Plotting fields. % Completed =".((100*i)/nt)

  
  set title "Time Step = #".(i+1)
  set multiplot layout 1,2 columnsfirst
  
  i = j
  
  set xlabel "x"
  set ylabel "rho"
  #set yrange[ne_min:ne_max]
  plot "fields_horcut.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "y"
  set ylabel "Bz"
  #set yrange[B_min:B_max]
  plot "fields_vercut.txt" every ::i*ny+1::ny+(i*ny)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange 
  unset multiplot
}

unset output




E_max = 8.6e-11


set terminal gif animate delay 10 size 1920,1080
set size square
set output 'E_vector.gif'
do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting fields. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_t=".i.".txt"
  set title "Time Step = #".j


  #stats filename using 7 nooutput
  #Jx_min = STATS_min 
  #Jx_max = STATS_max 
  #stats filename using 8 nooutput
  #Jy_min = STATS_min
  #Jy_max = STATS_max
  i = j

  set multiplot layout 1,2 columnsfirst

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  plot filename using 1:2:($5/(0.3*nx*E_max)):($6/(0.3*ny*E_max)) with vectors head filled lt 2 ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS" ,\
  "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange


  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  splot filename using 1:2:4 title "Bz"
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange
  
  unset title
  unset multiplot
}

unset output
reset


exit



set terminal png size 1920,1080
set size square
do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting fields. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_t=".i.".txt"
  filename2="Fields/frame_t=".i.".png"
  
  set output filename2

  #stats filename using 7 nooutput
  #Jx_min = STATS_min 
  #Jx_max = STATS_max 
  #stats filename using 8 nooutput
  #Jy_min = STATS_min
  #Jy_max = STATS_max
  i = j
 
  set title "Time Step = #".j

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  plot filename using 1:2:($5/(0.4*nx*E_max)):($6/(0.4*ny*E_max)) with vectors head filled lt 2 ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "red" title "ELECTRONS" ,\
  "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange
  

  unset title
  unset multiplot
  unset output
}


reset

exit







set terminal gif animate delay 10 size 1920,1080
set size square
set output 'E_vector.gif'
do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting fields. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_t=".i.".txt"
  set title "Time Step = #".j


  #stats filename using 7 nooutput
  #Jx_min = STATS_min 
  #Jx_max = STATS_max 
  #stats filename using 8 nooutput
  #Jy_min = STATS_min
  #Jy_max = STATS_max
  i = j

  set multiplot layout 1,3   

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  plot filename using 1:2:($5/(10*Lx*E_max)):($6/(10*Ly*E_max)) with vectors head filled lt 2 ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "red" title "ELECTRONS" ,\
  "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  splot filename using 1:2:5 title "Ex"
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange
  
  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  splot filename using 1:2:6 title "Ey"
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange
  

  unset title
  unset multiplot
}

unset output
reset

exit


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
reset

