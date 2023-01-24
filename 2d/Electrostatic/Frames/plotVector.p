set pm3d map
set palette color

nx = 32
ny = 32

set dgrid3d nx,ny

pi = 3.14159265359

nt=100
tskip=1
N=1#4*nx*ny

Lx = 2*pi
Ly = 2*pi
dx = Lx/nx
dy = Ly/ny
xmin = 0.5*dx
xmax = xmin+Lx
ymin = 0.5*dy
ymax = ymin+Ly

set key font ",12"

set autoscale


set terminal gif animate delay 30 size 2560,1440
set output 'E_vector.gif'
do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting fields. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_t=".i.".txt"
  set title "Time Step = #".j


  stats filename using 5 nooutput
  Ex_min = STATS_min 
  Ex_max = STATS_max 
  stats filename using 6 nooutput
  Ey_min = STATS_min
  Ey_max = STATS_max

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  plot filename using 1:2:($5/(Lx*Ex_max)):($6/(Ly*Ex_max)) with vectors head filled lt 2 ,\
  "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 3.0 lc rgb "blue" title "IONS" ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.5 lc rgb "red" title "ELECTRONS" 
  unset title
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

}

unset output
reset


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

