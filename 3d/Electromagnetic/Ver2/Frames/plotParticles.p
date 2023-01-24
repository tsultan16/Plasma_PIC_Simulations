reset
set pm3d map
set palette color

nx = 100
ny = 100
nz = 100

set dgrid3d nx,nz

pi = 3.14159265359

nt=200
tskip=5
N=1


Lx = nx
Ly = ny
Lz = nz
dx = Lx/nx
dy = Ly/ny
dz = Lz/nz

xmin = 0.5*dx
xmax = xmin+Lx
ymin = 0.5*dy
ymax = ymin+Ly
zmin = 0.5*dz
zmax = zmin+Lz

set key font ",12"

set autoscale



set terminal png size 1920,1080
set size square

E_max = 0.05

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_xy_t=".i.".txt"
  
  set output "Fields/fields_xy_t=".i.".png"
  set title "Time Step = #".j
  
  i = j

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  plot filename using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS" ,\
  "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  unset title
  unset output
}



set terminal gif animate delay 10 size 1920,1080
set size square
set output 'particles_xy_1.gif'


do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/fields_xy_t=".i.".txt"
  set title "Time Step = #".j
  
  i = j

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  plot filename using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS" ,\
  "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  unset title
}

unset output

reset
exit

