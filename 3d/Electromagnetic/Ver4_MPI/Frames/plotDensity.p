reset
set pm3d map
set palette color
set palette rgb 33,13,10

nx = 120
ny = 60
nz = 60

pi = 3.14159265359

nt=1400
tskip=2


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



set dgrid3d nx,nz
set terminal png size 1920,1080

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting charge density. % Completed =".((100*i)/(nt/tskip))
    
  i = j*tskip
  print "Time Steps Completed =".i


  filename="Snapshots/density_xz_t=".i.".dat"
  set output "Images/density_xz_t=".i.".png"

  set title "ne, ni, Time Step = #".i
  
  set multiplot layout 1,2 columnsfirst 

  set xlabel "x"
  set ylabel "z"
  set cbrange [0:5]
  set zrange [0:5]
  splot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:3
  unset title
  unset xlabel
  unset ylabel

  set xlabel "x"
  set ylabel "z"
  set cbrange [0:5]
  set zrange [0:5]
  splot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:4
  unset title
  unset xlabel
  unset ylabel
  
  unset multiplot
  unset output

}



