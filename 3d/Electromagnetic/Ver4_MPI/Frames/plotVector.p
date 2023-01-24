reset
set pm3d map
set palette color
set palette rgb 33,13,10

nx = 100
ny = 50
nz = 50

pi = 3.14159265359

nt=500
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


E_max = 10


############
# xy plane
############

set terminal gif animate delay 10 size 1920,1080
#set size square
set output 'current_xy.gif'

J_max = 0.01

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  
  filename="Snapshots/Jvector_xy_t=".i.".dat"
  set title "Time Step = #".j
  
  i = j

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(J_max)):($4/(J_max)) with vectors head filled lt 2 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  unset title
}

unset output






set terminal gif animate delay 10 size 1920,1080
set output 'Bvector_xz.gif'

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  
  filename="Snapshots/Bvector_xz_t=".i.".dat"


  set title "Time Step = #".j
  
  i = j

  set xlabel "x" 
  set ylabel "z"
  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  unset title
}

unset output

set terminal gif animate delay 10 size 1920,1080
set output 'Bvector_yz.gif'

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  
  filename="Snapshots/Bvector_yz_t=".i.".dat"


  set title "Time Step = #".j
  
  i = j

  set xlabel "y" 
  set ylabel "z"
  set xrange [ymin:ymax]
  set yrange [zmin:zmax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  unset title
}

unset output


#exit


exit


############
# yz plane
############


set terminal png size 1920,1080

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  
  filename="Snapshots/Bvector_yz_t=".i.".dat"
  filename_e="Snapshots/electrons_yz_t=".i.".dat"
  filename_i="Snapshots/ions_yz_t=".i.".dat"
  set output "Images/particles_yz_t=".i.".png"

  
  set title "Time Step = #".j
  
  i = j

  set xlabel "y" 
  set ylabel "z"
  set xrange [ymin:ymax]
  set yrange [zmin:zmax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 ,\
  filename_e binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.25 lc rgb "red" title "ELECTRONS" ,\
  filename_i binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  unset title
  unset output
}

exit



exit


