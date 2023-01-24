reset
set pm3d map
set palette color
set palette rgb 33,13,10

nx = 50
ny = 50
nz = 50

pi = 3.14159265359

nt=465
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


E_max = 1e-4


############
# xz plane
############


set terminal gif animate delay 10 size 1920,1080
set output "Images/particles_xz.gif"


do for [j=1:nt/tskip] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  
  i = j*tskip
  set title "Time Step = #".j

  print "Time Steps Completed =".i
  
  filename="Snapshots/Bvector_xz_t=".i.".dat"
  filename_e="Snapshots/electrons_xz_t=".i.".dat"
  filename_i="Snapshots/ions_xz_t=".i.".dat"  
  #set output "Images/particles_xz_t=".i.".png"


  i = j

  set xlabel "x" 
  set ylabel "z"
  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 #,\
  #filename_e binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 1.5 lc rgb "red" title "ELECTRONS" ,\
  #filename_i binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 2.5 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  unset title
  #unset output



 # i = j*tskip
 # set title "Time Step = #".j

  #print "Time Steps Completed =".i
  
  #filename="Snapshots/Bvector_xy_t=".i.".dat"
  #filename_e="Snapshots/electrons_xy_t=".i.".dat"
  #filename_i="Snapshots/ions_xy_t=".i.".dat"  
  #set output "Images/particles_xy_t=".i.".png"


  
  #i = j

  #set xlabel "x" 
  #set ylabel "y"
  #set xrange [xmin:xmax]
  #set yrange [ymin:ymax]
  
  #plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 ,\
  #filename_e binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 1.5 lc rgb "red" title "ELECTRONS" ,\
  #filename_i binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 2.5 lc rgb "blue" title "IONS" 
 
  #unset xlabel
  #unset ylabel
  #unset xrange
  #unset yrange

  unset title
  #unset output
  

}

unset output

exit

############
# xy plane
############

set terminal gif animate delay 10 size 1920,1080
#set size square
set output 'current_xy.gif'

J_max = 1e-3

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
set output 'Bvector_xy.gif'

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  
  filename="Snapshots/Bvector_xy_t=".i.".dat"


  set title "Time Step = #".j
  
  i = j

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  unset title
}

unset output



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





set terminal png size 1920,1080
#set size square
set output 'particles_xz.gif'


do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  
  set output "Images/particles_xz_t=".i.".png"
  filename="Snapshots2/Bvector_xz_t=".i.".dat"
  filename_e="Snapshots2/electrons_t=".i.".dat"
  filename_i="Snapshots2/ions_t=".i.".dat"
  set title "Time Step = #".j
  
  
  i = j

  set xlabel "x" 
  set ylabel "z"
  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 title "B-field",\
  filename_e binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.25 lc rgb "red" title "ELECTRONS" ,\
  filename_i binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  unset title
  unset output
}

unset output

exit




set dgrid3d nx,nz
set terminal png size 1920,1080

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting B field. % Completed =".((100*i)/(nt/tskip))
  
  
  i = j*tskip
  print "Time Steps Completed =".i

  set output "Images/Bfield_xz_t=".i.".png"

  set title "rho, ne, ni, |B field|, Time Step = #".i

  set multiplot layout 2,2 columnsfirst 

  filename1="Snapshots/density_xz_t=".i.".dat"
  filename2="Snapshots/B_xz_t=".i.".dat"

  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  set cbrange [-1e-6:1e-6]
  set zrange [-1e-6:1e-6]
  
  set xlabel "x"
  set ylabel "y"
  splot filename1 binary format="%*1int%double%double%double%double%double%*1int" using 1:2:4
  unset title
  unset xlabel
  unset ylabel
  unset cbrange
  unset zrange
  
  set cbrange [0:20]
  set zrange [0:20]
  
  set xlabel "x"
  set ylabel "y"
  splot filename1 binary format="%*1int%double%double%double%double%double%*1int" using 1:2:5
  unset title
  unset xlabel
  unset ylabel
  unset cbrange
  unset zrange
  
  
  set cbrange [0:20]
  set zrange [0:20]
  
  set xlabel "x"
  set ylabel "y"
  splot filename1 binary format="%*1int%double%double%double%double%double%*1int" using 1:2:3
  unset title
  unset xlabel
  unset ylabel
  
  
  unset cbrange
  unset zrange
  
  #set cbrange [0:3000]
  #set zrange [0:3000]
  
  set xlabel "x"
  set ylabel "y"
  splot filename2 binary format="%*1int%double%double%double%*1int" using 1:2:3
  unset title
  unset xlabel
  unset ylabel
  
  unset cbrange
  unset zrange
  unset xrange
  unset yrange
  
  
  unset multiplot
  unset output

}



reset
exit


set terminal png size 1920,1080
#set size square


do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  \
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/vector_xy_t=".i.".dat"
  filename2="Fields/fields_vercut_t=".i.".dat"
  set output "Fields/vector_xy_t=".i.".png"
  set title "Time Step = #".j
  
  i = j

  set multiplot layout 2,1 columnsfirst 

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 ,\
 "particles_e.dat" binary format="%*1int%double%double%*1int" every ::0*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS" ,\
 "particles_i.dat" binary format="%*1int%double%double%*1int" every ::0*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  set xlabel "z" 
  set ylabel "Bz"
  set xrange [zmin:zmax]
  plot filename2 binary format="%*1int%double%double%*1int" using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  unset title
  unset multiplot
  unset output
}




set terminal gif animate delay 10 size 1920,1080
#set size square
set output 'particles_xz.gif'


do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

  
  set title "Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i
  filename="Fields/vector_xy_t=".i.".dat"
  set title "Time Step = #".j
  
  i = j

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [zmin:zmax]

   plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 ,\
  "particles_e.dat" binary format="%*1int%double%double%*1int" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS" ,\
  "particles_i.dat" binary format="%*1int%double%double%*1int" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 1.0 lc rgb "blue" title "IONS" 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange

  
  unset title
}

unset output


set terminal gif animate delay 10 size 1920,1080
#set size square
set output 'current_xy.gif'
E_max = 1000
J_max = 1000

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))
  print "Time Steps Completed =".i
  
  set title "Time Step = #".(i+1)
  
  i = j*tskip

  filename="Fields/Jvector_xy_t=".i.".dat"
  filename2="Fields/Bvector_xz_t=".i.".dat"

  set title "Time Step = #".j
  
  i = j

  set multiplot layout 2,1 columnsfirst 

  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  
  plot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(J_max)):($4/(J_max)) with vectors head filled lt 2
  unset ylabel
  unset xrange
  unset yrange


  set xlabel "x" 
  set ylabel "z"
  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  plot filename2 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 
 
  unset xlabel
  unset ylabel
  unset xrange
  unset yrange


  unset multiplot
  unset title
}

unset output


exit


set dgrid3d nx,nz
set terminal png size 1920,1080

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting B field. % Completed =".((100*i)/(nt/tskip))
  
  set title "|B field|, Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i

  set output "Fields/Bfield_xz_t=".i.".png"

  set title "ne, ni, |B field|, Time Step = #".i

  set multiplot layout 3,1 columnsfirst 

  filename="Fields/fields_xz_t=".i.".dat"

  set xrange [xmin:xmax]
  set yrange [zmin:zmax]
  
  set cbrange [0:200]
  set zrange [0:200]
  
  set xlabel "x"
  set ylabel "y"
  splot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:3
  unset title
  unset xlabel
  unset ylabel

  set cbrange [0:200]
  set zrange [0:200]
  
  set xlabel "x"
  set ylabel "y"
  splot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:4
  unset title
  unset xlabel
  unset ylabel
  
  filename="Fields/B_xz_t=".i.".dat"

  set cbrange [-5600:5600]
  set zrange [-5600:5600]
  
  set xlabel "x"
  set ylabel "y"
  splot filename binary format="%*1int%double%double%double%double%*1int" using 1:2:4
  unset title
  unset xlabel
  unset ylabel
  
  unset cbrange
  unset zrange
  
  unset xrange
  unset yrange
  unset multiplot
  unset output

}

reset
exit

set dgrid3d nx,nz
set terminal png size 1920,1080

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting charge density. % Completed =".((100*i)/(nt/tskip))
    
  i = j*tskip
  print "Time Steps Completed =".i

  filename="Fields/fields_xz_t=".i.".dat"
  set output "Fields/rho_xz_t=".i.".png"

  set title "RHO, Time Step = #".i
  set xlabel "x"
  set ylabel "y"
  set cbrange [-1:1]
  set zrange [-1:1]
  splot filename binary format="%*1int%double%double%double%*1int" using 1:2:3
  unset title
  unset xlabel
  unset ylabel

  unset output

}


set dgrid3d nx,nz
set terminal png size 1920,1080

do for [j=0:nt/tskip-1] {

  i = j
  print "Plotting B field. % Completed =".((100*i)/(nt/tskip))
  
  set title "|B field|, Time Step = #".(i+1)
  
  i = j*tskip
  print "Time Steps Completed =".i

  filename="Fields/fields_xz_t=".i.".dat"
  set output "Fields/Bfield_xz_t=".i.".png"

  set title "|B field|, Time Step = #".i
  set xlabel "x"
  set ylabel "y"
  splot filename binary format="%*1int%double%double%double%*1int" using 1:2:3
  unset title
  unset xlabel
  unset ylabel

  unset output

}

