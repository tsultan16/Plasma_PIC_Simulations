reset
set pm3d map
set palette color
set palette rgb 33,13,10

nx = 100
ny = 200
nz = 1

pi = 3.14159265359

nt = 600
tskip= 10
tstart = tskip


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
nranks_x = 2
nranks_y = 1


ps = 0.05

set key font ",12"

E_max = 0.1

set terminal jpeg size 1920,1080
#set size ratio -1


############
# xy plane
############


do for [j=tstart/tskip:nt/tskip] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

 
  i = j*tskip
  set title "Time Step = #".j
  print "Time Steps Completed =".i
    

  set output "Images/particles_xy_t=".i.".jpg"

  filename_e1="Snapshots/electrons_xy_rank=00_t=".i.".dat"
  filename_e2="Snapshots/electrons_xy_rank=10_t=".i.".dat"  
  filename_e3="Snapshots/electrons_xy_rank=20_t=".i.".dat"
  filename_e4="Snapshots/electrons_xy_rank=30_t=".i.".dat"
  filename_e5="Snapshots/electrons_xy_rank=40_t=".i.".dat"
  filename_e6="Snapshots/electrons_xy_rank=50_t=".i.".dat"  
  filename_e7="Snapshots/electrons_xy_rank=60_t=".i.".dat"
  filename_e8="Snapshots/electrons_xy_rank=70_t=".i.".dat"
  
  filename_i1="Snapshots/ions_xy_rank=00_t=".i.".dat"
  filename_i2="Snapshots/ions_xy_rank=10_t=".i.".dat"   
  filename_i3="Snapshots/ions_xy_rank=20_t=".i.".dat"
  filename_i4="Snapshots/ions_xy_rank=30_t=".i.".dat"
  filename_i5="Snapshots/ions_xy_rank=40_t=".i.".dat"
  filename_i6="Snapshots/ions_xy_rank=50_t=".i.".dat"   
  filename_i7="Snapshots/ions_xy_rank=60_t=".i.".dat"
  filename_i8="Snapshots/ions_xy_rank=70_t=".i.".dat"
  
  filename_B1="Snapshots/Bvector_xy_rank=00_t=".i.".dat"
  filename_B2="Snapshots/Bvector_xy_rank=10_t=".i.".dat"  
  filename_B3="Snapshots/Bvector_xy_rank=20_t=".i.".dat"
  filename_B4="Snapshots/Bvector_xy_rank=30_t=".i.".dat"
  filename_B5="Snapshots/Bvector_xy_rank=40_t=".i.".dat"
  filename_B6="Snapshots/Bvector_xy_rank=50_t=".i.".dat"
  filename_B7="Snapshots/Bvector_xy_rank=60_t=".i.".dat"
  filename_B8="Snapshots/Bvector_xy_rank=70_t=".i.".dat"
  
  filename_Bz = "Snapshots/Bz_t=".i.".dat"
  
  
  
  set multiplot layout 1,2


  set xlabel "x" 
  set ylabel "y"
  set xrange [xmin:nranks_x*xmax]
  set yrange [ymin:nranks_y*ymax]
  set ytics 10.0 #lt 0 lw 1 lc rgb "#bbbbbb"
  set xtics 10.0 #lt 0 lw 1 lc rgb "#bbbbbb"
  set grid
  
  plot filename_B8 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B7 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B6 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B5 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B4 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B3 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B2 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  filename_B1 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2 notitle,\
  #filename_e1 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  #filename_e2 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  filename_e3 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  filename_e4 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  filename_e5 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  filename_e6 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  filename_e7 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  filename_e8 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "red" notitle,\
  #filename_i1 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  #filename_i2 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  filename_i3 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  filename_i4 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  filename_i5 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  filename_i6 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  filename_i7 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle,\
  filename_i8 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize ps lc rgb "blue" notitle

  unset xlabel
  unset ylabel
  unset xrange
  unset yrange
  
  
  set xlabel "x" 
  set ylabel "By"
  set xrange [xmin:nranks_x*xmax]
  set yrange [-9e-5:9e-5]
  
  plot filename_Bz binary format="%*1int%double%double%*1int" using 1:2 with linespoints

  unset xlabel
  unset ylabel
  unset xrange
  unset yrange  
  
  
  unset title
  unset multiplot
  unset output 
}
  

exit



  
############
# xz plane
############

#set terminal gif animate delay 10 size 1920,1080
set terminal jpeg size 1920,1080
#set output "ranks_xz.gif"
#set autoscale
set size ratio -1

do for [j=1:nt/tskip] {

  i = j
  print "Plotting particles. % Completed =".((100*i)/(nt/tskip))

 
  i = j*tskip
  set title "Time Step = #".j
  print "Time Steps Completed =".i
  
  
  set output "Images/particles_xz_t=".i.".jpg"
    
  filename_e1="Snapshots/electrons_xz_rank=00_t=".i.".dat"
  filename_e2="Snapshots/electrons_xz_rank=10_t=".i.".dat"  
  filename_e3="Snapshots/electrons_xz_rank=20_t=".i.".dat"
  filename_e4="Snapshots/electrons_xz_rank=30_t=".i.".dat"
  filename_i1="Snapshots/ions_xz_rank=00_t=".i.".dat"
  filename_i2="Snapshots/ions_xz_rank=10_t=".i.".dat"   
  filename_i3="Snapshots/ions_xz_rank=20_t=".i.".dat"
  filename_i4="Snapshots/ions_xz_rank=30_t=".i.".dat"
  
  filename_B1="Snapshots/Bvector_xz_rank=00_t=".i.".dat"
  filename_B2="Snapshots/Bvector_xz_rank=10_t=".i.".dat"  
  filename_B3="Snapshots/Bvector_xz_rank=20_t=".i.".dat"
  filename_B4="Snapshots/Bvector_xz_rank=30_t=".i.".dat"


  
  set xlabel "x" 
  set ylabel "z"
  set xrange [xmin:nranks_x*xmax]
  set yrange [zmin:zmax]
  set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
  set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"

  
  plot filename_B1 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2,\
  filename_B2 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2,\
  filename_B3 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2,\
  filename_B4 binary format="%*1int%double%double%double%double%*1int" using 1:2:($3/(E_max)):($4/(E_max)) with vectors head filled lt 2,\
  filename_e1 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "red" notitle,\
  filename_e2 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "red" notitle,\
  filename_e3 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "red" notitle,\
  filename_e4 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "red" notitle,\
  filename_i1 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "blue" notitle,\
  filename_i2 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "blue" notitle,\
  filename_i3 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "blue" notitle,\
  filename_i4 binary format="%*1int%double%double%*1int" using 1:2 with point pointtype 7 pointsize 0.35 lc rgb "blue" notitle


  unset xlabel
  unset ylabel
  unset xrange
  unset yrange
  unset title
  unset output 
  
}
  

#exit


