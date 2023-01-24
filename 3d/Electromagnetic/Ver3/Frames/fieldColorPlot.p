reset

set pm3d map
#set palette color
set palette rgb 33,13,10


nx = 200
ny = 100
nz = 100


set dgrid3d nx,nz

nt=1000
tskip=20

set key font ",12"
set autoscale

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

reset 
exit