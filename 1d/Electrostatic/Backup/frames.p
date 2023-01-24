
set key font ",12"

nx=32
nk=nx/2 -1
N=128
nt=100
tskip=1

option = 1 

#set ytics font "Verdana,6" 


do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/Fields/fields_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 3,1 
  
  set xlabel "x"
  set ylabel "RHO"
  #set yrange [-0.0009:0.0009]
  plot "fields.txt" every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "PHI"
  set yrange [-0.09:0.09]
  plot "fields.txt" every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "E"
  set yrange [-0.006:0.006]
  plot "fields.txt" every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
  unset output
    
  filename2="Frames/Fourier/fourier_t=".i.".png"
  set output filename2
  
  set title "Time Step = #".(i+1)
  set multiplot layout 3,1 
  
  set xlabel "k"
  set ylabel "RHOk"
  #set yrange [-0.0009:0.0009]
  plot "fieldsk.txt" every ::i*nk::nk+(i*nk)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "k"
  set ylabel "PHIk"
  #set yrange [-0.09:0.09]
  plot "fieldsk.txt" every ::i*nk::nk+(i*nk)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
 
  set xlabel "k"
  set ylabel "ESEk"
  #set yrange [-0.006:0.006]
  plot "fieldsk.txt" every ::i*nk::nk+(i*nk)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  unset multiplot
  unset output
   
  filename3="Frames/Particles/particles_t=".i.".png"
  set output filename3
  
  set title "Time Step = #".(i+1)
  #plot 1: fluid density
  #set label 1 'RHO' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "particles"
  set yrange [0:1]
  set xrange [-1:101]
  plot "particles.txt" every ::i*N::N+(i*N)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  unset xrange
  unset output
  
  
}

  set terminal png size 2560,1440

  filename="Frames/energy_momentum.png"
  set output filename
  
  set multiplot layout 4,1 
  
  set xlabel "t"
  set ylabel "KE"
  #set yrange [-0.0009:0.0009]
  plot "energy.txt" using 1:2 with lines lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "t"
  set ylabel "ESE"
  #set yrange [-0.09:0.09]
  plot "energy.txt" using 1:3 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
 
  set xlabel "t"
  set ylabel "KE+ESE"
  #set yrange [-0.006:0.006]
  plot "energy.txt" using 1:4 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "t"
  set ylabel "Px"
  set yrange [-0.005:0.005]
  plot "energy.txt" using 1:5 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  unset multiplot
  unset output

  reset