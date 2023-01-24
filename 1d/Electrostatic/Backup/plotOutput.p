
set key font ",12"

nx=32
nk=nx-1
N=128
nt=500

stats 'fields.txt' using 2
rho_min = STATS_min
rho_max = STATS_max 
stats 'fields.txt' using 3
phi_min = STATS_min
phi_max = STATS_max 
stats 'fields.txt' using 4
E_min = STATS_min
E_max = STATS_max 
stats 'fields.txt' using 5
ne_min = STATS_min
ne_max = STATS_max 
stats 'fields.txt' using 6
ni_min = STATS_min
ni_max = STATS_max 
stats 'fieldsk.txt' using 4
Ek_min = STATS_min
Ek_max = STATS_max 



set terminal gif animate delay 30 size 2560,1440
set output 'Frames/output_fields.gif'
do for [i=0:nt-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 3,2 columnsfirst
  
  set xlabel "x"
  set ylabel "RHO"
  set yrange [rho_min:rho_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "PHI"
  set yrange[phi_min:phi_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "E"
  set yrange[E_min:E_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  
  set xlabel "x"
  set ylabel "ne"
  set yrange[ne_min:ne_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "ni"
  set yrange[ni_min:ni_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:6 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "k"
  set ylabel "ESEk"
  set yrange[Ek_min:Ek_max]
  plot "fieldsk.txt" every ::i*nk::nk+(i*nk)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
}

unset output

set terminal gif animate delay 5 size 2560,1440
set output 'Frames/particles.gif'
do for [i=0:nt-1] {

  set title "Time Step = #".(i+1)
  set xlabel "x"
  set ylabel "particles"
  set yrange [0:1]
  set xrange [-1:101]
  plot "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with linespoint pointtype 7 pointsize 8.0 lc rgb "blue" title "IONS" ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with linespoint pointtype 7 pointsize 4.0 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  unset xrange
}
unset output


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
  unset yrange
  
  unset multiplot
  unset output

  filename="Frames/modes.png"
  set output filename
  
  set multiplot layout 7,1 
  
  set xlabel "t"
  set ylabel "ESE(k=1)"
  plot "modes.txt" using 1:2 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=8)"
  plot "modes.txt" using 1:3 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=16)"
  plot "modes.txt" using 1:4 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=20)"
  plot "modes.txt" using 1:5 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=25)"
  plot "modes.txt" using 1:6 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=27)"
  plot "modes.txt" using 1:7 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "t"
  set ylabel "ESE(k=30)"
  plot "modes.txt" using 1:8 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
  unset output



  filename="Frames/frequencies.png"
  set output filename
   
  set xlabel "k*dx"
  set ylabel "w(k)"
  set yrange [0:101]  
  plot "frequencies.txt" using 1:2 with points pointtype 7 lc rgb "blue" notitle ,\
  "frequencies.txt" using 1:3 with lines lc rgb "red" notitle
  unset xlabel
  unset ylabel
  unset yrange

  unset output
  
  
  reset