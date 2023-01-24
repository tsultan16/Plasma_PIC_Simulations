
set key font ",12"

nx=32
nk=nx-1
N=1024*4
vbins=1000

nt=4000/10
tskip = 1


stats 'energy.txt' using 2
KE_min = STATS_min
KE_max = STATS_max 


  set terminal png size 2560,1440
  filename="Frames/energy_momentum.png"
  set output filename
  
  set multiplot layout 4,2 columnsfirst 
  
  set xlabel "t"
  set ylabel "KE"
  set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:2 with lines lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE"
  set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:3 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "t"
  set ylabel "KE+ESE"
  set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:4 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
    
  set xlabel "t"
  set ylabel "Px"
  set yrange [-0.005:0.005]
  plot "energy.txt" using 1:5 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "t"
  set ylabel "vd1"
  #set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:6 with lines lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "t"
  set ylabel "vspread1"
  #set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:7 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
 
  set xlabel "t"
  set ylabel "vd2"
  #set yrange [0:1.1*KE_max]
  plot "energy.txt" using 1:8 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
    
  set xlabel "t"
  set ylabel "vspread2"
  #set yrange [-0.005:0.005]
  plot "energy.txt" using 1:9 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  unset multiplot
  unset output


  filename="Frames/ESE_log.png"
  set output filename
   
  set xlabel "t"
  set ylabel "LN(ESEk(k=1))"
  set grid ytics lc rgb "#bbbbbb" lw 1 lt 1
  set grid xtics lc rgb "#bbbbbb" lw 1 lt 1
  set grid mxtics mytics lw 1 lt 0
  plot "energy.txt" every ::nt/40::nt/15 using 1:(log($3)) with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  unset grid

 unset output
  
 
  filename="Frames/modes_log.png"
  set output filename
  set multiplot layout 3,1 
  
  set xlabel "t"
  set ylabel "LN(ESEk(k=1))"
  set grid ytics lc rgb "#bbbbbb" lw 1 lt 1
  set grid xtics lc rgb "#bbbbbb" lw 1 lt 1
  set grid mxtics mytics lw 1 lt 0
  plot "modes.txt" every ::nt/40::nt/15 using 1:(log($2)) with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  unset grid

  set xlabel "t"
  set ylabel "LN(ESEk(k=2))"
  set grid ytics lc rgb "#bbbbbb" lw 1 lt 1
  set grid xtics lc rgb "#bbbbbb" lw 1 lt 1
  set grid mxtics mytics lw 1 lt 0
  plot "modes.txt" every ::nt/40::nt/15 using 1:(log($3)) with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  unset grid


  set xlabel "t"
  set ylabel "LN(ESEk(k=4))"
  set grid ytics lc rgb "#bbbbbb" lw 1 lt 1
  set grid xtics lc rgb "#bbbbbb" lw 1 lt 1
  set grid mxtics mytics lw 1 lt 0
  plot "modes.txt" every ::nt/20::nt/15 using 1:(log($4)) with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  unset grid


  unset multiplot
  unset output

  exit
  ###################################################################################################

  filename="Frames/modes.png"
  set output filename
  
  set multiplot layout 7,1 
  
  set xlabel "t"
  set ylabel "ESE(k=1)"
  set logscale y
  plot "modes.txt" using 1:2 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  
  set xlabel "t"
  set ylabel "ESE(k=2)"
  plot "modes.txt" using 1:3 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=4)"
  plot "modes.txt" using 1:4 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=6)"
  plot "modes.txt" using 1:5 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=8)"
  plot "modes.txt" using 1:6 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "t"
  set ylabel "ESE(k=10)"
  plot "modes.txt" using 1:7 with lines lc rgb "blue" notitle
  set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]  
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "t"
  set ylabel "ESE(k=12)"
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
  set yrange [0:3]  
  plot "frequencies.txt" using 1:2 with points pointtype 7 lc rgb "blue" notitle ,\
  "frequencies.txt" using 1:3 with lines lc rgb "red" notitle
  unset xlabel
  unset ylabel
  unset yrange

  unset output
  
   
set terminal gif animate delay 20 size 2560,1440
set output 'Frames/fv.gif'
do for [j=0:nt/tskip-1] {

  i = j*tskip
  print "Plotting fv. % Completed =".((100*i)/nt)

  set title "Time Step = #".(i+1)
  set xlabel "vx"
  set ylabel "f(vx)"
  #set xrange [-0.1:6.3]
  #set yrange [-0.2:0.2]
  plot "fv.txt" every ::i*vbins::vbins+(i*vbins)-1 using 1:($2+$3) with boxes fill solid 0.2 lc rgb "blue" title "Species 1+ Species2"
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange
  
}
unset output

  
  
set terminal gif animate delay 20 size 2560,1440
set output 'Frames/particles.gif'
do for [j=0:nt/tskip - 1] {

  i = j*tskip
  print "Plotting particles. % Completed =".((100*i)/nt)


  set multiplot layout 2,1 columnsfirst

  set title "Time Step = #".(i+1)
  set xlabel "x"
  set ylabel "y"
  #set xrange [-0.1:6.3]
  set yrange [-0.2:0.2]
  plot "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 4.0 lc rgb "blue" title "IONS" ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:2 with point pointtype 7 pointsize 4.0 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  #unset xrange

  set xlabel "x"
  set ylabel "vx"
  #set xrange [-0.1:6.3]
  #set yrange [-0.2:0.2]
  plot "particles_i.txt" every ::i*N::N+(i*N)-1 using 1:3 with point pointtype 7 pointsize 0.5 lc rgb "blue" title "IONS" ,\
  "particles_e.txt" every ::i*N::N+(i*N)-1 using 1:3 with point pointtype 7 pointsize 0.5 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange
  
  unset multiplot
  
}
unset output


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
do for [j=0:nt/tskip-1] {

  i = j*tskip
  print "Plotting fields. % Completed =".((100*i)/nt)

  
  set title "Time Step = #".(i+1)
  set multiplot layout 3,2 columnsfirst
  
  set xlabel "x"
  set ylabel "RHO"
  #set yrange [rho_min:rho_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "x"
  set ylabel "PHI"
  #set yrange[phi_min:phi_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
 
  set xlabel "x"
  set ylabel "E"
  #set yrange[E_min:E_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  
  set xlabel "x"
  set ylabel "ne"
  #set yrange[ne_min:ne_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "x"
  set ylabel "ni"
  #set yrange[ni_min:ni_max]
  plot "fields.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:6 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  set xlabel "k"
  set ylabel "ESEk"
  #set yrange[Ek_min:Ek_max]
  plot "fieldsk.txt" every ::i*nk::nk+(i*nk)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #unset yrange
  
  unset multiplot
}

unset output

exit

set terminal png size 2560,1440
set output 'Frames/phase_space_i.png'

  set title "Phase Space"

  set multiplot layout 3,1 

  set xlabel "x"
  set ylabel "vx"
  #set yrange [0:1]
  #set xrange [-1:101]
  plot "particles_i.txt"every N using 1:3 with linespoint pointtype 7 pointsize 2.0 lc rgb "blue" title "IONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange
  
  
  set xlabel "vx"
  set ylabel "vy"
  #set yrange [0:1]
  #set xrange [-1:101]
  plot "particles_i.txt" every N using  3:4 with linespoint pointtype 7 pointsize 2.0 lc rgb "blue" title "IONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange

  set xlabel "x"
  set ylabel "y"
  #set yrange [0:1]
  #set xrange [-1:101]
  plot "particles_i.txt"every N using 1:2 with linespoint pointtype 7 pointsize 2.0 lc rgb "blue" title "IONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange

  unset multiplot
  
unset output

set terminal png size 2560,1440
set output 'Frames/phase_space_e.png'

  set title "Phase Space"

  set multiplot layout 3,1 

  set xlabel "x"
  set ylabel "vx"
  #set yrange [0:1]
  #set xrange [-1:101]
  plot "particles_e.txt" every N using  1:3 with linespoint pointtype 7 pointsize 2.0 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange

  set xlabel "vx"
  set ylabel "vy"
  #set yrange [0:1]
  #set xrange [-1:101]
  plot "particles_e.txt" every N using  3:4 with linespoint pointtype 7 pointsize 2.0 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange
  
  set xlabel "x"
  set ylabel "y"
  #set yrange [0:1]
  #set xrange [-1:101]
  plot "particles_e.txt" every N using  1:2 with linespoint pointtype 7 pointsize 2.0 lc rgb "red" title "ELECTRONS"  
  unset title
  unset xlabel
  unset ylabel
  #unset yrange
  #unset xrange
  
  unset multiplot

unset output



set terminal gif animate delay 5 size 2560,1440
set output 'Frames/single_particle.gif'
do for [i=0:nt-1] {

  set title "Time Step = #".(i+1)
  set xlabel "x"
  set ylabel "y"
  set xrange [2.9:3.4]
  set yrange [-0.25:0.45]
  plot "single_particle.txt" every ::0::i using 1:2 with linespoint pointtype 7 pointsize 2.0 lc rgb "red" title "ELECTRON" ,\
  "single_particle.txt" every ::0::i using 3:4 with linespoint pointtype 7 pointsize 2.0 lc rgb "blue" title "ION" 
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  unset xrange
}
unset output

  
reset