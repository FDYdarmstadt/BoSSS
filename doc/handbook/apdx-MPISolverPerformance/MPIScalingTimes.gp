set output 'MPIScalingTimes.tex'
set terminal cairolatex  pdf   size 17cm,17cm
set multiplot
set size 1,0.5
set size 0.98,0.49
set origin 0.01,0.505
unset lmargin 
unset rmargin 
unset tmargin 
set bmargin 1e00
unset logscale x
unset logscale y
unset logscale x2
unset logscale y2
set autoscale x
set autoscale y
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "Time [s]"
unset x2label
set y2label "SlvIter excl"
set title "Exclusive times"
unset key
set key at 2.5e00,1e-01 vertical maxrows 2 
set xtics 
set x2tics format " " 
set ytics 
set y2tics format " " 
plot "MPIScalingTimes_data_0.csv" title "Swz Kcycle w Coarse" with linespoints linecolor  5 dashtype 1 linewidth 1 pointtype 9 pointsize 0.5, "MPIScalingTimes_data_1.csv" title "Swz w Coarse" with linespoints linecolor  "blue" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5
set size 0.98,0.49
set origin 0.01,0.005
unset lmargin 
unset rmargin 
set tmargin 1e00
unset bmargin 
unset logscale x
unset logscale y
unset logscale x2
unset logscale y2
set autoscale x
set autoscale y
set autoscale x2
set autoscale y2
set xlabel "Processors"
set ylabel "Time [s]"
unset x2label
set y2label "SlvInit excl"
unset title 
set key off
set xtics 
set x2tics format " " 
set ytics 
set y2tics format " " 
plot "MPIScalingTimes_data_2.csv" title "Swz Kcycle w Coarse" with linespoints linecolor  5 dashtype 1 linewidth 1 pointtype 9 pointsize 0.5, "MPIScalingTimes_data_3.csv" title "Swz w Coarse" with linespoints linecolor  "blue" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5


exit
