set output 'Scaling_2.tex'
set terminal cairolatex  pdf   size 14cm,12cm
set multiplot
set size 1,0.5
set size 0.98,0.49
set origin 0.01,0.505
unset lmargin 
unset rmargin 
set tmargin 0e00
set bmargin 2e00
set logscale x 10
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [4e00:5.12e02]
set autoscale y
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "~wallclock time"
unset x2label
unset y2label
unset title 
unset key
set key font ",16"inside bottom right Left reverse 
set xtics format "$10^{%L}$" 
set xtics offset 0, 0-0.4 font "sans, 18" 
set xtics font "sans, 16" 
set x2tics format " " 
set ytics format "$10^{%L}$" 
set ytics font "sans, 16" 
set y2tics format " " 
set termoption dashed
plot "Scaling_2_data_0.csv" title "Kcycle w. add.-Schwarz DG2" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 9 pointsize 0.5, "Scaling_2_data_1.csv" title "2h limit" with lines linecolor  "black" dashtype 4 linewidth 2
set size 0.98,0.49
set origin 0.01,0.005
unset lmargin 
unset rmargin 
set tmargin 0e00
set bmargin 3e00
set logscale x 10
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [4e00:5.12e02]
set autoscale y
set autoscale x2
set autoscale y2
set xlabel "no of cores"
set ylabel "iterations"
unset x2label
unset y2label
unset title 
unset key
set key font ",16"inside bottom right Left reverse 
set xtics format "$10^{%L}$" 
set xtics offset 0, 0-0.4 font "sans, 18" 
set xtics font "sans, 16" 
set x2tics format " " 
set ytics format "$10^{%L}$" 
set ytics font "sans, 16" 
set y2tics format " " 
plot "Scaling_2_data_2.csv" title "Kcycle w. add.-Schwarz DG2" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 9 pointsize 0.5


exit
