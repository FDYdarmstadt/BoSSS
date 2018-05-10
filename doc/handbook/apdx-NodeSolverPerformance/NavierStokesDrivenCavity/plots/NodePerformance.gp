set output 'NodePerformance.tex'
set terminal cairolatex  pdf   size 18cm,10cm
set multiplot
set size 1,1
set size 0.98,0.98
set origin 0.01,0.01
unset lmargin 
unset rmargin 
set tmargin 4e00
set bmargin 4e00
set logscale x 10
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [1000:1000000]
set autoscale y
set autoscale x2
set autoscale y2
set xlabel "TotalDoFs"
set ylabel "Time [s]"
unset x2label
unset y2label
unset title 
unset key
set key outside right vertical maxrows 4 
set xtics format "$10^{%L}$" 
set x2tics format " " 
set ytics format "$10^{%L}$" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
plot "NodePerformance_data_0.csv" title "NewtonGmres Automatic" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "NodePerformance_data_1.csv" title "NewtonGmres Swz Kcycle w Coarse Overlap" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 1 pointsize 0.5, "NodePerformance_data_2.csv" title "Picard Mumps" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 10 pointsize 0.5, "NodePerformance_data_3.csv" title "linear" with lines linecolor  "black" dashtype 1 linewidth 1


exit
