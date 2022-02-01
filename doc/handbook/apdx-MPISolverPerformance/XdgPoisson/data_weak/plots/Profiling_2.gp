set output 'Profiling_2.tex'
set terminal cairolatex  pdf   size 14cm,12cm
set multiplot
set size 1,1
set size 0.98,0.98
set origin 0.01,0.01
unset lmargin 
unset rmargin 
unset tmargin 
set bmargin 3e00
set logscale x 10
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [4e00:5.12e02]
set yrange [1e00:1e04]
set autoscale x2
set autoscale y2
set xlabel "no of cores"
set ylabel "runtime of proc0"
unset x2label
unset y2label
unset title 
unset key
set key font ",16"Left reverse 
set xtics format "$10^{%L}$" 
set xtics offset 0, 0-0.4 font "sans, 18" 
set xtics font "sans, 16" 
set x2tics format " " 
set ytics format "$10^{%L}$" 
set ytics font "sans, 16" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
set termoption dashed
plot "Profiling_2_data_0.csv" title "Slv Iter" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "Profiling_2_data_1.csv" title "Slv Init" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 4 pointsize 0.5, "Profiling_2_data_2.csv" title "Agg Init" with linespoints linecolor  "black" dashtype 4 linewidth 3 pointtype 12 pointsize 0.5, "Profiling_2_data_3.csv" title "Mtx ass" with linespoints linecolor  "black" dashtype 5 linewidth 3 pointtype 10 pointsize 0.5


exit
