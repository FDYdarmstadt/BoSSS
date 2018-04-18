set output 'NodePerformance.tex'
set terminal cairolatex  pdf   size 17cm,17cm
set multiplot
set size 1,0.5
set size 0.98,0.49
set origin 0.01,0.505
unset lmargin 
unset rmargin 
set tmargin 4e00
set bmargin 1e00
unset logscale x
unset logscale y
unset logscale x2
unset logscale y2
set xrange [0:1000]
set autoscale y
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "Time [s]"
unset x2label
unset y2label
unset title 
unset key
set key outside right vertical maxrows 4 
set xtics 
set x2tics format " " 
set ytics 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
plot "NodePerformance_data_0.csv" title "Picard unknown MGLevels2" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 6 pointsize 0.5, "NodePerformance_data_1.csv" title "Picard SoftGMRES Swz w Coarse MGLevels2" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 6 pointsize 0.5, "NodePerformance_data_2.csv" title "NewtonGmres Swz w Coarse Overlap MGLevels2" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 4 pointsize 0.5, "NodePerformance_data_3.csv" title "NewtonGmres Swz Kcycle w Coarse Overlap MGLevels2" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 1 pointsize 0.5
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
set xrange [0:5000]
set autoscale y
set autoscale x2
set autoscale y2
set xlabel "Grid:NoOfCells"
set ylabel "Speedup"
unset x2label
unset y2label
unset title 
unset key
set key outside right vertical maxrows 5 
set xtics 
set x2tics format " " 
set ytics 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
plot "NodePerformance_data_4.csv" title "Picard unknown MGLevels2" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 6 pointsize 0.5, "NodePerformance_data_5.csv" title "Picard SoftGMRES Swz w Coarse MGLevels2" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 6 pointsize 0.5, "NodePerformance_data_6.csv" title "NewtonGmres Swz w Coarse Overlap MGLevels2" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 4 pointsize 0.5, "NodePerformance_data_7.csv" title "NewtonGmres Swz Kcycle w Coarse Overlap MGLevels2" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 1 pointsize 0.5, "NodePerformance_data_8.csv" title "Ideal" with lines linecolor  "black" dashtype 1 linewidth 3


exit
