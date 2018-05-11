set output 'Additional_1.tex'
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
set xrange [1:64]
set autoscale y
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "Time [s]"
unset x2label
set y2label "SolvIter: Precondition Solve"
set title "Exclusive times"
unset key
set key outside right vertical maxrows 2 
set xtics 
set x2tics format " " 
set ytics 
set y2tics format " " 
set termoption dashed
set termoption dashed
plot "Additional_1_data_0.csv" title "Automatic MGLevels3" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "Additional_1_data_1.csv" title "SoftGMRES Swz w Coarse Overlap MGLevels3" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 6 pointsize 0.5


exit
