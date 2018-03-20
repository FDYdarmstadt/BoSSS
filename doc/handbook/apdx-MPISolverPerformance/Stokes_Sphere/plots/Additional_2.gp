set output 'Additional_2.tex'
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
set xrange [1:16]
set autoscale y
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "Speedup"
unset x2label
set y2label "SolIter: Precondition Solve"
set title "Speedup"
unset key
set key outside right vertical maxrows 3 
set xtics 
set x2tics format " " 
set ytics 
set y2tics format " " 
set termoption dashed
set termoption dashed
plot "Additional_2_data_0.csv" title "Automatic MGLevels3" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "Additional_2_data_1.csv" title "SoftGMRES Swz w Coarse Overlap MGLevels3" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 6 pointsize 0.5, "Additional_2_data_2.csv" title "Ideal" with lines linecolor  "black" dashtype 1 linewidth 3


exit
