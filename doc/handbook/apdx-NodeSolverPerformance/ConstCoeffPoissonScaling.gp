set output 'ConstCoeffPoissonScaling.tex'
set terminal cairolatex  pdf   size 14cm,10.5cm
set multiplot
set size 0.5,0.333333333333333
set size 0.49,0.326666666666667
set origin 0.005,0.67
unset logscale x
set logscale y
set logscale x2
unset logscale y2
set autoscale x
set yrange [0.001:10000]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 2$"
unset x2label
unset y2label
unset title 
set key off
set xtics format " " 
set x2tics format "$10^{%T}$" 
set ytics format "$10^{%T}$" 
set y2tics format " " 
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 1
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_0.csv" title "Add Swz" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 9 pointsize 0.5, "ConstCoeffPoissonScaling_data_1.csv" title "Pardiso & Blk PC" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 4 pointsize 0.5, "ConstCoeffPoissonScaling_data_2.csv" title "Add Swz w Coarse" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_3.csv" title "Pardiso" axes x2y1  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_4.csv" title "CG" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5, "ConstCoeffPoissonScaling_data_5.csv" title "Mumps" axes x2y1  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.505,0.67
unset logscale x
unset logscale y
set logscale x2
set logscale y2
set autoscale x
set autoscale y
set autoscale x2
set y2range [0.001:10000]
unset xlabel
unset ylabel
unset x2label
set y2label "$k = 3$"
unset title 
set key off
set xtics format " " 
set x2tics format "$10^{%T}$" 
set ytics format " " 
set y2tics format "$10^{%T}$" 
set lmargin 1
set rmargin 10
set tmargin 1
set bmargin 1
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_6.csv" title "Add Swz w Coarse" axes x2y2  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_7.csv" title "Pardiso" axes x2y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_8.csv" title "CG" axes x2y2  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5, "ConstCoeffPoissonScaling_data_9.csv" title "Mumps" axes x2y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.005,0.336666666666667
set logscale x
set logscale y
unset logscale x2
unset logscale y2
set autoscale x
set yrange [0.001:10000]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 4$"
unset x2label
unset y2label
unset title 
set key off
set xtics format " " 
set x2tics format " " 
set ytics format "$10^{%T}$" 
set y2tics format " " 
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 1
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_10.csv" title "Add Swz w Coarse" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_11.csv" title "Pardiso" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_12.csv" title "CG" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5, "ConstCoeffPoissonScaling_data_13.csv" title "Mumps" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.505,0.336666666666667
set logscale x
unset logscale y
unset logscale x2
set logscale y2
set autoscale x
set autoscale y
set autoscale x2
set y2range [0.001:10000]
unset xlabel
unset ylabel
unset x2label
set y2label "$k = 5$"
unset title 
set key off
set xtics format "$10^{%T}$" 
set x2tics format " " 
set ytics format " " 
set y2tics format "$10^{%T}$" 
set lmargin 1
set rmargin 10
set tmargin 1
set bmargin 1
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_14.csv" title "Pardiso" axes x1y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_15.csv" title "Mumps" axes x1y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.005,0.00333333333333333
set logscale x
set logscale y
unset logscale x2
unset logscale y2
set autoscale x
set yrange [0.001:10000]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 6$"
unset x2label
unset y2label
unset title 
unset key
set key outside right vertical maxrows 2 
set xtics format "$10^{%T}$" 
set x2tics format " " 
set ytics format "$10^{%T}$" 
set y2tics format " " 
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 1
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_16.csv" title "Pardiso" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_17.csv" title "Mumps" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5


exit
