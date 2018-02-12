set output 'ConstCoeffPoissonScaling.tex'
set terminal cairolatex  pdf   size 17cm,17cm
set multiplot
set size 0.5,0.333333333333333
set size 0.49,0.326666666666667
set origin 0.005,0.67
set lmargin 1e01
set rmargin 1e00
set tmargin 1e00
set bmargin 1e00
unset logscale x
set logscale y
set logscale x2
unset logscale y2
set xrange [100:10000000]
set yrange [0.001:10000]
set x2range [100:10000000]
set y2range [0.001:10000]
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
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_0.csv" title "Mult Gr w Blk Jac" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 5 pointsize 0.5, "ConstCoeffPoissonScaling_data_1.csv" title "Add Swz w Coarse" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_2.csv" title "Add Swz" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 9 pointsize 0.5, "ConstCoeffPoissonScaling_data_3.csv" title "Pardiso w Blk PC" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 4 pointsize 0.5, "ConstCoeffPoissonScaling_data_4.csv" title "Pardiso" axes x2y1  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_5.csv" title "CG" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5, "ConstCoeffPoissonScaling_data_6.csv" title "Mumps" axes x2y1  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.505,0.67
set lmargin 1e00
set rmargin 1e01
set tmargin 1e00
set bmargin 1e00
unset logscale x
unset logscale y
set logscale x2
set logscale y2
set xrange [100:10000000]
set yrange [0.001:10000]
set x2range [100:10000000]
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
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_7.csv" title "Mult Gr w Blk Jac" axes x2y2  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 5 pointsize 0.5, "ConstCoeffPoissonScaling_data_8.csv" title "Add Swz w Coarse" axes x2y2  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_9.csv" title "Pardiso" axes x2y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_10.csv" title "CG" axes x2y2  with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5, "ConstCoeffPoissonScaling_data_11.csv" title "Mumps" axes x2y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.005,0.336666666666667
set lmargin 1e01
set rmargin 1e00
set tmargin 1e00
set bmargin 1e00
set logscale x
set logscale y
unset logscale x2
unset logscale y2
set xrange [100:10000000]
set yrange [0.001:10000]
set x2range [100:10000000]
set y2range [0.001:10000]
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
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_12.csv" title "Mult Gr w Blk Jac" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 5 pointsize 0.5, "ConstCoeffPoissonScaling_data_13.csv" title "Pardiso" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_14.csv" title "CG" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5, "ConstCoeffPoissonScaling_data_15.csv" title "Mumps" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.505,0.336666666666667
set lmargin 1e00
set rmargin 1e01
set tmargin 1e00
set bmargin 1e00
set logscale x
unset logscale y
unset logscale x2
set logscale y2
set xrange [100:10000000]
set yrange [0.001:10000]
set x2range [100:10000000]
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
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_16.csv" title "Pardiso" axes x1y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_17.csv" title "Mumps" axes x1y2  with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5
set size 0.49,0.326666666666667
set origin 0.005,0.00333333333333333
set lmargin 1e01
set rmargin 1e00
set tmargin 1e00
set bmargin 1e00
set logscale x
set logscale y
unset logscale x2
unset logscale y2
set xrange [100:10000000]
set yrange [0.001:10000]
set x2range [100:10000000]
set y2range [0.001:10000]
unset xlabel
set ylabel "$k = 6$"
unset x2label
unset y2label
unset title 
unset key
set key at 1e12,1e04 vertical maxrows 7 
set xtics format "$10^{%T}$" 
set x2tics format " " 
set ytics format "$10^{%T}$" 
set y2tics format " " 
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonScaling_data_18.csv" title "Pardiso" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoissonScaling_data_19.csv" title "Mumps" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 10 pointsize 0.5, "ConstCoeffPoissonScaling_data_20.csv" title "Mult Gr w Blk Jac" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 5 pointsize 0.5, "ConstCoeffPoissonScaling_data_21.csv" title "Add Swz w Coarse" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_22.csv" title "Add Swz" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 9 pointsize 0.5, "ConstCoeffPoissonScaling_data_23.csv" title "Pardiso w Blk PC" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 4 pointsize 0.5, "ConstCoeffPoissonScaling_data_24.csv" title "CG" with linespoints linecolor  "black" dashtype 1 linewidth 1 pointtype 3 pointsize 0.5


exit
