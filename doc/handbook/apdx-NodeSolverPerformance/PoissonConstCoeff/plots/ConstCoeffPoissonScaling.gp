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
set logscale y 10
set logscale x2 10
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
set x2tics format "$10^{%L}$" 
set ytics format "$10^{%L}$" 
set y2tics format " " 
plot "ConstCoeffPoissonScaling_data_0.csv" title "Mult Gr w Blk Jac" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 5 pointsize 0.5, "ConstCoeffPoissonScaling_data_1.csv" title "Add Swz w Coarse" axes x2y1  with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 11 pointsize 0.5, "ConstCoeffPoissonScaling_data_2.csv" title "linear" axes x2y1  with lines linecolor  "black" dashtype 1 linewidth 1
set size 0.49,0.326666666666667
set origin 0.505,0.67
set lmargin 1e00
set rmargin 1e01
set tmargin 1e00
set bmargin 1e00
unset logscale x
unset logscale y
set logscale x2 10
set logscale y2 10
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
set x2tics format "$10^{%L}$" 
set ytics format " " 
set y2tics format "$10^{%L}$" 
plot "ConstCoeffPoissonScaling_data_3.csv" title "linear" axes x2y2  with lines linecolor  "black" dashtype 1 linewidth 1
set size 0.49,0.326666666666667
set origin 0.005,0.336666666666667
set lmargin 1e01
set rmargin 1e00
set tmargin 1e00
set bmargin 1e00
set logscale x 10
set logscale y 10
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
set ytics format "$10^{%L}$" 
set y2tics format " " 
plot "ConstCoeffPoissonScaling_data_4.csv" title "linear" with lines linecolor  "black" dashtype 1 linewidth 1
set size 0.49,0.326666666666667
set origin 0.505,0.336666666666667
set lmargin 1e00
set rmargin 1e01
set tmargin 1e00
set bmargin 1e00
set logscale x 10
unset logscale y
unset logscale x2
set logscale y2 10
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
set xtics format "$10^{%L}$" 
set x2tics format " " 
set ytics format " " 
set y2tics format "$10^{%L}$" 
plot "ConstCoeffPoissonScaling_data_5.csv" title "linear" axes x1y2  with lines linecolor  "black" dashtype 1 linewidth 1
set size 0.49,0.326666666666667
set origin 0.005,0.00333333333333333
set lmargin 1e01
set rmargin 1e00
set tmargin 1e00
set bmargin 1e00
set logscale x 10
set logscale y 10
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
set key at 1e12,1e04 vertical maxrows 3 
set xtics format "$10^{%L}$" 
set x2tics format " " 
set ytics format "$10^{%L}$" 
set y2tics format " " 
plot "ConstCoeffPoissonScaling_data_6.csv" title "linear" with lines linecolor  "black" dashtype 1 linewidth 1, "ConstCoeffPoissonScaling_data_7.csv" title "Mult Gr w Blk Jac" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 5 pointsize 0.5, "ConstCoeffPoissonScaling_data_8.csv" title "Add Swz w Coarse" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 11 pointsize 0.5


exit
