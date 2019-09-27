set output 'XdgPoissonIterations.tex'
set terminal cairolatex  pdf   size 14cm,12cm
set multiplot
set size 1,0.333333333333333
set size 0.98,0.326666666666667
set origin 0.01,0.67
set lmargin 1e01
set rmargin 5e00
set tmargin 0e00
set bmargin 1e00
set logscale x 10
unset logscale y
unset logscale x2
unset logscale y2
set xrange [10:10000000]
set yrange [0:600]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 2$"
unset x2label
unset y2label
unset title 
unset key
set key font ",12"inside top left Left reverse 
set xtics format " " 
set x2tics format " " 
set ytics 
set y2tics format " " 
plot "XdgPoissonIterations_data_0.csv" title "Kcycle schwarz" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 11 pointsize 0.5, "XdgPoissonIterations_data_1.csv" title "gmres levelpmg" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 5 pointsize 0.5
set size 0.98,0.326666666666667
set origin 0.01,0.336666666666667
set lmargin 1e01
set rmargin 5e00
set tmargin 0e00
set bmargin 1e00
set logscale x 10
unset logscale y
unset logscale x2
unset logscale y2
set xrange [10:10000000]
set yrange [0:600]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 3$"
unset x2label
unset y2label
unset title 
set key off
set xtics format " " 
set x2tics format " " 
set ytics 
set y2tics format " " 
plot "XdgPoissonIterations_data_2.csv" title "Kcycle schwarz" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 11 pointsize 0.5, "XdgPoissonIterations_data_3.csv" title "gmres levelpmg" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 5 pointsize 0.5
set size 0.98,0.326666666666667
set origin 0.01,0.00333333333333333
set lmargin 1e01
set rmargin 5e00
set tmargin 0e00
set bmargin 1e00
set logscale x 10
unset logscale y
unset logscale x2
unset logscale y2
set xrange [10:10000000]
set yrange [0:600]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 5$"
unset x2label
unset y2label
unset title 
set key off
set xtics format "$10^{%L}$" 
set xtics offset 0, 0-0.4 font "sans, 18" 
set xtics font "sans, 16" 
set x2tics format " " 
set ytics 
set y2tics format " " 
plot "XdgPoissonIterations_data_4.csv" title "Kcycle schwarz" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 11 pointsize 0.5, "XdgPoissonIterations_data_5.csv" title "gmres levelpmg" with linespoints linecolor  "black" dashtype 1 linewidth 3 pointtype 5 pointsize 0.5


exit
