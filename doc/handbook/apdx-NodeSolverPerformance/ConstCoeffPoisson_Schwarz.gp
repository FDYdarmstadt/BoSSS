set output 'ConstCoeffPoisson_Schwarz.tex'
set terminal cairolatex  pdf   size 16cm,12cm
set multiplot
set size 1,0.5
set size 0.98,0.49
set origin 0.01,0.505
unset lmargin 
unset rmargin 
unset tmargin 
unset bmargin 
set logscale x
set logscale y
unset logscale x2
unset logscale y2
set autoscale x
set yrange [1:10000]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 2$"
unset x2label
unset y2label
unset title 
unset key
set key outside right vertical maxrows 4 
set xtics format "$10^{%T}$" 
set x2tics format " " 
set ytics format "$10^{%T}$" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
set termoption dashed
plot "ConstCoeffPoisson_Schwarz_data_0.csv" title "Slv Iter" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoisson_Schwarz_data_1.csv" title "Slv Init" with linespoints linecolor  "black" dashtype 2 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoisson_Schwarz_data_2.csv" title "Agg Init" with linespoints linecolor  "black" dashtype 4 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoisson_Schwarz_data_3.csv" title "Mtx ass" with linespoints linecolor  "black" dashtype 5 linewidth 1 pointtype 6 pointsize 0.5
set size 0.98,0.49
set origin 0.01,0.005
unset lmargin 
unset rmargin 
unset tmargin 
unset bmargin 
set logscale x
set logscale y
unset logscale x2
unset logscale y2
set autoscale x
set yrange [1:10000]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 3$"
unset x2label
unset y2label
unset title 
unset key
set key outside right vertical maxrows 4 
set xtics format "$10^{%T}$" 
set x2tics format " " 
set ytics format "$10^{%T}$" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
set termoption dashed
plot "ConstCoeffPoisson_Schwarz_data_4.csv" title "Slv Iter" with linespoints linecolor  "black" dashtype 3 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoisson_Schwarz_data_5.csv" title "Slv Init" with linespoints linecolor  "black" dashtype 2 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoisson_Schwarz_data_6.csv" title "Agg Init" with linespoints linecolor  "black" dashtype 4 linewidth 1 pointtype 6 pointsize 0.5, "ConstCoeffPoisson_Schwarz_data_7.csv" title "Mtx ass" with linespoints linecolor  "black" dashtype 5 linewidth 1 pointtype 6 pointsize 0.5


exit
