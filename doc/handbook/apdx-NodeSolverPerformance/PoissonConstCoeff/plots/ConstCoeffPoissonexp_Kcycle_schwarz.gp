set output 'ConstCoeffPoissonexp_Kcycle_schwarz.tex'
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
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [1e01:1e07]
set yrange [1e-02:1e04]
set autoscale x2
set autoscale y2
unset xlabel
set ylabel "$k = 2$"
unset x2label
unset y2label
unset title 
unset key
set key font ",16"inside top left Left reverse 
set xtics format " " 
set x2tics format " " 
set ytics format "$10^{%L}$" 
set ytics font "sans, 16" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonexp_Kcycle_schwarz_data_0.csv" title "Slv Iter" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_1.csv" title "Slv Init" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 4 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_2.csv" title "Agg Init" with linespoints linecolor  "black" dashtype 4 linewidth 3 pointtype 12 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_3.csv" title "Mtx ass" with linespoints linecolor  "black" dashtype 5 linewidth 3 pointtype 10 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_4.csv" title "linear" with lines linecolor  "black" dashtype 1 linewidth 1
set size 0.98,0.326666666666667
set origin 0.01,0.336666666666667
set lmargin 1e01
set rmargin 5e00
set tmargin 0e00
set bmargin 1e00
set logscale x 10
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [1e01:1e07]
set yrange [1e-02:1e04]
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
set ytics format "$10^{%L}$" 
set ytics font "sans, 16" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonexp_Kcycle_schwarz_data_5.csv" title "Slv Iter" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_6.csv" title "Slv Init" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 4 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_7.csv" title "Agg Init" with linespoints linecolor  "black" dashtype 4 linewidth 3 pointtype 12 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_8.csv" title "Mtx ass" with linespoints linecolor  "black" dashtype 5 linewidth 3 pointtype 10 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_9.csv" title "linear" with lines linecolor  "black" dashtype 1 linewidth 1
set size 0.98,0.326666666666667
set origin 0.01,0.00333333333333333
set lmargin 1e01
set rmargin 5e00
set tmargin 0e00
set bmargin 1e00
set logscale x 10
set logscale y 10
unset logscale x2
unset logscale y2
set xrange [1e01:1e07]
set yrange [1e-02:1e04]
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
set ytics format "$10^{%L}$" 
set ytics font "sans, 16" 
set y2tics format " " 
set termoption dashed
set termoption dashed
set termoption dashed
set termoption dashed
plot "ConstCoeffPoissonexp_Kcycle_schwarz_data_10.csv" title "Slv Iter" with linespoints linecolor  "black" dashtype 3 linewidth 3 pointtype 6 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_11.csv" title "Slv Init" with linespoints linecolor  "black" dashtype 2 linewidth 3 pointtype 4 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_12.csv" title "Agg Init" with linespoints linecolor  "black" dashtype 4 linewidth 3 pointtype 12 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_13.csv" title "Mtx ass" with linespoints linecolor  "black" dashtype 5 linewidth 3 pointtype 10 pointsize 0.5, "ConstCoeffPoissonexp_Kcycle_schwarz_data_14.csv" title "linear" with lines linecolor  "black" dashtype 1 linewidth 1


exit
