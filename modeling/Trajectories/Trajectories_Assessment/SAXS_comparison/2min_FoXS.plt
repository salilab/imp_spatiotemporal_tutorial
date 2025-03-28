set terminal png enhanced; set output "snapshot1_2min_2min_exp.png";
set lmargin 7; set rmargin 2;set multiplot
set origin 0,0;set size 1,0.3; set tmargin 0; set bmargin 3;set ylabel 'Residual';set format y '';set xtics nomirror;set xlabel 'q [Å^{-1}]';set ytics nomirror; set border 3
set style line 11 lc rgb '#808080' lt 1;set border 3 back ls 11;f(x)=1
plot f(x) notitle lc rgb '#333333', 'snapshot1_2min_2min_exp.fit' u 1:(($2-$4)/$3) notitle w lines lw 2.5 lc rgb '#e26261'
set origin 0,0.3;set size 1,0.69; set bmargin 0; set tmargin 1;set xlabel ''; set format x ''; set ylabel 'I(q) log-scale';set format y "10^{%L}"; set logscale y
plot 'snapshot1_2min_2min_exp.fit' u 1:2 t 'Experimental' lc rgb '#333333' pt 6 ps 0.8, 'snapshot1_2min_2min_exp.fit' u 1:4 t 'FoXS χ^2 = 5.49' w lines lw 2.5 lc rgb '#e26261'
unset multiplot
reset
set terminal postscript eps size 3.5,2.62 color enhanced solid linewidth 2.5 font 'Helvetica,22'; set output "snapshot1_2min_2min_exp.eps";
set lmargin 2; set rmargin 2;set multiplot
set origin 0,0;set size 1,0.3; set tmargin 0; set bmargin 3;set ylabel '';set format y '';
set xtics nomirror font 'Helvetica,18';
set ytics nomirror; set border 3
set style line 11 lc rgb '#808080' lt 1;
set border 3 back ls 11;f(x)=1
plot f(x) notitle lc rgb '#333333', 'snapshot1_2min_2min_exp.fit' u 1:(($2-$4)/$3) notitle w lines lw 2.5 lc rgb '#e26261'
set origin 0,0.3;set size 1,0.69; set bmargin 0; set tmargin 1;set xlabel ''; set format x ''; set ylabel 'I(q) log-scale';
set format y "10^{%L}"; set logscale y
plot 'snapshot1_2min_2min_exp.fit' u 1:2 t 'Experimental' lc rgb '#333333' pt 6 ps 0.8, 'snapshot1_2min_2min_exp.fit' u 1:4 t 'FoXS {/Symbol c}^2 = 5.49' w lines lw 2.5 lc rgb '#e26261'
unset multiplot
reset
