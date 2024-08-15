set terminal png enhanced;set output "snapshot3_1min.png";
set ylabel '';set format y '';set xtics nomirror;set ytics nomirror; set border 3
set style line 11 lc rgb '#808080' lt 1;set border 3 back ls 11;
plot 'snapshot3_1min.pdb.dat' u 1:(log($2)) t 'FoXS' w lines lw 2.5 lc rgb '#e26261'
reset
