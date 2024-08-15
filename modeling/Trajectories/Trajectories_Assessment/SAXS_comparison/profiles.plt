set terminal png enhanced;set output "profiles.png"
set ylabel '';set format y '';set xtics nomirror;set ytics nomirror; set border 3
set style line 11 lc rgb '#808080' lt 1;set border 3 back ls 11;
plot 'snapshot1_1min.pdb.dat' u 1:(log($2)) t "snapshot1_1min" w lines lw 2 lt 2,'snapshot2_1min.pdb.dat' u 1:(log($2)) t "snapshot2_1min" w lines lw 2 lt 3,'snapshot3_1min.pdb.dat' u 1:(log($2)) t "snapshot3_1min" w lines lw 2 lt 4
reset
