set term post eps color enh solid
set size square
set key left
set yrange[*:*]
set xlabel 'Distance, Angstroms'
set ylabel 'Magnetization at 200ms'
set title '600 rate_matrix_auto.eps'
set output 'eps/rate_matrix_auto.eps'
plot 'eps/outy.out' with points pt 7
