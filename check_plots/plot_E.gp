#!/usr/local/bin/gnuplot

fn = 'mg_hcp_C11z_E'

set terminal png font 'arial,14'
set output 'g_MgHCP_C11z.png'

set style line 1 lt 1 lw 4 lc rgb '#F42A2C'

set title "Structure energie for Mg HCP under C11z distortion"
set ylabel 'E [eV]'
set xlabel 'step number'

plot fn using 1:2 with lines ls 1 t ""
