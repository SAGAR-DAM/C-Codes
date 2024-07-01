set term wxt
set view map
set size ratio -1
set palette defined ( \
0 0 0 0, \
0.1 0 0 0.5, \
0.2 0 0 1, \
0.3 0 0.5 1, \
0.4 0 1 1, \
0.5 0.5 1 0.5, \
0.6 1 1 0, \
0.7 1 0.5 0, \
0.8 1 0 0, \
0.9 0.5 0 0, \
1 0.5 0 0 )
set xlabel 'x'
set ylabel 'y'
set cblabel 'z'
set xtics font ',10'
set ytics font ',10'
set xlabel 'x' font ',10'
set ylabel 'y' font ',10'
set cblabel 'z' font ',10'
set cbtics font ',10'
set xrange [-20:20]
set yrange [-20:20]
set title 'contour plot' font ',10'
unset key
splot 'data.txt' using 1:2:3 with image
