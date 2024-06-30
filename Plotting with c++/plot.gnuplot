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
set cblabel 'f(z)'
set xtics font ',10'
set ytics font ',10'
set xlabel 'x' font ',10'
set ylabel 'y' font ',10'
set cblabel 'z' font ',10'
set cbtics font ',10'
set xrange [-10:10]
set yrange [-10:10]
set title 'contour complex plot of f(z) = |ln(sin z)|; z=x+jy' font ',10'
unset key
splot 'data.txt' using 1:2:3 with image
