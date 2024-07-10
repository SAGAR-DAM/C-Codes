set term wxt size 1000,1000
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
set xlabel 'Re(z)'
set ylabel 'Im(z)'
set cblabel 'Root Index'
set xtics font ',10'
set ytics font ',10'
set xlabel 'Re(z)' font ',10'
set ylabel 'Im(z)' font ',10'
set cblabel 'Root Index' font ',10'
set cbtics font ',10'
set xrange [-1.1:1.1]
set yrange [-1.1:1.1]
set title 'Newton Fractal' font ',10'
unset key
splot 'data.txt' using 1:2:3 with image
