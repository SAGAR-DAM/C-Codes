set term wxt size 900,900
set view map
set size ratio -1
set palette defined (0 0 0 0, 0.1 0 0 0.5, 0.2 0 0 1, 0.3 0 0.5 1, 0.4 0 1 1, 0.5 0.5 1 0.5, 0.6 1 1 0, 0.7 1 0.5 0, 0.8 1 0 0, 0.9 0.5 0 0, 1 0.5 0 0)
set xlabel 'Re(z)'
set ylabel 'Im(z)'
set xrange [-1.52373:1.52373]
set yrange [-1.52373:1.52373]
set title 'Newton Fractal'
unset key
splot 'data.txt' using 1:2:3 with image, 'roots.txt' using 1:2:(0) with points pt 7 ps 2 lc rgb 'white'
