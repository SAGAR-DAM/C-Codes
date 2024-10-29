set title "Eigenvalues of the given matrix (roots of characteristic eq)\nx^{ 10 }+22x^{ 9 }+143x^{ 8 }-1852x^{ 7 }-13120x^{ 6 }+658278x^{ 5 }+7.57691e+06x^{ 4 }-2.15188e+07x^{ 3 }-8.4192e+08x^{ 2 }+5.80455e+08x+3.50507e+10 = 0"
set xlabel "Real"
set ylabel "Imaginary"
set grid
set xrange [-15.8944:11.1143]
set yrange [-14.764:14.764]
plot '-' with points pt 7 ps 1.5 lc rgb 'red' title 'Roots'
-7.6922 0
7.08803 2.48529
9.64939 8.31744
7.08803 -2.48529
-9.16129 8.48993
-7.73258 12.9673
-13.9949 0
-9.16129 -8.48993
-7.73258 -12.9673
9.64939 -8.31744
e
pause -1
