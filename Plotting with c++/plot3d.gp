# Gnuplot script to plot 3D line
set terminal pngcairo
set output '3d_plot.png'
set xlabel 'X Axis'
set ylabel 'Y Axis'
set zlabel 'Z Axis'
set view 60, 30
splot 'data.txt' using 1:2:3 with lines title '3D Line Plot'
pause -1
