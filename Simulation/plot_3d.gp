set title '3D Particle Trajectory'
set xlabel 'X Position (m)'
set ylabel 'Y Position (m)'
set zlabel 'Z Position (m)'
set grid
set terminal wxt size 800,600
splot '3d_particle_data.dat' using 1:2:3 with lines title 'Particle Path'
