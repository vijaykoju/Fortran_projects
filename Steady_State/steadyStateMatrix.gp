unset key
set tic scale 0
set palette color
set xrange [0:767]
set yrange [0:767]
set xlabel 'Row'
set ylabel 'Column'
set xtics out
set ytics out
unset border
set view map
set title 'Iteration 1800'
#set output 'steadyState_0.gif'
splot 'steadyState_2210.dat' matrix with image
