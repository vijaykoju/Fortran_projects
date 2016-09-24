filename = "steadyState_".i.".dat"
plotfile = "steadyState_".i.".gif"
iteration = "Iteration ".j.""
print filename." ".plotfile

unset key
set tic scale 0
set palette color
set xrange [0:767]
set yrange [0:767]
set xlabel 'Row'
set ylabel 'Column'
set xtics out
set ytics out
set view map
set title iteration
set output plotfile
splot filename matrix with image

j = j+3
i = i+3
if (i <= n) reread
