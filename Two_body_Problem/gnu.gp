#set terminal gif
#set output 'phaseDiagram.gif'
set xlabel 'X'
set ylabel 'Y'
set title 'Phase Diagram'
set key left
plot 'positionData.dat' u 1:2 w line title 'phase of object1', 'positionData.dat' u 3:4 w line title 'phase of object2'

#set output 'energyDiagram.gif'
#set xlabel 'Iterations'
set ylabel 'Energy'
set title 'Energy Diagram'
set key left
plot 'energyData.dat' u 1:2 w line title 'Potential energy', 'energyData.dat' u 1:3 w line title 'Kinetic energy', 'energyData.dat' u 1:4 w line title 'Total energy'
