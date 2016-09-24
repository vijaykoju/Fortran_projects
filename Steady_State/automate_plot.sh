#!/bin/bash
for FILE in *; do
    gnuplot <<- EOF
				set palette color
				set xrange [0:767]
				set yrange [0:767]
        set xlabel "Row"
        set ylabel "Column"
        set term png
				set view map
        set output "${FILE}.png"
        splot "${FILE}" matrix with image
    EOF
done
