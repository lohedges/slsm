#!/usr/bin/env gnuplot

# Loop over all frames.
do for [i=1:200] {
    reset

    # Colour definitions.
    set style line 1 lc rgb '#0060ad' lt 1 lw 5 pt 7 ps 2 # --- blue
    set style line 2 lc rgb '#ad0009' lt 1 lw 5 pt 7 ps 2 # --- red
    set style line 3 lc rgb '#000000' lt 1 lw 5 pt 7 ps 2 # --- black

    # Set the terminal style.
    set terminal pngcairo size 800, 400 enhanced font 'Verdana, 14'

    # Set the PNG file name.
    set output sprintf('frame%03.0f.png', i)

    # Two column plot.
    set multiplot layout 1,2

    # Set up left-hand plot pane.
    unset key
    unset tics
    set border 0
    set size square
    set xrange [5:95]
    set yrange [5:95]

    # Set up the boundary segment file names.
    filename1 = sprintf('boundary-segments_%04d.txt', 0)
    filename2 = sprintf('boundary-segments_%04d.txt', i)

    # Plot the level set zero contour.
    plot filename1 w l ls 1, \
         filename2 w l ls 2, \

    # Set up right-hand plot pane.
    set xrange [0:8000]
    set yrange [50:80]
    set xtics  0,2000,8000
    set ytics  50,10,80
    set border
    set xlabel 'Time'
    set ylabel 'Objective'

    # Set up the time series file name.
    filename = sprintf('dumbbell_0.0020.txt')

    # Plot the time series data.
    plot filename every ::1::i using 1:2 w l ls 3

    unset multiplot
}
