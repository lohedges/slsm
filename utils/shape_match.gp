#!/usr/bin/env gnuplot

# Loop over all frames.
do for [i=1:100] {
    reset

    # Colour definitions.
    set style line 1 lc rgb '#ad0009' lt 1 lw 5 pt 7 ps 2 # --- red
    set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 7 ps 2 # --- blue

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
    set xrange [5:395]
    set yrange [5:395]

    # Set up the boundary segment file names.
    filename1 = sprintf('boundary-segments_%04d.txt', i)

    # Plot the level set zero contour.
    plot filename1 w l ls 1

    # Set up right-hand plot pane.
    set xrange [0:100]
    set yrange [0:0.3]
    set xtics  0,20,100
    set ytics  0,0.1,0.3
    set border
    set xlabel 'Time'
    set ylabel 'Mismatch'

    # Set up the time series file name.
    filename2 = sprintf('shape-match_0.0000.txt')

    # Plot the time series data.
    plot filename2 every ::1::i using 1:2 w l ls 2

    unset multiplot
}
