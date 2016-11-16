#!/usr/bin/env gnuplot

# Loop over all frames.
do for [i=1:99] {
    reset

    # Colour definitions.
    set style line 1 lc rgb '#0060ad' lt 1 lw 3 pt 7 ps 2 # --- blue
    set style line 2 lc rgb '#ad0009' lt 1 lw 3 pt 7 ps 2 # --- red
    set style line 3 lc rgb '#000000' lt 1 lw 3 pt 7 ps 1 # --- black
    set style line 4 lc rgb '#ad0009' lt 1 lw 5 pt 7 ps 2 # --- red (thick)

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
    set xrange [5:195]
    set yrange [5:195]

    # Set up the boundary segment file names.
    filename1 = sprintf('boundary-segments_%04d.txt', i)

    # Plot the level set zero contour.
    plot filename1 w l ls 4

    # Set up right-hand plot pane.
    set xrange [0:3000]
    set yrange [0:0.08]
    set xtics  0,1000,3000
    set ytics  0,0.02,0.08
    set border
    set key left
    set xlabel 'Time'

    # Set up the time series file name.
    filename2 = sprintf('minimise_perimeter.txt')

    # Plot the time series data.
    plot filename2 every ::1::i using 1:3 w l ls 1 title 'velocity',  \
         filename2 every ::1::i using 1:4 w l ls 2 title 'curvature', \
         filename2 every ::1::i using 1:5 w l ls 3 title '1/R'

    unset multiplot
}
