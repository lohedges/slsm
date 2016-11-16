#!/usr/bin/env gnuplot

# Loop over all frames.
do for [i=1:100] {
    reset

    # Colour definitions.
    set style line 1 lc rgb '#ad0009' lt 1 lw 3 pt 7 ps 2 # --- red
    set style line 2 lc rgb '#0060ad' lt 1 lw 3 pt 7 ps 2 # --- blue

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
    plot filename1 w l ls 1

    # Set up right-hand plot pane.
    set xrange  [0:100]
    set yrange  [100:500]
    set y2range [0.5:1.0]
    set xtics   0,20,100
    set ytics   100,100,500
    set y2tics  0.5,0.1,1.0
    set ytics   nomirror
    set border
    set xlabel  'Time'
    set ylabel  'Perimeter' textcolor rgb '#0060ad'
    set y2label 'Area'      textcolor rgb '#ad0009'

    # Set up the time series file name.
    filename2 = sprintf('perimeter_0.0500.txt')

    # Plot the time series data.
    plot filename2 every ::1::i using 1:2 w l ls 2 axes x1y1, \
         filename2 every ::1::i using 1:3 w l ls 1 axes x1y2

    unset multiplot
}
