#!/usr/bin/env gnuplot

# Loop over all frames.
do for [i=1:1000:4] {
    reset

    # Colour definitions.
    set style line 1 lc rgb '#0060ad' lt 1 lw 5 pt 7 ps 2 # --- blue
    set style line 2 lc rgb '#ad0009' lt 1 lw 5 pt 7 ps 2 # --- red
    set style line 3 lc rgb '#000000' lt 1 lw 3 pt 7 ps 2 # --- black

    # Set the terminal style.
    set terminal pngcairo size 800, 400 enhanced font 'Verdana, 14'

    # Set the PNG file name.
    set output sprintf('frame%04.0f.png', i)

    # Two column plot.
    set multiplot layout 1,2

    # Set up left-hand plot pane.
    unset key
    set size ratio 0.86
    set xrange [-30:30]
    set yrange [-25:25]
    set xtics  -30,10,30
    set ytics  -20,10,20
    set xlabel 'X'
    set ylabel 'Y'

    # Set up the boundary segment file names.
    filename1 = sprintf('boundary-segments_%04d.txt', 0)
    filename2 = sprintf('boundary-segments_%04d.txt', i)

    # Plot the level set zero contour.
    plot filename1 using ($1-35):($2-30) w l ls 1, \
         filename2 using ($1-35):($2-30) w l ls 2, \

    # Set up right-hand plot pane.
    set xrange [0:75]
    set yrange [-20:20]
    set xtics  0,25,75
    set ytics  -20,10,20
    set border
    set xlabel 'Time (x 1e4)'
    set ylabel 'X Position'

    # Set up the time series file name.
    filename = sprintf('bimodal_0.0500.txt')

    # Plot the time series data.
    plot filename every ::1::i using ($1/10000):($4) w l ls 3

    unset multiplot
}
