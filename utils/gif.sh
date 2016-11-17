#!/usr/bin/env bash

# Create an animated gif from a bunch of PNG images.
# The first argument, $1, is the name of the output file.

# Convert the PNGs to GIFs.
printf "Converting PNG images ... "
for i in $(ls frame*.png)
do
    convert $i ${i%.png}.gif
done
printf "Done\n"

# Create the animated GIF.
printf "Creating animated GIF ... "
gifsicle --colors 16 --delay=1 --loop frame*.gif > $1.gif
printf "Done\n"

# Remove the PNG and GIF files.
printf "Cleaning files        ... "
rm -f frame*.{png,gif}
printf "Done\n"
