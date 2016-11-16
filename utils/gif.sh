#!/usr/bin/env bash

# Create an animated gif from a bunch of PNG images.

# Convert the PNGs to GIFs.
echo "Converting PNG images ..."
for i in $(ls frame*.png)                                                                                                             ⏎ master ◼jjjj
do
    convert $i ${i%.png}.gif
done

# Create the animated GIF.
echo "Creatine animated GIF ..."
gifsicle --colors 256 --delay=1 --loop frame*.gif > $1.gif

# Remove the PNG and GIF files.
echo "Cleaning files ..."
rm -f frame*.{png,gif}

echo "Done!"
