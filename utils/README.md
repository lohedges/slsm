# Utils

A set of utility scripts for processing output data and making animations.
For example, to produce an animated GIF of the
[shape matching](../demos/shape_match.cpp) demo run the following in your
terminal:

```bash
./demos/shape_match
./utils/shape_match.gp
./utils/gif.sh shape_match
```
(This assumes that you are in the top level directory of the repository.)

You will now have an animated GIF called `shape_match.gif` in your directory.

The scripts require [Gnuplot](http://gnuplot.sourceforge.net),
[ImageMagick](https://www.imagemagick.org/script/index.php)
and [Gifsicle](https://www.lcdf.org/gifsicle).

To clean all of the output from a demo code from your working directory, run:

```bash
./utils/clean.sh
```
