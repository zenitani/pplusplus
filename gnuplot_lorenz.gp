# This routine displays a particle orbit in the "data/lorenz.dat" file.
# To use, run the program in the following way.
#   $ ./sample_lorenz > data/lorenz.dat
# Then, load this routine from the gnuplot.
#   $ gnuplot
#   gnuplot> load "gnuplot_lorenz.gp"
# It works for gnuplot 5 or later.

unset key
set ticslevel 0
set view 60, 30

i=0
imax=6000
istep=6
interval=0.03
print "Press Ctrl-C to interrupt the animation."

file="data/lorenz.dat"
set xrange [-20:20]
set yrange [-30:30]
set zrange [0:50]

load "gnuplot_animation.gp"

set key

# end
