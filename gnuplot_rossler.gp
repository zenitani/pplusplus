# This routine displays a particle orbit in the "data/rossler.dat" file.
# To use, run the program in the following way.
#   $ ./sample_rossler > data/rossler.dat
# Then, load this routine from the gnuplot.
#   $ gnuplot
#   gnuplot> load "gnuplot_rossler.gp"
# It works for gnuplot 5 or later.

unset key
set ticslevel 0
set view 60, 30

i=0
imax=10000
istep=10
interval=0.03
print "Press Ctrl-C to interrupt the animation."

file="data/rossler.dat"
set xrange [-10:15]
set yrange [-15:10]
set zrange [0:30]

load "gnuplot_animation.gp"

set key
# end
