# This routine displays a particle orbit in the "data/ExB.dat" file.
# To use, run the program in the following way.
#   $ ./sample_ExB > data/ExB.dat
# Then, load this routine from the gnuplot.
#   $ gnuplot
#   gnuplot> load "gnuplot_ExB.gp"

unset key

file="data/ExB.dat"
plot [0:25][0:0.4] file us 1:2 w linesp

set key

# end
