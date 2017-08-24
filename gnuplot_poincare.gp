# This routine displays a particle orbit in the "data/poincare.dat" file.
# To use, run the program in the following way.
#   $ ./sample_poincare > data/poincare.dat
# Then, load this routine from the gnuplot.
#   $ gnuplot
#   gnuplot> load "gnuplot_poincare.gp"
#
# See also Appendix in Zenitani et al., Phys. Plasmas 20, 092120 (2013).

unset key

file="data/poincare.dat"

#set term pngcairo enhanced size 600, 500
set key left bottom
#set output "poi.png"
set size square
set pointsize 0.1

# X11
plot [-1:1][-1:1] file us 4:5 w points title "Zenitani+ 2013 Fig. 9"
# red (pngcairo)
#plot [-1:1][-1:1] file us 4:5 w points pt 7 lt 1
# blue (pngcairo)
#plot [-1:1][-1:1] file us 4:5 w points pt 7 lt 3
# gray (pngcairo)
#plot [-1:1][-1:1] file us 4:5 w points pt 7 lt 6

set key

# end
