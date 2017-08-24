# PNG movie
# set terminal pngcairo enhanced size 600, 500
# filename = sprintf ("movie/%04d.png", i)
# set output filename

# [-20:20][-30:30][0:50]

splot \
    file us 1:2:3 every ::0::i w l lt 4 lw 1 lc rgb "magenta" title "", \
    file us 1:2:3 every ::i::i     lt 1 pt 7 ps 4 lc rgb "red" title ""

if (i < imax) i=i+istep;pause interval;reread
exit
# end
