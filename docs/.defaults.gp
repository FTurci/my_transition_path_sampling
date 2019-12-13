# #+results:
reset
    
set term pngcairo font 'Helvetica, 9' #size 700,340
set ytics nomirror offset 0.3,0 scale 0.5
set xtics nomirror offset 0,0.3 scale 0.5

# Put borders in "the background", http://www.gnuplotting.org/tag/terminal
set style line 11 lc rgb '#808080' lt 1
set border 3 back # ls 11

# Colors
red="#D87791"
red_dark="#B11B2B"
red_dark="#B11B2B"  # from gnuplotting
blue="#94B3DA"
blue="#0060ad"  # from gnuplotting
blue_dark="#3360B1"
gray="gray60"
black="gray20"
green_dark="#5e9c36"  # gnuplotting 
green="#60A36E"
green_dark="#237634"
orange="#E48100"
orange_dark="#D74700"
 
pscale = 1
set style line 1 pt 7  lt 1 ps pscale*1.1 lw 1.5 lc rgb blue
set style line 2 pt 5  lt 1 ps pscale*1   lw 1.5 lc rgb green_dark
set style line 3 pt 13 lt 1 ps pscale*1.2 lw 1.5 lc rgb red_dark
set style line 4 pt 6  lt 1 ps pscale*1.1 lw 1.5 lc rgb orange
set style line 5 pt 4  lt 1 ps pscale*1   lw 1.5 lc rgb black
set style line 6 pt 12 lt 1 ps pscale*1.2 lw 1.5 lc rgb gray
set style line 7 pt 12 lt 1 ps pscale*1.2 lw 1.5 lc "purple"

# Trick to get these line styles as the default styles
# https://stackoverflow.com/questions/24800244/predefined-point-types-in-gnuplot
set style increment user
