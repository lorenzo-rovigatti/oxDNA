set xrange [-1:150]
set xlabel "base [id]"
set key above
set font ",30"
set size ratio 0.75

set term pdfcairo enhanced font ",18"

set ylabel "shear [deg]"
set output "shear.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "stretch [deg]"
set output "stretch.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "stagger [deg]"
set output "stagger.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "buckle [deg]"
set output "buckle.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "black" title "oxDNA"
set output

set yrange [-20:0]
set ylabel "propeller [deg]"
set output "propeller.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "black" title "oxDNA"
set output

set autoscale y
set ylabel "opening [deg]"
set output "opening.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "shift [A]"
set output "shift.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "slide [A]"
set output "slide.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "rise [A]"
set output "rise.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "tilt [deg]"
set output "tilt.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "roll [deg]"
set output "roll.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "black" title "oxDNA"
set output

set ylabel "twist [deg]"
set output "twist.pdf"
plot "av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "black" title "oxDNA"
set output
