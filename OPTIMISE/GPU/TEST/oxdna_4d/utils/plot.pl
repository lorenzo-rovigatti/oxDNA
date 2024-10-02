set xrange [-1:17]
set xtics ("G" 0, "C" 1, "G" 2, "G" 3, "G" 4, "A" 5, "T" 6, "T" 7, "A" 8, "C" 9, "G" 10, "C" 11, "A" 12, "G" 13, "C" 14, "G" 15, "C" 16)
set xlabel "base"
set ylabel "translation [oxdna units]"

set term pdfcairo enhanced
set output "intra_trans.pdf"
plot "av_int_coord.txt" u (column(0)):2 w l lw 3 lc rgb "red" title "shear", "av_int_coord.txt" u (column(0)):3 w l lw 3 lc rgb "blue" title "stretch", "av_int_coord.txt" u (column(0)):4 w l lw 3 lc rgb "black" title "stagger"
set output

set output "intra_rot.pdf"
set ylabel "rotations [degrees]"
plot "av_int_coord.txt" u (column(0)):5 w l lw 3 lc rgb "red" title "buckle", "av_int_coord.txt" u (column(0)):6 w l lw 3 lc rgb "blue" title "propeller", "av_int_coord.txt" u (column(0)):7 w l lw 3 lc rgb "black" title "opening"
set output

set output "inter_trans.pdf"
set ylabel "translations [oxdna units]"
plot "av_int_coord.txt" u (column(0)):8 w l lw 3 lc rgb "red" title "shift", "av_int_coord.txt" u (column(0)):9 w l lw 3 lc rgb "blue" title "slide", "av_int_coord.txt" u (column(0)):10 w l lw 3 lc rgb "black" title "rise"
set output

set output "inter_rot.pdf"
set ylabel "rotations [degrees]"
plot "av_int_coord.txt" u (column(0)):11 w l lw 3 lc rgb "red" title "tilt", "av_int_coord.txt" u (column(0)):12 w l lw 3 lc rgb "blue" title "roll", "av_int_coord.txt" u (column(0)):13 w l lw 3 lc rgb "black" title "twist"
set output
