set xrange [-1:17]
set xtics ("G" 0, "C" 1, "G" 2, "G" 3, "G" 4, "A" 5, "T" 6, "T" 7, "A" 8, "C" 9, "G" 10, "C" 11, "A" 12, "G" 13, "C" 14, "G" 15, "C" 16)
set xlabel "base"
set ylabel "roll [degrees]"
set key above
set font ",30"
set size ratio 0.75

set term pdfcairo enhanced font ",18"

set output "shear.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):2 w l lw 3 lc rgb "purple" title "Step 45"
set output


set ylabel "stretch [A]"
set output "stretch.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):3 w l lw 3 lc rgb "purple" title "Step 45"
set output

set yrange [-5:5]
set ylabel "buckle [degrees]"
set output "buckle.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):5 w l lw 3 lc rgb "purple" title "Step 45"
set output

set yrange [-25:0]
set ylabel "propeller [degrees]"
set output "propeller.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):6 w l lw 3 lc rgb "purple" title "Step 45"
set output

set yrange [-5:5]
set ylabel "opening [degrees]"
set output "opening.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):7 w l lw 3 lc rgb "purple" title "Step 45"
set output

set autoscale y
set ylabel "shift [A]"
set output "shift.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):8 w l lw 3 lc rgb "purple" title "Step 45"
set output

set ylabel "slide [A]"
set output "slide.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):9 w l lw 3 lc rgb "purple" title "Step 45"
set output

set ylabel "tilt [degrees]"
set output "tilt.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):11 w l lw 3 lc rgb "purple" title "Step 45"
set output

set ylabel "stagger [A]"
set output "stagger.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):4 w l lw 3 lc rgb "purple" title "Step 45"
set output

set ylabel "roll [degrees]"
set output "roll.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):12 w l lw 3 lc rgb "purple" title "Step 45"
set output

set ylabel "twist [degrees]"
set output "twist.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):13 w l lw 3 lc rgb "purple" title "Step 45"
set output

set ylabel "rise [A]"
set output "rise.pdf"
plot "Step0/Rep0/av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "black" title "Step 0", "Step5/Rep0/av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "blue" title "Step 5", "Step10/Rep0/av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "grey" title "Step 10", "Step15/Rep0/av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "red" title "Step 15", "Step20/Rep0/av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "gold" title "Step 20", "Step45/Rep0/av_int_coord_mapv4.txt" u (column(0)):10 w l lw 3 lc rgb "purple" title "Step 45"
set output
