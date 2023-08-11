set key above
set font ",30"
set size ratio 0.75
set xlabel "base"
set term pdfcairo enhanced font ",18"


#set xrange [-1:17]
#set xtics ("G" 0, "C" 1, "G" 2, "G" 3, "G" 4, "A" 5, "T" 6, "T" 7, "A" 8, "C" 9, "G" 10, "C" 11, "A" 12, "G" 13, "C" 14, "G" 15, "C" 16)

set ylabel "shear [degrees]"
set output "shear".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):2 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "stretch [A]"
set output "stretch".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):3 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "stagger [A]"
set output "stagger".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):4 w l lw 3 lc rgb colors[i] title titles[i]
set output

set yrange [-5:5]
set ylabel "buckle [degrees]"
set output "buckle".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):5 w l lw 3 lc rgb colors[i] title titles[i]
set output

set yrange [-25:0]
set ylabel "propeller [degrees]"
set output "propeller".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):6 w l lw 3 lc rgb colors[i] title titles[i]
set output

set yrange [-5:5]
set ylabel "opening [degrees]"
set output "opening".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):7 w l lw 3 lc rgb colors[i] title titles[i]
set output

set autoscale y
set ylabel "shift [A]"
set output "shift".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):8 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "slide [A]"
set output "slide".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):9 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "rise [A]"
set output "rise".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):10 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "tilt [degrees]"
set output "tilt".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):11 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "roll [degrees]"
set output "roll".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):12 w l lw 3 lc rgb colors[i] title titles[i]
set output

set ylabel "twist [degrees]"
set output "twist".output_file
plot for[i=1:Nfiles] files[i] u (column(0)):13 w l lw 3 lc rgb colors[i] title titles[i]
set output


