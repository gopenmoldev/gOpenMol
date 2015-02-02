#
# small utility routine to plot a truncated octahedron in
# an arbitrary place in space, and another to place one at (0,0,0)
# given the "side-to-side" cell dimension.
#

proc lulPlotOctaBox {xmin xmax ymin ymax zmin zmax {colour}} {

set x_centre [expr ($xmin + $xmax)/2.0];
set y_centre [expr ($ymin + $ymax)/2.0];
set z_centre [expr ($zmin + $zmax)/2.0];

set dim [expr ($xmax - $xmin)/4.0];
set x1 [expr $x_centre + $dim];
set _x1 [expr $x_centre - $dim];
set y1 [expr $y_centre + $dim];
set _y1 [expr $y_centre - $dim];
set z1 [expr $z_centre + $dim];
set _z1 [expr $z_centre - $dim];

# "top square"
plot line $_x1 $y_centre $zmin $x_centre $y1 $zmin "$colour"
plot line $x_centre $y1 $zmin $x1 $y_centre $zmin "$colour" append;
plot line $x1 $y_centre $zmin $x_centre $_y1 $zmin "$colour" append;
plot line $x_centre $_y1 $zmin $_x1 $y_centre $zmin "$colour" append;
# "bottom square"
plot line $_x1 $y_centre $zmax $x_centre $y1 $zmax "$colour" append;
plot line $x_centre $y1 $zmax $x1 $y_centre $zmax "$colour" append;
plot line $x1 $y_centre $zmax $x_centre $_y1 $zmax "$colour" append;
plot line $x_centre $_y1 $zmax $_x1 $y_centre $zmax "$colour" append;
# "front square"
plot line $xmin $y_centre $_z1 $xmin $y1 $z_centre "$colour" append;
plot line $xmin $y1 $z_centre $xmin $y_centre $z1 "$colour" append;
plot line $xmin $y_centre $z1 $xmin $_y1 $z_centre "$colour" append;
plot line $xmin $_y1 $z_centre $xmin $y_centre $_z1 "$colour" append;
# "back square"
plot line $xmax $y_centre $_z1 $xmax $y1 $z_centre "$colour" append;
plot line $xmax $y1 $z_centre $xmax $y_centre $z1 "$colour" append;
plot line $xmax $y_centre $z1 $xmax $_y1 $z_centre "$colour" append;
plot line $xmax $_y1 $z_centre $xmax $y_centre $_z1 "$colour" append;
# "left square"
plot line $x_centre $ymin $_z1 $_x1 $ymin $z_centre "$colour" append;
plot line $_x1 $ymin $z_centre $x_centre $ymin $z1 "$colour" append;
plot line $x_centre $ymin $z1 $x1 $ymin $z_centre "$colour" append;
plot line $x1 $ymin $z_centre $x_centre $ymin $_z1 "$colour" append;
# "right square"
plot line $x_centre $ymax $_z1 $_x1 $ymax $z_centre "$colour" append;
plot line $_x1 $ymax $z_centre $x_centre $ymax $z1 "$colour" append;
plot line $x_centre $ymax $z1 $x1 $ymax $z_centre "$colour" append;
plot line $x1 $ymax $z_centre $x_centre $ymax $_z1 "$colour" append;
# "bits in between":
# "top 4 joints"
plot line $x1 $y_centre $zmin $xmax $y_centre $_z1 "$colour" append;
plot line $_x1 $y_centre $zmin $xmin $y_centre $_z1 "$colour" append;
plot line $x_centre $y1 $zmin $x_centre $ymax $_z1 "$colour" append;
plot line $x_centre $_y1 $zmin $x_centre $ymin $_z1 "$colour" append;
# "middle 4 joints"
plot line $xmin $y1 $z_centre $_x1 $ymax $z_centre "$colour" append;
plot line $x1 $ymax $z_centre $xmax $y1 $z_centre "$colour" append;
plot line $xmin $_y1 $z_centre $_x1 $ymin $z_centre "$colour" append;
plot line $x1 $ymin $z_centre $xmax $_y1 $z_centre "$colour" append;
# "bottom 4 joints"
plot line $x1 $y_centre $zmax $xmax $y_centre $z1 "$colour" append;
plot line $_x1 $y_centre $zmax $xmin $y_centre $z1 "$colour" append;
plot line $x_centre $y1 $zmax $x_centre $ymax $z1 "$colour" append;
plot line $x_centre $_y1 $zmax $x_centre $ymin $z1 "$colour" append;

}

proc lulPlotOcta {side {colour}} {

set dim [expr $side / 2.0];
lulPlotOctaBox -$dim $dim -$dim $dim -$dim $dim "$colour";

}

