##########################################################################
# FlatBlack.tcl
#
#NOTE:  comment lines start with a # character
# simple demo tcl script 
#
    global env
#

#  I like the orthographic, rather than perspective view
define projection orthographic

#white backgroud
define bgcolour white
#improve display quality for publication
define cylinderquality 100
define spherequality 100

#make atoms black
atom  color * * * Black

# CPK display
atom cpk
atom scale cpk 0.1 * * *

#make matte finish
define material specular red 0.0 
define material specular blue 0.0 
define material specular green 0.0

