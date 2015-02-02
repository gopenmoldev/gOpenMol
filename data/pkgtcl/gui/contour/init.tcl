##############################################################################
#                       Copyright (c) 1994 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements 2003 - 2004 by: Eero Häkkinen
##############################################################################

package provide gom::gui::Contour 1.0
package require gom::PackageUtilities

gom::SourcePackageFiles gui/contour

#
# To achive fast display of contours the default display type is to save
# the data in memory and display the polygons from memory
# 2003-06-06/LUL
#
contour method save
