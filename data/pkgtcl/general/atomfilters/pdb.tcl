##############################################################################
#                        Copyright (c) 2002 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
#                      Confidential unpublished property of
#                              Leif Laaksonen
#                            All rights reserved
##############################################################################

namespace eval gom::FileFilters::AtomCoordinates::Pdb {

;
proc PostRead { fileName append errorValue } {
    if { 0 == $errorValue } { find ssbonds [show molstruct] }
}

}; # end of namespace gom::FileFilters::AtomCoordinates::Pdb
