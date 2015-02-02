/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "keys.h"
#include "gaussians.h"
#include "parser.h"

#ifndef impfchk_included

namespace Plugin {
namespace Fchk {
    
int ReadFCHK(int,const char**);
int FreeFCHKArrays(void);

} // namespace Fchk
} // namespace Plugin

#define impfchk_included 1
#endif
