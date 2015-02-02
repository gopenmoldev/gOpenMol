/*  

Copyright (c) 2004 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Coded by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gommath.h"

#include "stdafx.h"

#ifndef HAVE_NEARBYINT
/***********************************************************************/
double nearbyint(double x)
/***********************************************************************/
{
    return floor(x + 0.5);
}
/***********************************************************************/
#endif
