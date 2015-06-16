#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gommolecoord.h"
#  include "gomext/molecoord.h"
#endif

/***** PUBLIC GOMAPI BEGIN *****/

/** @defgroup gom_doc_atom Molecule Data
*** @ingroup  gom_doc
***/

/** @file
*** @ingroup  gom_doc_atom
*** @ingroup  gom_doc_atom_coord
***/

/** @brief Get the number of atoms in a molecule structure.
*** @ingroup gom_doc_atom
*** @param structure  Molecule structure index
***/
extern int gomp_GetNumAtomsInMolecStruct(
    GOM_ARG( int, structure ) );

/** @brief Get the total number of atoms in all molecule structures.
*** @ingroup gom_doc_atom
***/
extern int gomp_GetTotalNumberOfAtoms(void);


/** @defgroup gom_doc_atom_coord Atom Coordinates
*** @ingroup  gom_doc_atom
***/

/** @name Coordinate Functions
***/
/* @{ */
/** @brief Atom coordinate.
*** @ingroup  gom_doc_atom_coord
***/
extern       float  gomp_GetAtomXCoord(
    GOM_ARG( int  , structure ),
    GOM_ARG( int  , atom      ) );
extern       float  gomp_GetAtomYCoord(
    GOM_ARG( int  , structure ),
    GOM_ARG( int  , atom      ) );
extern       float  gomp_GetAtomZCoord(
    GOM_ARG( int  , structure ),
    GOM_ARG( int  , atom      ) );
extern       int    gomp_PutAtomXCoord(
    GOM_ARG( int  , structure ),
    GOM_ARG( float, x         ),
    GOM_ARG( int  , atom      ) );
extern       int    gomp_PutAtomYCoord(
    GOM_ARG( int  , structure ),
    GOM_ARG( float, y         ),
    GOM_ARG( int  , atom      ) );
extern       int    gomp_PutAtomZCoord(
    GOM_ARG( int  , structure ),
    GOM_ARG( float, z         ),
    GOM_ARG( int  , atom      ) );
extern const float *gomp_GetAtomXCoordPointer(
    GOM_ARG( int, structure ) );
extern const float *gomp_GetAtomYCoordPointer(
    GOM_ARG( int, structure ) );
extern const float *gomp_GetAtomZCoordPointer(
    GOM_ARG( int, structure ) );
extern       float *gomp_GetModifiableAtomXCoordPointer(
    GOM_ARG( int, structure ) );
extern       float *gomp_GetModifiableAtomYCoordPointer(
    GOM_ARG( int, structure ) );
extern       float *gomp_GetModifiableAtomZCoordPointer(
    GOM_ARG( int, structure ) );
extern int gomp_SaveAtomCoords(
    GOM_ARG( int, structure ) );
extern int gomp_GetSavedAtomCoords(
    GOM_ARG( int, structure ) );
extern int gomp_PushAtomCoordinates(void);
extern int gomp_PopAtomCoordinates(void);
/* @} */

/***** PUBLIC GOMAPI END *****/

#define IMPLEMENT_GET_MODIFIABLE_POINTER( \
    type,name,args,value,invalidate) \
    const type *gomp_Get##name args { return((type*)value); } \
    type       *gomp_##GetModifiable##name args { invalidate; return(value); }
