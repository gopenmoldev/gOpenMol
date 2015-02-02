/***** PUBLIC GOMAPI BEGIN *****/
#include "gomlistener.h"
/****** PUBLIC GOMAPI END ******/
#include "parser_types.h"

/* handles cases when there is no information available */
#define DEFAULT_SEGMENT_NAME   "omol"
#define DEFAULT_RESIDUE_NAME   "omol"

/* max length for segment, residue and atom names */
#define MAX_SEG_NAME_LEN 4 /* max four characters */
#define MAX_RES_NAME_LEN 4 /* max four characters */
#define MAX_ATM_NAME_LEN 4 /* max four characters */

#define COVdelt   0.4f   /* tolerance value when calculating connectivity */

/* conversion constant to convert between atomic units to Angstrom */
#define ATM2ANG   0.52917715f

#define BOHR_RADIUS 0.52917715  /* conversion constant */

/* atom location and name info structure  */

typedef struct {
    int     type;
    float   bndrad;
    float   vdwrad;
    float   plurad;
    char    global;
    float   emin;
    float   rmin;
    float   patom;
    char    hbond;
    char    atype[MAX_ATM_NAME_LEN];
    float   mass;
    int     ncharge;
    int     cnct;
} AtomData;

typedef struct {
    int params;              /* terms in the parameters list */
    AtomData *AtomParams;
} AtomTypes_t;

extern AtomTypes_t gomp_AtomTypes;

extern int SetCHARMStructureOrder(int         , int         , int);
/* reserve space for n numer of atoms
   into the Molec structure */
extern int gomp_GetSpaceForMolec(const char *,int);

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gommolecule.h"
#  include "gomext/molecule.h"
#else

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc_atom
***/

/** @file
*** @ingroup gom_doc_atom
*** @ingroup gom_doc_atom_id
*** @ingroup gom_doc_atom_param
*** @ingroup gom_doc_atom_level
*** @ingroup gom_doc_atom_disp
*** @ingroup gom_doc_atom_conn
*** @ingroup gom_doc_listener_atom_data
*** @ingroup gom_doc_listener_atom_coord
***/

/****************************/
/*                          */
/* Atom Identification Data */
/*                          */
/****************************/

/** @defgroup gom_doc_atom_id Atom Identification Data
*** @ingroup gom_doc_atom
***/

/** @name Identification Data Functions
***/
/** @brief Segment name.
*** @ingroup gom_doc_atom_id
***/
/* @{ */
extern const char  *gomp_GetAtomSegName(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomSegName(
    GOM_ARG( int         , structure   ),
    GOM_ARG( const char *, segmentName ),
    GOM_ARG( int         , atom        ) );
/* @} */
/** @brief Residue name.
*** @ingroup gom_doc_atom_id
***/
/* @{ */
extern const char  *gomp_GetAtomResName(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomResName(
    GOM_ARG( int         , structure      ),
    GOM_ARG( const char *, residueName    ),
    GOM_ARG( int         , atom           ) );
/* @} */
/** @brief Atom name.
*** @ingroup gom_doc_atom_id
***/
/* @{ */
extern const char  *gomp_GetAtomAtmName(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomAtmName(
    GOM_ARG( int         , structure      ),
    GOM_ARG( const char *, atomName       ),
    GOM_ARG( int         , atom           ) );
/* @} */
/** @brief Residue number 1.
*** @ingroup gom_doc_atom_id
***/
/* @{ */
extern       int    gomp_GetAtomResNum1(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomResNum1(
    GOM_ARG( int         , structure      ),
    GOM_ARG( int         , residueNumber1 ),
    GOM_ARG( int         , atom           ) );
extern const int   *gomp_GetAtomResNum1Pointer(
    GOM_ARG( int         , structure ) );
extern       int   *gomp_GetModifiableAtomResNum1Pointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Residue number 2.
*** @ingroup gom_doc_atom_id
***/
/* @{ */
extern       int    gomp_GetAtomResNum2(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomResNum2(
    GOM_ARG( int         , structure      ),
    GOM_ARG( int         , residueNumber2 ),
    GOM_ARG( int         , atom           ) );
extern const int   *gomp_GetAtomResNum2Pointer(
    GOM_ARG( int         , structure ) );
extern       int   *gomp_GetModifiableAtomResNum2Pointer(
    GOM_ARG( int         , structure ) );
/* @} */

/*******************/
/*                 */
/* Atom Parameters */
/*                 */
/*******************/

/** @defgroup gom_doc_atom_param Atom Parameters
*** @ingroup gom_doc_atom
***/

/** @name Parameter Functions
***/
/** @brief Atom partial charge.
*** @ingroup gom_doc_atom_param
***/
/* @{ */
extern       float  gomp_GetAtomCharge(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomCharge(
    GOM_ARG( int         , structure      ),
    GOM_ARG( float       , partialCharge  ),
    GOM_ARG( int         , atom           ) );
extern const float *gomp_GetAtomChargePointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomChargePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom nuclear charge.
*** @ingroup gom_doc_atom_param
***/
/* @{ */
extern       float  gomp_GetAtomNucCharge(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomNucCharge(
    GOM_ARG( int         , structure      ),
    GOM_ARG( float       , nuclearCharge  ),
    GOM_ARG( int         , atom           ) );
extern const float *gomp_GetAtomNucChargePointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomNucChargePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom covalent radius.
*** @ingroup gom_doc_atom_param
***/
/* @{ */
extern       float  gomp_GetAtomCovar(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomCovar(
    GOM_ARG( int         , structure      ),
    GOM_ARG( float       , covalentRadius ),
    GOM_ARG( int         , atom           ) );
extern const float *gomp_GetAtomCovarPointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomCovarPointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom B value.
*** @ingroup gom_doc_atom_param
***/
/* @{ */
extern       float  gomp_GetAtomBValue(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomBValue(
    GOM_ARG( int         , structure      ),
    GOM_ARG( float       , bValue         ),
    GOM_ARG( int         , atom           ) );
extern const float *gomp_GetAtomBValuePointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomBValuePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom Basis set tag.
*** @ingroup gom_doc_atom_param
***/
/* @{ */
extern const char  *gomp_GetAtomBasisSetTag(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomBasisSetTag(
    GOM_ARG( int         , structure      ),
    GOM_ARG( const char *, basisSetTag    ),
    GOM_ARG( int         , atom           ) );
/* @} */

/**************************/
/*                        */
/* Atom Level Information */
/*                        */
/**************************/

/** @defgroup gom_doc_atom_level Atom Level Information
*** @ingroup gom_doc_atom
***/

/** @name Parameter Functions
***/
/** @brief Atom type number (type).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       int    gomp_GetAtomType(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomType(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , type      ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Bonding radius (bndrad).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomBndRad(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomBndRad(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , bndrad    ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief van der Waals spheres radius (vdwrad).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomVdwRad(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomVdwRad(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , vdwrad    ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Sphere radius (plurad) in plots (not used in gOpenMol).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomPluRad(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomPluRad(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , plurad    ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Value (global) used in global search for bonds.
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       char   gomp_GetAtomGlobal(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomGlobal(
    GOM_ARG( int         , structure ),
    GOM_ARG( char        , global    ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Value (emin) used in the calculation of van der Waals and electrostatic energy.
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomEmin(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomEmin(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , emin      ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Value (rmin) used in the calculation of van der Waals and electrostatic energy.
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomRmin(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomRmin(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , rmin      ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Atom polarizabilities (patom).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomPatom(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomPatom(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , patom     ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Atom hydrogen bond identity.
*** @ingroup gom_doc_atom_level
***
*** The atom is either a hydrogen bond
*** acceptor (A),
*** donor (D or E) or
*** not hydrogen bonded (N).
***/
/* @{ */
extern       char   gomp_GetAtomHbond(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomHbond(
    GOM_ARG( int         , structure ),
    GOM_ARG( char        , hbond     ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief The CHARMm atom type (atype).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern const char  *gomp_GetAtomAtype(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomAtype(
    GOM_ARG( int         , structure ),
    GOM_ARG( const char *, atype     ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Atom mass (mass).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       float  gomp_GetAtomMass(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomMass(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , mass      ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Atom connectivity (not used at all).
*** @ingroup gom_doc_atom_level
***/
/* @{ */
extern       int    gomp_GetAtomCnct(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomCnct(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , cnct      ),
    GOM_ARG( int         , atom      ) );
/* @} */

/***************************/
/*                         */
/* Atom Display Properties */
/*                         */
/***************************/

/** @defgroup gom_doc_atom_disp Atom Display Properties
*** @ingroup gom_doc_atom
***/

/** @name Display Property Functions
***/
/** @brief Atom display state.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       char   gomp_GetAtomDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( char        , set       ),
    GOM_ARG( int         , atom      ) );
extern const char  *gomp_GetAtomDisplayStatePointer(
    GOM_ARG( int         , structure ) );
extern       char  *gomp_GetModifiableAtomDisplayStatePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom CPK display state.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       char   gomp_GetAtomCPKDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomCPKDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( char        , set       ),
    GOM_ARG( int         , atom      ) );
extern const char  *gomp_GetAtomCPKDisplayStatePointer(
    GOM_ARG( int         , structure ) );
extern       char  *gomp_GetModifiableAtomCPKDisplayStatePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom CPK scale.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       float  gomp_GetAtomCPKScale(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomCPKScale(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , scale     ),
    GOM_ARG( int         , atom      ) );
extern const float *gomp_GetAtomCPKScalePointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomCPKScalePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom licorice display state.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       char   gomp_GetAtomLicoDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomLicoDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( char        , set       ),
    GOM_ARG( int         , atom      ) );
extern const char  *gomp_GetAtomLicoDisplayStatePointer(
    GOM_ARG( int         , structure ) );
extern       char  *gomp_GetModifiableAtomLicoDisplayStatePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom licorice sphere radius.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       float  gomp_GetAtomLicoRadS(
    GOM_ARG( int         , structure    ) );
extern       int    gomp_PutAtomLicoRadS(
    GOM_ARG( int         , structure    ),
    GOM_ARG( float       , sphereRadius ) );
/* @} */
/** @brief Atom licorice cylinder radius.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       float  gomp_GetAtomLicoRadC(
    GOM_ARG( int         , structure      ) );
extern       int    gomp_PutAtomLicoRadC(
    GOM_ARG( int         , structure      ),
    GOM_ARG( float       , cylinderRadius ) );
/* @} */
/** @brief Atom label display state.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       char   gomp_GetAtomLabelDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomLabelDisplayState(
    GOM_ARG( int         , structure ),
    GOM_ARG( char        , set       ),
    GOM_ARG( int         , atom      ) );
extern const char  *gomp_GetAtomLabelDisplayStatePointer(
    GOM_ARG( int         , structure ) );
extern       char  *gomp_GetModifiableAtomLabelDisplayStatePointer(
    GOM_ARG( int         , structure ) );
/* @} */
/** @brief Atom colour.
*** @ingroup gom_doc_atom_disp
***/
/* @{ */
extern       int    gomp_GetAtomColour(
    GOM_ARG( int         , structure ),
    GOM_ARG( float      *, red       ),
    GOM_ARG( float      *, green     ),
    GOM_ARG( float      *, blue      ),
    GOM_ARG( int         , atom      ) );
extern       int    gomp_PutAtomColour(
    GOM_ARG( int         , structure ),
    GOM_ARG( float       , red       ),
    GOM_ARG( float       , green     ),
    GOM_ARG( float       , blue      ),
    GOM_ARG( int         , atom      ) );
extern const float *gomp_GetAtomColourRedPointer(
    GOM_ARG( int         , structure ) );
extern const float *gomp_GetAtomColourGreenPointer(
    GOM_ARG( int         , structure ) );
extern const float *gomp_GetAtomColourBluePointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomColourRedPointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomColourGreenPointer(
    GOM_ARG( int         , structure ) );
extern       float *gomp_GetModifiableAtomColourBluePointer(
    GOM_ARG( int         , structure ) );
/* @} */

/********************/
/*                  */
/* Atom Connections */
/*                  */
/********************/

/** @defgroup gom_doc_atom_conn Atom Connections
*** @ingroup gom_doc_atom
***/

/** @name Atom Connection Functions
***/
/** @brief Atom bonds.
*** @ingroup gom_doc_atom_conn
***/
/* @{ */
extern const int   *gomp_GetAtomConnection(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int   *gomp_GetModifiableAtomConnection(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
/* @} */
/** @brief Atom hydrogen bonds.
*** @ingroup gom_doc_atom_conn
***/
/* @{ */
extern const int   *gomp_GetAtomHydrogenBond(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
extern       int   *gomp_GetModifiableAtomHydrogenBond(
    GOM_ARG( int         , structure ),
    GOM_ARG( int         , atom      ) );
/* @} */

/********************/
/*                  */
/* Atom Listeners   */
/*                  */
/********************/

/** @weakgroup gom_doc_listener
***/

/** @defgroup gom_doc_listener_atom_data Atom Data Changed Listeners
*** @ingroup gom_doc_listener
*** @ingroup gom_doc_atom
***
*** Atom data changed listeners are called after data update is requested
*** if atom data has been changed.
***
*** @sa @ref gomp_UpdateData.
***
*** @{
***/

/** @brief Atom data changed listener handle type.
***/
GOM_DECLARE_HANDLE(gom_AtomDataChangedListener);

/** @brief Atom data changed listener callback type.
*** @param callbackData  Custom data pointer
*** @param mask          Structure mask
***
*** Note that the mask may have overflowed.
*** To check if a structure may have changed, use something similar to
*** @code
*** int               structure = ...; // structure index
*** unsigned long int bit       = ((unsigned long int)1) << structure;
*** if ( ( structureMask & bit ) || ! bit ) {
***     // The structure may have changed.
*** }
*** @endcode
***/
typedef int (*gom_AtomDataChangedListenerFunc)(
    GOM_ARG( void            *, callbackData  ),
    GOM_ARG( unsigned long int, structureMask ) );

/** @brief Register a new atom data changed listener callback.
***
*** @param callback      Listener callback function
*** @param callbackData  Pointer to be passed to @p listener
*** @return       New listener handle
*** @retval NULL  An error occured
***/
extern gom_AtomDataChangedListener* gomp_AddAtomDataChangedListener(
    GOM_ARG( gom_AtomDataChangedListenerFunc, callback    ),
    GOM_ARG( void                          *, callbackData ) );

/** @brief Unregister an atom data changed listener callback by using the handle.
***
*** @param listener  Listener handle
*** @return  Number of canceled listeners (always 1).
***/
extern int gomp_CancelAtomDataChangedListener(
    GOM_ARG( gom_AtomDataChangedListener *, listener ) );

/** @brief Unregister an atom data changed listener callback by using the callback.
***
*** @param callback  Listener callback
*** @return  Number of canceled listeners (always non-negative).
***/
extern int gomp_CancelAtomDataChangedListenersByFunc(
    GOM_ARG( gom_AtomDataChangedListenerFunc, callback ) );

/** @}
***/

/** weakgroup gom_doc_atom_coord
***/

/** @defgroup gom_doc_listener_atom_coord Atom Coordinate Changed Listeners
*** @ingroup gom_doc_listener
*** @ingroup gom_doc_atom_coord
***
*** Atom coordinate changed listeners are called after data update is requested
*** if an atom coordinate has been changed.
***
*** @sa @ref gomp_UpdateData.
***
*** @{
***/

/** @brief Atom coordinate changed listener handle type.
***/
GOM_DECLARE_HANDLE(gom_AtomCoordinateChangedListener);

/** @brief Atom coordinate changed listener callback type.
*** @param callbackData  Custom data pointer
*** @param mask          Structure mask
***
*** Note that the mask may have overflowed.
*** To check if a structure may have changed, use something similar to
*** @code
*** int               structure = ...; // structure index
*** unsigned long int bit       = ((unsigned long int)1) << structure;
*** if ( ( structureMask & bit ) || ! bit ) {
***     // The structure may have changed.
*** }
*** @endcode
***/
typedef int (*gom_AtomCoordinateChangedListenerFunc)(
    GOM_ARG( void            *, callbackData  ),
    GOM_ARG( unsigned long int, structureMask ) );

/** @brief Register a new atom coordinate changed listener callback.
***
*** @param callback      Listener callback function
*** @param callbackData  Pointer to be passed to @p listener
*** @return       New listener handle
*** @retval NULL  An error occured
***/
extern gom_AtomCoordinateChangedListener* gomp_AddAtomCoordinateChangedListener(
    GOM_ARG( gom_AtomCoordinateChangedListenerFunc, callback    ),
    GOM_ARG( void                                 *, callbackData ) );

/** @brief Unregister an atom coordinate changed listener callback by using the handle.
***
*** @param listener  Listener handle
*** @return  Number of canceled listeners (always 1).
***/
extern int gomp_CancelAtomCoordinateChangedListener(
    GOM_ARG( gom_AtomCoordinateChangedListener *, listener ) );

/** @brief Unregister an atom coordinate changed listener callback by using the callback.
***
*** @param callback  Listener callback
*** @return  Number of canceled listeners (always non-negative).
***/
extern int gomp_CancelAtomCoordinateChangedListenersByFunc(
    GOM_ARG( gom_AtomCoordinateChangedListenerFunc, callback ) );

/** @}
***/

/****** PUBLIC GOMAPI END ******/

#endif /* ! ENABLE_EXTENSIONS */

typedef char segment_name_t[MAX_SEG_NAME_LEN+1];
typedef char residue_name_t[MAX_RES_NAME_LEN+1];
typedef char atom_name_t   [MAX_ATM_NAME_LEN+1];

extern const segment_name_t *gomp_GetAtomSegNamePointer(
    GOM_ARG( int         , structure ) );
extern       segment_name_t *gomp_GetModifiableAtomSegNamePointer(
    GOM_ARG( int         , structure ) );
extern const residue_name_t *gomp_GetAtomResNamePointer(
    GOM_ARG( int         , structure ) );
extern       residue_name_t *gomp_GetModifiableAtomResNamePointer(
    GOM_ARG( int         , structure ) );
extern const atom_name_t *gomp_GetAtomAtmNamePointer(
    GOM_ARG( int         , structure ) );
extern       atom_name_t *gomp_GetModifiableAtomAtmNamePointer(
    GOM_ARG( int         , structure ) );

extern int gomp_AssignAtomProperties( int, int, int );
extern int gomp_AssignAtomColourProperties( int, int, const char * );
extern int gomp_AssignAtomBasisSet( int, int, const char * );
extern int gomp_AssignAtomInfoFromDictionary( int, int );

extern int gomp_IdentifyAtom( int, int, int * );
extern int gomp_IdentifyAtoms( int );
extern int gomp_IdentifyAtomColoursAllStructures( void );

extern int gomp_GetMaxAtomConnections( void );
extern int gomp_SetMaxAtomConnections( int );

extern int gomp_GetConnectionSearchWindow( void );

extern int gomp_ResetAtmConn( int );

extern int gomp_SetSearchWindow( int );
extern int gomp_GetSearchWindow( void );

extern int gomp_GetMinRes1Num( void );
extern int gomp_GetMaxRes1Num( void );

extern int gomp_ReadElementParams( const char * );
extern int gomp_ReadAtomParams( const char * );
extern int gomp_ReadAtomConversionTable( const char * );

extern int gomp_MatchAtom( const char * );

extern int         gomp_GetNumberOfAtomSymbols( void );
extern const char *gomp_GetAtomSymbol( int );
extern const char *gomp_Number2Name( int );

extern int gomp_ParseGetAtomInfo( const char * );
extern int gomp_ParseSetAtomCovar( const char *, const char * );

extern int gomp_UpdateNameStack( int, const char * );

extern int gomp_ParseDisplayList( int, const gom_SelectionList * );
extern int gomp_ControlSelectRoundAtoms(
    int, const gom_SelectionList *, float,
    const gom_SelectionList *, const gom_FloatColour * );
extern int gomp_ParseSelectionList(
    int, const char *, const char *, const char * );
extern int gomp_ParseTypeCPKList( int, const gom_SelectionList * );
extern int gomp_ParseTypeLicoList( int, const gom_SelectionList * );
extern int gomp_ParseCPKScaleList( float, const gom_SelectionList * );
extern int gomp_ParseColourList(
    const gom_SelectionList *, const gom_FloatColour * );
extern int gomp_ParseColourListByVector(
    const gom_SelectionList *, double *, double * );
extern int gomp_ParseColourListByCharge(
    const gom_SelectionList *, double *, double * );
extern int gomp_ParseColourListFourth( const gom_SelectionList * );
extern int gomp_ParseColourListResidueNumber( const gom_SelectionList * );
extern int gomp_ShowAtomNumber( int, const gom_SelectionList * );
