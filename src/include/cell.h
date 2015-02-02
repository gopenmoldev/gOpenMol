#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gomcell.h"
#  include "gomext/cell.h"
#endif

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc
***/

/** @defgroup gom_doc_cell Cell Box
*** @ingroup  gom_doc
***/

/** @file
*** @ingroup gom_doc_cell
***/

/** @name Property Functions
***/
/* @{ */
/** @brief Cell box dimension.
*** @ingroup gom_doc_cell
***/
extern float gomp_GetCellA(void);
extern float gomp_GetCellB(void);
extern float gomp_GetCellC(void);
extern int   gomp_SetCellA( GOM_ARG( float, a ) );
extern int   gomp_SetCellB( GOM_ARG( float, b ) );
extern int   gomp_SetCellC( GOM_ARG( float, c ) );
/* @} */
/** @brief Cell box angle.
*** @ingroup gom_doc_cell
***/
/* @{ */
extern float gomp_GetCellAlpha(void);
extern float gomp_GetCellBeta(void);
extern float gomp_GetCellGamma(void);
extern int   gomp_SetCellAlpha( GOM_ARG( float, alpha ) );
extern int   gomp_SetCellBeta( GOM_ARG( float, beta ) );
extern int   gomp_SetCellGamma( GOM_ARG( float, gamma ) );
/* @} */
/** @brief Cell box position.
*** @ingroup gom_doc_cell
***/
/* @{ */
extern float gomp_GetCellXtrans(void);
extern float gomp_GetCellYtrans(void);
extern float gomp_GetCellZtrans(void);
extern int   gomp_SetCellXtrans( GOM_ARG( float, x ) );
extern int   gomp_SetCellYtrans( GOM_ARG( float, y ) );
extern int   gomp_SetCellZtrans( GOM_ARG( float, z ) );
/* @} */
/** @brief Cell box line width.
*** @ingroup gom_doc_cell
***/
/* @{ */
extern int   gomp_GetCellLinewidth(void);
extern int   gomp_SetCellLinewidth( GOM_ARG( int, width ) );
/* @} */
/** @brief Cell box colour.
*** @ingroup gom_doc_cell
***/
/* @{ */
extern int   gomp_GetCellColour(
    GOM_ARG( float *, red   ),
    GOM_ARG( float *, green ),
    GOM_ARG( float *, blue  ) );
extern int   gomp_SetCellColour(
    GOM_ARG( float, red   ),
    GOM_ARG( float, green ),
    GOM_ARG( float, blue  ) );
/* @} */
/** @brief Cell box display state.
*** @ingroup gom_doc_cell
***/
/* @{ */
extern int   gomp_GetPlotCell(void);
extern int   gomp_SetPlotCell( GOM_ARG( int, set ) );
/* @} */

/****** PUBLIC GOMAPI END ******/

extern int   gomp_PlotCellBox(void*,int,int);
