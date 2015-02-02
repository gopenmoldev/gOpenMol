/*
  Copyright (c) 2001 - 2005 by:
  Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
  Confidential unpublished property of 
  Leif Laaksonen  
  All rights reserved

  Coded by: Eero Häkkinen

  Enhancements 2001 by:
  Leif Laaksonen

  Enhancements 2002 by:
  Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <ctype.h>
#include "gommath.h"
#include <stdlib.h>
#include <tcl.h>

#include "bond.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

typedef struct {
    /* d for donor, h for hydrogen, a for acceptor, */
    /* aa for atom bonded to acceptor               */
    struct {
        float da2_max;
        float da2_max_nohydrogens;
        float ha2_max;
        float dh_default;
    } length;
    struct {
        double dha_min;
        double da_aa_min;
        double ha_aa_min;
        double aromatic_max;
    } angle;
} HbondCriteria_t;

static const HbondCriteria_t DefaultHbondCriteria = {
    { 3.9f*3.9f , 3.5f*3.5f , 2.5f*2.5f , 1.0f },
    { 90*M_PI/180 , 90*M_PI/180 , 90*M_PI/180 , 20*M_PI/180 }
};

static HbondCriteria_t HbondCriteria;

static struct HbondData_t {
    int SetCount;
    struct {
        int   subset_length;
        int  *subset_list;
        int   do_search;
        char *search_list;
    } *Sets;
    int   SearchHydrogens;
    struct {
        float red;
        float green;
        float blue;
    } Colour;
} HbondData = { 0 , NULL , 1 , { 0.0f , 1.0f , 1.0f } };

typedef struct {
    int from;
    int to;
} HbondSearchData_t;

/************************************************************************/
/*                                                                      */
/* The hydrogen bond test used by these functions are the one described */
/* in HBPLUS userguide (http://www.csb.yale.edu/userguides/             */
/* datamanip/hbplus/hbplus_descrip.html).                               */
/*                                                                      */
/* There are five criteria for all atoms:                               */
/* -distance between donor and acceptor is less than 3.9 Angstroms      */
/* -distance between hydrogen and acceptor is less than 2.5 Angstroms   */
/* -angle donor-hydrogen-acceptor is greater than 90 degrees            */
/* -angle donor-acceptor-any_atom_bonded_to_acceptor                    */
/*          is greater than 90 degrees                                  */
/* -angle hydrogen-acceptor-any_atom_bonded_to_acceptor                 */
/*          is greater than 90 degrees                                  */
/*                                                                      */
/* The length criteria are tested by function TestHbondLengthCriteria   */
/* and the angle criteria are tested by function TestHbondAngleCriteria.*/
/*                                                                      */
/* If the acceptor is at an amino-aromatic ring or is a side atom of an */
/* aromatic ring the potential hydrogen bond must pass two tests more   */
/* (function IsAtomAtAminoAromaticRing tests if atom is at an amino-    */
/* aromatic ring). The extra test are:                                  */
/* -angle between the aromatic plane and the line donor-acceptor must   */
/*  be smaller than 20 degree.                                          */
/* -angle between the aromatic plane and the line hydrogen-acceptor     */
/*  must be smaller than 20 degree.                                     */
/*                                                                      */
/* These function require the atom connection table (returned by        */
/* function gomp_GetAtomConnection).                                          */
/*                                                                      */
/* These functions can also be used to find hydrogen bond in a molecule */
/* without hydrogens (the result isn't perfect then). In that case      */
/* the distance and angle test which need hydrogens are omitted.        */
/* The distance between the donor and the hydrogen is supposed to be    */
/* 1.0 Angstroms and therefore the distance between the donor and the   */
/* acceptor must be less than 3.5 Angstroms.                            */
/* If the donor is a nitrogen at an aromatic ring the hydrogen is       */
/* supposed to at the same plane as the nitrogen and the two atoms      */
/* bonded to the nitrogen and all the angle in that plane are supposed  */
/* to be 120 degree.                                                    */
/*                                                                      */
/************************************************************************/

static int TestHbondLengthCriteria(int i,int j,float max_length2,
                                   const float *X_p,
                                   const float *Y_p,
                                   const float *Z_p)
{
    register float diff,r2;

    diff = X_p[i] - X_p[j];
    r2 = diff * diff;
    if(r2 > max_length2) return(1);

    diff = Y_p[i] - Y_p[j];
    r2 += diff * diff;
    if(r2 > max_length2) return(1);

    diff = Z_p[i] - Z_p[j];
    r2 += diff * diff;
    if(r2 > max_length2) return(1);

    return(0);
}

typedef struct {float x,y,z;} vector_t;

static int TestAminoAromaticHbondLengthCriteria(int donor,int acceptor,
                                                const int * ring_atoms,int ring_atom_count,
                                                const float *X_p,
                                                const float *Y_p,
                                                const float *Z_p)
{
    float s,w;
    vector_t r;

    /* Hydrogen atom bonded to an aromatic nitrogen atom is */
    /* straight away from the centre of the aromatic ring.  */
    /* Find the position of the hydrogen atom.              */
    r.x = X_p[ring_atoms[0]] - ( X_p[ring_atoms[1]] + X_p[ring_atoms[ring_atom_count-1]] )/2;
    r.y = Y_p[ring_atoms[0]] - ( Y_p[ring_atoms[1]] + Y_p[ring_atoms[ring_atom_count-1]] )/2;
    r.z = Z_p[ring_atoms[0]] - ( Z_p[ring_atoms[1]] + Z_p[ring_atoms[ring_atom_count-1]] )/2;

    s  = r.x*r.x;
    s += r.y*r.y;
    s += r.z*r.z;
/*  s  = (float)(default_hb_length_dh / sqrt( s )); */
    s  = (float)(HbondCriteria.length.dh_default / sqrt( s ));

    r.x *= s;
    r.y *= s;
    r.z *= s;

    /* Test the length of the hydrogen bond.    */
    w  = X_p[donor]+r.x-X_p[acceptor];
    s  = w*w;
    w  = Y_p[donor]+r.y-Y_p[acceptor];
    s += w*w;
    w  = Z_p[donor]+r.z-Z_p[acceptor];
    s += w*w;

/*  if( s > max_hb_length_ha2 ) */
    if( s > HbondCriteria.length.ha2_max )
        return(1);
    return(0);
}

static int TestHbondAngleCriteria(int i,int j,int k,
                                  const float *X_p,
                                  const float *Y_p,
                                  const float *Z_p,
                                  float angle_min)
{
    register float a,s;
    vector_t ji,jk;

    /* calculate unit vector from point j to point i */
    ji.x = X_p[i] - X_p[j];
    ji.y = Y_p[i] - Y_p[j];
    ji.z = Z_p[i] - Z_p[j];

    s  = ji.x*ji.x;
    s += ji.y*ji.y;
    s += ji.z*ji.z;
    s  = (float)(1.0 / sqrt( s ));

    ji.x *= s;
    ji.y *= s;
    ji.z *= s;

    /* calculate unit vector from point j to point k */
    jk.x = X_p[k] - X_p[j];
    jk.y = Y_p[k] - Y_p[j];
    jk.z = Z_p[k] - Z_p[j];

    s  = jk.x*jk.x;
    s += jk.y*jk.y;
    s += jk.z*jk.z;
    s  = (float)(1.0 / sqrt( s ));

    jk.x *= s;
    jk.y *= s;
    jk.z *= s;

    a = ji.x*jk.x + ji.y*jk.y + ji.z*jk.z;
    /* A is cosine of the angle between vectors ji and jk.  */
    if( a < cos(angle_min) )
        return(0);
    return(1);
}

static int TestAminoAromaticCriteria(int donor,int hydrogen,int acceptor,
                                     const int *ring_atoms,int ring_atom_count,
                                     const float *X_p,
                                     const float *Y_p,
                                     const float *Z_p)
{
    vector_t tangent1,tangent2,normal,position;
    float a,s;
    int i;

    /* Calculate two tangent vector for the aromatic plane. */
    tangent1.x = X_p[ring_atoms[0]] - ( X_p[ring_atoms[1]] + X_p[ring_atoms[ring_atom_count-1]] )/2;
    tangent1.y = Y_p[ring_atoms[0]] - ( Y_p[ring_atoms[1]] + Y_p[ring_atoms[ring_atom_count-1]] )/2;
    tangent1.z = Z_p[ring_atoms[0]] - ( Z_p[ring_atoms[1]] + Z_p[ring_atoms[ring_atom_count-1]] )/2;

    tangent2.x = X_p[ring_atoms[1]] - X_p[ring_atoms[ring_atom_count-1]];
    tangent2.y = Y_p[ring_atoms[1]] - Y_p[ring_atoms[ring_atom_count-1]];
    tangent2.z = Z_p[ring_atoms[1]] - Z_p[ring_atoms[ring_atom_count-1]];

    /* Normal vector will be an unit normal vector of the aromatic plane    */
    normal.x = tangent1.y*tangent2.z - tangent1.z*tangent2.y;
    normal.y = tangent1.z*tangent2.x - tangent1.x*tangent2.z;
    normal.z = tangent1.x*tangent2.y - tangent1.y*tangent2.x;

    s  = normal.x*normal.x;
    s += normal.y*normal.y;
    s += normal.z*normal.z;
    s  = (float)(1.0 / sqrt( s ));

    normal.x *= s;
    normal.y *= s;
    normal.z *= s;

    /* On the first run test the angle between the line     */
    /* donor-acceptor and the the aromatic plane.           */
    /* On the second run test the angle between the line    */
    /* hydrogen-acceptor and the aromatic plane.            */
    for(i=0;i<2 && donor>=0;++i,donor=hydrogen)
    {
        position.x = X_p[donor] - X_p[acceptor];
        position.y = Y_p[donor] - Y_p[acceptor];
        position.z = Z_p[donor] - Z_p[acceptor];

        s  = position.x*position.x;
        s += position.y*position.y;
        s += position.z*position.z;
        s  = (float)(1.0 / sqrt( s ));

        position.x *= s;
        position.y *= s;
        position.z *= s;

        a = position.x*normal.x + position.y*normal.y + position.z*normal.z;
        /* A is cosine of the angle between aromatic plane normal and   */
        /* the position vector. We want to test the angle between       */
        /* the aromatic plane and the position vector.                  */
        if( fabs(a) > sin(HbondCriteria.angle.aromatic_max) )
            /* cos(M_PI/2-max_aromatic_angle) == sin(HbondCriteria.angle.aromatic_max)    */
            return(1);
    }
    return(0);
}

/* Count should be zero if this function is called outside this function.   */
static int IsAtomAtAminoAromaticRing(int Wstr,int atom,int count,
                                     int first_atom,int *ring_atoms)
{
    int i;
    int ring_atom_count=-1;
    const int *AtmConn;
    int atom_type=-1;   /* -1 for carbon, +1 for nitrogen   */
    const char * atom_name;

    for(i=0; i<count; ++i)
    {
        if( ring_atoms[i] == atom )
            /* We have already found this atom. */
            return(0);
    }

    atom_name = gomp_GetAtomAtype(Wstr,atom);

    if( strcmp(atom_name,"C") )
    {
        if( strcmp(atom_name,"N") )
            return(0);
        else
            atom_type=1;
    }

    ring_atoms[count] = atom;

    AtmConn = gomp_GetAtomConnection(Wstr, atom);
    for(i=1; i<=AtmConn[0]; ++i) {
        if( count == 5 )
        {
            if( AtmConn[i] == first_atom )
                /* There are six carbon and nitrogen atoms in a ring.   */
                return(atom_type);
        }
        else if( (count == 4) && (AtmConn[i] == first_atom) )
            /* There are five carbon and nitrogen atoms in a ring.   */
            return(atom_type);
        else
        {
            ring_atom_count = IsAtomAtAminoAromaticRing(
                Wstr,AtmConn[i],count+1,first_atom,ring_atoms);
            if( ring_atom_count > 0 )
                /* We have found a ring of five or six atoms and    */
                /* at least one of them is a nitrogen.              */
                return(ring_atom_count+1);
            else if( ring_atom_count < 0 )
                /* We have found a ring of five or six atoms but    */
                /* we aren't sure if any of them is a nitrogen.     */
                /* Return positive number only if this atom is a    */
                /* nitrogen.                                        */
                return(atom_type*(ring_atom_count-1));
        }
    }
    return(0);
}

static int TestHbondCriteria3(int Wstr,int donor,int hydrogen,int acceptor,
                              const float *X_p,
                              const float *Y_p,
                              const float *Z_p)
{
    int i;
    int last_bonded_atom=-1,bonded_atom_count=0;
    const int *AtmConn;
    int ring_atoms[6];
    int ring_atom_count = 0;
    const char * atom_name;

    /* Test the distance between the hydrogen and the acceptor. */
    if(TestHbondLengthCriteria(hydrogen,acceptor,
                               HbondCriteria.length.ha2_max,X_p,Y_p,Z_p))
/*      max_hb_length_ha2,X_p,Y_p,Z_p)) */
        return(1);

    /* Test the angle of donor-hydrogen-acceptor.   */
    if(TestHbondAngleCriteria(donor,hydrogen,acceptor,X_p,Y_p,Z_p,
                              HbondCriteria.angle.dha_min))
        return(1);
    AtmConn = gomp_GetAtomConnection(Wstr, acceptor);
    for(i=1; i<=AtmConn[0]; ++i) {
        /* Test the angle of donor-acceptor-atom_bonded_to_acceptor.    */
        if(TestHbondAngleCriteria(donor,acceptor,AtmConn[i],X_p,Y_p,Z_p,
                                  HbondCriteria.angle.da_aa_min))
            return(1);
        /* Test the angle of hydrogen-acceptor-atom_bonded_to_acceptor. */
        if(TestHbondAngleCriteria(hydrogen,acceptor,AtmConn[i],X_p,Y_p,Z_p,
                                  HbondCriteria.angle.ha_aa_min))
            return(1);

        atom_name = gomp_GetAtomAtype(Wstr,AtmConn[i]);
        if( !strcmp(atom_name,"C") )
        {
            last_bonded_atom = AtmConn[i];
            ++bonded_atom_count;
        }
    }
    if( bonded_atom_count == 1 )
        /* The acceptor may be a side atom of an amino aromatic ring. */
        ring_atom_count = IsAtomAtAminoAromaticRing(
            Wstr,last_bonded_atom,0,last_bonded_atom,ring_atoms);
    else
        /* The acceptor may be in an amino aromatic ring. */
        ring_atom_count = IsAtomAtAminoAromaticRing(
            Wstr,acceptor,0,acceptor,ring_atoms);
    if( ring_atom_count > 0 )
        return TestAminoAromaticCriteria(donor,hydrogen,acceptor,
                                         ring_atoms,ring_atom_count,X_p,Y_p,Z_p);
    return(0);
}

static int TestHbondCriteria2(int Wstr,int donor,int acceptor,
                              const float *X_p,
                              const float *Y_p,
                              const float *Z_p,
                              int no_extra_test)
{
    int i;
    int last_bonded_atom=-1,bonded_atom_count=0;
    const int *AtmConn;
    int ring_atoms[6];
    int ring_atom_count=0;
    const char * atom_name;

    /* Test the distance between the hydrogen and the acceptor. */
    /* We don't have a hydrogen atom so we suppose it to be     */
    /* between the donor and the acceptor.                      */
    if(TestHbondLengthCriteria(donor,acceptor,
/*      max_hb_length_dha2,X_p,Y_p,Z_p)) */
                               HbondCriteria.length.da2_max_nohydrogens,X_p,Y_p,Z_p))
        return(1);

    AtmConn = gomp_GetAtomConnection(Wstr, acceptor);
    for(i=1; i<=AtmConn[0]; ++i) {
        /* Test the angle of donor-acceptor-atom_bonded_to_acceptor.    */
        if(TestHbondAngleCriteria(donor,acceptor,AtmConn[i],X_p,Y_p,Z_p,
                                  HbondCriteria.angle.da_aa_min))
            return(1);

        atom_name = gomp_GetAtomAtype(Wstr,AtmConn[i]);
        if( !strcmp(atom_name,"C") )
        {
            last_bonded_atom = AtmConn[i];
            ++bonded_atom_count;
        }
    }
    if( bonded_atom_count == 1 )
        /* The acceptor may be a side atom of an amino aromatic ring. */
        ring_atom_count = IsAtomAtAminoAromaticRing(
            Wstr,last_bonded_atom,0,last_bonded_atom,ring_atoms);
    else if( bonded_atom_count > 1 )
        /* The acceptor may be in an amino aromatic ring. */
        ring_atom_count = IsAtomAtAminoAromaticRing(
            Wstr,acceptor,0,acceptor,ring_atoms);
    if( ring_atom_count > 0) {
        if(TestAminoAromaticCriteria(donor,-1,acceptor,
                                     ring_atoms,ring_atom_count,X_p,Y_p,Z_p))
            return(1);

        if( !no_extra_test && (bonded_atom_count > 1) ) {
            atom_name = gomp_GetAtomAtype(Wstr,acceptor);
            if( !strcmp(atom_name,"N") &&
                TestHbondCriteria2(Wstr,acceptor,donor,X_p,Y_p,Z_p,1) ) {
                /* Our acceptor is a nitrogen at an aromatic ring but we    */
                /* can swap the acceptor and the donor. Test if this would  */
                /* be a hydrogen bond after swapping and hydrogen           */
                /* localizating, too.                                       */
                if(TestAminoAromaticHbondLengthCriteria(
                    acceptor,donor,ring_atoms,ring_atom_count,X_p,Y_p,Z_p))
                    return(1);
            }
        }
    }
    if( !no_extra_test ) {
        atom_name = gomp_GetAtomAtype(Wstr,donor);
        if( !strcmp(atom_name,"N") ) {
            /* Our donor is a nitrogen. Test if it is at an aromatic ring.  */
            ring_atom_count = IsAtomAtAminoAromaticRing(
                Wstr,donor,0,donor,ring_atoms);
            if( ring_atom_count > 0) {
                /* Test if this would be a hydrogen bond after hydrogen     */
                /* localizating, too.                                       */
                if(TestAminoAromaticHbondLengthCriteria(
                    donor,acceptor,ring_atoms,ring_atom_count,X_p,Y_p,Z_p))
                    return(1);
            }
        }
    }
    return(0);
}

static int FindHbondAcceptor(int Wstr,int donor,int hydrogen,
                             HbondSearchData_t *pSearchData)
{
    const float *X_p,*Y_p,*Z_p;
    int i,j,acceptor,hydrogen_or_donor;
    int *AtmHbond;
    const char *search_list;
    const int *subset_list;
    const char *atom_name;

    search_list = HbondData.Sets[Wstr].search_list;
    subset_list = HbondData.Sets[Wstr].subset_list;

    if( hydrogen >= 0 )
        hydrogen_or_donor = hydrogen;
    else
        hydrogen_or_donor = donor;

    X_p = gomp_GetAtomXCoordPointer(Wstr);
    Y_p = gomp_GetAtomYCoordPointer(Wstr);
    Z_p = gomp_GetAtomZCoordPointer(Wstr);

    /* Search the acceptor from the search window */
    for( i=pSearchData->from; i<pSearchData->to; i++ ) {

        if( subset_list )
            acceptor = subset_list[i];
        else
            acceptor = i;

        if( search_list &&
            !search_list[donor] &&
            !search_list[acceptor] )
            continue;

        if( acceptor==donor || acceptor==hydrogen ) continue;

        atom_name = gomp_GetAtomAtype(Wstr,acceptor);

        if( strcmp(atom_name,"O") && strcmp(atom_name,"N") )
            /* We don't recognize this atom to be an acceptor.   */
            continue;

        /* Test the distance between the donor and the acceptor.    */
        if(TestHbondLengthCriteria(donor,acceptor,
                                   HbondCriteria.length.da2_max,X_p,Y_p,Z_p))
            continue;

        if( hydrogen >= 0 ) {
            if(TestHbondCriteria3(Wstr,donor,hydrogen,acceptor,X_p,Y_p,Z_p))
                continue;
        }
        else {
            /* Test that there aren't already a hydrogen bond   */
            /* between these atoms. We don't have a hydrogen    */
            /* so can't absolutely sure which atom is a donor   */
            /* and which one is an acceptor.                    */
            const int *AtmHbond = gomp_GetAtomHydrogenBond(Wstr, acceptor);
            for( j=1; j<=AtmHbond[0]; j++ )
            {
                if( AtmHbond[j] == donor )
                    break;
            }
            if( j <= AtmHbond[0] )
                /* There is already a hydrogen bond between */
                /* the acceptor and the donor.              */
                continue;

            if(TestHbondCriteria2(Wstr,donor,acceptor,X_p,Y_p,Z_p,0))
                continue;
        }

        AtmHbond = gomp_GetModifiableAtomHydrogenBond(Wstr, acceptor);
        if( AtmHbond[0]+1 >= gomp_GetMaxAtomConnections() ) {
            gomp_PrintMessage("*** ERROR . Max atom hydrogen bond exceeded");
            gomp_PrintERROR("Problems with your coordinates.\nYou have reached max number of currently allowed bonds.\nPlease check internal format.\nNew value can be defined with the command:\n'define atom maxco Ivalue'.");
            return(1);
        }

        AtmHbond[++AtmHbond[0]] = hydrogen_or_donor;

        AtmHbond = gomp_GetModifiableAtomHydrogenBond(Wstr, hydrogen_or_donor);
        if( AtmHbond[0]+1 >= gomp_GetMaxAtomConnections() ) {
            gomp_PrintMessage("*** ERROR . Max atom hydrogen bond exceeded");
            gomp_PrintERROR("Problems with your coordinates.\nYou have reached max number of currently allowed bonds.\nPlease check internal format.\nNew value can be defined with the command:\n'define atom maxco Ivalue'.");
            return(1);
        }

        AtmHbond[++AtmHbond[0]] = acceptor;
    }
    return(0);
}

static int GetHydrogenBondingDistanceParam(float *target,const char *name,int Square)
{
    const char *Text;
    char   Temp[BUFF_LEN];
    float  Value;

    sprintf(Temp,"gomHydrogenBondingParams(%s)",name);
    Text = Tcl_GetVar(gomp_GetTclInterp(),Temp,TCL_GLOBAL_ONLY);

    if(!Text)
        return(1);

    Value = 0.0f;
    sscanf(Text,"%f",&Value);
    if(Value <= 0.0)
        return(1);

    if(Square)
        *target = Value * Value;
    else
        *target = Value;

    return(0);
}

static int GetHydrogenBondingAngleParam(double  *target,const char *name)
{
    const char *Text;
    char    Temp[BUFF_LEN];
    double  Value;

    sprintf(Temp,"gomHydrogenBondingParams(%s)",name);
    Text = Tcl_GetVar(gomp_GetTclInterp(),Temp,TCL_GLOBAL_ONLY);

    if(!Text)
        return(1);

    Value = 0.0;
    sscanf(Text,"%lf",&Value);
    if(Value <= 0.0)
        return(1);

    *target = M_PI*Value/180;

    return(0);
}

static int ReserveSpaceForSet(int Wstr)
{
    int  i;
    void *NewTable;

    if( Wstr >= HbondData.SetCount ) {
        /* Reserve space for new set. */
        if( HbondData.SetCount )
            NewTable = gomp_ReallocateVoidVector(
                HbondData.Sets,(Wstr+1)*sizeof(*HbondData.Sets));
        else
            NewTable = gomp_AllocateCharVector((Wstr+1)*sizeof(*HbondData.Sets));

        if( !NewTable )
            return(1);

        HbondData.Sets = NewTable;

        /* Init newly created subsets. */
        for( i=HbondData.SetCount; i<=Wstr; i++) {
            HbondData.Sets[i].subset_length = 0;
            HbondData.Sets[i].subset_list   = NULL;
            HbondData.Sets[i].do_search     = 0;
            HbondData.Sets[i].search_list   = NULL;
        }

        HbondData.SetCount = Wstr + 1;
    }

    return(0);
}

/***********************************************************************/
int gomp_SetHbondSubset(int Wstr,int slong,const int *sel_list)
/***********************************************************************/
{
    int  NAtoms,*new_list;

    NAtoms = gomp_GetNumAtomsInMolecStruct(Wstr);

    if(ReserveSpaceForSet(Wstr))
        return(1);

    if( slong >= NAtoms )
        /* List contains all atoms. So it is needless. */
        new_list = NULL;
    else {
        new_list = gomp_AllocateIntVector( (slong>0) ? slong : 1 );
        if( !new_list )
            return(1);
    }

    /* Overwrite the old subset. */
    if( HbondData.Sets[Wstr].subset_length )
        free(HbondData.Sets[Wstr].subset_list);

    /* Copy subset. */
    if( new_list && slong > 0 )
        memcpy(new_list, sel_list, slong*sizeof(*new_list));

    HbondData.Sets[Wstr].subset_length = slong;
    HbondData.Sets[Wstr].subset_list   = new_list;
    HbondData.Sets[Wstr].do_search     = 1;

    return(0);
}

/***********************************************************************/
int gomp_RecalculateHbonds(int Wstr)
/***********************************************************************/
{
    int   i,j,donor;
    int *AtmHbond;
    const int *AtmConn;
    const int *subset_list;
    const char *search_list;
    const char *atom_name;
    HbondSearchData_t searchData;

    if( Wstr >= HbondData.SetCount ||
        !HbondData.Sets[Wstr].do_search )
        /* Nothing to be done. */
        return(0);

    subset_list     = HbondData.Sets[Wstr].subset_list;
    search_list     = HbondData.Sets[Wstr].search_list;

    if( subset_list ) {
        searchData.from = 0;
        searchData.to   = HbondData.Sets[Wstr].subset_length;

        for( i=searchData.from; i<searchData.to; ++i ) {
            if(search_list && !search_list[subset_list[i]])
                continue;
            /* turn off everything previously defined */
            AtmHbond = gomp_GetModifiableAtomHydrogenBond(Wstr, subset_list[i]);
            for(j = AtmHbond[0]; j > 0; j--)
                gomp_BreakBond(1, Wstr, subset_list[i], AtmHbond[j]);
            AtmHbond[0] = 0;
        }
    }
    else {
        searchData.from = 0;
        searchData.to   = gomp_GetNumAtomsInMolecStruct(Wstr);

        for( i=searchData.from; i<searchData.to; ++i ) {
            if(search_list && !search_list[i])
                continue;
            /* turn off everything previously defined */
            AtmHbond = gomp_GetModifiableAtomHydrogenBond(Wstr, i);
            for(j = AtmHbond[0]; j > 0; j--)
                gomp_BreakBond(1, Wstr, i, AtmHbond[j]);
            AtmHbond[0] = 0;
        }
    }

    for( i=searchData.from; i<searchData.to; ++i ) {

        if( subset_list )
            donor = subset_list[i];
        else
            donor = i;

        atom_name   = gomp_GetAtomAtype(Wstr,donor);

        if( strcmp(atom_name,"O") && strcmp(atom_name,"N") )
            /* We don't recognize this atom to be a donor.   */
            continue;

        if( HbondData.SearchHydrogens ) {
            /* Find a hydrogen bonded to the donor.  */
            AtmConn = gomp_GetAtomConnection(Wstr, donor);
            for(j=1; j<=AtmConn[0]; ++j) {
                
                if( strcmp(gomp_GetAtomAtype(Wstr,AtmConn[j]),"H") )
                    continue;
                
                FindHbondAcceptor(Wstr,donor,AtmConn[j],&searchData);
            }
        }
        else
            FindHbondAcceptor(Wstr,donor,-1,&searchData);
    }

    return(0);
}

/***********************************************************************/
int gomp_FindHbonds(int Wstr, int SearchHydrogens, int slong, const int *sel_list)
/***********************************************************************/
{
    int   NAtoms,i,from,to;
    char  Text[BUFF_LEN];

    HbondData.SearchHydrogens = SearchHydrogens;
    HbondCriteria             = DefaultHbondCriteria;

    /* Get customized parametres. */
    GetHydrogenBondingDistanceParam(&HbondCriteria.length.da2_max,"d-a",1);
    GetHydrogenBondingDistanceParam(&HbondCriteria.length.da2_max_nohydrogens,"d-a_nohydrogens",1);
    GetHydrogenBondingDistanceParam(&HbondCriteria.length.ha2_max,"h-a",1);
    GetHydrogenBondingDistanceParam(&HbondCriteria.length.dh_default,"d-h",0);

    GetHydrogenBondingAngleParam(&HbondCriteria.angle.dha_min,"d-h-a");
    GetHydrogenBondingAngleParam(&HbondCriteria.angle.da_aa_min,"d-a-aa");
    GetHydrogenBondingAngleParam(&HbondCriteria.angle.ha_aa_min,"h-a-aa");
    GetHydrogenBondingAngleParam(&HbondCriteria.angle.aromatic_max,"aromatic");

    gomp_PrintMessage("\n******* Hydrogen bond calculation *******");
    /* Print distance criteria. */
    if(HbondData.SearchHydrogens) {
        sprintf(Text,"Distance between donor and acceptor is less than %f Angstroms",
                sqrt(HbondCriteria.length.da2_max));
        gomp_PrintMessage(Text);
        sprintf(Text,"Distance between hydrogen and acceptor is less than %f Angstroms",
                sqrt(HbondCriteria.length.ha2_max));
        gomp_PrintMessage(Text);
    }
    else {
        sprintf(Text,"Distance between donor and acceptor is less than %f Angstroms",
                sqrt(HbondCriteria.length.da2_max_nohydrogens));
        gomp_PrintMessage(Text);
    }
    sprintf(Text,"Default hydrogen bond length %f Angstroms",
            HbondCriteria.length.dh_default);
    gomp_PrintMessage(Text);
    /* Print angle criteria. */
    if(HbondData.SearchHydrogens) {
        sprintf(Text,"Angle donor-hydrogen-acceptor is greater than %f degrees",
                180.0 * HbondCriteria.angle.dha_min / M_PI);
        gomp_PrintMessage(Text);
        sprintf(Text,"Angle hydrogen-acceptor-atom_bonded_to_acceptor is greater than %f degrees",
                180.0 * HbondCriteria.angle.ha_aa_min / M_PI);
        gomp_PrintMessage(Text);
    }
    sprintf(Text,"Angle donor-acceptor-atom_bonded_to_acceptor is greater than %f degrees",
            180.0 * HbondCriteria.angle.da_aa_min / M_PI);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max aromatic angle %f degrees",
            180.0 * HbondCriteria.angle.aromatic_max / M_PI);
    gomp_PrintMessage(Text);

    /* Create a new search list. */
    if( ReserveSpaceForSet(Wstr) )
        return(1);

    NAtoms = gomp_GetNumAtomsInMolecStruct(Wstr);

    if( slong <= 0 || slong >= NAtoms ) {
        free(HbondData.Sets[Wstr].search_list);
        HbondData.Sets[Wstr].do_search   = (slong > 0);
        HbondData.Sets[Wstr].search_list = NULL;
    }
    else {
        if( !HbondData.Sets[Wstr].search_list ) {
            HbondData.Sets[Wstr].search_list = gomp_AllocateCharVector(NAtoms);

            if( !HbondData.Sets[Wstr].search_list )
                return(1);
        }

        HbondData.Sets[Wstr].do_search = 1;
    }

    from = 0;
    to   = gomp_GetNumAtomsInMolecStruct(Wstr);

    if( HbondData.Sets[Wstr].search_list ) {
        for( i=from; i<to; i++ )
            HbondData.Sets[Wstr].search_list[i] = 0;
        /* Save from the selection list */
        for(i = 0 ; i < slong ; i++)
            HbondData.Sets[Wstr].search_list[sel_list[i]] = 1;
    }

    if( slong > 0 )
        (void)gomp_RecalculateHbonds(Wstr);

    gomp_PrintMessage("Done!");

    return(0);
}

/***********************************************************************/
int  gomp_DeleteAllHbonds()
/***********************************************************************/
{
    int   i;
    int   Wstr;
    int   from,to;
    int *AtmHbond;

    for(Wstr = 0 ; Wstr < gomp_GetNumMolecStructs() ; Wstr++) {

        if(!gomp_GetSelectedStructure(Wstr)) continue;

        from = 0;
        to   = gomp_GetNumAtomsInMolecStruct(Wstr);

        for(i=from; i<to; ++i) {
/* turn off everything previously defined */
            AtmHbond    = gomp_GetModifiableAtomHydrogenBond(Wstr, i);
            AtmHbond[0] = 0;
        }
    }

    gomp_DeleteAllHbondSearchData();
    
    return(0);
}
/***********************************************************************/
int  gomp_DeleteAllHbondSearchData()
/***********************************************************************/
{
    int   i,still_in_use;

    still_in_use = 0;

    /* Free search lists. */
    for( i=0; i<HbondData.SetCount; i++ ) {
        free(HbondData.Sets[i].search_list);
        HbondData.Sets[i].do_search   = 0;
        HbondData.Sets[i].search_list = NULL;

        if( HbondData.Sets[i].subset_list )
            still_in_use = 1;
    }

    if( !still_in_use ) {
        free(HbondData.Sets);
        HbondData.SetCount = 0;
        HbondData.Sets     = NULL;
    }

    gomp_DeleteAllHbondSubsets();
    
    return(0);
}
/***********************************************************************/
int  gomp_DeleteAllHbondSubsets()
/***********************************************************************/
{
    int i,still_in_use;

    still_in_use = 0;

    /* Free subsets. */
    for( i=0; i<HbondData.SetCount; i++ ) {
        free(HbondData.Sets[i].subset_list);
        HbondData.Sets[i].subset_length = 0;
        HbondData.Sets[i].subset_list   = NULL;

        if( HbondData.Sets[i].do_search )
            /* There is still need for the sets. */
            still_in_use = 1;
    }

    if( !still_in_use ) {
        free(HbondData.Sets);
        HbondData.SetCount = 0;
        HbondData.Sets     = NULL;
    }

    return(0);
}

/***********************************************************************/
int  gomp_SetHBondColour(float red , float green , float blue)
/***********************************************************************/
{
    HbondData.Colour.red    = red;
    HbondData.Colour.green  = green;
    HbondData.Colour.blue   = blue;

    return(0);
}
/***********************************************************************/
int  gomp_GetHBondColour(float *red , float *green , float *blue)
/***********************************************************************/
{
    *red   = HbondData.Colour.red;
    *green = HbondData.Colour.green;
    *blue  = HbondData.Colour.blue;

    return(0);
}
