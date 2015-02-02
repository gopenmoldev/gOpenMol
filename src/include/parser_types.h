#ifndef INC_GOPENMOL_PARSER_TYPES
#define INC_GOPENMOL_PARSER_TYPES

typedef struct {
    const char *Segment;
    const char *Residue;
    const char *Atom;
} gom_SelectionList;

typedef struct {
    float red;
    float green;
    float blue;
} gom_FloatColour;

#endif /* INC_GOPENMOL_PARSER_TYPES */
