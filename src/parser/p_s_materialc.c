#include "maindefs.h"

#include "parser.h"
#include "light_model.h"

#include "stdafx.h"

/* show material specular {red|green|blue} */
#define SHOW_MATERIAL_SPECULAR_ENTRY \
GOM_PARSER_CMD_PART("specular",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowMaterialSpecular,\
                    GOM_PARSER_UNUSED_VALUE)
static const FloatFunc parseShowMaterialSpecularRed   = {
    gomp_GetMaterialSpecularRed };
static const FloatFunc parseShowMaterialSpecularGreen = {
    gomp_GetMaterialSpecularGreen };
static const FloatFunc parseShowMaterialSpecularBlue  = {
    gomp_GetMaterialSpecularBlue };
static const gom_ParserArgumentList parseShowMaterialSpecular[] = {
    GOM_PARSER_FINAL_FLOAT_CMD("red"  ,&parseShowMaterialSpecularRed),
    GOM_PARSER_FINAL_FLOAT_CMD("green",&parseShowMaterialSpecularGreen),
    GOM_PARSER_FINAL_FLOAT_CMD("blue" ,&parseShowMaterialSpecularBlue),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show material shininess */
#define SHOW_MATERIAL_SHININESS_ENTRY \
GOM_PARSER_FINAL_FLOAT_CMD("shininess",&parseShowMaterialShininess)
static const FloatFunc parseShowMaterialShininess = {
    gomp_GetMaterialShininess };

const gom_ParserArgumentList gomp_ShowMaterialCommand[] = {
    SHOW_MATERIAL_SPECULAR_ENTRY,
    SHOW_MATERIAL_SHININESS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
