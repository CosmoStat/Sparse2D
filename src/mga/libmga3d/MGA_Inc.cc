
#include "MGA_Inc.h"

// strings associated to enumerations
const char * StringBlockType (type_3dcurvelet_block type)
{
    switch (type)
    {
        case BL3D_CONST:
	      return ("Constant block size"); 
        case BL3D_UP: 
              return ("Block size is doubled at each resolution"); 
        case BL3D_UP2: 
              return ("Block size is doubled at each two resolutions"); 
        case BL3D_DOWN:
	      return ("Block size is divided by 2 at each resolution");
        case BL3D_DOWN2:
	      return ("Block size is divided by 2 at each two resolutions");
        default: 
	      return ("Undefined block type");
    }
}

/****************************************************************************/
 

// String manipulation
char *add_fits(char * NameStep)
{
    char Sreturn[256];
    char Name[256]; char* name;
    char *ValRet;

    strcpy(Name, NameStep); name = &Name[0];
    if (strstr(name+max((int)strlen(name)-5,0), ".fits") != NULL)  strcpy(Sreturn, Name);
    else sprintf(Sreturn, "%s.%s", Name, "fits");
    ValRet = strdup(Sreturn);
    return (ValRet);
}

char *add_mr(char * NameStep)
{
    char Sreturn[256];
    char Name[256]; char* name;
    char *ValRet;

    strcpy(Name, NameStep); name = &Name[0];
    if (strstr(name+max((int)strlen(name)-3,0), ".mr") != NULL) strcpy(Sreturn, Name);
    else sprintf(Sreturn, "%s.%s", Name, "mr");
    ValRet = strdup(Sreturn);
    return (ValRet);
}

