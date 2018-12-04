#include <string.h>
#include "xxx_LogStatus.h"
#include "xxx_Types.h"
#include "gxx_proj.h"

/****************************************************************************/
/**
 * @file gxx_get_units.c
 * @brief Get the unit code from the name of the units.  The unit code
 * returned from this function can be entered into gxx_projtran.
 * @ingroup gpsGeo
 */
/**
 * @brief Get the unit code from the name of the units.  The unit code
 * returned from this function can be entered into gxx_projtran.
 */
/****************************************************************************/


int gxx_get_units
(
    const char *unit_name,    //!<[in] Units name 
    int *unit_num       //!<[out] Units number 
)
{
    /* ensure that the unit name is upper case, compare names and assign the
       code number appropriatly
       ---------------------------------------------------------------------*/
    if(strcasecmp(unit_name, "RADIANS") == 0)
        *unit_num = 0;
    else if(strcasecmp(unit_name, "FEET") == 0) 
        *unit_num = 1;
    else if(strcasecmp(unit_name, "METERS") == 0) 
        *unit_num = 2;
    else if(strcasecmp(unit_name, "SECONDS") == 0) 
        *unit_num = 3;
    else if(strcasecmp(unit_name, "DEGREES") == 0) 
        *unit_num = 4;
    else if(strcasecmp(unit_name, "DMS") == 0) 
        *unit_num = 5;
    else
    {
        xxx_LogStatus("GPSLIB",__FILE__,__LINE__,
                      "Illegal projection units for gctp\n"
                      "Valid units include RADIANS, FEET, METERS, SECONDS, "
                      "DEGREES, and DMS");
        return ERROR;
    }

    return SUCCESS;
}
