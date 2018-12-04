#include <stdlib.h>
#include <string.h>
#include "xxx_LogStatus.h"
#include "xxx_Types.h"
#include "gxx_coord_conv.h"

#define MAXMIN 60060.0
#define MAXSEC 60.0
#define MINMIN -60060.0
#define MINSEC -60.0

/***************************************************************************/
/**
 * @file gxx_dmsdeg.c
 * @brief To convert packed degress, minutes, seconds
 * to total degrees, total minutes, or total seconds.
 * @ingroup gpsGeo
 */
/**
 * @brief To convert packed degress, minutes, seconds
 * to total degrees, total minutes, or total seconds.
 *
 * ALGORITHM DESCRIPTION:
 * Receive an angle in seconds, minutes, or degrees
 * Convert it to DMS.
 * The angle is then checked to be sure it is within the limits
 * of its use (LAT, LON, or DEGREES).
 *
 * @return SUCCESS: Successful completion
 * @return ERROR: Operation failed
 */
/****************************************************************************/

int gxx_dmsdeg
(
    double dms,         /* //!<[in] Angle in DMS (DDDMMMSSS) format */
    double *deg,        /* //!<[out] Angle in decimal degrees */
    char *check         /* //!<[in] Angle usage type (LAT, LON, or DEGREES) */
)
{
    float second;     /* Second of DMS angle value */
    float upper;      /* Upper bound of the angle value for its use */
    float lower;      /* Lower bound of the angle value for its use */

    int  degree;      /* Degree of DMS angle value */
    int  minute;      /* Minute of DMS angle value */
    short sign;       /* Sign of the angle */
    char ErrMsg[STRLEN];    /* Error Message */

    /* Find the upper and lower bound limits for the angle based on its
       usage type.  */
    if (strncmp (check,"LAT",3) == 0)
    { /* LAT */
        upper = 90.0;
        lower = -90.0;
    } /* LAT */
    else if (strncmp (check,"LON",3) == 0)
    { /* LON */
        upper = 180.0;
        lower = -180.0;
    } /* LON */
    else
    { /* DEGREES */
        upper = 360.0;
        lower = 0.0;
    }

   /* Convert the angle to total degrees based on DMS.  */
    if (dms < 0)
    {
        sign = -1;
        dms = dms * -1;
    }
    else
        sign = 1;

    degree = (int) (dms / 1000000);
    dms = dms - (degree * 1000000);
    minute = (int) (dms / 1000);
    second = dms - (minute * 1000);
    *deg = sign * (degree + (minute/60.0) + (second/3600.0));

    /* Check to make sure the angle is within the limits of its use (LAT, LON,
       or DEGREES) */
    if ((*deg > upper) || (*deg < lower))
    {
        sprintf (ErrMsg, "Illegal coordinate value in decdeg");
        xxx_LogStatus("GPSLIB", __FILE__, __LINE__, ErrMsg);
        sprintf (ErrMsg, "Calculated angle of %f outside bounds of %f to %f",
            *deg, lower, upper);
        xxx_LogStatus("GPSLIB", __FILE__, __LINE__, ErrMsg);
        return ERROR;
    }

    return SUCCESS;
}
