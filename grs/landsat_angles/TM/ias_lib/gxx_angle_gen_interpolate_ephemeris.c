/* IAS Library Includes */
#include "xxx_LogStatus.h"
#include "gxx_geo_math.h"
#include "gxx_angle_gen_private.h"

/* Local Defines */
#define NUM_INTERP  4 /* Number of points to use in Lagrange interpolation */

/******************************************************************************/
/**
 * @file gxx_angle_gen_interpolate_ephemeris.c
 * @brief Interpolate the provided ephemeris structure points at the provided
 *  input time.
 * @ingroup gpsAngleCoef
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
int gxx_angle_gen_interpolate_ephemeris
(
    const gxx_angle_gen_ephemeris_TYPE *ephemeris, /* I: Metadata ephemeris 
                                                         points */
    unsigned int ephem_count,    /* I: Number of entries in ephem and time */
    double in_time,              /* I: Time from epoch to interpolate */
    VECTOR *vector               /* O: Output interpolated vector */
)  
{
    int point;                   /* Loop index */
    unsigned int index;          /* Subindex at which to start interpolation */
    unsigned int used_count;     /* Number of points used */
    VECTOR ephem[NUM_INTERP];    /* Local array of vectors for interpolation */
    double time[NUM_INTERP];     /* Local array of sample times */
    char msg[STRLEN];            /* status messge */

    if (ephem_count < 1)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
                      "Not enough points provided to interpolate");
        return ERROR;
    }

    /* If not enough points load all of them */
    if (ephem_count <= NUM_INTERP)
    {
        for (point = 0; point < (int)ephem_count; point++)
        {
            ephem[point] = ephemeris[point].position;
            time[point] = ephemeris[point].sample_time;
        }
        used_count = ephem_count;
    }
    else
    {
        /* Find the samples that straddle the desired time */
        for (point = 0; point < (int)ephem_count && ephemeris[point].sample_time 
             < in_time; point++);
      
        point -= (NUM_INTERP + 1) / 2;

        if (point < 0)
        { 
            point = 0;
        }
        if (point > (int)(ephem_count - NUM_INTERP))
        { 
            point = ephem_count - NUM_INTERP;
        }

        /* Load the local array */
        for (index = 0; index < NUM_INTERP; index++)
        {
            ephem[index] = ephemeris[point + index].position;
            time[index] = ephemeris[point + index].sample_time;
        }
        used_count = NUM_INTERP;
    }

    /* Run the interpolation on the data window */
    gxx_vector_lagrange(ephem, time, used_count, in_time, vector, msg);

    return SUCCESS;
}
