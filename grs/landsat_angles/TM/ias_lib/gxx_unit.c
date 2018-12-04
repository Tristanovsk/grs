#include "xxx_LogStatus.h"
#include "gxx_geo_math.h"

/****************************************************************************/
/**
 * @file gxx_unit.c
 * @brief Calculate a vector's unit vector
 * @ingroup gpsMath
 */
/**
 * @brief Calculate a vector's unit vector
 *
 * @return int (SUCCESS or ERROR)
 */
/****************************************************************************/

int gxx_unit
(
    const VECTOR *vec, //!<[in] Input vector 
    VECTOR *uni        //!<[out] Unit vector of the input vector 
)
{
    double  mag;            /* Magnitude of the input vector */

    
    /*
      Find the magnitude of the input vector
    */
    mag = gxx_norm(vec);

    if (mag == 0.0)
    {
        xxx_LogStatus ("GPSLIB", __FILE__, __LINE__, 
                       "Error: input vector magnitude is zero");
        return  ERROR;
    }

    /*
      Calculate the unit vector
    */
    uni->x = vec->x / mag;
    uni->y = vec->y / mag;
    uni->z = vec->z / mag;
      
    return SUCCESS;
}
