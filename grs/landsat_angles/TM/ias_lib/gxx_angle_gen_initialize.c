/****************************************************************************/
/**
 * @file gxx_angle_gen_initialize.c
 * @brief Initialize the angle generation metadata structure.
 * @ingroup gpsAngleCoef
 */
/****************************************************************************/

/* Standard Library Includes */
#include <stdlib.h>

/* IAS Library Includes */
#include "xxx_LogStatus.h"
#include "gxx_angle_gen_private.h"

/******************************************************************************/
/**
 * @brief Initialize the angle generation metadata structure. 
 * Note: ephem_count must be initialized to the size of the ephemeris and 
 * solar vector structs.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/****************************************************************************/
int gxx_angle_gen_initialize
(
    gxx_angle_gen_metadata_TYPE *metadata //!<[in/out] Angle metadata struct
)           
{
    /* Check the number of points requested */
    if (metadata->ephem_count < 1)
    {
        xxx_LogStatus(PROGRAM, __FILE__,__LINE__,
                      "Invalid ephemeris count");
        return ERROR;
    }

    /* Allocate the ephemeris */
    metadata->ephemeris = (gxx_angle_gen_ephemeris_TYPE *)malloc(
        metadata->ephem_count * sizeof(gxx_angle_gen_ephemeris_TYPE));
    if (!metadata->ephemeris)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
                      "Error allocating ephemeris");
        return ERROR;
    }

    /* Allocate the solar vector */
    metadata->solar_vector = (gxx_angle_gen_ephemeris_TYPE *)malloc(
        metadata->ephem_count * sizeof(gxx_angle_gen_ephemeris_TYPE));
    if (!metadata->solar_vector)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
                      "Error allocating solar vector");
        free(metadata->ephemeris);
        metadata->ephemeris = NULL;
        return ERROR;
    }

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Free up the angle generation metadata structure.
 *
 */
/****************************************************************************/
void gxx_angle_gen_free
(
    gxx_angle_gen_metadata_TYPE *metadata //!<[in] Metadata structure
)           
{
    free(metadata->ephemeris);
    metadata->ephemeris = NULL;

    free(metadata->solar_vector);
    metadata->solar_vector = NULL;
}
