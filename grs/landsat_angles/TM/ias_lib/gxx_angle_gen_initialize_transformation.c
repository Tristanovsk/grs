/****************************************************************************/
/**
 * @file gxx_angle_gen_initialize_transformation.c
 * @brief Initialize the angle generation metadata transformation
 * @ingroup gpsAngleCoef
 */
/****************************************************************************/

/* Standard Library Includes */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* IAS Library Includes */
#include "xxx_LogStatus.h"
#include "gxx_proj.h"
#include "gxx_angle_gen_private.h"

/******************************************************************************/
/**
 * @brief Use projtran to transform a set of coordinates to a new projection.
 *
 * @return integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int gxx_angle_gen_transform_projection
(
    gxx_projection_TYPE in_projection,  /* I: input projection parms */
    gxx_projection_TYPE out_projection, /* I: output projection parms */
    double in_x,                        /* I: x coordinate */
    double in_y,                        /* I: y coordinate */
    double *out_x,                      /* O: x coordinate */
    double *out_y                       /* O: y coordinate */
)
{
    int status;          /* Projection package return code */
    double temp_in_x = in_x;
    double temp_in_y = in_y;
    int projection_units_in;
    int projection_units_out;
    int i;
    char msg[STRLEN];

    for (i = 0; i<2; i++)
    {
        if (i)
        {
            /* Set the input projection units. */
            if (gxx_get_units(in_projection.units, &projection_units_in) 
                != SUCCESS)
            {
                sprintf(msg, "Getting the projection units for %s", 
                        in_projection.units);
                xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
                return ERROR;
            }
        }
        else
        {
            /* Set the output projection units. */
            if (gxx_get_units(out_projection.units, &projection_units_out) 
                != SUCCESS)
            {
                sprintf(msg, "Getting the projection units for %s", 
                        out_projection.units);
                xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
                return ERROR;
            }
        }
    }

    /* Apply transformation */
    status = gxx_projtran(&in_projection.code, &projection_units_in, 
        &in_projection.zone, in_projection.projprms, &in_projection.spheroid,
        &out_projection.code, &projection_units_out, &out_projection.zone,
        out_projection.projprms, &out_projection.spheroid, &temp_in_x, 
        &temp_in_y, out_x, out_y);
    return status;

}

/****************************************************************************/
/**
 * @brief A routine to set the members into an IAS projection.
 */
/****************************************************************************/
static void gxx_angle_gen_set_projection
(
    int proj_code,            /* I: input projection code */
    int zone,                 /* I: input zone */
    char *units,              /* I: input units */
    int spheroid,             /* I: input spheroid */
    const double *parms,      /* I: input projection parameters */
    gxx_projection_TYPE *proj /* I: target projection structure */
)
{
    int i;

    /* Set the projection code, zone, units, and spheroid to the GCTP proj */
    proj->code = proj_code;
    proj->zone = zone;
    proj->spheroid = spheroid;
    strcpy(proj->units, units);

    /* Set the projection parameters */
    for (i = 0; i < PROJPRMS_SIZE; i++)
        proj->projprms[i] = parms[i];
}



/******************************************************************************/
/**
 * @brief Initialize the angle generation metadata transformation
 *
 * @return integer (SUCCESS or ERROR)
 */
/******************************************************************************/
int gxx_angle_gen_initialize_transformation
(
    gxx_angle_gen_metadata_TYPE *metadata,  /* I/O: Angle metadata struct */
    double in_x,                            /* I: X coordinate */
    double in_y,                            /* I: y coordinate */
    double *out_x,                          /* O: x coordinate */
    double *out_y                           /* O: y coordinate */
)           
{
    int index;                              /* Loop variable */
    double parameters[PROJPRMS_SIZE];       /* Geo parameter array */
    gxx_projection_TYPE in_projection;      /* Source projection */
    gxx_projection_TYPE out_projection;     /* Destination projection */
    char out_units[STRLEN] = "RADIANS";

    /* Initialize the geographic projection parameters */
    for (index = 0; index < PROJPRMS_SIZE; index++)
    { 
        parameters[index] = 0.0;
    }

    /* Load the projection structures */
    gxx_angle_gen_set_projection(metadata->projection.code, 
        metadata->projection.zone, metadata->units, 
        metadata->projection.spheroid, metadata->projection.projprms, 
        &in_projection);
    gxx_angle_gen_set_projection(GEO, NULLZONE, out_units, 
        metadata->projection.spheroid,
        parameters, &out_projection);

    /* Initialize the projection transformation */
    if (gxx_angle_gen_transform_projection(in_projection, out_projection,
        in_x, in_y, out_x, out_y) != SUCCESS)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__,
                      "Error calculating the projection transformation");
        return ERROR;
    }

    return SUCCESS;
}
