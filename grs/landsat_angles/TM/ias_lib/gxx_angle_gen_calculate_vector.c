/****************************************************************************/
/**
 * @file gxx_angle_gen_calculate_vector.c
 * @brief Calculate the satellite viewing and solar illumination vectors
 * @ingroup gpsAngleCoef
 */
/****************************************************************************/

/* Standard Library Includes */
#include <stdio.h>
#include <math.h>

/* IAS Library Includes */
#include "xxx_LogStatus.h"
#include "gxx_angle_gen_geo_utilities.h"
#include "gxx_angle_gen_private.h"

/******************************************************************************/
/**
 * @brief Converts the L1R line to image time.
 *
 * @return: Image time
 */
/******************************************************************************/
static double gxx_angle_gen_l1r_to_time
(
    const gxx_angle_gen_band_TYPE *band_metadata,//!<[in] Band-specific metadata
    gxx_angle_gen_scan_time_TYPE *scan_time, //!<[in] scan time polynomials
    int dir,                                 //!<[in] Scan direction
    double l1r_line,                         //!<[in] Level 1R line coordinate
    double l1r_samp                          //!<[in] Level 1R sample coordinate
)      
{
    unsigned int count;                 /* Integer counter */
    double time = 0.0;                  /* Time associated with L1R location */
    double time_in_scan = 0.0;          /* Time within scan */
    double l1r_line_sub;
    double l1r_samp_sub;
    double l1r_scale = band_metadata->lines_per_scan / 16.0;

    l1r_line_sub = l1r_line / l1r_scale;
    l1r_samp_sub = l1r_samp / l1r_scale;
    for( count=0 ; count<(unsigned int)scan_time->ncoeff ; count++ )
        time += pow(l1r_line_sub, (double)count) 
        * scan_time->scan_time_poly[dir][count];

    time_in_scan = (l1r_samp_sub / scan_time->mean_eol[dir]) 
        * scan_time->mean_activescan[dir];

    time += time_in_scan;

    return time;
}

/******************************************************************************/
/**
 * @brief Uses the scene corners and pixel size to compute the map 
 * projection X/Y coordinates for the input L1T line/sample
 *
 * @return SUCCESS: Successful completion
 * @return ERROR: Operation failed
 */
/******************************************************************************/
static void gxx_angle_gen_l1t_to_proj
(
    const gxx_angle_gen_metadata_TYPE *metadata, //!<[in] Metadata info
    unsigned int band_index,    //!<[in] Band index
    double l1t_line,            //!<[in] L1T line number
    double l1t_samp,            //!<[in] L1T sample number
    double *projection_x,       //!<[out] Projection X coordinate
    double *projection_y        //!<[out] Projection Y coordinate
)           
{
    /* Convert L1T line/sample to projection X/Y */
    /* Note, this only works for north-up images */
    *projection_x = metadata->corners.upleft.x 
        + metadata->band_metadata[band_index].pixel_size * l1t_samp;
    *projection_y = metadata->corners.upleft.y 
        - metadata->band_metadata[band_index].pixel_size * l1t_line;
}

/******************************************************************************/
/**
 * @brief Calculate the satellite viewing and solar illumination vectors at a 
 * specified L1T line/sample/height location
 *
 * @return SUCCESS: Successful completion
 * @return ERROR: Operation failed
 */
/******************************************************************************/
int gxx_angle_gen_calculate_vector
(
    gxx_angle_gen_metadata_TYPE *metadata,//!<[in] Metadata structure
    double l1t_line,         //!<[in] Current L1T line number
    double l1t_samp,         //!<[in] Current L1T sample number
    double l1r_line,         //!<[in] Current L1R line number
    double l1r_samp,         //!<[in] Current L1R sample number
    double height,           //!<[in] Current L1T height
    unsigned int band_index, //!<[in] Current band index
    unsigned int dir,        //!<[in] scan direction
    gxx_angle_gen_TYPE sat_or_sun_type, //!<[in] Angle type
    VECTOR *view         //!<[out] View vector
)   
{
    double projection_x;               /* Projection X coordinate */
    double projection_y;               /* Projection Y coordinate */
    double latitude;                   /* Latitude */
    double longitude;                  /* Longitude */
    double sample_time;                /* Sample time from epoch */
    VECTOR ecef_vertical[3];       /* Vertical ECEF basis vectors */
    VECTOR ecef_vector;            /* ECEF vector */
    gxx_angle_gen_ephemeris_TYPE *ephemeris;/* Ephemeris pointer */
    char err_msg[STRLEN];

    /* Convert the line/sample to projection coordinates */
    gxx_angle_gen_l1t_to_proj(metadata, band_index, l1t_line, l1t_samp, 
        &projection_x, &projection_y);

    /* Convert the projection coordinates to latitude/longitude */
    if (gxx_angle_gen_initialize_transformation(metadata, projection_x,
        projection_y, &longitude, &latitude) != SUCCESS)
    {
        sprintf(err_msg, "Error converting projection X/Y to lat/long for l1t "
            "line %lf and sample %lf", l1t_line, l1t_samp);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, err_msg);
        return ERROR;
    }

    /* Calculate local vertical basis vectors */
    ecef_vertical[0].z = 0.0;
    ecef_vertical[0].x = -sin(longitude);
    ecef_vertical[0].y = cos(longitude);
    ecef_vertical[1].z = cos(latitude);
    ecef_vertical[2].z = sin(latitude);
    ecef_vertical[1].x = -ecef_vertical[2].z * ecef_vertical[0].y;
    ecef_vertical[1].y = -ecef_vertical[2].z * -ecef_vertical[0].x;
    ecef_vertical[2].x = ecef_vertical[1].z * ecef_vertical[0].y;
    ecef_vertical[2].y = ecef_vertical[1].z * -ecef_vertical[0].x;

    /* Calculate observation time from L1R line number */
    sample_time = gxx_angle_gen_l1r_to_time(
        &metadata->band_metadata[band_index], 
        &metadata->scan_time, dir, l1r_line, l1r_samp);
    
    /* Determine the correct ephemeris */
    ephemeris = metadata->solar_vector;
    if (sat_or_sun_type == GXX_ANGLE_GEN_SATELLITE)
    {
        ephemeris = metadata->ephemeris;
    }

    /* Interpolate the ecef vector */
    if (gxx_angle_gen_interpolate_ephemeris(ephemeris, metadata->ephem_count, 
        sample_time, &ecef_vector) != SUCCESS)
    {
        xxx_LogStatus(PROGRAM, __FILE__,__LINE__, 
                      "Error interpolating the ephemeris");
        return ERROR;
    }  

    VECTOR ground_position;        /* Ground ECEF vector */
    /* Adjust the ecef vector for the satellite */
    if (sat_or_sun_type == GXX_ANGLE_GEN_SATELLITE)
    {
        /* Calculate the ground position vector */
        gxx_angle_gen_geo_to_ecef(metadata, latitude, longitude, height, 
            &ground_position);

        ecef_vector.x -= ground_position.x;
        ecef_vector.y -= ground_position.y;
        ecef_vector.z -= ground_position.z;
    } 

    /* Calculate the ECEF vector */
    if (gxx_unit(&ecef_vector, &ground_position) != SUCCESS)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__,
                      "Error normalizing the ECEF vector");
        return ERROR;
    }

    /* Calculate the view vector */
    view->x = gxx_dot(&ground_position, &ecef_vertical[0]);
    view->y = gxx_dot(&ground_position, &ecef_vertical[1]);
    view->z = gxx_dot(&ground_position, &ecef_vertical[2]);

    return SUCCESS;
}
