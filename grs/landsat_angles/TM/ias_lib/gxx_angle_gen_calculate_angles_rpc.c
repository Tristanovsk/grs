/******************************************************************************/
/**
 * @file gxx_angle_gen_calculate_angle_rpc.c
 * @brief Calculates the satellite viewing and solar angles using rational
 * polynomial coefficients
 * @ingroup gpsAngleCoef
 */
/******************************************************************************/

/* Standard Library Includes */
#include <math.h>

/* IAS Library Includes */
#include "xxx_LogStatus.h"
#include "gxx_angle_gen_distro.h"
#include "gxx_angle_gen_private.h"

/******************************************************************************/
/**
 * @brief Calculates the individual vector value used to evaluate the rational
 * polynomial coefficient
 */
/*****************************************************************************/
static void calculate_rpc_vector_value
(
    double l1t_line,           //!<[in] L1T line coordinate
    double l1t_samp,           //!<[in] L1T sample coordinate 
    double l1r_line,           //!<[in] L1R line coordinate 
    double l1r_samp,           //!<[in] L1R sample coordinate 
    double mean_offset,        //!<[in] Vector mean offset 
    double height,             //!<[in] Input height 
    const double *numerator,   //!<[in] Numerator values pointer 
    const double *denominator, //!<[in] Denominator values pointer 
    double *output_value       //!<[out] Output vector value 
)
{
    double equation_num;    /* Equation numerator */
    double equation_den;    /* Equation denominator */

    /* Calculate the numerator */
    equation_num = numerator[0] + numerator[1] * l1t_line
        + numerator[2] * l1t_samp
        + numerator[3] * height
        + numerator[4] * l1r_line
        + numerator[5] * l1t_line * l1t_line
        + numerator[6] * l1t_samp * l1t_line 
        + numerator[7] * l1t_samp * l1t_samp
        + numerator[8] * l1r_samp * l1r_line * l1r_line
        + numerator[9] * l1r_line * l1r_line * l1r_line;
       
    /* Calculate the denominator */
    equation_den = 1.0 + denominator[0] * l1t_line
        + denominator[1] * l1t_samp
        + denominator[2] * height
        + denominator[3] * l1r_line 
        + denominator[4] * l1t_line * l1t_line 
        + denominator[5] * l1t_line * l1t_samp
        + denominator[6] * l1t_samp * l1t_samp
        + denominator[7] * l1r_samp * l1r_line * l1r_line
        + denominator[8] * l1r_line * l1r_line * l1r_line;
   
    *output_value = mean_offset + (equation_num / equation_den);
}

/******************************************************************************/
/**
 * @brief Calculates the satellite viewing and solar angles using rational
 * polynomial coefficients
 *
 * Calculates the satellite viewing and solar illumination zenith and 
 * azimuth angles for a specified L1T line/sample and height (from DEM)
 * using the rational polynomial coefficients for the current band. 
 *
 * @return ERROR: Failed to calculate angles
 * @return SUCCESS: Successfully calculated angles
 */
/******************************************************************************/
int gxx_angle_gen_calculate_angles_rpc
(
    const gxx_angle_gen_metadata_TYPE *metadata, //!<[in] Metadata structure 
    double l1t_line,        //!<[in] Output space line coordinate 
    double l1t_samp,        //!<[in] Output space sample coordinate 
    const double *elev,     //!<[in] Pointer to input elevation or NULL if mean
                            // scene height should be used
    int band_index,         //!<[in] Current band index
    double scan_buffer,     //!<[in] Scan buffer
    int subsamp,            //!<[in] Sub sample factor
    gxx_angle_gen_TYPE sat_or_sun_type,     //!<[in] Angle calculation type 
    int *outside_image_flag,//!<[out] Flag indicating return was outside image 
    double *angle          //!<[out] Array containing zenith and azimuth angles
)       
{
    double height;      /* Model height */
    double l1r_line[2]; /* Input space (L1R) line coordinate declared size 2
                           so it can support 2 SCAs */
    double l1r_samp[2]; /* Input space (L1R) sample coordinate declared size
                           2 so it can support 2 SCAs */
    const gxx_angle_gen_band_TYPE *band_ptr;    /* Pointer to current band */
    const gxx_angle_gen_ang_rpc_TYPE *data_ptr; /* Solar or satellite pointer */
    int dir;            /* Scan direction */
    gxx_scan_direction_TYPE scan_dir; /*Scan direction */
    int num_dir;

    /* Initialize the output zenith and azimuth*/
    angle[IAS_ANGLE_GEN_ZENITH_INDEX] = 0.0;
    angle[IAS_ANGLE_GEN_AZIMUTH_INDEX] = 0.0;
    *outside_image_flag = 0;

    /* Setup the band pointer */
    band_ptr = &metadata->band_metadata[band_index];

    /* Set the height to use */
    height = band_ptr->satellite.mean_height;
    if (elev) 
        height = *elev;

    /* Get the scan direction the point falls in. */
    dir = gxx_angle_gen_find_dir(l1t_line, l1t_samp, height, scan_buffer, 
                                 subsamp,
                                 &(metadata->band_metadata[band_index]),
                                 l1r_line, l1r_samp, &num_dir, &scan_dir);

    if (dir == SUCCESS && scan_dir != no_scan_direction)
    {
        dir = 0;
    }
    
    if( num_dir > 1 )
    {
        if( scan_dir == first_scan_direction )
            dir = 0;
        else if( scan_dir == second_scan_direction )
            dir = 1;
        else
            dir = 0;
    }
    else
            dir = 0;

    /* Offset the output space coordinates */
    l1t_line -= band_ptr->satellite.line_terms.l1t_mean_offset;
    l1t_samp -= band_ptr->satellite.samp_terms.l1t_mean_offset;
    height -= band_ptr->satellite.mean_height;

    data_ptr = &band_ptr->solar;
    if (sat_or_sun_type == GXX_ANGLE_GEN_SATELLITE)
    {
        data_ptr = &band_ptr->satellite;
    }

    VECTOR vector;          /* Viewing vector */

    /* Determine the line and sample location using the L1R offset */
    l1r_line[dir] -= data_ptr->line_terms.l1r_mean_offset;
    l1r_samp[dir] -= data_ptr->samp_terms.l1r_mean_offset;

    /* Calculate the rpc vector */ 
    calculate_rpc_vector_value(l1t_line, l1t_samp, l1r_line[dir], 
        l1r_samp[dir], data_ptr->mean_offset.x, height, 
        data_ptr->x_terms.numerator, data_ptr->x_terms.denominator, 
        &vector.x);

    calculate_rpc_vector_value(l1t_line, l1t_samp, l1r_line[dir], 
        l1r_samp[dir], data_ptr->mean_offset.y, height, 
        data_ptr->y_terms.numerator, data_ptr->y_terms.denominator, 
        &vector.y);

    calculate_rpc_vector_value(l1t_line, l1t_samp, l1r_line[dir], 
        l1r_samp[dir], data_ptr->mean_offset.z, height, 
        data_ptr->z_terms.numerator, data_ptr->z_terms.denominator, 
        &vector.z);

    /* Calculate zenith and azimuth angles */
    angle[IAS_ANGLE_GEN_ZENITH_INDEX] += acos(vector.z);
    angle[IAS_ANGLE_GEN_AZIMUTH_INDEX] += atan2(vector.x, vector.y);

    return SUCCESS;
}
