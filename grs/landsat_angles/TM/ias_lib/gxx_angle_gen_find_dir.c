/******************************************************************************/
/**
 * @file gxx_angle_gen_find_dir.c
 * @brief Finds the scan direction with input height as an optional factor.
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
 * @brief Checkes to see if the line is in a buffered scan gap.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int gxx_angle_gen_check_buffered_scan_gap
(
    double scan_buffer,     //!<[in] Scan buffer
    double l1r_l,           //!<[in] Current line location
    const gxx_angle_gen_band_TYPE *eband //!<[in] Enhanced metadata for 
                                         // current band
)
{
    int iscan_number;      /* Scan number associated with L1r line integer */
    double fscan_number;   /* Scan number associated with L1r line fractional */

    if (l1r_l <= scan_buffer || l1r_l >= (eband->l1r_lines-scan_buffer))
        return ERROR;

    fscan_number = l1r_l / (double)eband->lines_per_scan;
    iscan_number = (int)ceil(fscan_number);
    fscan_number = fabs(iscan_number - fscan_number);

    if (fscan_number < scan_buffer) 
        return SUCCESS;

    fscan_number = l1r_l / (double)eband->lines_per_scan;
    iscan_number = (int)floor(fscan_number);
    fscan_number = fabs(fscan_number - iscan_number);

    if (fscan_number < scan_buffer) 
        return SUCCESS;

    return ERROR;
}

/******************************************************************************/
/**
 * @brief Finds the scan direction with input height as an optional factor.
 *
 * @ returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
int gxx_angle_gen_find_dir
(
    double l1t_line,                //!<[in] L1T line
    double l1t_samp,                //!<[in] L1T sample
    double height,                  //!<[in] height
    double scan_buffer,             //!<[in] scan buffer
    int subsamp,                    //!<[in] sub sample factor
    const gxx_angle_gen_band_TYPE *eband, //!<[in/out] metadata current band
    double *l1r_line,               //!<[out] Array of output L1R line numbers
    double *l1r_samp,               //!<[out] Array of output L1R sample 
                                    // numbers
    int *num_dir_found,             //!<[out] Number of directions found
    gxx_scan_direction_TYPE *scan_dir //!<[out] Scan direction found
)
{
    double  l1t_l;                  /* Offset value of L1T line */
    double  l1t_s;                  /* Offset value of L1T sample */
    double  hgt;                    /* Offset value of height */
    double  l1r_l;                  /* Local L1R line */
    double  l1r_s;                  /* Local L1R sample */
    int     dir;                    /* Scan direction */
    int     num_dir = 0;            /* Number of scan directions found */
    int     scan_number;            /* Scan number associated with L1r line */
    int     l1r_dir;                /* Direction associated with L1r line */

    *scan_dir = no_scan_direction;
    *num_dir_found = 0;

    l1r_line[0] = l1r_line[1] = 0.0;
    l1r_samp[0] = l1r_samp[1] = 0.0;

    for (dir=0; dir<eband->number_scan_dirs; dir++)
    {
        l1t_l = l1t_line - eband->scan_metadata[dir].line_terms.l1t_mean_offset;
        l1t_s = l1t_samp - eband->scan_metadata[dir].samp_terms.l1t_mean_offset;

        if (height != 0)  /* Factor in height. */
        {
            hgt = height - eband->scan_metadata[dir].mean_height;
            l1r_l   = (eband->scan_metadata[dir].line_terms.numerator[0] 
                    + eband->scan_metadata[dir].line_terms.numerator[1]
                    * l1t_l 
                    + eband->scan_metadata[dir].line_terms.numerator[2] * l1t_s 
                    + eband->scan_metadata[dir].line_terms.numerator[3] * hgt 
                    + eband->scan_metadata[dir].line_terms.numerator[4] * l1t_l
                    * l1t_s) / (1.0
                    + eband->scan_metadata[dir].line_terms.denominator[0]
                    * l1t_l 
                    + eband->scan_metadata[dir].line_terms.denominator[1]
                    * l1t_s
                    + eband->scan_metadata[dir].line_terms.denominator[2] * hgt 
                    + eband->scan_metadata[dir].line_terms.denominator[3]
                    * l1t_l * l1t_s)
                    + eband->scan_metadata[dir].line_terms.l1r_mean_offset;
            l1r_s   = (eband->scan_metadata[dir].samp_terms.numerator[0] 
                    + eband->scan_metadata[dir].samp_terms.numerator[1] 
                    * l1t_l + eband->scan_metadata[dir].samp_terms.numerator[2]
                    * l1t_s 
                    + eband->scan_metadata[dir].samp_terms.numerator[3]*hgt 
                    + eband->scan_metadata[dir].samp_terms.numerator[4] 
                    * l1t_l * l1t_s) / (1.0 
                    + eband->scan_metadata[dir].samp_terms.denominator[0]
                    * l1t_l 
                    + eband->scan_metadata[dir].samp_terms.denominator[1]
                    * l1t_s 
                    + eband->scan_metadata[dir].samp_terms.denominator[2] * hgt 
                    + eband->scan_metadata[dir].samp_terms.denominator[3]
                    * l1t_l * l1t_s)
                    + eband->scan_metadata[dir].samp_terms.l1r_mean_offset;
        }
        else  /* Does not factor in height. */
        {
            l1r_l   = (eband->scan_metadata[dir].line_terms.numerator[0] 
                    + eband->scan_metadata[dir].line_terms.numerator[1] 
                    * l1t_l + (eband->scan_metadata[dir].line_terms.numerator[2] 
                    + eband->scan_metadata[dir].line_terms.numerator[4] * l1t_l)
                    * l1t_s)
                    / (1.0 + eband->scan_metadata[dir].line_terms.denominator[0]
                    * l1t_l 
                    + (eband->scan_metadata[dir].line_terms.denominator[1]
                    + eband->scan_metadata[dir].line_terms.denominator[3]
                    * l1t_l) * l1t_s) 
                    + eband->scan_metadata[dir].line_terms.l1r_mean_offset;
            l1r_s   = (eband->scan_metadata[dir].samp_terms.numerator[0] 
                    + eband->scan_metadata[dir].samp_terms.numerator[1] 
                    * l1t_l + (eband->scan_metadata[dir].samp_terms.numerator[2] 
                    + eband->scan_metadata[dir].samp_terms.numerator[4]
                    * l1t_l) * l1t_s)
                    / (1.0 + eband->scan_metadata[dir].samp_terms.denominator[0]
                    * l1t_l 
                    + (eband->scan_metadata[dir].samp_terms.denominator[1]
                    + eband->scan_metadata[dir].samp_terms.denominator[3]
                    * l1t_l) * l1t_s) 
                    + eband->scan_metadata[dir].samp_terms.l1r_mean_offset;
        }

        if (l1r_l > 0 && l1r_l < eband->l1r_lines && l1r_s > 0 
            && l1r_s < eband->l1r_samps)
        {
            scan_number = (int)( l1r_l / eband->lines_per_scan );
            if (scan_number % 2 == 0) 
                l1r_dir = 0;
            else
                l1r_dir = 1;

            if (l1r_dir == dir)
            {
                l1r_line[num_dir] = l1r_l;
                l1r_samp[num_dir] = l1r_s;
                num_dir++;
                /* This might be dangerous, assumes counter matches 
                   scan_direction_type enumerated type. */
                *scan_dir = dir;
            }
            else if (subsamp && *scan_dir == no_scan_direction)
            {
                if (gxx_angle_gen_check_buffered_scan_gap(scan_buffer, l1r_l, 
                    eband) == SUCCESS)
                {
                    l1r_line[num_dir] = l1r_l;
                    l1r_samp[num_dir] = l1r_s;
                    num_dir++;
                    /* This might be dangerous, assumes counter matches 
                       scan_direction_type enumerated type. */
                    *scan_dir = dir;
                }
            }
        }
    }

    *num_dir_found = num_dir;
    return(SUCCESS);
}
