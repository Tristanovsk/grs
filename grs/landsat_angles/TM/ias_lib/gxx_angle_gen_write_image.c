/****************************************************************************/
/**
 * @file gxx_angle_gen_write_image.c
 * @brief Write angle images as a binary file.
 * @ingroup gpsAngleCoef
 */
/****************************************************************************/

/* Standard Library Includes */
#include <stdio.h>
#include <string.h>
#include <limits.h>

/* IAS Library Includes */
#include "xxx_LogStatus.h"       
#include "gxx_angle_gen_distro.h"

/******************************************************************************/
/**
 * @brief Write out the header file that envi needs.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int write_envi_header
(
    char *ang_filename,                    /* I: Envi file name */
    int num_lines,                         /* I: Number of image lines */
    int num_samps,                         /* I: Number of image samples */
    DBL_XY ul_corner,                      /* I: Image upper left corner */
    double pixel_size,                     /* I: Image pixel size */
    const gxx_projection_TYPE *projection  /* I: Image framing information */
)
{
    char msg[STRLEN];
    FILE *ofp;               /* Output file pointer. */

    /* Create the output header file name */
    strcat(ang_filename, ".hdr");
    printf("Writing view angle band header file %s.\n", ang_filename);

    /* Open the output header file */
    if ((ofp = fopen(ang_filename, "w")) == NULL)
    {
        sprintf(msg,"Opening output ENVI header file %s.", ang_filename);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Write the output header file */
    fprintf(ofp, "ENVI\n");
    fprintf(ofp, "description = {\n");
    fprintf(ofp, "  View Angle Band File for %s.}\n", ang_filename);
    fprintf(ofp, "samples = %d\n", num_samps);
    fprintf(ofp, "lines   = %d\n", num_lines);
    fprintf(ofp, "bands   = 2\n");
    fprintf(ofp, "header offset = 0\n");
    fprintf(ofp, "file type = ENVI Standard\n");
    fprintf(ofp, "data type = 2\n");
    fprintf(ofp, "interleave = bsq\n");
    fprintf(ofp, "sensor type = Unknown\n");
    fprintf(ofp, "byte order = 0\n");
    if (projection->code == 1)
    {
        fprintf(ofp, "map info = {UTM, 1.500, 1.500, %6.3lf, %6.3lf, %6.3lf,"
                     " %6.3lf, %d, North, WGS-84, units=Meters}\n",
                     ul_corner.x, ul_corner.y, pixel_size, pixel_size, 
                     projection->zone);
    }
    else if (projection->code == 6 || projection->code == 3)
    {
        double semi_major = projection->projprms[0];
        double semi_minor = projection->projprms[1];
        double false_easting = projection->projprms[6];
        double false_northing = projection->projprms[7];
        double latitude_origin;
        double central_meridian_longitude;
        double standard_parallel_1_lat;
        double standard_parallel_2_lat;

        if (gxx_dmsdeg(projection->projprms[5], &latitude_origin, "LAT")
            != SUCCESS)
        {
            strcpy(msg, "Error converting the latitude of the projection "
                "origin from DMS to degrees");
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
            return ERROR;
        }

        if (gxx_dmsdeg(projection->projprms[4], &central_meridian_longitude, 
            "LON") != SUCCESS)
        {
            strcpy(msg, "Error converting the longitude of the central "
                "meridian from DMS to degrees");
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
            return ERROR;
        }

        if (gxx_dmsdeg(projection->projprms[2], &standard_parallel_1_lat, 
            "LAT") != SUCCESS)
        {
            strcpy(msg, "Error converting the latitude of the first standard "
                "parallel from DMS to degree");
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
            return ERROR;
        }

        if (gxx_dmsdeg(projection->projprms[3], &standard_parallel_2_lat, 
            "LAT") != SUCCESS)
        {
            strcpy(msg, "Error converting the latitude of the second standard "
                "parallel from DMS to degree");
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
            return ERROR;
        }

        if (projection->code == 6)
        {
            fprintf(ofp, "map info = {Polar Stereographic, 1.500, 1.500, "
                    "%6.3lf, %6.3lf, %6.3lf, %6.3lf, WGS-84, units=Meters}\n",
                     ul_corner.x, ul_corner.y, pixel_size, pixel_size);
            fprintf(ofp, "projection info = {31, %lf, %lf, %lf, "
                     "%lf, %lf, %lf, %lf, %lf, WGS-84, " 
                     "Polar Stereographic, units=Meters}\n",
                     semi_major, semi_minor, latitude_origin,
                     central_meridian_longitude, false_easting, false_northing,
                     standard_parallel_1_lat, standard_parallel_2_lat);
        }
        else
        {
            fprintf(ofp, "map info = {Albers Conical Equal Area, 1.00, 1.00, "
                     "%6.3lf, %6.3lf, %6.3lf, %6.3lf, WGS-84, units=Meters}\n",
                     ul_corner.x, ul_corner.y, pixel_size, pixel_size);

            fprintf(ofp, "projection info = {9, %lf, %lf, %lf, "
                     "%lf, %lf, %lf, %lf, %lf, WGS-84, "
                     "Albers Conical Equal Area, units=Meters}\n",
                     semi_major, semi_minor, latitude_origin,
                     central_meridian_longitude, false_easting, false_northing,
                     standard_parallel_1_lat, standard_parallel_2_lat);
        }
    }
    else
    {
        strcpy(msg,"Invalid map projection, not UTM, LIMA PS, or Albers");
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }
    fprintf(ofp, "wavelength units = Unknown\n");

    /* Close the output header file */
    fclose(ofp);

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Write out the solar illumination or the satellite viewing 
 * angle images as a binary file.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/*****************************************************************************/
int gxx_angle_gen_write_image
(
    const char *image_filename, /* I: Image file name */
    const short *azimuth,       /* I: Array of azimuth angles */
    const short *zenith,        /* I: Array of zenith angles */
    gxx_angle_gen_TYPE sat_or_sun_type, /* I: Image type to write */
    int band_index,             /* I: Output band index */
    int band_number,            /* I: Name of the band */
    int num_lines,              /* I: Number of image lines */
    int num_samps,              /* I: Number of image samples */
    DBL_XY ul_corner,           /* I: Image upper left corner */
    double pixel_size,          /* I: Image pixel size */
    const gxx_projection_TYPE *projection /* I: Image framing information */
)        
{
    FILE *output_file;           /* Output file pointer */
    char ang_filename[PATH_MAX]; /* Output angle file name */
    int count;                   /* Total number of samples */
    int status;                  /* Status placeholder */
    char msg[STRLEN];

    /* Check the input band index */
    if (band_index < 0 || band_index > MAXBANDS)
    {
        sprintf(msg, "Invalid band index %d", band_index);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Calculate the total number of samples */
    count = num_lines * num_samps;
    
    if (sat_or_sun_type == GXX_ANGLE_GEN_SATELLITE)
    {
        /* Construct satellite view angle output file name */
        status = snprintf(ang_filename, sizeof(ang_filename), 
            "%s_sensor_B%02d.img", image_filename, band_number);
    }
    else
    {
        /* Construct solar angle output file name */
        status = snprintf(ang_filename, sizeof(ang_filename), 
            "%s_solar_B%02d.img", image_filename, band_number);
    }

    if (status < 0 || (unsigned int)status >= sizeof(ang_filename))
    {
        sprintf(msg, "Error creating the image filename from filename %s", 
            image_filename);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Open output file */
    output_file = fopen(ang_filename, "wb");
    if (!output_file)
    {
        sprintf(msg, "Error opening output view angle band file %s", 
            ang_filename);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Write the Zenith Angle layer */
    if (fwrite(zenith, sizeof(short), count, output_file) 
        != (unsigned int)count)
    {
        sprintf(msg, "Error writing zenith angle layer for band %d to file %s", 
            band_number, ang_filename);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        fclose(output_file);
        return ERROR;
    }

    /* Write the Azimuth Angle layer */
    if (fwrite(azimuth, sizeof(short), count, output_file) 
        != (unsigned int)count)
    {
        sprintf(msg, "Error writing azimuth angle layer for band %d to file %s",
                band_number, ang_filename);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        fclose(output_file);
        return ERROR;
    }

    /* Close the output image file */
    if (fclose(output_file) != 0)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
            "Error closing the image file");
        return ERROR;
    }
   
    /* Write the envi header */
    if (write_envi_header(ang_filename, num_lines, num_samps, ul_corner,
                          pixel_size, projection) != SUCCESS)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__,
            "Error creating the image header");
        return ERROR;
    }

    return SUCCESS;
}
