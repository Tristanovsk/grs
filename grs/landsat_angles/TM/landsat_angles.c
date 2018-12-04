#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <string.h>
#include <math.h>

#include "gxx_angle_gen_distro.h"
#include "xxx_Band.h"
#include "xxx_LogStatus.h"
#include "xxx_Sensor.h"
#include "gxx_sensor.h"
#include "landsat_angles.h"

/* Local defines */
#define SCALED_R2D 4500.0 / atan(1.0)

/**************************************************************************
NAME: landsat_angles (main)

PURPOSE: Calculates satellite viewing and solar illumination zenith and 
         azimuth angles for each L1T pixel (with optional sub-sampling), 
         and writes the angles out to binary flat image files with 
         ENVI headers.      

RETURNS: EXIT_SUCCESS or EXIT_FAILURE
***************************************************************************/
int main
(
    int argc,
    char *argv[]
)
{
    int subsamp = 1;                    /* Sub-sampling factor */
    double scan_buffer = 0.0;           /* Scan buffering */
    gxx_angle_gen_metadata_TYPE metadata; /* Angle metadata structure */
    char sat_ang_filename[STRLEN];      /* Satellite output filename */
    char sun_ang_filename[STRLEN];      /* Solar output filename */
    char acf_name[STRLEN];              /* Angle coefficient filename */
    char fullpath[STRLEN];               /* Absolute path of angle file */
    char sensor_type[STRLEN];           /* Sensor Type */
    char msg[STRLEN];                   /* Error messages */
    int index;                              /* Loop variable */

    /* If there are no arguments, generate a usage message */
    if (argc < 2)
    {
        printf("Usage: landsat_angles <angle_coefficient_file> "
            "<optional parameters>\n\n"
            "Optional parameters:\n\t"
                "-s <subsample factor>\n\t"
                "-b <scan buffer factor>\n");
        exit(EXIT_FAILURE);
    }

    /* Get the angle coefficient file name. */
    strcpy(acf_name, argv[1]);
    printf("angle coefficient filename is %s.\n", acf_name);

    /* change harmel, jan 2017, to keep full path name */
    strcpy(fullpath, acf_name);
    strcpy(fullpath, dirname(fullpath));
    strcat(fullpath,"/angle");

    /* Read the optional input parameters.  */
    for (index = 2; index < argc; index++)
    {
        if (strcmp(argv[index], "-s") == 0)
        {
            subsamp = atoi(argv[index + 1]);
            printf("subsamp = %d\n", subsamp);
            if (subsamp < 1)
            {
                xxx_LogStatus(PROGRAM, __FILE__, __LINE__,
                              "Error:  Subsamp must be greater than 0.");
                exit(EXIT_FAILURE);
            }
        }
        if (strcmp(argv[index], "-b") == 0)
        {
            scan_buffer = atoi(argv[index+1]);
            printf("scan_buffer = %f\n", scan_buffer);
        }
    }

    /* Read the angle coefficient file. */
    if (gxx_angle_gen_read_ang(acf_name, &metadata))
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__,
            "Error reading the angle coefficient file.");
	/*exit(EXIT_FAILURE);*/
    }

    /* Get the sat/sensor from the ACF. */
    strcpy(sensor_type, metadata.spacecraft_id);
    if (xxx_initialize_sensor_type(sensor_type) == IAS_SENSOR_UNKNOWN)
    {
        sprintf(msg, "Invalid SENSOR_TYPE string: %s", sensor_type);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        exit(EXIT_FAILURE);
    }

    /* Process the angles for each band */
    for (index = 0; index < (int)metadata.num_bands; index++)
    {

     
        /* change harmel, jan 2017, to skip temperature bands */
        /*if (index == 5 || index == 6) continue;*/

        int nlines, nsamps;             /* Number of lines and samples */
        short *sat_zn = NULL;           /* Satellite zenith angle */
        short *sat_az = NULL;           /* Satellite azimuth angle */
        short *sun_zn = NULL;           /* Solar zenith angle */
        short *sun_az = NULL;           /* Solar azimuth angle */
        double satang[2];               /* Satellite viewing angles */
        double sunang[2];               /* Solar viewing angles */
        int count;                      /* Output sample counter */
        int line;                       /* Line index */
        int samp;                       /* Sample index */
        int outside_image;              /* Return was outside image */
        double *height = NULL;          /* Height to evaluate */
        angle_frame_TYPE frame;         /* Output image frame info. */
        short max_sat_zn, min_sat_zn;
        short max_sat_az, min_sat_az;
        size_t angle_size;              /* angle array size */

        /* Calculate size of subsampled output image */
        nlines = (metadata.band_metadata[index].l1t_lines - 1) / subsamp + 1;
        nsamps = (metadata.band_metadata[index].l1t_samps - 1) / subsamp + 1;

       /* Allocate the output buffers */
        angle_size = nlines*nsamps*sizeof(short);
        sat_zn = malloc(angle_size);
        sat_az = malloc(angle_size);
        sun_zn = malloc(angle_size);
        sun_az = malloc(angle_size);
        if (sat_zn == NULL || sat_az == NULL || sun_zn == NULL ||
            sun_az == NULL)
        {
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, "Error allocating "
                          "satellite and solar angle arrays.");
            free(sat_zn);
            free(sat_az);
            free(sun_zn);
            free(sun_az);
            gxx_angle_gen_free(&metadata);
            exit (EXIT_FAILURE);
        }
        count = 0;

        /* Establish the output image file frame */
        frame.band_number =
                 xxx_get_user_band(metadata.band_metadata[index].band_number);
        frame.num_lines = nlines;
        frame.num_samps = nsamps;
        frame.projection.code = metadata.projection.code;
        frame.projection.zone = metadata.projection.zone;
        frame.pixel_size = metadata.band_metadata[index].pixel_size;
        frame.ul_corner.x = metadata.corners.upleft.x;
        frame.ul_corner.y = metadata.corners.upleft.y;

        /* Loop through the L1T lines and samples */
        for (line = 0; line < metadata.band_metadata[index].l1t_lines; 
            line += subsamp)
        {
            for (samp = 0; samp < metadata.band_metadata[index].l1t_samps; 
                samp += subsamp)
            {
                if (gxx_angle_gen_calculate_angles_rpc(&metadata,
                                                       (double)line,
                                                       (double)samp, height,
                                                       index, scan_buffer,
                                                       subsamp,
                                                       GXX_ANGLE_GEN_SATELLITE,
                                                       &outside_image, satang)
                    != SUCCESS)
                {
                    sprintf(msg,"Error evaluating view angles in band %d.",
                            metadata.band_metadata[index].band_number);
                    xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
                }
                if (gxx_angle_gen_calculate_angles_rpc(&metadata,
                                                       (double)line,
                                                       (double)samp, height,
                                                       index, scan_buffer,
                                                       subsamp,
                                                       GXX_ANGLE_GEN_SOLAR,
                                                       &outside_image, sunang)
                    != SUCCESS)
                {
                    sprintf(msg,"Error evaluating solar angles in band %d.",
                            metadata.band_metadata[index].band_number);
                    xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
                }

                /* Calculate satellite angles from vector */
                sat_zn[count] = (short)floor(SCALED_R2D*satang[0] + 0.5);
                sat_az[count] = (short)floor(SCALED_R2D*satang[1] + 0.5);

                /* Calculate solar angles from vector */
                sun_zn[count] = (short)floor(SCALED_R2D*sunang[0] + 0.5);
                sun_az[count] = (short)floor(SCALED_R2D*sunang[1] + 0.5);
                if( count == 0 )
                {
                   max_sat_az = sat_az[count];
                   min_sat_az = sat_az[count];
                   max_sat_zn = sat_zn[count];
                   min_sat_zn = sat_zn[count];
                }
                else
                {
                   if( sat_az[count] > max_sat_az ) max_sat_az = sat_az[count];
                   if( sat_az[count] < min_sat_az ) min_sat_az = sat_az[count];
                   if( sat_zn[count] > max_sat_zn ) max_sat_zn = sat_zn[count];
                   if( sat_zn[count] < min_sat_zn ) min_sat_zn = sat_zn[count];
                }
                count++;
            }
        }

        /* change harmel, jan 2017, to keep full path name */
        strcpy(sat_ang_filename, fullpath);
        strcpy(sun_ang_filename, fullpath);

        /* Write angles to output file */
        if (gxx_angle_gen_write_image(sat_ang_filename, sat_az, sat_zn,
            GXX_ANGLE_GEN_SATELLITE, index, frame.band_number, frame.num_lines,
            frame.num_samps, frame.ul_corner, frame.pixel_size, 
            &frame.projection) != SUCCESS)
        {
            sprintf(msg,"Error writing view angle band for band %d.", 
                metadata.band_metadata[index].band_number);
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        }

        if (gxx_angle_gen_write_image(sun_ang_filename, sun_az, sun_zn,
            GXX_ANGLE_GEN_SOLAR, index, frame.band_number, frame.num_lines,
            frame.num_samps, frame.ul_corner, frame.pixel_size,
            &frame.projection) != SUCCESS)
        {
            sprintf(msg,"Eror writing solar angle band for band %d.",
                metadata.band_metadata[index].band_number);
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        }

        /* Fmaree the angle arrays */
        free(sat_zn);
        free(sat_az);
        free(sun_zn);
        free(sun_az);
    }

    /* Free the ephemeris structure */
    gxx_angle_gen_free(&metadata);

    exit (SUCCESS);
}

