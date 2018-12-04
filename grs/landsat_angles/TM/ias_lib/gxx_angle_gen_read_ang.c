/****************************************************************************/
/**
 * @file gxx_angle_gen_read_ang.c
 * @brief Read the ANG metadata file.
 * @ingroup gpsAngleCoef
 */
/****************************************************************************/

/* Standard Library Includes */
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* IAS Library Includes */
#include "xxx_LogStatus.h"
#include "xxx_Band.h"
#include "xxx_string_utils.h"
#include "xxx_ODL.h"      
#include "gxx_proj.h"
#include "gxx_angle_gen_distro.h"
#include "gxx_angle_gen_private.h"

/* Local Defines */
#define MAX_FIELD_NAME  32 /* Maximum number of characters in an ODL field */
#define SECS_TOLERANCE .0000015 /* Epoch seconds tolerance */ 

/******************************************************************************/
/**
 * @brief Retrieves the ODL field and does error checking on that field.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int get_odl_field 
(
    int memory_size,       /* I: Size of memory to retrieve */
    xxx_ValType_TYPE type, /* I: Type of ODL field to retreive */
    OBJDESC *tree,         /* I: ODL tree to retrieve the field */
    char *group,           /* I: ODL group where field is located */
    char *field,           /* I: ODL field name */
    int expected_count,    /* I: Expected values to retrieve */
    void *address          /* O: Pointer to storage address */
)
{
    int count;  /* Number of values retrieved */
    char msg[STRLEN];
    char return_msg[STRLEN];

    /* Retrieve the ODL field */
    if (xxx_GetODLField(address, memory_size, type, tree, group, field, 
        return_msg, &count) != SUCCESS)
    {
        sprintf(msg, "Error reading from the ODL tree: %s", return_msg);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }
    
    /* Compare the actual and expected count */
    if (count != expected_count)
    {
        sprintf(msg, "Retrieved %d values expected values is %d for %s",
                     count, expected_count, field);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Constructs the ODL field name and then retrieves the odl field.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int construct_and_get_odl_field
(
    int memory_size,           /* I: Size of memory to retrieve */
    xxx_ValType_TYPE type,     /* I: Type of ODL field to retreive */
    OBJDESC *tree,             /* I: ODL tree to retrieve the field */
    char *group,               /* I: ODL group where field is located */
    char *field,               /* I: ODL field name */
    char *field_prepend,       /* I: Value prepended to field */
    int expected_count,        /* I: Expected values to retrieve */
    void *address              /* O: Pointer to storage address */
)
{
    char field_name[MAX_FIELD_NAME]; /* Field name */
    int status;
    char msg[STRLEN]; 

    /* Read the number of L1T samples for this band */
    status = snprintf(field_name, sizeof(field_name), "%s_%s", field_prepend, 
        field);
    if (status < 0 || (unsigned int)status >= sizeof(field_name))
    {
        sprintf(msg, "Constructing the odl field for %s_%s", field_prepend,
            field);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Retrieve the ODL field */
    if (get_odl_field(memory_size, type, tree, group, field_name, 
        expected_count, address) != SUCCESS)
    {
        return ERROR;
    }

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Reads the header group from the ANG metadata file.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int gxx_angle_gen_read_ang_header
(
    OBJDESC *odl_data,                    /* I: Metadata ODL object */
    gxx_angle_gen_metadata_TYPE *metadata /* O: Metadata structure to load */
)       
{
    unsigned int index;             /* Loop index */
    BandNumber band_list[MAXBANDS]; /* Local array for loading band  numbers */
    unsigned int temp_bands[10];    /* Temporary array to hold band list */

    /* Read landsat scene id */
    if (get_odl_field(sizeof(metadata->landsat_scene_id), XXX_String, 
        odl_data, "FILE_HEADER", "LANDSAT_SCENE_ID", 1,
        &metadata->landsat_scene_id) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the WRS Path. */
    if (get_odl_field(sizeof(metadata->wrs_path), XXX_Int, odl_data,
        "FILE_HEADER", "WRS_PATH", 1, &metadata->wrs_path) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the WRS row. */
    if (get_odl_field(sizeof(metadata->wrs_row), XXX_Int, odl_data,
        "FILE_HEADER", "WRS_ROW", 1, &metadata->wrs_row) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the satellite mode. */
    if (get_odl_field(sizeof(metadata->mode), XXX_String, odl_data,
        "FILE_HEADER", "MODE", 1, &metadata->mode) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the first scan direction. */
    if (get_odl_field(sizeof(metadata->first_scan_direction)+1, XXX_String, 
        odl_data, "FILE_HEADER", "FIRST_SCAN_DIRECTION", 1, 
        &metadata->first_scan_direction) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the satellite name */
    if (get_odl_field(sizeof(metadata->spacecraft_id), XXX_String, odl_data,
        "FILE_HEADER", "SPACECRAFT_ID", 1, &metadata->spacecraft_id) != SUCCESS)
    {
        return ERROR;
    }

    /* Read number of bands */
    if (get_odl_field(sizeof(metadata->num_bands), XXX_Int, odl_data, 
        "FILE_HEADER", "NUMBER_OF_BANDS", 1, &metadata->num_bands) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the band list */
    if (get_odl_field(metadata->num_bands * sizeof(metadata->num_bands), 
        XXX_Int, odl_data, "FILE_HEADER", "BAND_LIST", metadata->num_bands, 
        temp_bands) != SUCCESS)
    {
        return ERROR;
    }

    /* Store the band list numbers in the data structure */
    for (index = 0; index < metadata->num_bands; index++) 
    {
        /* Convert the human readable band list into an internal band list. */
        switch (temp_bands[index])
        {
            case 1:
                band_list[index] = 0;
                break;
            case 2:
                band_list[index] = 1;
                break;
            case 3:
                band_list[index] = 2;
                break;
            case 4:
                band_list[index] = 3;
                break;
            case 5:
                band_list[index] = 4;
                break;
            case 6:
                band_list[index] = 9;
                break;
            case 61:
                band_list[index] = 5;
                break;
            case 62:
                band_list[index] = 6;
                break;
            case 7:
                band_list[index] = 7;
                break;
            case 8:
                band_list[index] = 8;
                break;
            default:
                return ERROR;
        }            

        metadata->band_metadata[index].band_number = band_list[index];
        metadata->band_present[index] = TRUE;
    }

    /* Read the number of scan time poly coefficients */
    if (get_odl_field(sizeof(metadata->scan_time.ncoeff), XXX_Int, odl_data,
        "SCAN_TIME_POLY", "SCAN_TIME_POLY_NCOEFF", 1, 
        &metadata->scan_time.ncoeff) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of scan time poly scan directions */
    if (get_odl_field(sizeof(metadata->scan_time.number_scan_dirs), XXX_Int,
        odl_data, "SCAN_TIME_POLY", "SCAN_TIME_POLY_NUMBER_DIRECTIONS", 1,
        &metadata->scan_time.number_scan_dirs) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the per-scan information */
    for (index = 0; index < metadata->scan_time.number_scan_dirs; index++)
    {
        char parm[STRLEN];
        sprintf(parm, "SCAN_TIME%02d_MEAN_ACTIVESCAN", index);

        /* Read the per-scan mean active scan time */
        if (get_odl_field(sizeof(metadata->scan_time.mean_activescan[index]),
            XXX_Double, odl_data, "SCAN_TIME_POLY", parm, 1,
            &metadata->scan_time.mean_activescan[index]) != SUCCESS)
        {
            return ERROR;
        }

        sprintf(parm, "SCAN_TIME%02d_MEAN_EOL", index);
        /* Read the per-scan mean end of line length */
        if (get_odl_field(sizeof(metadata->scan_time.mean_eol[index]),
            XXX_Double, odl_data, "SCAN_TIME_POLY", parm, 1,
            &metadata->scan_time.mean_eol[index]) != SUCCESS)
        {
            return ERROR;
        }

        sprintf(parm, "SCAN_TIME%02d_POLY_COEFF", index);
        /* Read the per-scan polynomial coefficients */
        if (get_odl_field(32,
            XXX_Double, odl_data, "SCAN_TIME_POLY", parm, 
            metadata->scan_time.ncoeff, 
            &metadata->scan_time.scan_time_poly[index]) != SUCCESS)
        {
            return ERROR;
        }
    }

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Reads the projection group from the ANG metadata file.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int gxx_angle_gen_read_ang_projection
(
    OBJDESC *odl_data,                    /* I: Metadata ODL object */
    gxx_angle_gen_metadata_TYPE *metadata /* O: Metadata structure to load */
)        
{
    double coords[2];      /* Local array for loading coordinate pairs  */
    char map_projection[4];/* Local map projection */    
    char ellipsoid[6];     /* Local ellipsoid */
    char msg[STRLEN];

    /* Read the ellipsoid axes*/
    if (get_odl_field(sizeof(coords), XXX_Double, odl_data, "PROJECTION",
        "ELLIPSOID_AXES", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->wgs84_major_axis = coords[0];
    metadata->wgs84_minor_axis = coords[1];

    /* Read the map projection */
    if (get_odl_field(sizeof(map_projection), XXX_String, odl_data, 
        "PROJECTION", "MAP_PROJECTION", 1, map_projection) != SUCCESS)
    {
        return ERROR;
    }

    /* Process the map projection */
    xxx_strtoupper(map_projection);
    if (strcmp(map_projection, "UTM") == 0)
    {
        metadata->projection.code = UTM;
    }
    else if (strcmp(map_projection, "PS") == 0)
    {
        metadata->projection.code = PS;
    }
    else if (strcmp(map_projection, "AEA") == 0)
    {
        metadata->projection.code = AEA;
    }
    else
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
            "MAP_PROJECTION parameter is invalid must be "
            "UTM, PS, or AEA");
        return ERROR;
    }

    /* Read the projection units */
    if (get_odl_field(sizeof(metadata->units), XXX_String, odl_data, 
        "PROJECTION", "PROJECTION_UNITS", 1, metadata->units) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the projection datum */
    if (get_odl_field(sizeof(metadata->datum), XXX_String, odl_data, 
        "PROJECTION", "DATUM", 1, metadata->datum) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the projection spheroid */
    if (get_odl_field(sizeof(ellipsoid), XXX_String, odl_data, "PROJECTION",
        "ELLIPSOID", 1, ellipsoid) != SUCCESS)
    {
        return ERROR;
    }
    
    metadata->projection.spheroid = WGS84_SPHEROID;
    xxx_strtoupper(ellipsoid);
    if (strcmp(ellipsoid, "WGS84") != 0)
    {
        sprintf(msg, "Unexpected ELLIPSOID: %s", ellipsoid);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Read the projection zone only if projection is UTM */
    if (metadata->projection.code == UTM)
    {
        if (get_odl_field(sizeof(&metadata->projection.zone), XXX_Int, 
            odl_data, "PROJECTION", "UTM_ZONE", 1, &metadata->projection.zone) 
            != SUCCESS)
        {
            return ERROR;
        }

        /* Validate the UTM zone */
        if (metadata->projection.zone < 1 || metadata->projection.zone > 60)
        {
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
                          "UTM_ZONE must be (1 - 60)");
            return ERROR;
        }
    }

    /* Read in all the projection parameters */
    if (get_odl_field(sizeof(metadata->projection.projprms), XXX_Double,
        odl_data, "PROJECTION", "PROJECTION_PARAMETERS", PROJPRMS_SIZE,
        metadata->projection.projprms) != SUCCESS)
    {
        return ERROR;
    }

    /* Read in the upper left corner */
    if (get_odl_field(sizeof(coords), XXX_Double, odl_data, "PROJECTION", 
        "UL_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.upleft.x = coords[0];
    metadata->corners.upleft.y = coords[1];

    /* Read in the upper right corner */
    if (get_odl_field(sizeof(coords), XXX_Double, odl_data, "PROJECTION", 
        "UR_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.upright.x = coords[0];
    metadata->corners.upright.y = coords[1];

    /* Read in the lower left corner */
    if (get_odl_field(sizeof(coords), XXX_Double, odl_data,
        "PROJECTION", "LL_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.loleft.x = coords[0];
    metadata->corners.loleft.y = coords[1];

    /* Read in the lower right corner */
    if (get_odl_field(sizeof(coords), XXX_Double, odl_data,
        "PROJECTION", "LR_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.loright.x = coords[0];
    metadata->corners.loright.y = coords[1];

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Reads the ephemeris group from the ANG metadata file.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
static int gxx_angle_gen_read_ang_ephemeris
(
    OBJDESC *odl_data,                    /* I: Metadata ODL object */
    gxx_angle_gen_metadata_TYPE *metadata /* O: Metadata structure to load */
)        
{
    unsigned int index; /* Loop index */
    int year;           /* Year of acquisition */
    int day;            /* Day of acquisition */
    double seconds;     /* Seconds of acquisition */
    int solar_points;   /* Number of solar points */
    size_t dbuffer_size;/* Size of dbuffer array */
    double *dbuffer;    /* Local array for loading values */
    char msg[STRLEN];

    /* Read the epoch year */
    if (get_odl_field(sizeof(year), XXX_Int, odl_data, "EPHEMERIS",
        "EPHEMERIS_EPOCH_YEAR", 1, &year) != SUCCESS)
    {
        return ERROR;
    }
    metadata->ephem_epoch_time.year = year;
    
    /* Ephemeris epoch time - day of year */
    if (get_odl_field(sizeof(day), XXX_Int, odl_data, "EPHEMERIS",
        "EPHEMERIS_EPOCH_DAY", 1, &day) != SUCCESS)
    {
        return ERROR;
    }
    metadata->ephem_epoch_time.day = day;

    /* Ephemeris epoch time - second of day */
    if (get_odl_field(sizeof(metadata->ephem_epoch_time.seconds), 
        XXX_Double, odl_data, "EPHEMERIS", "EPHEMERIS_EPOCH_SECONDS", 1,
        &metadata->ephem_epoch_time.seconds) != SUCCESS)
    {
        return ERROR;
    }

    /* Number of ephemeris points */
    if (get_odl_field(sizeof(metadata->ephem_count), XXX_Int, odl_data, 
        "EPHEMERIS", "NUMBER_OF_POINTS", 1, &metadata->ephem_count) != SUCCESS)
    {
        return ERROR;
    }

    /* Allocate space for the input buffer */
    dbuffer_size = metadata->ephem_count * sizeof(double);
    dbuffer = (double *)malloc(dbuffer_size);
    if (!dbuffer)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
                      "Allocating ephemeris buffer");
        return ERROR;
    }
    
    /* Initialize the metadata structure */
    if (gxx_angle_gen_initialize(metadata) != SUCCESS)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, 
                      "Allocating ephemeris structures");
        free(dbuffer);
        return ERROR;
    }

    /* Ephemeris sample time offsets from epoch */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_TIME", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    { 
        metadata->ephemeris[index].sample_time = dbuffer[index];
    }

    /* Ephemeris ECEF X coordinates */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_ECEF_X", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->ephemeris[index].position.x = dbuffer[index];
    }

    /* Ephemeris ECEF Y coordinates */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_ECEF_Y", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    {
        metadata->ephemeris[index].position.y = dbuffer[index];
    }

    /* Ephemeris ECEF Z coordinates */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_ECEF_Z", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->ephemeris[index].position.z = dbuffer[index];
    }

    /* Solar vector data group */
    /* Read the Solar Vector epoch year */
    if (get_odl_field(sizeof(year), XXX_Int, odl_data, "SOLAR_VECTOR", 
        "SOLAR_EPOCH_YEAR", 1, &year) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    if ((int)metadata->ephem_epoch_time.year != year)
    {
        sprintf(msg, "Solar epoch year %d and Ephemeris epoch year %d are not"
            "equal", year, (int)metadata->ephem_epoch_time.year);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        free(dbuffer);
        return ERROR;
    }
    
    /* Solar vector epoch time - day of year */
    if (get_odl_field(sizeof(day), XXX_Int, odl_data, "SOLAR_VECTOR",
        "SOLAR_EPOCH_DAY", 1, &day) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }
    
    if ((int)metadata->ephem_epoch_time.day != day)
    {
        sprintf(msg, "Solar epoch day %d and Ephemeris epoch day %d are not"
            "equal", day, (int)metadata->ephem_epoch_time.day);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        free(dbuffer);
        return ERROR;
    }

    /* Solar vector epoch time - seconds of day */
    if (get_odl_field(sizeof(seconds), XXX_Double, odl_data, "SOLAR_VECTOR",
        "SOLAR_EPOCH_SECONDS", 1, &seconds) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    if (fabs(metadata->ephem_epoch_time.seconds - seconds) > SECS_TOLERANCE)
    {
        sprintf(msg, "Solar epoch seconds %lf and Ephemeris epoch seconds "
            "%lf difference is higher than %lf tolerance", seconds, 
            metadata->ephem_epoch_time.seconds, SECS_TOLERANCE);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        free(dbuffer);
        return ERROR;
    }

    /* Earth to Sun distance */
    if (get_odl_field(sizeof(metadata->earth_sun_distance), XXX_Double, 
        odl_data, "SOLAR_VECTOR", "EARTH_SUN_DISTANCE", 1, 
        &metadata->earth_sun_distance) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    /* Number of solar vector points */
    if (get_odl_field(sizeof(solar_points), XXX_Int, odl_data, 
        "SOLAR_VECTOR", "NUMBER_OF_POINTS", 1, &solar_points) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    if ((int)metadata->ephem_count != solar_points)
    {
        sprintf(msg, "Solar points %d and Ephemeris points %d are not"
            "equal", solar_points, metadata->ephem_count);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        free(dbuffer);
        return ERROR;
    }

    /* Solar vector sample time offsets from epoch */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "SOLAR_VECTOR", 
        "SAMPLE_TIME", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->solar_vector[index].sample_time = dbuffer[index];
    }

    /* Solar vector ECEF X coordinates */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "SOLAR_VECTOR", 
        "SOLAR_ECEF_X", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    { 
        metadata->solar_vector[index].position.x = dbuffer[index];
    }

    /* Solar vector ECEF Y coordinates */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "SOLAR_VECTOR",
        "SOLAR_ECEF_Y", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    {
        metadata->solar_vector[index].position.y = dbuffer[index];
    }

    /* Solar vector ECEF Z coordinates */
    if (get_odl_field(dbuffer_size, XXX_Double, odl_data, "SOLAR_VECTOR",
        "SOLAR_ECEF_Z", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        gxx_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->solar_vector[index].position.z = dbuffer[index];
    }

    /* Free the input buffer */
    free(dbuffer);

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Reads the band group from the ANG metadata file.
 *
 * @returns interger (SUCCESS or ERROR)
 */
/******************************************************************************/
static int gxx_angle_gen_read_ang_band
(
    OBJDESC *odl_data,                     /* I: ODL object */
    gxx_angle_gen_band_TYPE *band_metadata /* O: Metadata band structure to 
                                                 load */
)       
{
    int status;                     /* Status placeholder */
    int dir;                        /* Direction based loop variable */
    char band_group[11];            /* Band group name */
    double offsets[2];              /* Local offsets */
    char band_string[STRLEN];            /* Current band string */
    size_t num_size;                /* Numerator data size */
    size_t den_size;                /* Denominator data size */

    /* Initialize the numerator and denominator data sizes */
    num_size = IAS_ANGLE_GEN_ANG_RPC_COEF * sizeof(double);
    den_size = (IAS_ANGLE_GEN_ANG_RPC_COEF - 1) * sizeof(double);

    /* Construct the band string */
    sprintf(band_string, "BAND%02d",
            xxx_get_user_band(band_metadata->band_number));

    /* Read data */
    /* Band information */
    /* Construct the current band's group name */
    status = snprintf(band_group, sizeof(band_group), "RPC_%s", band_string);
    if (status < 0 || (unsigned int)status >= sizeof(band_group))
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__,
                      "Constructing the band group");
        return ERROR;
    }

    if (construct_and_get_odl_field(sizeof(band_metadata->l1t_lines),
        XXX_Int, odl_data, band_group, "NUM_L1T_LINES", band_string, 1,
        &band_metadata->l1t_lines) != SUCCESS)
    {
        return ERROR;
    }
    
    /* Read the number of L1T samples for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->l1t_samps),
        XXX_Int, odl_data, band_group, "NUM_L1T_SAMPS", band_string, 1,
        &band_metadata->l1t_samps) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of L1R lines for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->l1r_lines),
        XXX_Int, odl_data, band_group, "NUM_L1R_LINES", band_string, 1,
        &band_metadata->l1r_lines) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of L1R samples for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->l1r_samps),
        XXX_Int, odl_data, band_group, "NUM_L1R_SAMPS", band_string, 1,
        &band_metadata->l1r_samps) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the pixel size for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->pixel_size),
        XXX_Double, odl_data, band_group, "PIXEL_SIZE", band_string, 1,
        &band_metadata->pixel_size) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of scan directions for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->number_scan_dirs),
        XXX_Int, odl_data, band_group, "NUMBER_OF_DIRECTIONS", band_string, 1,
        &band_metadata->number_scan_dirs) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of lines per scan */
    if (construct_and_get_odl_field(sizeof(band_metadata->lines_per_scan),
        XXX_Int, odl_data, band_group, "LINES_PER_SCAN", band_string, 1,
        &band_metadata->lines_per_scan) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the image start time for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->image_start_time),
        XXX_Double, odl_data, band_group, "START_TIME", band_string, 1,
        &band_metadata->image_start_time) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the line time for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->seconds_per_line),
        XXX_Double, odl_data, band_group, "LINE_TIME", band_string, 1,
        &band_metadata->seconds_per_line) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite mean height */
    if (construct_and_get_odl_field(sizeof(
        band_metadata->satellite.mean_height), XXX_Double, odl_data,
        band_group, "MEAN_HEIGHT", band_string, 1, 
        &band_metadata->satellite.mean_height) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite mean L1r line/samp offsets */
    if (construct_and_get_odl_field(sizeof(offsets), XXX_Double, odl_data, 
        band_group, "MEAN_L1R_LINE_SAMP", band_string, 2, offsets) != SUCCESS)
    {
        return ERROR;
    }
    band_metadata->satellite.line_terms.l1r_mean_offset = offsets[0];
    band_metadata->satellite.samp_terms.l1r_mean_offset = offsets[1];

    /* Satellite mean L1T line/samp offsets */
    if (construct_and_get_odl_field(sizeof(offsets), XXX_Double, odl_data,
        band_group, "MEAN_L1T_LINE_SAMP", band_string, 2, offsets) != SUCCESS)
    {
        return ERROR;
    }
    band_metadata->satellite.line_terms.l1t_mean_offset = offsets[0];
    band_metadata->satellite.samp_terms.l1t_mean_offset = offsets[1];

    /* Satellite mean vector */
    if (construct_and_get_odl_field(sizeof(
        band_metadata->satellite.mean_offset), XXX_Double, odl_data, 
        band_group, "MEAN_SAT_VECTOR", band_string, 3, 
        &band_metadata->satellite.mean_offset) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite x numerator coefficients */
    if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
        band_group, "SAT_X_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF, 
        band_metadata->satellite.x_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite x denominator coefficients */
    if (construct_and_get_odl_field(den_size, XXX_Double, odl_data,
        band_group, "SAT_X_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->satellite.x_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite y numerator coefficients */
    if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
        band_group, "SAT_Y_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->satellite.y_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite y denominator coefficients */
    if (construct_and_get_odl_field(den_size, XXX_Double, odl_data,
        band_group, "SAT_Y_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->satellite.y_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite z numerator coefficients */
    if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
        band_group, "SAT_Z_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->satellite.z_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite z denominator coefficients */
    if (construct_and_get_odl_field(den_size, XXX_Double, odl_data,
        band_group, "SAT_Z_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->satellite.z_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar vector RPC model */
    band_metadata->solar.mean_height = band_metadata->satellite.mean_height;
    band_metadata->solar.line_terms.l1r_mean_offset 
        = band_metadata->satellite.line_terms.l1r_mean_offset;
    band_metadata->solar.samp_terms.l1r_mean_offset 
        = band_metadata->satellite.samp_terms.l1r_mean_offset;
    band_metadata->solar.line_terms.l1t_mean_offset 
        = band_metadata->satellite.line_terms.l1t_mean_offset;
    band_metadata->solar.samp_terms.l1t_mean_offset 
        = band_metadata->satellite.samp_terms.l1t_mean_offset;

    /* Solar mean vector */
    if (construct_and_get_odl_field(sizeof(band_metadata->solar.mean_offset),
        XXX_Double, odl_data, band_group, "MEAN_SUN_VECTOR", band_string, 3,
        &band_metadata->solar.mean_offset) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar x numerator coefficients */
    if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
        band_group, "SUN_X_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->solar.x_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar x denominator coefficients */
    if (construct_and_get_odl_field(den_size, XXX_Double, odl_data,
        band_group, "SUN_X_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->solar.x_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar y numerator coefficients */
    if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
        band_group, "SUN_Y_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->solar.y_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar y denominator coefficients */
    if (construct_and_get_odl_field(den_size, XXX_Double, odl_data,
        band_group, "SUN_Y_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->solar.y_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar z numerator coefficients */
    if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
        band_group, "SUN_Z_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->solar.z_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar z denominator coefficients */
    if (construct_and_get_odl_field(den_size, XXX_Double, odl_data,
        band_group, "SUN_Z_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->solar.z_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Read in the scan based RPC coefficients. */
    for (dir = 0; dir < band_metadata->number_scan_dirs; dir++)
    {
        /* reconstruct the band string including scan direction */
        sprintf(band_string, "BAND%02d_DIR%02d", 
            xxx_get_user_band(band_metadata->band_number), dir);

        /* Satellite mean height */
        if (construct_and_get_odl_field(sizeof(
            band_metadata->scan_metadata[dir].mean_height), XXX_Double,
            odl_data, band_group, "MEAN_HEIGHT", band_string, 1,
            &band_metadata->scan_metadata[dir].mean_height) != SUCCESS)
        {
            return ERROR;
        }

        /* Satellite mean L1r line/samp offsets */
        if (construct_and_get_odl_field(sizeof(offsets), XXX_Double, odl_data,
            band_group, "MEAN_L1R_LINE_SAMP", band_string, 2, offsets) != SUCCESS)
        {
            return ERROR;
        }
        band_metadata->scan_metadata[dir].line_terms.l1r_mean_offset = offsets[0];
        band_metadata->scan_metadata[dir].samp_terms.l1r_mean_offset = offsets[1];
    
        /* Satellite mean L1T line/samp offsets */
        if (construct_and_get_odl_field(sizeof(offsets), XXX_Double, odl_data,
            band_group, "MEAN_L1T_LINE_SAMP", band_string, 2, offsets) != SUCCESS)
        {
            return ERROR;
        }
        band_metadata->scan_metadata[dir].line_terms.l1t_mean_offset = offsets[0];
        band_metadata->scan_metadata[dir].samp_terms.l1t_mean_offset = offsets[1];

        /* Directional Line numerator coefficients */
        if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
            band_group, "LINE_NUM_COEF", band_string, 5,
            band_metadata->scan_metadata[dir].line_terms.numerator) != SUCCESS)
        {
            return ERROR;
        }

        /* Directional Line Denominator coefficients */
        if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
            band_group, "LINE_DEN_COEF", band_string, 4,
            band_metadata->scan_metadata[dir].line_terms.denominator) != SUCCESS)
        {
            return ERROR;
        }

        /* Directional Sample numerator coefficients */
        if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
            band_group, "SAMP_NUM_COEF", band_string, 5,
            band_metadata->scan_metadata[dir].samp_terms.numerator) != SUCCESS)
        {
            return ERROR;
        }

        /* Directional Sample Denominator coefficients */
        if (construct_and_get_odl_field(num_size, XXX_Double, odl_data,
            band_group, "LINE_DEN_COEF", band_string, 4,
            band_metadata->scan_metadata[dir].samp_terms.denominator) != SUCCESS)
        {
            return ERROR;
        }
    }

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Read the ANG metadata file.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
int gxx_angle_gen_read_ang
(
    char *ang_filename,                   /* I: Angle file name to read */
    gxx_angle_gen_metadata_TYPE *metadata /* O: Metadata structure to load */
)       
{
    OBJDESC *odl_data;      /* Metadata ODL object */
    unsigned int index;     /* Loop index */
    char msg[STRLEN];

    /* Open file */
    
    odl_data = xxx_OpenODL(ang_filename, xxx_NoEDCMeta, msg);
    if (!odl_data)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Read data */
    /* File information group */
    if (gxx_angle_gen_read_ang_header(odl_data, metadata) != SUCCESS)
    {
        xxx_CloseODL(odl_data, msg);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Projection information group */
    if (gxx_angle_gen_read_ang_projection(odl_data, metadata) 
        != SUCCESS)
    {
        xxx_CloseODL(odl_data, msg);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Ephemeris data group and solar vector data group */
    if (gxx_angle_gen_read_ang_ephemeris(odl_data, metadata) != SUCCESS)
    {
        xxx_CloseODL(odl_data, msg);
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    /* Band information */
    for (index = 0; index < metadata->num_bands; index++)
    {
        if (!metadata->band_present[index]) continue;

        if (gxx_angle_gen_read_ang_band(odl_data, 
            &metadata->band_metadata[index]) != SUCCESS)
        {
            xxx_CloseODL(odl_data, msg);
            xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
            gxx_angle_gen_free(metadata);
            return ERROR;
        }
    }

    /* Release the ODL structure */
    if (xxx_CloseODL(odl_data, msg) != SUCCESS)
    {
        xxx_LogStatus(PROGRAM, __FILE__, __LINE__, msg);
        return ERROR;
    }

    return SUCCESS;
}
