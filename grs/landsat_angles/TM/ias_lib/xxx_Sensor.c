/****************************************************************************/
/**
 * @file xxx_Sensor.c
 * @brief To handle multi-sensor initialization and data access.
 * @ingroup iasLib
 */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>
#include <limits.h>
#include "xxx_Const.h"
#include "xxx_Sensor.h"
#include "xxx_LogStatus.h"
#include "xxx_LogError.h"

static xxx_FuncPtrs xxx;
static IAS_SENSOR_TYPE SensorType = IAS_SENSOR_UNKNOWN;
static IAS_SATELLITE_TYPE SatelliteType = IAS_SAT_UNKNOWN;
static IAS_FORMAT_TYPE format_type = IAS_INVALID;
static IAS_MSS_X_TYPE mss_x_type = IAS_MSS_X_INVALID;
static char Sensor_str[MAX_SSS] = "UNK"; /* sensor type string */
static char Satellite_str[MAX_SSS] = "UNK"; /* satellite string ("L5" "L7") */


/******************************************************************************/
/**
 * @brief converts the input string to an enumeration value.
 */
/******************************************************************************/
IAS_SENSOR_TYPE xxx_initialize_sensor_type
(
    const char * SensorTypeStr   //!<[in] sensor type string 
)
{
    /* If TM setup */
    if (strcasecmp(SensorTypeStr,"TM") == 0 ||
        strcasecmp(SensorTypeStr,"L5_TM") == 0 ||
        strcasecmp(SensorTypeStr,"L4_TM") == 0)
    {
        SensorType = IAS_SENSOR_TM;
        strcpy(Sensor_str, "TM");
        if (strcasecmp(SensorTypeStr,"L5_TM") == 0)
        {
            SatelliteType = IAS_SAT_L5;
            strcpy(Satellite_str, "L5");
        }
        else if (strcasecmp(SensorTypeStr,"L4_TM") == 0)
        {
            SatelliteType = IAS_SAT_L4;
            strcpy(Satellite_str, "L4");
        }
    }
    /* If ETM setup */
    else if (strcasecmp(SensorTypeStr,"ETM") == 0 ||
             strcasecmp(SensorTypeStr,"L7_ETM") == 0)
    {
        SensorType = IAS_SENSOR_ETM;
        SatelliteType = IAS_SAT_L7;
        strcpy(Sensor_str, "ETM");
        strcpy(Satellite_str, "L7");
    }
    /* MSS */
    else if (strcasecmp(SensorTypeStr,"MSS") == 0 ||
        strcasecmp(SensorTypeStr,"L1_MSS") == 0 ||
        strcasecmp(SensorTypeStr,"L2_MSS") == 0 ||
        strcasecmp(SensorTypeStr,"L3_MSS") == 0 ||
        strcasecmp(SensorTypeStr,"L4_MSS") == 0 ||
        strcasecmp(SensorTypeStr,"L5_MSS") == 0)
    {
        SensorType = IAS_SENSOR_MSS;
        strcpy(Sensor_str, "MSS");
        if (strcasecmp(SensorTypeStr,"L1_MSS") == 0)
        {
            SatelliteType = IAS_SAT_L1;
            strcpy(Satellite_str, "L1");
        }
        else if (strcasecmp(SensorTypeStr,"L2_MSS") == 0)
        {
            SatelliteType = IAS_SAT_L2;
            strcpy(Satellite_str, "L2");
        }
        else if (strcasecmp(SensorTypeStr,"L3_MSS") == 0)
        {
            SatelliteType = IAS_SAT_L3;
            strcpy(Satellite_str, "L3");
        }
        else if (strcasecmp(SensorTypeStr,"L4_MSS") == 0)
        {
            SatelliteType = IAS_SAT_L4;
            strcpy(Satellite_str, "L4");
        }
        else if (strcasecmp(SensorTypeStr,"L5_MSS") == 0)
        {
            SatelliteType = IAS_SAT_L5;
            strcpy(Satellite_str, "L5");
        }
    }
    /* Unknown sensor type */
    else
    {
        SensorType = IAS_SENSOR_UNKNOWN;
        SatelliteType = IAS_SAT_UNKNOWN;
        strcpy(Sensor_str, "UNK");
        strcpy(Satellite_str, "UNK");
        xxx_LogError(__FILE__, __LINE__,"Setting sensor to UNKNOWN");
    }

    return SensorType;
}


/******************************************************************************/
/** 
 * @brief Determine the sensor type by examining the first few characters
 * of the file name in PathName.  
 *
 * @return A character string suitable for passing to 
 * xxx_initialize_sensor_type()
 */
/******************************************************************************/
const char *xxx_get_sensor_from_fname
(
    const char *PathName         //!<[in] full path filename 
)
{
    char pathcopy[PATH_MAX];      /* Copy of PathName; passed to basename() */
    char *filename;               /* The file name derived from basename()  */
    char error_msg[PATH_MAX + STRLEN];

    strncpy(pathcopy, PathName, PATH_MAX - 1);
    pathcopy[PATH_MAX - 1] = '\0';     /* Make sure path is null-terminated */

    filename = basename(pathcopy);  

    if (strncmp("L1", filename, 2) == 0)
        return "L1_MSS";
    else if (strncmp("L2", filename, 2) == 0)
        return "L2_MSS";
    else if (strncmp("L3", filename, 2) == 0)
        return "L3_MSS";
    else if (strncmp("L4", filename, 2) == 0)
    {
        /* Discern between TM and MSS by the "processor number":
           0 = TM, 1 = MSS. */
        if (filename[7] == '0')
        {
            /* Set the data format. */
            if (filename[6] == '1')
                format_type = IAS_TM_R;
            else
                format_type = IAS_TM_A;
            return "L4_TM";
        }
        else
            return "L4_MSS";
    }
    else if (strncmp("L5", filename, 2) == 0)
    {
        /* Discern between TM and MSS by the "processor number":
           0 = TM, 1 = MSS. */
        if (filename[7] == '0')
        {
            /* Set the data format. */
            if (filename[6] == '1')
                format_type = IAS_TM_R;
            else
                format_type = IAS_TM_A;
            return "L5_TM";
        }
        else
            return "L5_MSS";
    }
    else if (strncmp("tmodel", filename, 6) == 0)
        return "TM";
    else if (strncmp("mmodel", filename, 6) == 0)
        return "MSS";
    if (strncmp("L7", filename, 2) == 0 ||
        strncmp("etmodel", filename, 7) == 0)
        return "ETM";
    else
    {
        sprintf(error_msg, 
                "Cannot identify sensor from filename: %s", filename);
        xxx_LogError(__FILE__, __LINE__, error_msg);
        return "UNK"; /* Unknown sensor type */
    }
}


/******************************************************************************/
/**
 * @brief Returns the sensor type
 */
/******************************************************************************/
IAS_SENSOR_TYPE xxx_get_sensor_type(void)
{
    return SensorType;
}


/******************************************************************************/
/**
 * @brief returns the sensor type string
 */
/******************************************************************************/
char *xxx_get_sensor_type_str(void)
{
    return Sensor_str;
}


/******************************************************************************/
/**
 * @brief returns the satellite type 
 */
/******************************************************************************/
IAS_SATELLITE_TYPE xxx_get_satellite_type(void)
{
    return SatelliteType;
}


/******************************************************************************/
/**
 * @brief returns the satellite string
 */
/******************************************************************************/
char *xxx_get_satellite_str(void)
{
    return Satellite_str;
}


/******************************************************************************/
/**
 * @brief returns format type
 */
/******************************************************************************/
IAS_FORMAT_TYPE xxx_set_format
(
    const char *format           //!<[in] TM or MSS format string 
)
{
    if (format[0] == 'R')
        format_type = IAS_TM_R;
    else if (format[0] == 'X')
        format_type = IAS_MSS_X;
    else if (format[0] == 'P')
        format_type = IAS_MSS_P;
    else if (format[0] == 'A')
    {
        if (SensorType == IAS_SENSOR_TM)
            format_type = IAS_TM_A;
        else
            format_type = IAS_MSS_A;
    }
    
    return format_type;
}


/******************************************************************************/
/**
 * @brief returns format type
 */
/******************************************************************************/
IAS_FORMAT_TYPE xxx_get_format(void)
{
    return format_type;
}


/******************************************************************************/
/**
 * @brief sets the mss x type
 */
/******************************************************************************/
IAS_MSS_X_TYPE xxx_set_mss_x_type
(
    const char *type             //!<[in] MSS-X data type string 
)
{
    if (strncmp(type, "WBV", 3) == 0)
        mss_x_type = IAS_MSS_X_WBV;
    else if (strncmp(type, "CCT", 3) == 0)
        mss_x_type = IAS_MSS_X_CCT;
    else if (strncmp(type, "ORF", 3) == 0)
        mss_x_type = IAS_MSS_X_ORF;

    return mss_x_type;
}


/******************************************************************************/
/**
 * @brief Returns the mss x type
 */
/******************************************************************************/
IAS_MSS_X_TYPE xxx_get_mss_x_type(void)
{
    return mss_x_type;
}


/******************************************************************************/
/**
 * @brief Returns the current function pointers 
 */
/******************************************************************************/
xxx_FuncPtrs *xxx_get_funcs( void )
{
    return ( &xxx );
}

/******************************************************************************/
/**
 * @brief Retrieve the database connection string ("username/password") for the
 * relevant sensor schema. 
 */
/******************************************************************************/
char *xxx_get_dbsensor(void)
{
    char *dbuser = NULL;

    
    switch (SensorType)
    {
        case IAS_SENSOR_ETM:
            dbuser = getenv("IAS_DB_L7_ETM");
            if (dbuser == NULL)
                xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                              "Error: IAS_DB_L7_ETM environment variable "
                              "not set.");
            break;

        case IAS_SENSOR_TM:
            switch (SatelliteType)
            {
                case IAS_SAT_L5:
                    dbuser = getenv("IAS_DB_L5_TM");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L5_TM environment "
                                      "variable not set.");
                    break;

                case IAS_SAT_L4:
                    dbuser = getenv("IAS_DB_L4_TM");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L4_TM environment "
                                      "variable not set.");
                    break;

                default:
                case IAS_SAT_UNKNOWN:
                    printf("Satellite type %d with TM sensor not "
                           "supported.\n", SatelliteType);
            }
            break;

        case IAS_SENSOR_MSS:
            switch (SatelliteType)
            {
                case IAS_SAT_L1:
                    dbuser = getenv("IAS_DB_L1_MSS");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L1_MSS environment "
                                      "variable not set.");
                    break;

                case IAS_SAT_L2:
                    dbuser = getenv("IAS_DB_L2_MSS");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L2_MSS environment "
                                      "variable not set.");
                    break;

                case IAS_SAT_L3:
                    dbuser = getenv("IAS_DB_L3_MSS");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L3_MSS environment "
                                      "variable not set.");
                    break;

                case IAS_SAT_L4:
                    dbuser = getenv("IAS_DB_L4_MSS");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L4_MSS environment "
                                      "variable not set.");
                    break;

                case IAS_SAT_L5:
                    dbuser = getenv("IAS_DB_L5_MSS");
                    if (dbuser == NULL)
                        xxx_LogStatus("xxx_get_dbsensor", __FILE__, __LINE__,
                                      "Error: IAS_DB_L5_MSS environment "
                                      "variable not set.");
                    break;

                default:
                case IAS_SAT_UNKNOWN:
                    printf("Satellite type %d with MSS sensor not "
                           "supported.\n", SatelliteType);
            }
            break;

        default:
        case IAS_SENSOR_UNKNOWN:
            printf("Sensor type %d is not supported. \n", SensorType);
            break;
    }
    
    return dbuser;
}


/******************************************************************************/
/** 
 * @brief Map internal MSS band number to "public" value.  
 * 
 * The mapping only applies to L4 and L5, and is as follows:
 *
 *   public  internal
 *   ------  --------
 *      1        4
 *      2        5
 *      3        6
 *      4        7
 *
 *  The newbandstr string must of equal or greater length than bandstr.
 *  
 * @return public band number  
 * @return On error, returns 0. 
 */
/*******************************************************************************/
int xxx_map_band
(
    const int band,       //!<[in] Internal band number 
    const char *bandstr,  //!<[in] Internal band string; may be NULL 
    char *newbandstr      //!<[out] Public band string 
)
{
    int newband = band - 3; /* new band number */
    char bstr[2];           /* string representation of band number */
    char newbstr[2];        /* string representation of new band number */
    char *ptr;              /* pointer to band number in band string */

    /* Copy the band string. */
    if (bandstr != NULL)
        strcpy(newbandstr, bandstr);
    
    if (SensorType == IAS_SENSOR_MSS &&
        (SatelliteType == IAS_SAT_L4 || SatelliteType == IAS_SAT_L5))
    {
        /* Validate the band number. */
        if (newband < 1 || band > 9)
        {
            xxx_LogStatus("xxx_map_band", __FILE__, __LINE__,
                          "Band number out of range.");
            return 0;
        }

        /* If the band string was specified, create the new band string. */
        if (bandstr != NULL)
        {
            /* Find the band number in the string, and replace it. */
            sprintf(bstr, "%1d", band);
            sprintf(newbstr, "%1d", newband);
            ptr = strstr(newbandstr, bstr);
            if (ptr != NULL)
                *ptr = newbstr[0];
        }
        return newband;
    }
    else
        return band;
}
