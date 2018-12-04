/****************************************************************************/
/**
 * @file xxx_Band.c
 * @ingroup iasLib
 *
 * @brief This contains band-related functions shared by all sensors.
 *
 *  These functions manage the band list or individual band nodes.  
 *  These lists hold various information about multiple bands.  A band
 *  node holds information about a single band.
 */
/****************************************************************************/

/* System headers */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers */

#include <xxx_LogStatus.h>

#include <xxx_Band.h>
#include "xxx_Const.h"
#include <xxx_Types.h>
#include <xxx_LogStatus.h>
#include <xxx_Sensor.h>
#include <xdb_Defines.h>       /* This is just to define SensorTypeStr.   */

static BAND_HDR_TYPE *band_hdr = NULL;
static int MaxBands;           /* This is the physical bands.             */
static int MaxLogicalBands;    /* This includes bands like B6H.           */
static int NumReflectiveBands; /* This includes 1 2 3 4 5 7 and 8.        */
static int NumResolutions;     /* Sensor's number of distinct resolutions */


/****************************************************************************/
/**
 * @brief Build a band node, which is a member of a band linked list, 
 * given all sorts of initialization values. 
 */
/****************************************************************************/

int xxx_BuildBandNode_full
(
    BAND_HDR_TYPE *band_list,//!<[in/out] Band list to put the band node in
    BAND_OPTION_TYPE option, //!<[in] Standard band list, or customized one
    int band_number,         //!<[in] physical band number            
    int band_long_number,    //!<[in] long version of the band number   
    BandNumber band_id,      //!<[in] enumerated type (eg: b1, etc) used in GPS
    BAND_TYPE band_id_rps,   //!<[in] enumerated type (eg: B1, etc) used in RPS
    char *band_short_name,   //!<[in] short version of the band's name  
    char *band_name,         //!<[in] band's name for display         
    char *sds_band_name,     //!<[in] SDS band filename section       
    char *band_name_ext,     //!<[in] band name with a "_"            
    char *band_name_ext_81,  //!<[in] band name with a "_" (and "81") 
    int num_lines,           //!<[in] number of lines
    int num_samples,         //!<[in] number of samples
    float band_resolution,   //!<[in] pixel size (in meters)            
    band_mode_type band_mode,//!<[in] enumerated type (THERMAL...)    
    temp_corr_type temp_corr,//!<[in] enumerated type (NO_TEMP_CORRECTION...)
    int num_detectors,       //!<[in] number of detectors for the band 
    long int maxImgPixels,   //!<[in] from rxx_Constants.h   
    long int maxCalPixels,   //!<[in] from rxx_Constants.h   
    char band_id_msg[5],   //!<[in] from rxx_BAND_TYPE_MSG     
    char IMGnames[5],      //!<[in] from RXX_BAND1, etc...     
    char CALnames[5],      //!<[in] from RXX_CAL_BAND1, etc... 
    char SLOvdata[6],      //!<[in] from XXX_SLO_BAND1, etc... 
    int joffset,           //!<[in] offset - used in IC Trend write          
    int active             //!<[in] indicates whether the band is in use     
)
{
    BAND_NODE_TYPE *band_node = NULL;
    BAND_NODE_TYPE *band_tmp  = NULL;

    if (band_list == NULL)         /* The band header does not exist. */
    {
        xxx_LogStatus(BAND_PROG_NAME, __FILE__, __LINE__,
                      "The band header does not exist.");
        return ERROR;
    }

    /* Create the band node. */

    band_node = malloc (BAND_NODE_TYPE_SIZE);

    if (band_node == NULL)
    {
        xxx_LogStatus(BAND_PROG_NAME, __FILE__, __LINE__,
                      "malloc error for band node in xxx_BuildBandNode");
        return ERROR; 
    }

    band_node->Prev = NULL;
    band_node->Next = NULL;

    band_node->band_number       = band_number;
    band_node->band_long_number  = band_long_number; 
    band_node->band_id           = band_id; 
    band_node->band_id_rps       = band_id_rps; 

    band_node->band_short_name = malloc (strlen (band_short_name) + 1);
    strcpy (band_node->band_short_name, band_short_name);

    band_node->band_name = malloc (strlen (band_name) + 1);
    strcpy (band_node->band_name, band_name); 

    band_node->sds_band_name = malloc (strlen (sds_band_name) + 1);
    strcpy (band_node->sds_band_name, sds_band_name); 

    band_node->band_name_ext = malloc (strlen(band_name_ext) + 1);
    strcpy (band_node->band_name_ext, band_name_ext); 

    band_node->band_name_ext_81 = malloc (strlen(band_name_ext_81) + 1);
    strcpy (band_node->band_name_ext_81, band_name_ext_81); 

    band_node->num_lines = num_lines;
    band_node->num_samples = num_samples;
    band_node->band_resolution     = band_resolution; 
    band_node->band_mode           = band_mode; 
    band_node->temp_corr           = temp_corr;
    band_node->num_detectors       = num_detectors; 
    band_node->maxImgPixels        = maxImgPixels;
    band_node->maxCalPixels        = maxCalPixels;

    strcpy (band_node->band_id_msg, band_id_msg);
    strcpy (band_node->IMGnames,    IMGnames);
    strcpy (band_node->CALnames,    CALnames);
    strcpy (band_node->SLOvdata,    SLOvdata);

    band_node->joffset = joffset;

    band_node->active = active; 

    /* Insert the band node into the band linked list. */

    band_list->num_bands++;

    if (band_list->Head == NULL)  /* The linked list is empty. */
    {
        band_list->Head = band_node;
        band_list->End  = band_node;
    }
    else  /* Insert the new band node at the end of the list.    */
    {     /* This is to keep the list in band order while making */
          /* the insertion code more readable (in order).        */

        band_tmp = band_list->End;
        band_list->End = band_node;
        band_node->Prev = band_tmp;
        band_tmp->Next = band_node;
    }

    if (option == STANDARD)
    {
         band_hdr = band_list;
    }

    return SUCCESS;
}


/****************************************************************************/
/**
 * @brief Build a band node, which is a member of a band linked list, given a
 * typical set of initialization values.
 **/
/****************************************************************************/

int xxx_BuildBandNode 
(
    BAND_HDR_TYPE *band_list,//!<[in/out] Band list to put the band node in 
    BAND_OPTION_TYPE option, //!<[in] Standard band list, or customized one
    int band_number,         //!<[in] Physical band number
    char *band_name,         //!<[in] Band's name for display
    int num_lines,           //!<[in] number of lines
    int num_samples,         //!<[in] number of samples
    float band_resolution    //!<[in] Pixel size (in meters)
)
{
    char band_name_ext[11];
    BandNumber bnum;
    BAND_TYPE btype;

    snprintf(band_name_ext, 11, "_%s", band_name);

    if (band_number < 6)
    {
        bnum = band_number - 1;
        btype = band_number - 1;
    }
    else if (band_number == 6)
    {
        bnum = b6;
        btype = B6;
    }
    else if (band_number == 61)
    {
        bnum = b6l;
        btype = B6L;
    }
    else if (band_number == 62)
    {
        bnum = b6h;
        btype = B6H;
    }
    else
    {
        bnum = band_number;
        btype = band_number;
    }

    return xxx_BuildBandNode_full(band_list, option, band_number, band_number,
                                  bnum, btype, band_name, band_name, band_name,
                                  band_name_ext, band_name_ext,
                                  num_lines, num_samples, band_resolution,
                                  REFLECTIVE, TEMP_NO_CORRECTION, 1,
                                  100000, 100000, band_name, band_name,
                                  band_name, band_name, 0, 0);
}


/****************************************************************************/
/** 
 * @brief Display the Header of the band linked list.
 *
 * These display functions are inspired by the LMask display functions. 
 * They are meant to be linked list debugging/dump tools.               
 */
/****************************************************************************/

void xxx_DisplayBandHdr 
(
    BAND_HDR_TYPE *band_list //!<[in] header of the band list being displayed 
)
{
    if (band_list == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "The band header does not exist.");
        return;
    }

    printf ("---------------- BAND HEADER ------------------\n");
    printf ("POINTERS---------------------------------------\n");
    printf ("band_hdr->Head      = 0x%p\n", band_list->Head);
    printf ("band_hdr->End       = 0x%p\n", band_list->End);
    printf ("band_hdr->Current   = 0x%p\n", band_list->Current);
    printf ("band_hdr->ref_band  = 0x%p\n", band_list->ref_band);
    printf ("STATIC VALUES----------------------------------\n");
    printf ("band_hdr->num_bands = %d\n",   band_list->num_bands);
    printf ("band_hdr->last_band = %d\n",   band_list->last_band);
    printf ("-----------------------------------------------\n");
}


/****************************************************************************/
/**
 * @brief Display a single Band Node in the band linked list.                      
 */ 
/****************************************************************************/

void xxx_DisplayBandNode 
(
    BAND_NODE_TYPE *band_node //!<[in] a band node being displayed
)
{
    if (band_node == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "The band node does not exist.");
        return;
    }

    printf ("----------------- BAND NODE -------------------\n");
    printf ("POINTERS---------------------------------------\n");
    printf ("band_node                    = 0x%p\n", band_node);
    printf ("band_node->Next              = 0x%p\n", band_node->Next);
    printf ("band_node->Prev              = 0x%p\n", band_node->Prev);
    printf ("STATIC VALUES----------------------------------\n");
    printf ("band_node->band_number       = %d\n", band_node->band_number);
    printf ("band_node->band_long_number  = %d\n", band_node->band_long_number);
    printf ("band_node->band_id           = %d\n", band_node->band_id);
    printf ("band_node->band_id_rps       = %d\n", band_node->band_id_rps);
    printf ("band_node->band_short_name   = %s\n", band_node->band_short_name);
    printf ("band_node->band_name         = %s\n", band_node->band_name);
    printf ("band_node->sds_band_name     = %s\n", band_node->sds_band_name);
    printf ("band_node->band_name_ext     = %s\n", band_node->band_name_ext);
    printf ("band_node->band_name_ext_81  = %s\n", band_node->band_name_ext_81);
    printf ("band_node->num_lines         = %d\n", band_node->num_lines);
    printf ("band_node->num_samples       = %d\n", band_node->num_samples);
    printf ("band_node->band_resolution   = %f\n", band_node->band_resolution);
    printf ("band_node->band_mode         = %d\n", band_node->band_mode);
    printf ("band_node->temp_corr         = %d\n", band_node->temp_corr);
    printf ("band_node->num_detectors     = %d\n", band_node->num_detectors);
    printf ("band_node->maxImgPixels      = %ld\n",band_node->maxImgPixels);
    printf ("band_node->maxCalPixels      = %ld\n",band_node->maxCalPixels);
    printf ("band_node->CalDataLen        = %ld\n",band_node->CalDataLen);
    printf ("band_node->CalDkRegLen       = %ld\n",band_node->CalDkRegLen);
    printf ("band_node->CalDkRegLength[0] = %d\n", 
                                            band_node->CalDkRegLength[0]);
    printf ("band_node->CalDkRegLength[1] = %d\n", 
                                            band_node->CalDkRegLength[1]);
    printf ("band_node->CalDkRegOffset[0] = %d\n", 
                                            band_node->CalDkRegOffset[0]);
    printf ("band_node->CalDkRegOffset[1] = %d\n", 
                                            band_node->CalDkRegOffset[1]);
    printf ("band_node->IC_Usable_Length[0] = %d\n", 
                                            band_node->IC_Usable_Length[0]);
    printf ("band_node->IC_Usable_Length[1] = %d\n", 
                                            band_node->IC_Usable_Length[1]);
    printf ("band_node->band_id_msg       = %s\n", band_node->band_id_msg);
    printf ("band_node->IMGnames          = %s\n", band_node->IMGnames);
    printf ("band_node->CALnames          = %s\n", band_node->CALnames);
    printf ("band_node->SLOvdata          = %s\n", band_node->SLOvdata);
    printf ("band_node->joffset           = %d\n", band_node->joffset);
    printf ("DYNAMIC VALUES---------------------------------\n");
    printf ("band_node->active            = %d\n", band_node->active);
    printf ("-----------------------------------------------\n");
    printf ("\n");
}


/****************************************************************************/
/**
 * @brief Display a whole band linked list, including the header and all nodes.  
 */
/****************************************************************************/

void xxx_DisplayBandList 
(
    BAND_HDR_TYPE *band_list  //!<[in] the band list that is being displayed 
)
{
    BAND_NODE_TYPE *band_node;

    if (band_list == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "The band header does not exist.");
        return;
    }

    printf ("-------------------------------------------------\n");
    printf ("--------------- ALL BAND NODES ------------------\n");
    printf ("-------------------------------------------------\n");

    if (band_list->Head == NULL)
    {
        xxx_DisplayBandHdr (band_list);
        printf ("There are no band nodes to display.\n");
    }
    else
    {
        xxx_DisplayBandHdr (band_list);
        band_node = band_list->Head;

        while (band_node != NULL)
        {
            xxx_DisplayBandNode (band_node);
            band_node = band_node->Next;
        }
    }
}


/**************************************************************/
/** 
 * @brief Given the header of a band linked list, and a band number, 
 * return the band node that matches the band number.         
 *                                                            
 * This will work for matches with various fields. Currently  
 * it is restricted to integer fields.  It is not C++ after   
 * all.  (Note: I'm also using it for enumerated type         
 * fields, which are integers, so I don't have to duplicate   
 * the function for them).
 */
/**************************************************************/
                                  
BAND_NODE_TYPE * xxx_GetBandNode 
(
    BAND_HDR_TYPE *band_hdr,  //!<[in] header of the band list being searched
    int search_value,         //!<[in] the value we are searching for         
    BAND_FIELD band_field     //!<[in] where we are looking                  
)   
{
    BAND_NODE_TYPE *band_node;
    char msg[STRLEN];   /* message buffer */

    if (band_hdr == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "The band header does not exist.");
        return NULL;
    }

    if (band_hdr->Head == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "There are no band nodes to match.");
        return NULL;
    }

    band_node = band_hdr->Head;

    switch (band_field)
    {
        case BAND_NUMBER_TYPE:

            while ((band_node != NULL) && 
                   (band_node->band_number != search_value))
            {
                band_node = band_node->Next;
            }

            break;

        case BAND_LONG_NUMBER_TYPE:

            while ((band_node != NULL) && 
                   (band_node->band_long_number != search_value))
            {
                band_node = band_node->Next;
            }

            break;

        case BAND_ID_TYPE:

            while ((band_node != NULL) && 
                   (band_node->band_id != search_value))
            {
                band_node = band_node->Next;
            }

            break;

        case BAND_ID_RPS_TYPE:

            while ((band_node != NULL) && 
                   (band_node->band_id_rps != (BAND_TYPE)search_value))
            {
                band_node = band_node->Next;
            }

            break;

        case BAND_RESOLUTION_TYPE:
    
            while ((band_node != NULL) && 
                   (band_node->band_resolution != search_value))
            {
                band_node = band_node->Next;
            }

            break;

        default:

            sprintf(msg, "Unsupported band field type %d.", band_field);
            xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, msg);
            return NULL;

            break;
    }

    return band_node;
}


/**************************************************************/
/**
 * @brief Given the header of a band linked list, and a band name,
 * return the band node that matches the band name.           
 *                                                            
 * This will work for matches with various fields. Currently  
 * it is restricted to string fields.  Otherwise it is just   
 * like xxx_GetBandNode.
 */
/**************************************************************/
                                  
BAND_NODE_TYPE * xxx_GetBandNodeStr 
(
    BAND_HDR_TYPE *band_hdr,  //!<[in] header of the band list being searched 
    char *search_value,       //!<[in] the value we are searching for         
    BAND_FIELD band_field     //!<[in] where we are looking                   
)
{
    BAND_NODE_TYPE *band_node;
    char msg[STRLEN];    /* message buffer */

    if (band_hdr == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "The band header does not exist.");
        return NULL;
    }

    if (band_hdr->Head == NULL)
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "There are no band nodes to match.");
        return NULL;
    }

    band_node = band_hdr->Head;

    switch (band_field)
    {
        case BAND_SHORT_NAME_TYPE:

            while ((band_node != NULL) && 
                   (strcmp (band_node->band_short_name, search_value) != 0))
            {
                band_node = band_node->Next;
            }

            break;

        case BAND_NAME_TYPE:

            while ((band_node != NULL) && 
                   (strcmp (band_node->band_name, search_value) != 0))
            {
                band_node = band_node->Next;
            }

            break;

        case BAND_SDS_NAME_TYPE:

            while ((band_node != NULL) && 
                   (strcmp (band_node->sds_band_name, search_value) != 0))
            {
                band_node = band_node->Next;
            }

            break;

        default:

            sprintf(msg, "Unsupported band field type %d.", band_field);
            xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, msg);
            return NULL;

            break;
    }

    if (band_node == NULL)
    {
        sprintf(msg, "No band matched the band search value %s.", search_value);
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, msg);
    }

    return band_node;
}


/**************************************************************/
/** 
 * @brief To retrieve a band structure
 */
/**************************************************************/

BAND_HDR_TYPE *xxx_GetBands (void)
{
    return band_hdr;
}

void xxx_SetBands(BAND_HDR_TYPE *bh)
{
    band_hdr = bh;
}

/******************************************************************/
/** 
 * @brief To free the band header and all nodes in the band linked list.
 */
/******************************************************************/

int xxx_FreeBands
(
    BAND_HDR_TYPE *band_hdr  //!<[in] Pointer to band header structure 
)
{
    BAND_NODE_TYPE *band_node, *pDel;

    if (band_hdr == NULL)         /* band header does not exist */
    {
        xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, 
            "band header does not exist during call to xxx_FreeBands");
        return ERROR;
    }

    band_node = band_hdr->Head;    /* point to start of band node list */

    while (band_node != NULL)
    {
        pDel = band_node;
        band_node = band_node->Next;
        free(pDel);
    }

    free(band_hdr);

    return SUCCESS;
}


/***************************************************************/
/** 
 * @brief Create a blank header for a band list.
 */
/***************************************************************/

BAND_HDR_TYPE *xxx_CreateBandHdr 
(
    BAND_OPTION_TYPE option //!<[in] create the standard or customized band list
)
{
    BAND_HDR_TYPE *band_list;

    band_list = malloc (BAND_HDR_TYPE_SIZE);

    if (band_list == NULL)
    {
        xxx_LogStatus(BAND_PROG_NAME, __FILE__, __LINE__, 
            "malloc error for band_list in xxx_CreateBandHdr");
        return NULL;
    }

    /* Initialize the band header. */

    band_list->last_band = 0;   /* These need to be computed later. */
    band_list->num_bands = 0;   /* These need to be computed later. */

    band_list->Head    = NULL;
    band_list->End     = NULL;
    band_list->Current = NULL;

    if (option == STANDARD) /* building the standard band list */
        band_hdr = band_list;

    return band_list;
}


/****************************************************************************/
/** 
 * @brief Initialize the number of physical bands for the sensor.
 */ 
/****************************************************************************/

int xxx_initialize_max_bands(void)
{
    IAS_SENSOR_TYPE SensorType = xxx_get_sensor_type();
    char msg[STRLEN];        /* message buffer      */

    switch( SensorType )
    {
        case IAS_SENSOR_ETM:
            MaxBands = MAX_BANDS_ETM;
            MaxLogicalBands = MAX_LOGICAL_BANDS_ETM;
            NumReflectiveBands = NUM_REFLECTIVE_BANDS_ETM;
            return SUCCESS;

        case IAS_SENSOR_TM:
            MaxBands = MAX_BANDS_TM;
            MaxLogicalBands = MAX_LOGICAL_BANDS_TM;
            NumReflectiveBands = NUM_REFLECTIVE_BANDS_TM;
            return SUCCESS;

        case IAS_SENSOR_MSS:
            MaxBands = MAX_BANDS_MSS;
            MaxLogicalBands = MAX_LOGICAL_BANDS_MSS;
            NumReflectiveBands = NUM_REFLECTIVE_BANDS_MSS;
            return SUCCESS;

        /* Something bad happened */
        default:
            sprintf(msg, "Sensor type %d is not supported.", SensorType);
            xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, msg);
            return ERROR;
    }
}


/****************************************************************************/
/** 
 * @brief Initialize the number of band resolutions for the sensor.
 */ 
/****************************************************************************/

int xxx_initialize_num_resolutions(void)
{
    IAS_SENSOR_TYPE SensorType = xxx_get_sensor_type();
    char msg[STRLEN];        /* message buffer      */

    switch( SensorType )
    {
        case IAS_SENSOR_ETM:
            NumResolutions = NUM_RESOLUTIONS_ETM;
            return SUCCESS;

        case IAS_SENSOR_TM:
            NumResolutions = NUM_RESOLUTIONS_TM;
            return SUCCESS;

        case IAS_SENSOR_MSS:
            NumResolutions = NUM_RESOLUTIONS_MSS;
            return SUCCESS;

        /* Something bad happened */
        default:
            sprintf(msg, "Sensor type %d is not supported.", SensorType);
            xxx_LogStatus (BAND_PROG_NAME, __FILE__, __LINE__, msg);
            return ERROR;
    }
}


/****************************************************************************/
/** 
 * @brief Get the number of physical bands for the sensor.
 */ 
/****************************************************************************/

int xxx_get_max_bands(void)
{
    return MaxBands;
}


/****************************************************************************/
/** 
 * @brief Get the number of physical bands for the sensor.
 */ 
/****************************************************************************/

int xxx_get_max_logical_bands(void)
{
    return MaxLogicalBands;
}


/****************************************************************************/
/** 
 * @brief Get the number of reflective bands for the sensor.
 */ 
/****************************************************************************/

int xxx_get_num_reflective_bands(void)
{
    return NumReflectiveBands;
}


/****************************************************************************/
/** 
 * @brief Get the number of band resolutions for the sensor.
 */ 
/****************************************************************************/

int xxx_get_num_resolutions(void)
{
    return NumResolutions;
}


/****************************************************************************/
/**
 * @brief Build/update a list of pixel resolutions for a given band list.
 * Resolutions are in descending order.
 *
 * @return SUCCESS: if successful completion
 * @return ERROR: if operation failed
 *
 *
 */
/****************************************************************************/
int xxx_get_res_list
(
    BAND_HDR_TYPE *band_hdr,     //!<[in] band list header 
    int array_size,              //!<[in] size of res and numbands_in_res arrays
    int *numres,                 //!<[in/out] number of resolutions in band list
    float res[],           //!<[in/out] different pixel resolutions in band list
    int numbands_in_res[],       //!<[in/out] number of bands in each resolution 
    int bandflag[]             //!<[in/out] array of flags indicating bands that
                               //! have been accounted for; must be
                               //! large enough for all bands in list 
)
{
    int i, j;    /* loop counters */
    int newres=1;/* flag identifying new resolution to be added to the list */
    BAND_NODE_TYPE *band_node;    /* band node in list */

    
    band_node = band_hdr->Head;    
    while (band_node != NULL)
    {
        /* Determine whether this band's resolution is new to the list. */
        for (i=0; i<*numres; i++)
            if (band_node->band_resolution >= res[i])
            {
                if (band_node->band_resolution == res[i])
                    newres = 0;
                else
                    newres = 1;
                break;
            }
        if (i == *numres)
            newres = 1;
        
        if (!newres)
        {
            /* Resolution already in list. */
            if (!bandflag[band_node->band_id])
            {
                bandflag[band_node->band_id] = 1;
                numbands_in_res[i]++;
            }
            band_node = band_node->Next;
            continue;
        }
        else if (*numres == array_size)
        {
            xxx_LogStatus(BAND_PROG_NAME, __FILE__, __LINE__,
                        "Too many different pixel resolutions in band list.");
            return ERROR;
        }

        /* Add new resolution to list (in descending order). */
        for (j=*numres; j>i; j--)
        {
            res[j] = res[j-1];
            numbands_in_res[j] = numbands_in_res[j-1];
        }
        res[i] = band_node->band_resolution;
        bandflag[band_node->band_id] = 1;
        numbands_in_res[i] = 1;
        (*numres)++;
        band_node = band_node->Next;
    }

    return SUCCESS;
}

/****************************************************************************/
/**
 * @brief Return the enumerated BandNumber type to the number entered by
 *        the user
 *
 * @return user_band_number: user band number for specified band
 */
/****************************************************************************/
int xxx_get_user_band
(
    BandNumber band_id   //!<[in] enumerated band number
)
{
    static int user_band_number[] =
    {
        [b1]  = 1,
        [b2]  = 2,
        [b3]  = 3,
        [b4]  = 4,
        [b5]  = 5,
        [b6l] = 61,
        [b6h] = 62,
        [b7]  = 7,
        [b8]  = 8,
        [b6]  = 6
    };
    /* lookup table to convert the enumerated number to the user band
       number */

    /* If the enumerated type is in the legal range, return the
       user band number, otherwise return berror -- NOTE that
       b6 comes at the end of the enumerated type (i.e. it has a
       larger value than b8) */
    if (band_id >= b1 && band_id <= b6)
        return user_band_number[band_id];
    else
        return berror;
}


/****************************************************************************/
/**
 * @brief Convert between user band numbers and the internal enumerated
 *        band number.
 *
 * @return Bandnumber enum type
 */
/****************************************************************************/
BandNumber xxx_parse_user_band
(
    int user_band_number   //!<[in] integer band number
)
{
    switch(user_band_number)
    {
        case 1:
            return b1;
        case 2:
            return b2;
        case 3:
            return b3;
        case 4:
            return b4;
        case 5:
            return b5;
        case 6:
            return b6;
        case 61:
            return b6l;
        case 62:
            return b6h;
        case 7:
            return b7;
        case 8:
            return b8;
        default:
            return berror;
    }
}
