
#include <stdio.h>
#include <string.h>
#include <lablib3.h>                 /* prototypes for the ODL functions */
#include <xxx_ODL.h>                 /* prototype for the xxx ODL functions */

/*****************************************************************************/
/**
 * @file xxx_ConvertString.c
 * @brief To convert the string to the new type using the memory address given
 * by caller
 * @ingroup iasLib
 */
/**
 * @brief To convert the string to the new type using the memory address given
 * by caller
 *
 * @return SUCCESS: Converted the String in to requested type
 * @return ERROR: Allocation error or unknown parm_type
 * @return XXX_E_ODL_INVALID_DATA_TYPE: Invalid data type
 */
/*****************************************************************************/

int xxx_ConvertString
(
    void *p_destination,        //!<[out] Location where to copy converted
                                // values
    xxx_ValType_TYPE parm_type, //!<[in] What type the p_destination is, enum
                                // field
    char *kvalue,               //!<[in] String to be converted
    char *p_ErrorMessage        //!<[out] Error Message
)
{
    char *p_endptr;
    char *p_a;

    
    switch(parm_type)
    {
        case XXX_Long:
            (*(long *)p_destination) = strtol(kvalue,&p_endptr,10);
            if (*p_endptr != '\0')
                return XXX_E_ODL_INVALID_DATA_TYPE;
            break;

        case XXX_Int:
            (*(int *)p_destination) = (int)strtol(kvalue,&p_endptr,10);
            if (*p_endptr != '\0')
                return XXX_E_ODL_INVALID_DATA_TYPE;
            break;

        case XXX_Float:
            (*(float *)p_destination) = (float)strtod(kvalue,&p_endptr);
            if (*p_endptr != '\0')
                return XXX_E_ODL_INVALID_DATA_TYPE;
            break;

        case XXX_Sci_Not:
        case XXX_Double:
            (*(double *)p_destination) = strtod(kvalue,&p_endptr);
            if (*p_endptr != '\0')
                return XXX_E_ODL_INVALID_DATA_TYPE;
            break;

        case XXX_ArrayOfString:
            p_a = malloc (strlen(kvalue)+1);
            if (p_a == NULL)
            {
                strcpy(p_ErrorMessage, "Error : Allocating Memory ... ");
                return (ERROR);
            }
            strcpy(p_a,kvalue);
            (*(char **)p_destination) = p_a;
            break;

        default:
            strcpy(p_ErrorMessage, "Invalid Conversion type ");
            strcat(p_ErrorMessage, kvalue);
            strcat(p_ErrorMessage, " ...");
            return ERROR;
    }

    return SUCCESS;
}
