#include <ctype.h>            /* toupper prototype           */
#include "xxx_string_utils.h" /* prototype for this function */

/****************************************************************************/
/**
 * @file xxx_strtoupper.c
 * @brief xxx_strtoupper converts all the lower case characters in a string 
 * to upper case.
 * @ingroup iasLib
 */
/**
 * @brief xxx_strtoupper converts all the lower case characters in a string 
 * to upper case.
 *
 * The original string is converted in place.
 *
 * @return pointer: Returns a pointer to the string so it can immediately 
 * be used as a parameter to another function (such as strcmp)
 */
/****************************************************************************/


char *xxx_strtoupper 
(
    char *string_ptr  //!<[in/out] pointer to string to convert   
)
{
    char *c_ptr = string_ptr;

    
    while (*c_ptr != '\0')
    {
        *c_ptr = toupper(*c_ptr);
        c_ptr++;
    }

    return string_ptr;
}
