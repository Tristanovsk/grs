/****************************************************************************/
/**
 * @file xxx_GetTime.c
 * @brief Three different functions dealing with creating time stamps, 
 * returning the system time, and converting seconds since mission epoch
 * to system time
 * @ingroup iasLib
 */
/****************************************************************************/

#include <time.h>

#include "xxx_Types.h"
#include "xxx_GetTime.h"

/****************************************************************************/
/**
 * @brief To create a time and/or date stamp.
 *
 * The input p_format will contain the format specification of the 
 * time stamp.  Format specifications can be found in the UNIX manual
 * page of strftime(3).
 *
 * Examples:
 *             "%H:%M:%S"     14:55:20
 *             "%X on %x"     4:25 PM on 7/13/93
 *             "%c"           Sun Apr 19 16:10:20 1993
 *
 * The input p_stamp must be large enough to hold the result of the
 * format specification.
 *
 * @return SUCCESS: Successfully created time stamp.
 * @return ERROR: Error creating time stamp.
 */
/****************************************************************************/

int xxx_GetTime
(
    char *p_stamp,        //!<[in/out] Pointer to string to place time/date
    int stampsize,        //!<[in] Size of p_stamp
    const char *p_format  //!<[in] Pointer to format string
)
{
    time_t   ptime;               /*  Time in seconds  */

    
    ptime = time((time_t *) 0);

    if (strftime(p_stamp, stampsize, p_format, localtime(&ptime)) <= 0)
        return ERROR;

    return SUCCESS;
}

/****************************************************************************/
/**
 * @brief Returns system time (GMT) in the format "MM/DD/YYYY HH:MM:SS"
 */
/****************************************************************************/

void xxx_GetTime_system
(
    char *s  //!<[out] formatted time string 
)
{
    time_t curtime = time(NULL);
    strftime( s, TIMELEN, "%m/%d/%Y %H:%M:%S", gmtime(&curtime) );
}

/*****************************************************************************/
/**
 * @brief Converts seconds since mission epoch to system time (GMT) in the
 * format "MM/DD/YYYY HH:MM:SS"
 */
/****************************************************************************/

void xxx_psfTime
(
    int  sat,               //!<[in] Satellite
    double sec,             //!<[in] Seconds since mission epoch 
    char s[TIMELEN]         //!<[out] formatted time string: 
                            // MM/DD/YYYY HH:MM:SS
)
{
    time_t sec70;
    
    switch (sat)
    {
        case 7: /* epoch = 00:00 Jan 1, 1993 */
            sec70 = sec + 725846400.0;
            break;
        case 5: /* epoch = 00:00 Jan 6, 1980 */
        case 4:
            sec70 = sec + 315964800.0;
            break;
        case 1: /* epoch = 00:00 Jan 1, 1970 */
        case 2:
        case 3:
            sec70 = sec;
            break;
        default: /* shouldn't be here */
            sec70 = 0.0;
    }

    strftime(s, TIMELEN, "%m/%d/%Y %H:%M:%S", gmtime(&sec70));
}

