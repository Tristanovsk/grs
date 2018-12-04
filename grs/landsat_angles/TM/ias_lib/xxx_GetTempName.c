#undef _POSIX_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include "xxx_Const.h"
#include "xxx_Types.h"

/****************************************************************************/
/**
 * @file xxx_GetTempName.c
 * @brief Create a name for a temporary file using the non-POSIX tempnam()
 * system call.
 * @ingroup iasLib
 */
/**
 * @brief Create a name for a temporary file using the non-POSIX tempnam()
 * system call.
 *
 * @return NULL: Fail to generate name            
 * @return NON-NULL: Generated name
 */
/****************************************************************************/

char *xxx_GetTempName
(
    const char *p_dir,  //!<[in] Name of the directory in which the file is to
                        // be created
    const char *p_pfx   //!<[in] Favorite initial letter sequences 
                        // (up to 5 characters)
)
{
    struct timeval tv;
    static char TempName[STRLEN] = "";
    char Temp[STRLEN] = "";

    gettimeofday(&tv, NULL);   /* initalizes the structure */
 
    if (p_pfx != NULL && p_dir != NULL)
    {
        strcpy(TempName, p_dir);
        strcat(TempName, "/");
        strcat(TempName, p_pfx);
        sprintf(Temp, "%ld%ld%ld", (long)tv.tv_sec, (long)tv.tv_usec, 
                                   (long)getpid());
        strcat(TempName, Temp);
    } 
    else
    {
        return NULL;   
    }
 
    return TempName;
}
