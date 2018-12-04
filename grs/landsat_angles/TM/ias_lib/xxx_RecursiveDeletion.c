#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <sys/wait.h>

#include <xxx_Types.h>
#include <xxx_RecursiveDeletion.h>
#include <xxx_SysCommands.h>

/****************************************************************************/
/**
 * @file xxx_RecursiveDeletion.c
 * @brief Delete the Pathname using UNIX rm -rf command.
 * @ingroup iasLib
 */
/**
 * @brief Delete the Pathname using UNIX rm -rf command.
 *
 * @return SUCCESS: Pathname deleted
 * @return ERROR: I/O error deleting pathname
 */
/****************************************************************************/


int xxx_RecursiveDeletion
(
    char *p_Pathname,    //!<[in] Path name 
    char *p_ErrorMessage //!<[out] Error message 
)
{
    char DeleteCommand[MAXPATHLEN];
    int  status;

    
    sprintf(DeleteCommand, "%s -fr %s 2>/dev/null 1>/dev/null", RM_COMMAND,
            p_Pathname);
    if ((WEXITSTATUS(status = system(DeleteCommand))) != 0)
    {
        sprintf(p_ErrorMessage, "Delete failed for %s, status was: %d\n",
                p_Pathname, status);
        return ERROR;
    }

    return SUCCESS;
}
