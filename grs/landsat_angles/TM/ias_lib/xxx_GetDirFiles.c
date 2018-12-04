#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include <xxx_Types.h>

#include <xdb_Defines.h>
#include <xxx_OpenMap.h>
#include <xxx_CloseUnmap.h>
#include <xxx_GetTempName.h>
#include <xxx_GetDirFiles.h>
#include <xxx_RecursiveDeletion.h>
#include <xxx_SysCommands.h>

/****************************************************************************/
/**
 * @file xxx_GetDirFiles.c
 * @brief Returns a directory listing of the target directory.
 * @ingroup iasLib
 */
/**
 * @brief Returns a directory listing of the target directory.
 *
 * @return positive: Pointer to the xxx_DirList_TYPE
 * @return NULL: Failure
 */
/****************************************************************************/


xxx_DirList_TYPE *xxx_GetDirFiles
(
    int *p_Count,        //!<[out] Number of rows in p_PartionList
    char *p_PathName,    //!<[in] Directory path name
    char *p_Key,         //!<[in] List options like "*"
    char *p_ErrorMessage //!<[out] Error message
)
{
    char    Command[400];
    int     status;
    int     FileId;                 /* file descriptor for CPF file */
    long    FileLength;             /* size of CPF file */

    char    *p_TempName;
    xxx_DirList_TYPE *p_Listing = NULL;
    xxx_DirList_TYPE *p_TempPtr;
    char    ErrorMessage[XDB_ORAERRMSGLEN];
    char    line[MAXPATHLEN];
    char    *p_Ptr;

    
    *p_Count = 0;
    /* jan2018 harmel: ./ instead of /usr/tmp/ */
    if ((p_TempName = xxx_GetTempName("./", "lsTmp")) == NULL)
    {
        strcpy(p_ErrorMessage, "Can't generate temporary path name\n");
        return(NULL);
    }

    if (p_PathName[strlen(p_PathName)-1] != '/')
        sprintf(Command, "%s %s/%s 1>%s 2>/dev/null", LS_COMMAND, p_PathName,
                p_Key, p_TempName);
    else
        sprintf(Command, "%s %s%s 1>%s 2>/dev/null", LS_COMMAND, p_PathName,
                p_Key, p_TempName);

    if ((WEXITSTATUS(status = system(Command))) != 0)
    {
        sprintf(p_ErrorMessage, "%s: bad command, status: %d, quitting...\n",
                Command, status);
        goto ERROR_PROC;
    }

    if (xxx_OpenMap(p_TempName, &FileId, O_RDONLY, (void *)NULL, -1, 
                    &FileLength, (long *)NULL, FALSE, FALSE, p_ErrorMessage))
        goto ERROR_PROC;

    p_Ptr = line;
    while(read(FileId, p_Ptr, 1) == 1)
    {
        if (*p_Ptr == '\0' || *p_Ptr == '\n')
        {
            *p_Ptr = '\0';

            (*p_Count)++;
            if ((p_TempPtr =
                 (xxx_DirList_TYPE *)realloc(
                     p_Listing, (*p_Count)*sizeof(xxx_DirList_TYPE))) == NULL)
            {
                strcpy(p_ErrorMessage, "Can't realloc memory");
                free(p_Listing);
                goto ERROR_PROC;
            }
            p_Listing = p_TempPtr;

            strcpy(p_Listing[(*p_Count)-1].FilePathName, line);

            p_Ptr = line;
        }
        else
            p_Ptr++;
    }

    xxx_CloseUnmap(FileId, (void *)NULL, FileLength, FALSE, FALSE,
                   p_ErrorMessage);

    if (xxx_RecursiveDeletion(p_TempName, p_ErrorMessage))
    {
        free(p_Listing);
        goto ERROR_PROC2;
    }

    return p_Listing;

  ERROR_PROC:
    xxx_RecursiveDeletion(p_TempName, ErrorMessage);
    
  ERROR_PROC2:
    *p_Count = 0;
    return NULL;
}
