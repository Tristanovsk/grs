/* Feature test switches */
#undef _POSIX_SOURCE

#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <xxx_FileLock.h>
#include <xxx_Types.h>
#include <xdb_Defines.h>

#include <xxx_CloseUnmap.h>
#include <xxx_OpenMap.h>

/****************************************************************************/
/**
 * @file xxx_OpenMap.c
 * @brief To open a file, get its length, and optionally memory map the file.
 * @ingroup iasLib
 */
/**
 * @brief To open a file, get its length, and optionally memory map the file.
 *
 * Assume the file_start and p_file_len values are correct if the user is not
 * defaulting.
 *
 * @return XXX_OPEN_ERR: File open error
 * @return XXX_STAT_ERR: File status error
 * @return XXX_MMAP_ERR: Memory map error
 * @return XXX_MMAP_ZERO_ERR: Indicates zero file size for memory map error
 * @return SUCCESS: Terminated successfully
 */
/****************************************************************************/

int xxx_OpenMap
(
    char *p_filename,   //!<[in] Name of file to open 
    int *p_fd,          //!<[out] File descriptor 
    int file_options,   //!<[in] File options for opening 
    void **p_pa,        //!<[out] Memory map address 
    off_t file_start,   //!<[out] Start offset of file or -1 for map enter file 
    long *p_file_len,   //!<[out] File length 
    long *p_remainder,  //!<[out] Byte offset into the page, used only if not
                        // defaulting 
    boolean mapflag,    //!<[in] Mapping indicator 
    boolean lockflag,   //!<[in] locking indicator
    char *p_ErrorMessage//!<[out] Error Message 
)
{
    struct  stat buf;
    int flags;
    int org_file_options = file_options;
    char    ErrorMessage[XDB_ORAERRMSGLEN];

    
    /* open file and return error if one occurred */

    if (lockflag)
    {
        if (!((file_options & O_RDWR) == O_RDWR))
        {
            if ((file_options & O_WRONLY) == O_WRONLY)
                file_options = (file_options ^ O_WRONLY); /* clear out
                                                             O_WRONLY option*/
            file_options = file_options | O_RDWR; /* set to O_RDWR */
        }
    }

    if ((*p_fd = open(p_filename, file_options, 0664)) < 0)
    {
        strcpy(p_ErrorMessage, strerror(errno));
        strcat(p_ErrorMessage, ": ");
        strcat(p_ErrorMessage, p_filename);
        return(XXX_OPEN_ERR);
    }

    if (lockflag)
        XXX_WRITE_WAIT_LOCK(*p_fd, 0, SEEK_SET, 0);

    /* call fstat to get file status for file size. Return error if one
       occurred */

    if (fstat(*p_fd,&buf) < 0)
    {
        strcpy(p_ErrorMessage, strerror(errno));
        strcat(p_ErrorMessage, ": ");
        strcat(p_ErrorMessage, p_filename);
        xxx_CloseUnmap(*p_fd, (void *)NULL, (long)NULL, FALSE, lockflag,
                       ErrorMessage);
        *p_fd = -1;
        return(XXX_STAT_ERR);
    }

    /* memory map the file if requested. Return error if one occurred */

    if (mapflag)
    {
        if (buf.st_size == 0)
        {
            strcpy(p_ErrorMessage,
                   "Can't map file because its file size is zero: ");
            strcat(p_ErrorMessage, p_filename);
            xxx_CloseUnmap(*p_fd, (void *)NULL, (long)NULL, FALSE, lockflag,
                           ErrorMessage);
            *p_fd = -1;
            return(XXX_MMAP_ZERO_ERR);
        }

        if (file_start == -1)  /* map enter file */
        {
            file_start = 0;
            *p_file_len = buf.st_size;
            if (p_remainder != NULL)
                *p_remainder = 0;
        }
        else
        {
            /* turn into multiple of page size */
            *p_remainder =  file_start % (off_t) sysconf(_SC_PAGESIZE);
            *p_file_len += *p_remainder;
            file_start = (file_start / (off_t) sysconf(_SC_PAGESIZE)) * 
                (off_t) sysconf(_SC_PAGESIZE);
        }

        /* Memory map file as it was opened */
        if ((org_file_options & O_RDWR) == O_RDWR)
            flags = PROT_READ|PROT_WRITE;
        else if ((org_file_options & O_WRONLY) == O_WRONLY)
            flags = PROT_WRITE;
        else
            flags = PROT_READ;

        if ((*p_pa=mmap((void *)NULL, (size_t)*p_file_len, flags,
                        MAP_SHARED, *p_fd, file_start)) == MAP_FAILED)
        {
            strcpy(p_ErrorMessage, strerror(errno));
            strcat(p_ErrorMessage, ": ");
            strcat(p_ErrorMessage, p_filename);
            xxx_CloseUnmap(*p_fd, (void *)NULL, (long)NULL, FALSE, lockflag,
                           ErrorMessage);
            *p_fd = -1;
            return(XXX_MMAP_ERR);
        }
    }
    else
        *p_file_len = buf.st_size;

    return SUCCESS;
}
