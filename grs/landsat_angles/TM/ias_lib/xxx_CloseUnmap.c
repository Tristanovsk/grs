
#undef _POSIX_SOURCE

#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>

#include <xxx_FileLock.h>
#include <xxx_Types.h>

#include <xxx_CloseUnmap.h>

/****************************************************************************/
/**
 * @file xxx_CloseUnmap.c
 * @brief To close a file and optionally unmemory map the file.
 * @ingroup iasLib
 */
/**
 * @brief To close a file and optionally unmemory map the file.
 *
 * @return XXX_UNMAP_ERR: Indicates unmap error
 * @return XXX_CLOSE_ERR: Indicates file close error
 * @return SUCCESS: Terminated successfully
 */
/****************************************************************************/

int xxx_CloseUnmap
(
    int fd,              /* I: File descriptor */
    void *p_pa,          /* I: Memory map address */
    long file_len,       /* I: File length */
    boolean mapflag,     /* I: Mapping indicator */
    boolean lockflag,    /* I: locking indicator */
    char *p_ErrorMessage /* O: Error Message */
)
{
    int status = SUCCESS;   /* return status */

    /* memory map the file if requested. Return error if one occurred */

    if (mapflag)
    {
        if(munmap(p_pa, file_len) == -1)
        {
            strcpy(p_ErrorMessage, strerror(errno));
            status = XXX_UNMAP_ERR;
        }
    }

    if (lockflag)
        XXX_UNLOCK(fd, 0, SEEK_SET, 0);

    /* close file and return error if one occurred */

    if(close(fd) < 0)
    {
        strcpy(p_ErrorMessage, strerror(errno));
        status = XXX_CLOSE_ERR;
    }

    return status;
}
