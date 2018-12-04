
#include <fcntl.h>

#include <xxx_Types.h>

#include <xxx_FileLock.h>

/****************************************************************************/
/**
 * @file xxx_FileLock.c
 * @brief To lock/unlock files.
 * @ingroup iasLib
 */
/**
 * @brief To lock/unlock files.
 *
 *
 * On a no-wait error, check for EACCES or EAGAIN to see if file is locked 
 * already.
 *
 * The read-lock will prevent another process from write-locking the file.
 * This does not stop another process from reading from the file.
 *
 * There is a possibility of deadlocking if you lock a file and then
 * try to lock another file in the wait mode that is already lock.  So be
 * careful.
 *
 * @return >= 0: Successful
 * @return -1: On failure, fcntl() returns -1 and sets errno to indicate the
 * error.
 */
/****************************************************************************/

int xxx_FileLock
(
    int fd,       //!<[in] File descriptor
    int command,  //!<[in] fcntl commands: F_SETLK or F_SETLK
    short type,   //!<[in] Type of lock desired: F_RDLCK, FWRLCK, or F_UNLOCK
    off_t offset, //!<[in] Byte offset, relative to l_whence
    short whence, //!<[in] SEEK_SET, SEEK_CUR, SEEK_END
    off_t length  //!<[in] Number of bytes (0 means to EOF)
)
{
    struct flock lock;

    
    lock.l_type = type;     /* F_RDLCK, FWRLCK, F_UNLOCK */
    lock.l_start = offset;  /* byte offset, relative to l_whence */
    lock.l_whence = whence; /* SEEK_SET, SEEK_CUR, SEEK_END */
    lock.l_len = length;    /* number of bytes (0 means to EOF */
    
    return (fcntl(fd, command, &lock));
}   
