/****************************************************************************/
/**
 * @file xxx_LogError.c
 * @brief Write error notification to the IAS software log.  
 * @ingroup iasLib
 */
/****************************************************************************/

#include <errno.h>
#include <fcntl.h>
#include <malloc.h>
#include <netdb.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <syslog.h>
#include <stdlib.h>
#include <sys/param.h>
#include <sys/stat.h>

#include <xxx_Types.h>

/* Macros */
#define NAMESIZE 80
#define SYSLOGFORMAT "%3s   %-50s    %6d    %s\n" /*syslog handles time & pid*/

#ifdef LPGS
/*LPGS log specifics*/
#define LOG_NAME  "LPGS-Log"

#define Header1 "                                                LPGS local software log\n\n"
#define Header2 "                      Task Process id                  Unit Name                         Line Number  Error Message\n"
#define Header3 "                      ---- ---------- -------------------------------------------------- -----------  -------------\n"
#define Dummy   "10-Jun-2000 34:34:34  PSC   1234567   12345678911234567892123456789312345678941234567895    123456    " 

#define TIMESTAMPSIZE 23
#define TIMEFORMAT "%d-%b-%Y %T  "
#define LOGFORMAT "%22s%3s   %-7ld   %-50s    %6d    %s\n"

#else
/*IAS log specifics*/
#define LOG_NAME  "IAS-Log"
    
#define Header1 "                                              IAS local software log\n\n"
#define Header2 "                   Task Process id                  Unit Name                         Line Number  Error Message\n"
#define Header3 "                   ---- ---------- -------------------------------------------------- -----------  -------------\n"
#define Dummy   "1997 192 34:34:34  PSC   1234567   12345678911234567892123456789312345678941234567895    123456    " 
    
#define TIMESTAMPSIZE 20
#define TIMEFORMAT "%Y %j %H:%M:%S  "
#define LOGFORMAT "%19s%3s   %-7ld   %-50s    %6d    %s\n"
#endif

/* Global variables */
char xxx_LogErrorPath[MAXPATHLEN];
char xxx_TaskName[4];

/* External functions */
#include <xxx_FileLock.h>
#include <xxx_GetTime.h>
#include <xxx_LogError.h>

/****************************************************************************/
/**
 * @brief Write error notification to the IAS software log.  
 * 
 * This function is meant to be called by Process Control (PCS) and 
 * Data Management (DMS) subsystems and NOT Radiometry (RPS) and Geometry 
 * (GPS) subsystems.
 *
 * The task must set the xxx_TaskName string to the task executable, e.g. 
 * "DID".
 *
 * The local IAS software log file is created in the local directory and is
 * IAS-log.
 *
 * The local LPGS software log file is created in the local directory and is
 * LPGS-log.
 *
 * If any open, malloc, or write error occurs, the unit will return immediately 
 * with a ERROR status.
 *
 * @return SUCCESS: Wrote error message to log
 * @return ERROR: Failed to write error message to log
 */
/****************************************************************************/

int xxx_LogError
(
    char *p_UnitName,    //!<[in] Unit name 
    int LineNumber,      //!<[in] Line number 
    char *p_ErrorMessage //!<[out] Error message 
)
{
    char TimeStamp[TIMESTAMPSIZE];

    /* Determine the logging style. 
     * For IAS, default is File logging.
     * For LPGS, default is syslog logging.
     * Unless overridden by the LEVEL1_LOGGING environment variable. */
    enum { xxx_Logging_SYSLOG, xxx_Logging_FILE } logging_style;
#ifdef LPGS
    logging_style = xxx_Logging_SYSLOG;
#else
    logging_style = xxx_Logging_FILE;
#endif
    char *logging = getenv("LEVEL1_LOGGING");
    if (logging != NULL)
    {
        if (strcasecmp(logging, "FILE") == 0) 
            logging_style = xxx_Logging_FILE;
        if (strcasecmp(logging, "SYSLOG") == 0) 
            logging_style = xxx_Logging_SYSLOG;
    }

    if (logging_style == xxx_Logging_SYSLOG)
    {
        /* The LOG_LOCAL1 facility must match the facility selector 
         * in the rsyslog.conf file */
        openlog(getenv("USER"), LOG_PID, LOG_LOCAL1);
        if (p_ErrorMessage != (char *)NULL)
        {
            syslog(LOG_LOCAL1|LOG_INFO, SYSLOGFORMAT, xxx_TaskName, 
                p_UnitName, LineNumber, p_ErrorMessage);
        }

    }
    else /* logging_style == xxx_Logging_FILE */
    {
        int fdLog;
        char LocalLogFilename[NAMESIZE];

        struct stat info;
        char * Buffer;
        long    status;
        char Cantmalic[150];
        int Not_Lock;

        strcpy(LocalLogFilename, xxx_LogErrorPath);
        strcat(LocalLogFilename, LOG_NAME);

        if (stat(LocalLogFilename, &info) == -1)
        {
            if ((fdLog = open(LocalLogFilename, O_CREAT | O_WRONLY,
                        S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH |
                        S_IWOTH  )) == -1)
                return ERROR;

            if ((Buffer = (char *)malloc(strlen(Header1)+strlen(Header2)+
                        strlen(Header3)+1)) == (char *)NULL)
            {
                xxx_GetTime(TimeStamp, sizeof(TimeStamp), TIMEFORMAT);
                sprintf(Cantmalic, LOGFORMAT, TimeStamp, xxx_TaskName,
                    (long)getpid(), __FILE__, __LINE__,
                    "Can't malloc memory for write");
                while ((status = write(fdLog, Cantmalic, strlen(Cantmalic))) == -1
                    && errno == EINTR);
                close(fdLog);
                return ERROR;
            }

            strcpy(Buffer, Header1);
            strcat(Buffer, Header2);
            strcat(Buffer, Header3);

            while ((status = write(fdLog, Buffer, strlen(Buffer))) == -1 &&
                errno == EINTR);

            free(Buffer);
            if (status < 0) 
            {
                close(fdLog);
                return ERROR;
            }
        }
        else
        {
            if ((fdLog = open(LocalLogFilename, O_WRONLY | O_APPEND)) == -1)
                return ERROR;
        }

        Not_Lock = SUCCESS;
        while (Not_Lock)
        {   
            if ((Not_Lock = XXX_WRITE_LOCK(fdLog, 0, SEEK_SET, 0)))
            {
                close(fdLog);
                sleep(1);
                if ((fdLog = open(LocalLogFilename, O_WRONLY | O_APPEND)) == -1)
                    return ERROR;
            }
        }

        if (p_ErrorMessage != (char *)NULL)
        {
            /* the extra 100 bytes are just in case of an overflow in one of
               the other fields */
            if ((Buffer =
                    (char *)malloc(strlen(Dummy)+strlen(p_ErrorMessage)+2+100))
                == (char *)NULL)
            {
                xxx_GetTime(TimeStamp, sizeof(TimeStamp), TIMEFORMAT);
                sprintf(Cantmalic, LOGFORMAT, TimeStamp, xxx_TaskName,
                    (long)getpid(), __FILE__, __LINE__,
                    "Can't malloc memory for write");
                while ((status = write(fdLog, Cantmalic, strlen(Cantmalic))) == -1
                    && errno == EINTR);
                XXX_UNLOCK(fdLog, 0, SEEK_SET, 0);
                close(fdLog);
                return ERROR;
            }

            xxx_GetTime(TimeStamp, sizeof(TimeStamp), TIMEFORMAT);
            sprintf(Buffer, LOGFORMAT, TimeStamp, xxx_TaskName, (long)getpid(),
                p_UnitName, LineNumber, p_ErrorMessage);

            while ((status = write(fdLog, Buffer, strlen(Buffer))) == -1 &&
                errno == EINTR);

            free(Buffer);
            if (status < 0)
            {
                XXX_UNLOCK(fdLog, 0, SEEK_SET, 0);
                close(fdLog);
                return ERROR;
            }
        }

        XXX_UNLOCK(fdLog, 0, SEEK_SET, 0);
        close(fdLog);
    }

    return SUCCESS;
}
