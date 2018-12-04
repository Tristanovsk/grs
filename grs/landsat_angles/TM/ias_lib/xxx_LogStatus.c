#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/param.h>
#include <sys/stat.h>

#include "xxx_GetTime.h"
#include "xxx_Types.h"
#include "xxx_LogStatus.h"

#ifdef LPS
void lps_LogMessage(int, char *, int, char *);
#endif

#ifdef LPGS
/*LPGS specifics*/
#define TIMESTAMPSIZE 23
#define TIMEFORMAT "%d-%b-%Y %T  "
#define LOGFORMAT "%22s%3s   %-7ld   %-20s    %6d    %s\n"
#else
/*IAS specifics*/
#define TIMESTAMPSIZE 20
#define TIMEFORMAT "%Y %j %H:%M:%S  "
#define LOGFORMAT "%19s%3s   %-7ld   %-20s    %6d    %s\n"
#endif

/****************************************************************************/
/**
 * @file xxx_LogStatus.c
 * @brief Echo an application's status message to standard out.  
 * @ingroup iasLib
 */
/**
 * @brief Echo an application's status message to standard out.  
 *
 * This function is intended to be called by Radiometry (RPS) and Geometry(GPS)
 * functions only and NOT the Process Control (PCS) and Data Management (DMS) 
 * subsystems.
 */        
/****************************************************************************/

void xxx_LogStatus
(
    const char *p_program,  //!<[in] Program name 
    const char *p_function, //!<[in] Filename of calling unit 
                            // (see C macro __FILE__) 
    int LineNumber,         //!<[in] Line number of calling unit 
                            // (see C macro __LINE__) 
    const char *p_msg       //!<[in] Error message to write to standard out 
)
{
    char TimeStamp[TIMESTAMPSIZE];

#ifdef LPS    
    lps_LogMessage(6,p_function, LineNumber, p_msg);
#else
    xxx_GetTime(TimeStamp, sizeof(TimeStamp), TIMEFORMAT);
    printf(LOGFORMAT, TimeStamp, p_program, (long)getpid(), p_function,
           LineNumber, p_msg);
#endif
}
