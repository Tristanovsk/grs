#include <stdio.h>
#include <string.h>
#include <sys/param.h>
#include <libgen.h>

#include <lablib3.h>                 /* prototypes for the ODL functions */
#include "xxx_Errno.h"
#include "xxx_GetDirFiles.h"
#include "xxx_ODL.h"                 /* prototype for the xxx ODL functions */
#include "xxx_Types.h"
#include "xxx_GetTempName.h"
#include "xxx_RecursiveDeletion.h"

extern char ODLErrorMessage[];       /* External Variables */

/****************************************************************************/
/**
 * @file xxx_OpenODL.c
 * @brief Open and parse the ODL file.
 * @ingroup iasLib
 */
/**
 * @brief Open and parse the ODL file.
 *
 * The MTA files are shared by L1R and L0R products.  If the user requested a
 * LPS MTA 1 or 2 file, the code first sees if the L0R version is there (Only
 * the eas readMetaL1 dumper reads L1R products).  If a match is not found, the
 * L1R version is tried.  To make sure the L0R naming convertion is not
 * confused with the L1R naming convertion, the L0R pattern match makes sure 
 * the 10 character in the file name is NOT an '_' character as it is in the
 * L1R naming convertion, but a 0->9 number.
 *
 * @return positive: Pointer to ODL tree 
 * @return NULL: Failure
 */
/****************************************************************************/

OBJDESC *xxx_OpenODL
(
    char *p_ODLFile,         //!<[in] If a normal odl file is to be opened, the
                             // path is to the ODL file otherwise the path
                             // is to the L0R directory file
    xxx_EDCOdl_TYPE ODLType, //!<[in] Indicates if a ODL file to read or 
                             // a normal odl file (ie.xxx_NoEDCMeta) 
    char *p_ErrorMessage     //!<[out] Error message 
)
{
    OBJDESC *p_lp= NULL;
    char ODLPathName[MAXPATHLEN];
    char ODLErrorFile[MAXPATHLEN];
    char *p_TempName;
    int Count;
    int i;
    xxx_DirList_TYPE *p_DirectoryList;
    char        tmpDirname[PATH_MAX];/* used so dirname does not modify path*/

    typedef struct xxx_MetaIndex_TYPE
    {
        xxx_EDCOdl_TYPE  Want;
        char ProductMeta[24];
    } xxx_MetaIndex_TYPE;


    /* The following array has been modified to make it more generic.  This
       will work fine for Landsat products but will most likely break for
       other sensors.  'L7' was changed in each of the strings to 'L?'.
       --------------------------------------------------------------------*/
    xxx_MetaIndex_TYPE ProductMetaData[] = {
        {xxx_EDCMeta,  "L?*MTP*"},
        {xxx_LPSMeta0, "L?????0??[0-9]*_MTA.*"}, /* Make sure the 10 character
                                                    is not '_' as in L1R
                                                    product */
        {xxx_LPSMeta0, "L?0*_MTA.L1R"},
        {xxx_LPSMeta1, "L?????1??[0-9]*_MTA.*"}, /* Make sure the 10 character
                                                    is not '_' as in L1R
                                                    product */
        {xxx_LPSMeta1, "L?1*_MTA.L1R"},
        {xxx_LPSMeta2, "L?????2??[0-9]*_MTA.*"}, /* Make sure the 10 character
                                                    is not '_' as in L1R
                                                    product */
        {xxx_LPSMeta2, "L?2*_MTA.L1R"},
        {xxx_CPFMeta,  "L?*CPF*"},
        {xxx_LPGSMeta, "L?*MTL*"},
    };

    printf("harmel %s",p_ODLFile);    
    ODLErrorMessage[0] = '\0';

    if(ODLType != xxx_NoEDCMeta)
    {
        /* Use tmpdirname to ensure that path is not modified by dirname */
        strcpy(tmpDirname, p_ODLFile);
        strcpy(ODLPathName,dirname(tmpDirname));
        printf("harmel %s",p_ODLFile);
        for (i=0;  i<(int)XXX_ARRAY_SIZE(ProductMetaData) &&
                 ProductMetaData[i].Want != ODLType;  i++);
        p_DirectoryList = xxx_GetDirFiles(&Count, ODLPathName,
                                          ProductMetaData[i].ProductMeta,
                                          p_ErrorMessage);

        /* Check for other product type */
        if (p_DirectoryList == NULL &&
            (ODLType == xxx_LPSMeta1 || ODLType == xxx_LPSMeta2))
        {
            i++;
            p_DirectoryList = xxx_GetDirFiles(&Count, ODLPathName,
                                              ProductMetaData[i].ProductMeta,
                                              p_ErrorMessage);

            /* If we're still failing, the data type could be TM-A, which
               has a '0' instead of a '1' for the format in the base
               filename.  Give it a shot. */
            if (p_DirectoryList == NULL)
            {
                for (i=0;  i<(int)XXX_ARRAY_SIZE(ProductMetaData) &&
                         ProductMetaData[i].Want != xxx_LPSMeta0;  i++);
                p_DirectoryList = xxx_GetDirFiles(&Count, ODLPathName,
                                               ProductMetaData[i].ProductMeta,
                                                  p_ErrorMessage);
                if (p_DirectoryList == NULL)
                {
                    i++;
                    p_DirectoryList = xxx_GetDirFiles(&Count, ODLPathName,
                                               ProductMetaData[i].ProductMeta,
                                                      p_ErrorMessage);
                }
            }
        }

        if (p_DirectoryList == NULL)
        {
            xxx_errno = XXX_E_FILE_NOT_FOUND;
            return NULL;
        }
        strcpy(ODLPathName, p_DirectoryList[0].FilePathName);
        free (p_DirectoryList);
    }
    else
        strcpy(ODLPathName, p_ODLFile);

    /* Build ODL error file name 
     * jan2018 harmel: ./ instead of /usr/tmp/ */
    if ( (p_TempName =  xxx_GetTempName("./", "ODLEr")) == NULL)
    {
        strcpy(p_ErrorMessage, "Can't generate temporary path name\n");
        xxx_errno = XXX_E_NO_MEMORY;
        return(NULL);
    }
    sprintf(ODLErrorFile,"%s.dat", p_TempName);

    /*   open and parse ODL file - error messages to ODLPathName */
    if ((p_lp = OdlParseLabelFile(ODLPathName, ODLErrorFile, ODL_NOEXPAND,
                                  FALSE)) == NULL)
    {
        sprintf(p_ErrorMessage, "ODL Error on File: %s:  %s", ODLPathName,
                ODLErrorMessage);
        xxx_errno = XXX_E_ODL_SYNTAX;
        return NULL;
    }

    /* Check of ODL syntax errors */
    if (ODLErrorMessage[0] != '\0')
    {
        sprintf(p_ErrorMessage,"%s: check out %s", ODLErrorMessage,
                ODLErrorFile);
        OdlFreeTree(p_lp);
        xxx_errno = XXX_E_ODL_SYNTAX;
        return NULL;
    }

    /* Delete the ODL error log file */
    xxx_RecursiveDeletion(ODLErrorFile, p_ErrorMessage);

    return p_lp;
}
