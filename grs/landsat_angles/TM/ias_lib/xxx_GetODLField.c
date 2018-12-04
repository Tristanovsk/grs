#include <stdio.h>
#include <string.h>
#include <toolbox.h>
#include "xxx_ODL.h"                 /* prototype for the xxx ODL functions */

extern char ODLErrorMessage[];       //!< External Variables

/****************************************************************************/
/**
 * @file xxx_GetODLField.c
 * @brief Get the requested ODL field and convert it, if necessary.
 * @ingroup iasLib
 */
/**
 * @brief Get the requested ODL field and convert it, if necessary.
 *
 * Caller must pass in a pointer to a memory location large enough to contain
 * the converted field value(s).  The unit will then convert the ascii values
 * to the specified value types.
 * 
 * If the caller is expecting an array of strings, the p_Memory array is
 * an array of char *.  The caller will have to free the malloc memory.
 * 
 * If any of the ODL calls fails the system will build the error in the
 * p_ErrorMessage, free memory (if any), and return error to the calling
 * function.
 *
 * @return SUCCESS: Found and converted field(s) into requested type
 * @return XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED: Not enough memory passed in
 * @return XXX_E_ODL_NOT_FOUND: Group/label not found
 * @return XXX_E_ODL_INVALID_DATA_TYPE: Data type mismatch
 * @return ERROR: Fatal error                       
 */
/****************************************************************************/


int xxx_GetODLField
(
    void *p_MemoryAddr,         //!<[in] p_Memory to copy fields into
    int MemorySize,             //!<[in] Size in bytes of p_Memory 
    xxx_ValType_TYPE ValueType, //!<[in] What type the P_FieldName is, enum
                                // field 
    OBJDESC *p_ODLTree,         //!<[in] ODL tree
    char *p_ClassName,          //!<[in] Group/Object name, optional
    char *p_LabelName,          //!<[in] Field to retrieve
    char *p_ErrorMessage,       //!<[out] Error Message
    int *p_Count                //!<[out] Count the number of values in a array 
)
{

    OBJDESC *p_lp;         /* Object Descriptor */
    KEYWORD *p_kw;         /* Keyword Name */
    char *p_kwv;           /* Keyword Value */
    char *p_keyword;       /* Copy of the Keyword Value */
    int i;
    char *p_word;
    char string[TB_MAXLINE];
    int ret_code = 0;

    
    *p_Count = 0;

    if (p_LabelName == NULL || strlen(p_LabelName) == 0)
    {
        strcpy(p_ErrorMessage,"Error : Invalid Keyword");
        return ERROR;
    }

    if ((p_lp = OdlFindObjDesc(p_ODLTree, p_ClassName, p_LabelName, NULL, 1,
                               ODL_RECURSIVE_DOWN)) == NULL)
    {
        strcpy(p_ErrorMessage, ODLErrorMessage);
        if ((long)strlen(p_ErrorMessage) <= 1)
            sprintf(p_ErrorMessage,
                    "Error : %s Object not found with %s keyword",
                    p_ClassName, p_LabelName);
        return XXX_E_ODL_NOT_FOUND;
    }

    if ((p_kw = OdlFindKwd(p_lp, p_LabelName, NULL, 1,
                           ODL_RECURSIVE_DOWN)) == NULL)
    {
        strcpy(p_ErrorMessage, ODLErrorMessage);
        if ((long)strlen(ODLErrorMessage) <= 1 )
            sprintf(p_ErrorMessage, " Error : %s Keyword not found",
                    p_LabelName);
        return XXX_E_ODL_NOT_FOUND;
    }

    if ((p_kwv = OdlGetKwdValue(p_kw)) == NULL)
    {
        strcpy(p_ErrorMessage, ODLErrorMessage);
        if ((long)strlen(p_ErrorMessage) <= 1 )
            sprintf(p_ErrorMessage, " Error : %s Keyword not found",
                    p_LabelName);
        return XXX_E_ODL_NOT_FOUND;
    }
    if ((p_keyword= malloc(strlen(p_kwv)+1)) == NULL)
    {
        strcpy(p_ErrorMessage,"Malloc error");
        return (ERROR);
    }
    strcpy(p_keyword,p_kwv);

    if (OdlGetKwdValueType (p_kw) == ODL_SET)
    {
        p_word = strtok(p_keyword,"(,) \"\n");
        while(p_word != NULL)
        {
            if (strlen(p_word))
            {
                switch (ValueType)
                {
                    case XXX_Long :
                        MemorySize -= sizeof(long);
                        if (MemorySize < 0 )
                        {
                            strcpy(p_ErrorMessage, " Input ODL value "
                                   "overflows allocated space ");
                            free(p_keyword);
                            return XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
                        }
                        if ((ret_code = xxx_ConvertString(p_MemoryAddr,
                                                          ValueType, p_word,
                                                          p_ErrorMessage))
                            != SUCCESS )
                        {
                            sprintf(p_ErrorMessage, "Error : Converting %s "
                                    "keyword's value %s", p_LabelName,
                                    p_word);
                            free(p_keyword);
                            return ret_code;
                        }
                        p_MemoryAddr = (long *)p_MemoryAddr + 1;
                        break;

                    case XXX_Int :
                        MemorySize -= sizeof(int);
                        if (MemorySize < 0 )
                        {
                            strcpy(p_ErrorMessage, " Input ODL value "
                                   "overflows allocated space ");
                            free(p_keyword);
                            return XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
                        }
                        if ((ret_code = xxx_ConvertString(p_MemoryAddr,
                                                          ValueType, p_word,
                                                          p_ErrorMessage))
                            != SUCCESS)
                        {
                            sprintf(p_ErrorMessage, "Error : Converting %s "
                                    "keyword's value %s",p_LabelName, p_word);
                            free(p_keyword);
                            return ret_code;
                        }
                        p_MemoryAddr = (int *)p_MemoryAddr + 1;
                        break;

                    case XXX_Float :
                        MemorySize -= sizeof(float);
                        if (MemorySize < 0 )
                        {
                            strcpy(p_ErrorMessage, " Input ODL value "
                                   "overflows allocated space ");
                            free(p_keyword);
                            return XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
                        }
                        if ((ret_code = xxx_ConvertString(p_MemoryAddr,
                                                          ValueType, p_word,
                                                          p_ErrorMessage))
                            != SUCCESS)
                        {
                            sprintf(p_ErrorMessage, "Error : Converting %s "
                                    "keyword's value %s", p_LabelName,
                                    p_word);
                            free(p_keyword);
                            return ret_code;
                        }
                        p_MemoryAddr = (float *)p_MemoryAddr + 1;
                        break;

                    case XXX_Double :
                    case XXX_Sci_Not :
                        MemorySize -= sizeof(double);
                        if (MemorySize < 0 )
                        {
                            strcpy(p_ErrorMessage, " Input ODL value "
                                   "overflows allocated space ");
                            free(p_keyword);
                            return XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
                        }
                        if ((ret_code = xxx_ConvertString(p_MemoryAddr,
                                                          ValueType, p_word,
                                                          p_ErrorMessage))
                            != SUCCESS)
                        {
                            sprintf(p_ErrorMessage,"Error : Converting %s "
                                    "keyword's value %s", p_LabelName,
                                    p_word);
                            free(p_keyword);
                            return ret_code;
                        }
                        p_MemoryAddr = (double *)p_MemoryAddr + 1;
                        break;

                    case XXX_ArrayOfString :
                        MemorySize -= sizeof(char *);
                        if (MemorySize < 0 )
                        {
                            strcpy(p_ErrorMessage, " Input ODL value "
                                   "overflows allocated space ");
                            free(p_keyword);
                            return  XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
                        }
                        if ((ret_code = xxx_ConvertString(p_MemoryAddr,
                                                          ValueType, p_word,
                                                          p_ErrorMessage))
                            != SUCCESS)
                        {
                            for( i=0;i<*p_Count;i++)
                                free((char *)p_MemoryAddr + i);
                            sprintf(p_ErrorMessage, "Error : Converting %s "
                                    "keyword's value %s", p_LabelName,
                                    p_word);
                            free(p_keyword);
                            return ret_code;
                        }
                        p_MemoryAddr = (char *)p_MemoryAddr + sizeof(char *);
                        break;

                    default:
                        strncpy( p_MemoryAddr, p_word, MemorySize);
                        sprintf(p_ErrorMessage, "Type XXX_String is not "
                                "valid for arrays of strings");
                        free(p_keyword);
                        return ERROR;
                }
                *p_Count += 1;
                p_word = strtok(NULL,",() \"\n");
            }
        }
    }

    else if (OdlGetKwdValueType (p_kw) == ODL_SEQUENCE)
    {
        p_word = strtok(p_keyword,"{,} \"\n");
        while(p_word != NULL)
        {
            if (strlen(p_word))
            {
                switch (ValueType)
                {
                    case XXX_ArrayOfString :
                        MemorySize -= sizeof(char *);
                        if (MemorySize < 0 )
                        {
                            strcpy(p_ErrorMessage, " Input ODL value "
                                   "overflows allocated space ");
                            free(p_keyword);
                            return  XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
                        }
                        if ((ret_code = xxx_ConvertString(p_MemoryAddr,
                                                          ValueType, p_word,
                                                          p_ErrorMessage))
                            != SUCCESS)
                        {
                            for( i=0;i<*p_Count;i++)
                                free((char *)p_MemoryAddr + i);
                            sprintf(p_ErrorMessage, "Error : Converting %s "
                                    "keyword's value %s", p_LabelName,
                                    p_word);
                            free(p_keyword);
                            return ret_code;
                        }
                        p_MemoryAddr = (char *)p_MemoryAddr + sizeof(char *);
                        break;

                    default:
                        strncpy( p_MemoryAddr, p_word, MemorySize);
                        sprintf(p_ErrorMessage, "An ODL sequence is of "
                                "type XXX_ArrayOfString only");
                        free(p_keyword);
                        return ERROR;
                }
                *p_Count += 1;
                p_word = strtok(NULL,",{} \"\n");
            }
        }
    }
    else
    {
        if (ValueType == XXX_String)
        {
            if (p_kwv[0] == '\"')
            {
                strcpy(string,  &p_kwv[1]);
                string[strlen(string)-1] = '\0';
            }
            else /* date */
                strcpy(string, p_kwv);

            if ((int)(MemorySize - (strlen(string) + 1)) < 0 )
            {
                strcpy(p_ErrorMessage,
                       " Input ODL value overflows allocated space ");
                strncpy(p_MemoryAddr, string, MemorySize);
                ((char *)p_MemoryAddr)[MemorySize-1] = '\0';
                free(p_keyword);
                return XXX_E_NOT_ENOUGH_MEMORY_SUPPLIED;
            }
            else
                strcpy(p_MemoryAddr, string);
        }
        else
        {
            switch(ValueType)
            {
                case XXX_Long:
                    MemorySize -= sizeof(long);
                    break;

                case XXX_Int:
                    MemorySize -= sizeof(int);
                    break;

                case XXX_Float:
                    MemorySize -= sizeof(float);
                    break;

                case XXX_Double:
                case XXX_Sci_Not:
                    MemorySize -= sizeof(double);
                    break;

                default:
                    sprintf(p_ErrorMessage,"Invalid Type specified");
                    free(p_keyword);
                    return ERROR;
            }

            if (MemorySize < 0 )
            {
                strcpy(p_ErrorMessage,
                       "Input ODL value overflows allocated space ");
                free(p_keyword);
                return ERROR;
            }

            if ((ret_code = xxx_ConvertString(p_MemoryAddr, ValueType, p_kwv,
                                              p_ErrorMessage)) != SUCCESS)
            {
                sprintf(p_ErrorMessage,
                        "Error : Converting %s keyword's value %s",
                        p_LabelName, p_kwv);
                free(p_keyword);
                return ret_code;
            }
        }
        *p_Count += 1;
    }
    free(p_keyword);
    return SUCCESS;
}
