
#include <stdio.h>
#include <lablib3.h>                 /* prototypes for the ODL functions */
#include <xxx_ODL.h>                 /* prototype for the xxx ODL functions */

extern char ODLErrorMessage[];       /* External Variables */

/*****************************************************************************/
/**
 * @file xxx_CloseODL.c
 * @brief Close the ODL tree.
 * @ingroup iasLib
 */
/**
 * @brief Close the ODL tree.
 *
 * @return SUCCESS: ODL tree freed
 * @return ERROR: Error freeing ODL tree
 */
/*****************************************************************************/

int xxx_CloseODL
(
    OBJDESC *p_lp,       //!<[in] ODL tree
    char *p_ErrorMessage //!<[out] Error message
)
{
    p_lp = OdlFreeTree(p_lp);
    if (p_lp != NULL)
    {
        strcpy(p_ErrorMessage, "Can't Close ODL file:  ");
        strcat(p_ErrorMessage, ODLErrorMessage);
        return ERROR;
    }
    return SUCCESS;
}
