
#include <stdio.h>

#include <xxx_Errno.h>
#include <xxx_Types.h>

/****************************************************************************/
/**
 * @file xxx_Errno.c
 * @brief Definition of error return values and strings.
 * @ingroup iasLib
 */
/**
 * @brief Definition of error return values and strings.
 *
 */
/*****************************************************************************/

char *xxx_errlist[] = {
	/* XXX_E_NO_MEMORY */				"no memory available",
	/* XXX_E_NULL_INPUT */				"NULL input to function",
	/* XXX_E_NO_MORE_OBJECTS */			"no more objects",
	/* XXX_E_OBJECT_NOT_FOUND */			"object not found",
	/* XXX_E_OUT_OF_RANGE */			"input out of range",
	/* XXX_E_INVALID_DATA */			"invalid data",
	/* XXX_E_TIMEOUT */				"timeout",
	/* XXX_E_FAILURE */				"unknown failure",
	/* XXX_E_INSUFFICIENT_INPUT */			"insufficient input",
	/* XXX_E_FILE_NOT_FOUND */			"File not found",
	/* XXX_E_ODL_SYNTAX */			       	"ODL syntax",
	NULL
};

long xxx_nerr = XXX_ARRAY_SIZE(xxx_errlist) - 1;

long xxx_errno;
