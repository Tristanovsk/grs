#include "gxx_geo_math.h"

/****************************************************************************/
/**
 * @file gxx_dot.c
 * @brief Multiply two vectors of type VECTOR (Dot product)
 * @ingroup gpsMath
 */
/**
 * @brief Multiply two vectors of type VECTOR (Dot product)
 *
 * @return Value of computed dot product
 */
/****************************************************************************/


double gxx_dot
(
    const VECTOR *vec1, //!<[in] Vector one to be multiplied     
    const VECTOR *vec2  //!<[in] Vector two to be multiplied     
)
{
    double dot_product = 0.0;
    

    dot_product = vec1->x * vec2->x + vec1->y * vec2->y + vec1->z * vec2->z;

    return dot_product;
}
