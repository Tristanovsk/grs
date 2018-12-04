#include <math.h>
#include "gxx_coord_conv.h"

/****************************************************************************/
/**
 * @file gxx_geod2cart.c
 * @brief Convert geodetic coordinate (lat, lon, height) into Cartesian
 * coordinates (x, y, z).
 * @ingroup gpsGeo
 */
/**
 * @brief Convert geodetic coordinate (lat, lon, height) into Cartesian
 * coordinates (x, y, z).
 *
 * Input lat & lon are in radians, height, semi-major axis, and 
 * output Cartesian position vector are in meters
 * flattening is a unitless number.
 *
 * ALGORITHM REFERENCES:
 * Geodesy: The Concept, by Vanicek and Karkivsky, 1986;
 * or any other geodesy textbook.
 */
/****************************************************************************/


void gxx_geod2cart
(
    double latitude,    //!<[in] Lat of geodetic coordinates in radians 
    double longitude,   //!<[in] Long of geodetic coordinates in radians
    double height,      //!<[in] Height (elevation) of geodetic coord in meters
    double semimajor,   //!<[in] Reference ellipsoid semi-major axis in meters 
    double flattening,  //!<[in] Flattening of the ellipsoid 
                        // (semimajor-semiminor)/semimajor    
    VECTOR *cart        //!<[out] Cartesian vector for the coord    
)
{
    double prime_vertical_radius,  /* radius of prime vertical     */
        ecc2,                      /* square of eccentricity       */
        coslat, sinlat;            /* cosine and sine of lat       */

    /*
     * Calculate the radius of prime vertical of the ellipsoid
     */ 
    ecc2 = flattening * (2.0 - flattening);
    coslat = cos(latitude);
    sinlat = sin(latitude);

    /*
      Since sinlat is less than or equal to 1 and ecc2 is less than 1 no
      test is performed for negative square root or division by zero
    */
    prime_vertical_radius = semimajor / sqrt(1.0 - ecc2 * sinlat * sinlat);

    /*
     * Calculate the cartesian coordinates
     */
    cart->x = (prime_vertical_radius + height)*coslat*cos(longitude);
    cart->y = (prime_vertical_radius + height)*coslat*sin(longitude);
    cart->z = (prime_vertical_radius*(1.0 - flattening)*(1.0 - flattening)
               + height)*sinlat;
}
