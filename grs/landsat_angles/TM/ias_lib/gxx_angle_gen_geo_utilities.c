/****************************************************************************/
/**
 * @file gxx_angle_gen_geo_utilities.c
 * @brief A set of utilities for performing geometric transformations
 * @ingroup gpsAngleCoef
 */
/****************************************************************************/

#include "xxx_fulljd.h"
#include "gxx_coord_conv.h"
#include "gxx_proj.h"
#include "gxx_angle_gen_geo_utilities.h"
#include "gxx_angle_gen_leap_seconds.h"

/******************************************************************************/
/**
 * @brief Transformation of coordinates through rotation around the x-axis and 
 * specified angle.
 */
/****************************************************************************/
static void  geo_rotate_around_z
(
    const VECTOR *r_old, //!<[in] coordinates (x, y, z) in the old system
    double angle,        //!<[in] angle to be rotated, in radian, 
                         // positive anticlockwise
    VECTOR *r_new        //!<[out] coordinates in the new system
)
{
    double ca, sa;       /* cosine and sine of angle */

    /* Compute the sine and cosine of the angle. */
    ca = cos(angle);
    sa = sin(angle);

    /* Calculate the new coords using the following:
     * new x coord = cosine * old x coord + sine * old y coord
     * new y coord = cosine * old y coord - sine * old x coord
     * new z coord = old z coord */
    r_new->x = ca * r_old->x + sa * r_old->y;
    r_new->y = ca * r_old->y - sa * r_old->x;
    r_new->z = r_old->z;
}

#ifndef IAS_NOVAS_NOT_AVAILABLE
#include "gxx_angle_gen_novas_wrapper.h"

/***************************************************************************/
/**
 * @brief Calculate the Greenwich apparent sidereal time (GAST) for given UT1
 * and Terrestrial (TT) times.
 *
 * @returns interger (SUCCESS or ERROR)
 */
/****************************************************************************/
static int geo_get_sidereal_time
(
    double jd_ut1,        //!<[in] UT1 Julian date of ephemeris time
    double jd_tt,         //!<[in] TT Julian date of ephemeris time
    double *gast          //!<[out] Greenwich apparent sidereal time, in rad
)
{
    int status;
    double gast_hrs;   /* Greenwich apparent sidereal time (GAST), in hours */
    double delta_t;    /* Seconds difference between TT and UT1 times */
    double pi = 4.0 * atan(1.0);
    double radians_per_degree = pi / 180.0;
    double hours2radians = 15.0 * radians_per_degree;
        /* The conversion factor to convert the gast from hours to radians.
           There are 15 degrees per hour. */

    /* Calculate the GAST at the ephemeris time and convert it to radians */
    delta_t = (jd_tt - jd_ut1) * SEC_PER_DAY;
    status = NOVAS_SIDEREAL_TIME(jd_ut1, 0.0, delta_t, NOVAS_APPARENT_GAST,
        NOVAS_EQUINOX_METHOD, NOVAS_FULL_ACCURACY, &gast_hrs);
    if (status != SUCCESS)
    {
        xxx_LogStatus("GEO_GET_SIDEREAL_TIME", __FILE__, __LINE__,
                      "NOVAS sidereal_time routine returned error");
        return ERROR;
    }
    *gast = hours2radians * gast_hrs;

    return SUCCESS;
}

/***************************************************************************/
/**
 * @brief Convert UTC time to Universal Time (UT1), Terrestrial Time (TT), 
 * and Barycentric Dynamical Time (TDB).
 *  TT  = UTC + Leap Seconds + TAI
 *  UT1 = UTC + (UT1_UTC)
 *  TT differs from TDB by peroidic variations
 *
 * NOTES:
 * Reference: United States Naval Observatory Circular 180.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/***************************************************************************/
static int geo_convert_utc2times
(
    double ut1_utc,     //!<[in] UT1-UTC, in seconds, due to variation of 
                        // Earth's spin rate
    const double ephem_time[3],//!<[in] Time year, day and seconds of day (UTC)
    double *ut1,        //!<[out] Univeral Time (UT1) Time
    double *tdb,        //!<[out] Barycentric Dynamical Time (TDB)
    double *tt          //!<[out] Terrestrial Time (TT)
)
{
    double j2000_seconds;/* Input time converted to J2000 seconds */
    double jd_tt;        /* Julian date (TT) for ephemeris time */
    double jd_tdb;       /* Julian date (TDB) for ephemeris time */
    double secdiff;      /* Difference between TDB and TT in seconds. */
    double jd_utc;       /* Julian date (UTC) at start of year */
    double dummy;

    /* Convert the spacecraft ephemeris time to seconds since J2000 */
    if (math_convert_year_doy_sod_to_j2000_seconds(ephem_time,
        &j2000_seconds) != SUCCESS)
    {
        xxx_LogStatus("GEO_CONVERT_UTC2TIMES", __FILE__, __LINE__,
                      "Failed converting time to J2000 seconds");
        return ERROR;
    }

    /* Add the Terrestrial Time at noon of the J2000 epoch to the seconds
       since J2000 to get the total Terrestrial Time */
    jd_tt = EPOCH_2000 + j2000_seconds / SEC_PER_DAY;

    /* Use NOVAS routine to calculate small delta between Terrestrial and
       Barycentric times.  This difference is a small periodic difference that
       must be solved for.  Note that the TT time is used as an input even 
       though the NOVAS code considers it the TDB time, but they're close
       enough that the difference is valid. */
    jd_tdb = jd_tt;
    NOVAS_TDB2TT(jd_tdb, &dummy, &secdiff);

    /* Calculate the TDB time */
    jd_tdb = jd_tt + secdiff / SEC_PER_DAY;

    /* Calculate full julian day for current year.  Note this will contain leap
       seconds up to the start of the current year. */
    jd_utc = xxx_fulljd((int)ephem_time[0], 1, 1, 0);

    /* Note ephem_time[1] and ephem_time[2] is UTC and will contain leap
       seconds */
    *ut1 = jd_utc + (ephem_time[1] - 1.0)
        + (ut1_utc + ephem_time[2]) / SEC_PER_DAY;

    *tt  = jd_tt;
    *tdb = jd_tdb;

    return SUCCESS;
}

/******************************************************************************/
/**
 * @brief Transformation from the true-of-date system to mean-of-date system 
 * through nutation angles at a specified Barycentric Dynamical Time (TDB).
 */
/******************************************************************************/
static void geo_transform_nutation_tod2mod
(
    const VECTOR *r_old,    //!<[in] coordinates (x, y, z) in the true-of-date 
                            // system
    double jd_tdb,          //!<[in] Julian date (Barycentric) for conversion
    VECTOR *r_new           //!<[out] coordinates in the mean-of-date system
)
{
    double in_vec[3];       /* input vector as an array for NOVAS routine */
    double out_vec[3];      /* output vector as an array for NOVAS routine */

    /* Copy the input vector to the array */
    in_vec[0] = r_old->x;
    in_vec[1] = r_old->y;
    in_vec[2] = r_old->z;

    /* Do the nutation using NOVAS */
    NOVAS_NUTATION(jd_tdb, NOVAS_TRUE_TO_MEAN_DIRECTION, NOVAS_FULL_ACCURACY,
        in_vec, out_vec);

    /* Copy the output to the output vector */
    r_new->x = out_vec[0];
    r_new->y = out_vec[1];
    r_new->z = out_vec[2];
}

/******************************************************************************/
/**
 * @brief Transform from Earth's true spin axis (old system) to the mean 
 * pole (CIO) (new system).
 */
/******************************************************************************/
void gxx_angle_gen_geo_transform_polar_motion_true_pole_to_mean
(
    const VECTOR *r_old, //!<[in] coordinates (x, y, z) in the old system
    double xp,    //!<[in] true pole position in the mean pole coords system, 
                  // x-axis pointing along Greenwich meridian; in arc seconds
    double yp,    //!<[in] true pole position in the mean pole coords system, 
                  // y-axis pointing along west 90 degree meridian; in arc 
                  // seconds
    double jd_tdb,//!<[in] Julian date (Barycentric)
    VECTOR *r_new //!<[out] coordinates in the new system
)
{
    double in_vec[3];   /* Intermediate vector for NOVAS */
    double out_vec[3];  /* Intermediate vector for NOVAS */

    in_vec[0] = r_old->x;
    in_vec[1] = r_old->y;
    in_vec[2] = r_old->z;

    /* Calculate the wobble using NOVAS */
    NOVAS_WOBBLE(jd_tdb, NOVAS_TRUE_TO_MEAN_DIRECTION, xp, yp, in_vec, out_vec);

    r_new->x = out_vec[0];
    r_new->y = out_vec[1];
    r_new->z = out_vec[2];
}

/******************************************************************************/
/**
 * @brief Transformation from the mean-of-date system to true-of-date system 
 * through nutation angles at a specified Barycentric Dynamical Time (TDB).
 */
/******************************************************************************/
void gxx_angle_gen_geo_transform_nutation_mod2tod
(
    const VECTOR *r_old,//!<[in] coordinates (x, y, z) in the mean-of-date 
                        // system
    double jd_tdb,      //!<[in] Julian date (Barycentric) for conversion
    VECTOR *r_new       //!<[out] coordinates in the true-of-date equator and  
                        // equinox sys.
)
{
    double in_vec[3];       /* input vector as an array for NOVAS routine */
    double out_vec[3];      /* output vector as an array for NOVAS routine */

    /* Copy the input vector to the array */
    in_vec[0] = r_old->x;
    in_vec[1] = r_old->y;
    in_vec[2] = r_old->z;

    /* Do the nutation using NOVAS */
    NOVAS_NUTATION(jd_tdb, NOVAS_MEAN_TO_TRUE_DIRECTION, NOVAS_FULL_ACCURACY,
        in_vec, out_vec);

    /* Copy the output to the output vector */
    r_new->x = out_vec[0];
    r_new->y = out_vec[1];
    r_new->z = out_vec[2];
}

/*****************************************************************************/
/**
 * @brief Transformation from J2000.0 system to mean-of-date system through 
 * precession angles at a specified Barycentric Dynamical Time (TDB).
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
int gxx_angle_gen_geo_transform_precession_j2k2mod
(
    const VECTOR *r_old,//!<[in] coordinates (x, y, z) in the J2000.0 system
    double jd_tdb,      //!<[in] Julian date (Barycentric) for conversion
    VECTOR *r_new       //!<[out] coordinates in the mean-of-date equator and 
                        // equinox sys.
)
{
    double in_vec[3];       /* input vector as an array for NOVAS routine */
    double out_vec[3];      /* output vector as an array for NOVAS routine */
    int status;

    /* Copy the input vector to the array */
    in_vec[0] = r_old->x;
    in_vec[1] = r_old->y;
    in_vec[2] = r_old->z;

    /* Do the precession using NOVAS */
    status = NOVAS_PRECESSION(EPOCH_2000, in_vec, jd_tdb, out_vec);
    if (status != SUCCESS)
    {
        /* This error is very unlikely to happen since it can only happen
           if EPOCH_2000 isn't passed to the NOVAS precession routine as
           one of the parameters */
        xxx_LogStatus("GEO_TRANSFORM_PRECESSION_J2K2MOD", __FILE__, __LINE__,
                      "NOVAS precession routine returned error");
        return ERROR;
    }

    /* Copy the output to the output vector */
    r_new->x = out_vec[0];
    r_new->y = out_vec[1];
    r_new->z = out_vec[2];

    return SUCCESS;
}

/*******************************************************************************/
/**
 * @brief Transformation from mean-of-date system to J2000.0 system through
 * precession angles at a specified Barycentric Dynamical Time (TDB).
 *
 * @returns integer (SUCCESS or ERROR)
 */
/*******************************************************************************/
int gxx_angle_gen_geo_transform_precession_mod2j2k
(
    const VECTOR *r_old,    /* I: coordinates (x, y, z) in the mean-of-date
                                  system */
    double jd_tdb,          /* I: Julian date (Barycentric) for conversion */
    VECTOR *r_new           /* O: coordinates in the J2000.0 system */
)
{
    double in_vec[3];       /* input vector as an array for NOVAS routine */
    double out_vec[3];      /* output vector as an array for NOVAS routine */
    int status;

    /* Copy the input vector to the array */
    in_vec[0] = r_old->x;
    in_vec[1] = r_old->y;
    in_vec[2] = r_old->z;

    /* Do the precession using NOVAS */
    status = NOVAS_PRECESSION(jd_tdb, in_vec, EPOCH_2000, out_vec);
    if (status != SUCCESS)
    {
        /* This error is very unlikely to happen since it can only happen
           if IAS_EPOCH_2000 isn't passed to the NOVAS precession routine as
           one of the parameters */
        xxx_LogStatus("GEO_TRANSFORM_PRECESSION_MOD2J2K", __FILE__, __LINE__,
                      "NOVAS precession routine returned error");
        return ERROR;
    }

    /* Copy the output to the output vector */
    r_new->x = out_vec[0];
    r_new->y = out_vec[1];
    r_new->z = out_vec[2];

    return SUCCESS;
}

/***************************************************************************/
/**
 * @brief Transform Earth-centered inertial Cartesian coordinates 
 * (ECI/true-of-date) to Earth-centered, Earth-fixed Cartesian 
 * coordinates (ECEF) at the specified GMT (UTC) time.
 *
 * @returns integer (SUCCESS or ERROR)
 */
/******************************************************************************/
int gxx_angle_gen_geo_eci2ecef
(
    double xp, //!<[in] Earth's true pole offset from mean pole, in arc second
    double yp, //!<[in] Earth's true pole offset from mean pole, in arc second
    double ut1_utc, //!<[in] UT1-UTC, in seconds, due to variation of Earth's 
                    // spin rate
    const VECTOR *craft_pos, //!<[in] Satellite position in ECI
    const VECTOR *craft_vel, //!<[in] Satellite velocity in ECI
    const double ephem_time[3],  //!<[in] UTC Ephemeris time (year, doy and sod)
    VECTOR *fe_satpos, //!<[out] Satellite position in ECEF
    VECTOR *fe_satvel  //!<[out] Satellite velocity in ECEF
)
{
    double gast;    /* Greenwich apparent sidereal time, in rad */
    VECTOR eci_vec; /* vector in ECI system */
    VECTOR pre_vec; /* vector after precession transformation */
    VECTOR nut_vec; /* vector after nutation transformation */
    VECTOR mid_vec; /* vector in intermediate step */
    VECTOR ecr_vec; /* vector in ECEF system */
    double jd_tdb;  /* Barycentric (TDB) Julian date of ephem_time */
    double jd_tt;   /* Terrestrial Julian date of ephem_time */
    double jd_ut1;  /* UT1 Julian date of ephem_time */

    /* Convert the input time into the different time standards needed */
    if (geo_convert_utc2times(ut1_utc, ephem_time, &jd_ut1, &jd_tdb,
        &jd_tt) != SUCCESS)
    {
        xxx_LogStatus("GEO_ECI2ECEF", __FILE__, __LINE__,
                       "Unable to convert UTC time to other time standards");
        return ERROR;
    }

    /* Get the Greenwich apparent sidereal time */
    if (geo_get_sidereal_time(jd_ut1, jd_tt, &gast) != SUCCESS)
    {
        xxx_LogStatus("GEO_ECI2ECEF", __FILE__, __LINE__,
                       "Unable to get Greenwich sidereal time");
        return ERROR;
    }

    /* Convert the satellite position vector from ECI to ECEF */
    eci_vec.x = craft_pos->x; /* copy from the input vector */
    eci_vec.y = craft_pos->y;
    eci_vec.z = craft_pos->z;

    /* ACS ephemeris data is in the J2000.0 system, precession and nutation
       transformation should be added here before the rotation transformation
       around Earth's spin axis */
    if (gxx_angle_gen_geo_transform_precession_j2k2mod(&eci_vec, jd_tdb, 
        &pre_vec) != SUCCESS)
    {
        xxx_LogStatus("GEO_ECI2ECEF", __FILE__, __LINE__,
                       "Failed performing the precession tranformation");
        return ERROR;
    }

    gxx_angle_gen_geo_transform_nutation_mod2tod(&pre_vec, jd_tdb, &nut_vec);
    eci_vec.x = nut_vec.x;
    eci_vec.y = nut_vec.y;
    eci_vec.z = nut_vec.z;

    /* rotate around true pole for gast */
    geo_rotate_around_z(&eci_vec, gast, &mid_vec);
    gxx_angle_gen_geo_transform_polar_motion_true_pole_to_mean(&mid_vec, xp, 
         yp, jd_tdb, &ecr_vec);
    fe_satpos->x = ecr_vec.x;
    fe_satpos->y = ecr_vec.y;
    fe_satpos->z = ecr_vec.z;

    /* Convert the satellite velocity vector from ECI to ECEF, no magnitude 
       change */
    eci_vec.x = craft_vel->x; /* copy vector from input */
    eci_vec.y = craft_vel->y;
    eci_vec.z = craft_vel->z;

    /* ACS ephemeris data is in the J2000.0 system, precession and nutation
       transformation should be added here before the rotation transformation
       around Earth's spin axis */
    if (gxx_angle_gen_geo_transform_precession_j2k2mod(&eci_vec, jd_tdb, 
        &pre_vec) != SUCCESS)
    {
        xxx_LogStatus("GEO_ECI2ECEF", __FILE__, __LINE__,
                      "Failed performing the precession tranformation");
        return ERROR;
    }

    gxx_angle_gen_geo_transform_nutation_mod2tod(&pre_vec, jd_tdb, &nut_vec);
    eci_vec.x = nut_vec.x;
    eci_vec.y = nut_vec.y;
    eci_vec.z = nut_vec.z;

    /* rotation around true pole   */
    geo_rotate_around_z(&eci_vec, gast, &mid_vec);

    /* from true pole to mean pole */
    gxx_angle_gen_geo_transform_polar_motion_true_pole_to_mean(&mid_vec, xp, 
         yp, jd_tdb, &ecr_vec);
    fe_satvel->x = ecr_vec.x;
    fe_satvel->y = ecr_vec.y;
    fe_satvel->z = ecr_vec.z;

    return SUCCESS;
}

/***************************************************************************/
/**
 * @brief Transform true-of-date inertial cartesian coordinates (ECITOD) to 
 * inertial cartesian coordinates (ECI of epoch J2000) at the specified GMT
 * (UTC) time.
 *
 * @returns interger (SUCCESS or ERROR)
 */
/******************************************************************************/
int geo_transform_tod2j2k
(
    double ut1_utc,              //!<[in] UT1-UTC, in seconds, due to variation
                                 // of Earth's spin rate
    const VECTOR *ecitod_pos,    //!<[in] Satellite position in ECITOD
    const double ephem_time[3],  //!<[in] UTC Ephemeris time (year, doy and sod)
    VECTOR *ecij2k_pos           //!<[out] Satellite position in ECIJ2K
)
{
    double jd_tdb;          /* TDB Julian date of ephemeris time */
    double jd_tt;           /* TT Julian date of ephemeris time */
    double jd_ut1;          /* UT1 Julian date of ephemeris time */
    VECTOR ecimod_pos;      /* ECI mean of date position after precession */

    /* Convert the input time into the different time standards needed */
    if (geo_convert_utc2times(ut1_utc, ephem_time, &jd_ut1, &jd_tdb,
        &jd_tt) != SUCCESS)
    {
        xxx_LogStatus("GEO_TRANSFORM_TOD2J2K", __FILE__, __LINE__,
                      "Unable to convert UTC time to other time standards");
        return ERROR;
    }

    /* Convert the satellite position vector from ECITOD to ECIJ2K */

    /* Nutation transformation converts ECI true of date to ECI mean of date */
    geo_transform_nutation_tod2mod(ecitod_pos, jd_tdb, &ecimod_pos);

    /* Precession transformation converts ECI mean of date to ECI of epoch
       J2000 */
    if (gxx_angle_gen_geo_transform_precession_mod2j2k(&ecimod_pos, 
            jd_tdb, ecij2k_pos) != SUCCESS)
    {
        xxx_LogStatus("GEO_TRANSFORM_TOD2J2K", __FILE__, __LINE__, 
                      "Error returned from transform precession MOD2J2K");
        return ERROR;
    }

    return SUCCESS;
}
#endif /* NOVAS library available */

/******************************************************************************/
/**
 * @brief Calculates the geocentric vector corresponding to the input latitude, 
 * longitude, and height.
 */
/******************************************************************************/
void gxx_angle_gen_geo_to_ecef
(
    const gxx_angle_gen_metadata_TYPE *metadata, //!<[in] Projection information
    double latitude,                        //!<[in] Latitude in radians
    double longitude,                       //!<[in] Longitude in radians
    double height,                          //!<[in] Height in meters
    VECTOR *ecef_vector                 //!<[out] ECEF vector
)
{
    double wgs84_flattening; /* Flattening value calculated from axis */

    /* Calculate the WGS84 flattening */
    wgs84_flattening = (metadata->wgs84_major_axis
        - metadata->wgs84_minor_axis) / metadata->wgs84_major_axis;

    /* Use the GEO library routine to perform the computations */
    gxx_geod2cart(latitude, longitude, height,
        metadata->wgs84_major_axis, wgs84_flattening, ecef_vector);
}
