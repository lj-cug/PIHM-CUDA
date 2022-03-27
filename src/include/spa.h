/////////////////////////////////////////////
//          HEADER FILE for SPA.C          //
//                                         //
//      Solar Position Algorithm (SPA)     //
//                   for                   //
//        Solar Radiation Application      //
//                                         //
//               May 12, 2003              //
//                                         //
//   Filename: SPA.H                       //
//                                         //
//   Afshin Michael Andreas                //
//   afshin_andreas@nrel.gov (303)384-6383 //
//                                         //
//   Measurement & Instrumentation Team    //
//   Solar Radiation Research Laboratory   //
//   National Renewable Energy Laboratory  //
//   1617 Cole Blvd, Golden, CO 80401      //
/////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//                                                                    //
// Usage:                                                             //
//                                                                    //
//   1) In calling program, include this header file,                 //
//      by adding this line to the top of file:                       //
//           #include "spa.h"                                         //
//                                                                    //
//   2) In calling program, declare the SPA structure:                //
//           spa_data spa;                                            //
//                                                                    //
//   3) Enter the required input values into SPA structure            //
//      (input values listed in comments below)                       //
//                                                                    //
//   4) Call the SPA calculate function and pass the SPA structure    //
//      (prototype is declared at the end of this header file):       //
//           spa_calculate(&spa);                                     //
//                                                                    //
//   Selected output values (listed in comments below) will be        //
//   computed and returned in the passed SPA structure.  Output       //
//   will based on function code selected from enumeration below.     //
//                                                                    //
//   Note: A non-zero return code from spa_calculate() indicates that //
//         one of the input values did not pass simple bounds tests.  //
//         The valid input ranges and return error codes are also     //
//         listed below.                                              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef __solar_position_algorithm_header
#define __solar_position_algorithm_header


//enumeration for function codes to select desired final outputs from SPA
enum
{
    SPA_ZA,                     //calculate zenith and azimuth
    SPA_ZA_INC,                 //calculate zenith, azimuth, and incidence
    SPA_ZA_RTS,                 //calculate zenith, azimuth, and sun rise/transit/set values
    SPA_ALL,                    //calculate all SPA output values
};

typedef struct
{
    //----------------------INPUT VALUES------------------------

    int             year;       // 4-digit year,    valid range: -2000 to 6000, error code: 1
    int             month;      // 2-digit month,         valid range: 1 to 12, error code: 2
    int             day;        // 2-digit day,           valid range: 1 to 31, error code: 3
    int             hour;       // Observer local hour,   valid range: 0 to 24, error code: 4
    int             minute;     // Observer local minute, valid range: 0 to 59, error code: 5
    int             second;     // Observer local second, valid range: 0 to 59, error code: 6

    realtype          delta_ut1;  // Fractional second difference between UTC and UT which is used
    // to adjust UTC for earth's irregular rotation rate and is derived
    // from observation only and is reported in this bulletin:
    // http://maia.usno.navy.mil/ser7/ser7.dat,
    // where delta_ut1 = DUT1
    // valid range: -1 to 1 second (exclusive), error code 17

    realtype          delta_t;    // Difference between earth rotation time and terrestrial time
    // It is derived from observation only and is reported in this
    // bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
    // where delta_t = 32.184 + (TAI-UTC) - DUT1
    // valid range: -8000 to 8000 seconds, error code: 7

    realtype          timezone;   // Observer time zone (negative west of Greenwich)
    // valid range: -18   to   18 hours,   error code: 8

    realtype          longitude;  // Observer longitude (negative west of Greenwich)
    // valid range: -180  to  180 degrees, error code: 9

    realtype          latitude;   // Observer latitude (negative south of equator)
    // valid range: -90   to   90 degrees, error code: 10

    realtype          elevation;  // Observer elevation [meters]
    // valid range: -6500000 or higher meters,    error code: 11

    realtype          pressure;   // Annual average local pressure [millibars]
    // valid range:    0 to 5000 millibars,       error code: 12

    realtype          temperature;    // Annual average local temperature [degrees Celsius]
    // valid range: -273 to 6000 degrees Celsius, error code; 13

    realtype          slope;      // Surface slope (measured from the horizontal plane)
    // valid range: -360 to 360 degrees, error code: 14

    realtype          azm_rotation;   // Surface azimuth rotation (measured from south to projection of
    //     surface normal on horizontal plane, negative west)
    // valid range: -360 to 360 degrees, error code: 15

    realtype          atmos_refract;  // Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    // valid range: -5   to   5 degrees, error code: 16

    int             function;   // Switch to choose functions for desired output (from enumeration)

    //-----------------Intermediate OUTPUT VALUES--------------------

    realtype          jd;         //Julian day
    realtype          jc;         //Julian century

    realtype          jde;        //Julian ephemeris day
    realtype          jce;        //Julian ephemeris century
    realtype          jme;        //Julian ephemeris millennium

    realtype          l;          //earth heliocentric longitude [degrees]
    realtype          b;          //earth heliocentric latitude [degrees]
    realtype          r;          //earth radius vector [Astronomical Units, AU]

    realtype          theta;      //geocentric longitude [degrees]
    realtype          beta;       //geocentric latitude [degrees]

    realtype          x0;         //mean elongation (moon-sun) [degrees]
    realtype          x1;         //mean anomaly (sun) [degrees]
    realtype          x2;         //mean anomaly (moon) [degrees]
    realtype          x3;         //argument latitude (moon) [degrees]
    realtype          x4;         //ascending longitude (moon) [degrees]

    realtype          del_psi;    //nutation longitude [degrees]
    realtype          del_epsilon; //nutation obliquity [degrees]
    realtype          epsilon0;   //ecliptic mean obliquity [arc seconds]
    realtype          epsilon;    //ecliptic true obliquity  [degrees]

    realtype          del_tau;    //aberration correction [degrees]
    realtype          lamda;      //apparent sun longitude [degrees]
    realtype          nu0;        //Greenwich mean sidereal time [degrees]
    realtype          nu;         //Greenwich sidereal time [degrees]

    realtype          alpha;      //geocentric sun right ascension [degrees]
    realtype          delta;      //geocentric sun declination [degrees]

    realtype          h;          //observer hour angle [degrees]
    realtype          xi;         //sun equatorial horizontal parallax [degrees]
    realtype          del_alpha;  //sun right ascension parallax [degrees]
    realtype          delta_prime;    //topocentric sun declination [degrees]
    realtype          alpha_prime;    //topocentric sun right ascension [degrees]
    realtype          h_prime;    //topocentric local hour angle [degrees]

    realtype          e0;         //topocentric elevation angle (uncorrected) [degrees]
    realtype          del_e;      //atmospheric refraction correction [degrees]
    realtype          e;          //topocentric elevation angle (corrected) [degrees]

    realtype          eot;        //equation of time [minutes]
    realtype          srha;       //sunrise hour angle [degrees]
    realtype          ssha;       //sunset hour angle [degrees]
    realtype          sta;        //sun transit altitude [degrees]

    //---------------------Final OUTPUT VALUES------------------------

    realtype          zenith;     //topocentric zenith angle [degrees]
    realtype          azimuth180; //topocentric azimuth angle (westward from south) [-180 to 180 degrees]
    realtype          azimuth;    //topocentric azimuth angle (eastward from north) [   0 to 360 degrees]
    realtype          incidence;  //surface incidence angle [degrees]

    realtype          suntransit; //local sun transit time (or solar noon) [fractional hour]
    realtype          sunrise;    //local sunrise time (+/- 30 seconds) [fractional hour]
    realtype          sunset;     //local sunset time (+/- 30 seconds) [fractional hour]

} spa_data;

//------ Utility functions for other applications (such as NREL's SAMPA) -----
realtype          deg2rad(realtype degrees);
realtype          rad2deg(realtype radians);
realtype          limit_degrees(realtype degrees);
realtype          third_order_polynomial(realtype a, realtype b, realtype c,
    realtype d, realtype x);
realtype          geocentric_right_ascension(realtype lamda, realtype epsilon,
    realtype beta);
realtype          geocentric_declination(realtype beta, realtype epsilon,
    realtype lamda);
realtype          observer_hour_angle(realtype nu, realtype longitude,
    realtype alpha_deg);
void            right_ascension_parallax_and_topocentric_dec(realtype latitude,
    realtype elevation, realtype xi, realtype h, realtype delta, realtype *delta_alpha,
    realtype *delta_prime);
realtype          topocentric_right_ascension(realtype alpha_deg,
    realtype delta_alpha);
realtype          topocentric_local_hour_angle(realtype h, realtype delta_alpha);
realtype          topocentric_elevation_angle(realtype latitude,
    realtype delta_prime, realtype h_prime);
realtype          atmospheric_refraction_correction(realtype pressure,
    realtype temperature, realtype atmos_refract, realtype e0);
realtype          topocentric_elevation_angle_corrected(realtype e0,
    realtype delta_e);
realtype          topocentric_zenith_angle(realtype e);
realtype          topocentric_azimuth_angle_neg180_180(realtype h_prime,
    realtype latitude, realtype delta_prime);
realtype          topocentric_azimuth_angle_zero_360(realtype azimuth180);


//Calculate SPA output values (in structure) based on input values passed in structure
int             spa_calculate(spa_data * spa);

#endif
