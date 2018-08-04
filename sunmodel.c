
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "sunmodel.h"

#define RAD2DEG(X) (180.0 * X / M_PI)
#define DEG2RAD(X) (M_PI * X / 180.0)

static double sunmodel_time_julian_cent(double jd) {
    return (jd - 2451545.0) / 36525.0;
}

static double sunmodel_jd_from_julian_cent(double t) {
    return t * 36525.0 + 2451545.0;
}

static int sunmodel_is_leap_year(int yr) {
    return (yr % 4 == 0 && yr % 100 != 0) || (yr % 400 == 0);
}

static double sunmodel_doy_from_jd(double jd) {
    double z, f, alpha;
    double a, b, c, d, e, k;
    double day, month, year;
    
    z = floor(jd + 0.5);
    f = (jd + 0.5) - z;
  
    if (z < 2299161.0) {
        a = z;
    } else {
        alpha = floor((z - 1867216.25) / 36524.25);
        a = z + 1.0 + alpha - floor(alpha * 0.25);
    }
    
    b = a + 1524.0;
    c = floor((b - 122.1) / 365.25);
    d = floor(365.25 * c);
    e = floor((b - d) / 30.6001);
  
    day = b - d - floor(30.6001 * e) + f;
    month = (e < 14.0) ? e - 1.0 : e - 13.0;
    year = (month > 2.0) ? c - 4716.0 : c - 4715.0;

    k = sunmodel_is_leap_year(year) ? 1.0 : 2.0;
    
    return floor((275.0 * month) / 9.0) 
            - k * floor((month + 9.0) / 12.0) + day 
            - 30.0;
}

static double sunmodel_geom_mean_long_sun(double t) {
    double lon;
    
    lon = 280.46646 + t * (36000.76983 + t * 0.0003032);
    
    while(lon > 360.0) {
      lon -= 360.0;
    }
    
    while(lon < 0.0) {
      lon += 360.0;
    }
    
    return lon; // in degrees
}

static double sunmodel_geom_mean_anomaly_sun(double t) {
    double m;
    m = 357.52911 + t * (35999.05029 - 0.0001537 * t);
    return m; // in degrees
}

static double sunmodel_eccentricity_earth_orbit(double t) {
    double e;
    e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t);
    return e; // unitless
}

static double sunmodel_sun_eq_of_center(double t) {
    double m, mrad, sinm, sin2m, sin3m, c;
    
    m = sunmodel_geom_mean_anomaly_sun(t);
    mrad = DEG2RAD(m);
    sinm = sin(mrad);
    sin2m = sin(mrad + mrad);
    sin3m = sin(mrad + mrad + mrad);
    c = sinm * (1.914602 - t * (0.004817 + 0.000014 * t))
            + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289;
    
    return c; // in degrees
}

static double sunmodel_sun_true_long(double t) {
    double l0, c, o;
    l0 = sunmodel_geom_mean_long_sun(t);
    c = sunmodel_sun_eq_of_center(t);
    o = l0 + c;
    return o; // in degrees
}

static double sunmodel_sun_true_anomaly(double t) {
    double m, c, v;
    m = sunmodel_geom_mean_anomaly_sun(t);
    c = sunmodel_sun_eq_of_center(t);
    v = m + c;
    return v; // in degrees
}

static double sunmodel_sun_rad_vector(double t) {
    double v, e, r;
    v = sunmodel_sun_true_anomaly(t);
    e = sunmodel_eccentricity_earth_orbit(t);
    r = (1.000001018 * (1.0 - e * e)) / (1.0 + e * cos(DEG2RAD(v)));
    return r; // in AUs
}

static double sunmodel_sun_apparent_long(double t) {
    double o, omega, lambda;
    o = sunmodel_sun_true_long(t);
    omega = 125.04 - 1934.136 * t;
    lambda = o - 0.00569 - 0.00478 * sin(DEG2RAD(omega));
    return lambda; // in degrees
}

static double sunmodel_mean_obliquity_of_ecliptic(double t) {
    double seconds, e0;
    seconds = 21.448 - t * (46.8150 + t * (0.00059 - t * 0.001813));
    e0 = 23.0 + (26.0 + seconds / 60.0) / 60.0;
    return e0; // in degrees
}

static double sunmodel_obliquity_correction(double t) {
    double e0, omega, e;
    e0 = sunmodel_mean_obliquity_of_ecliptic(t);
    omega = 125.04 - 1934.136 * t;
    e = e0 + 0.00256 * cos(DEG2RAD(omega));
    return e; // in degrees
}

//function calcSunRtAscension(t)
//{
//  var e = calcObliquityCorrection(t);
//  var lambda = calcSunApparentLong(t);
//  var tananum = (cos(DEG2RAD(e)) * sin(DEG2RAD(lambda)));
//  var tanadenom = (cos(DEG2RAD(lambda)));
//  var alpha = RAD2DEG(atan2(tananum, tanadenom));
//  return alpha;		// in degrees
//}
//
static double sunmodel_sun_declination(double t) {
    double e, lambda, sint, theta;
    e = sunmodel_obliquity_correction(t);
    lambda = sunmodel_sun_apparent_long(t);

    sint = sin(DEG2RAD(e)) * sin(DEG2RAD(lambda));
    theta = RAD2DEG(asin(sint));
    return theta; // in degrees
}

static double sunmodel_equation_of_time(double t) {
    double epsilon, l0, e, m, y, sin2l0, sinm, cos2l0, sin4l0, sin2m, etime;
    
    epsilon = sunmodel_obliquity_correction(t);
    l0 = sunmodel_geom_mean_long_sun(t);
    e = sunmodel_eccentricity_earth_orbit(t);
    m = sunmodel_geom_mean_anomaly_sun(t);

    y = tan(DEG2RAD(epsilon) * 0.5);
    y *= y;

    sin2l0 = sin(2.0 * DEG2RAD(l0));
    sinm   = sin(DEG2RAD(m));
    cos2l0 = cos(2.0 * DEG2RAD(l0));
    sin4l0 = sin(4.0 * DEG2RAD(l0));
    sin2m  = sin(2.0 * DEG2RAD(m));

    etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 
            - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;
    return RAD2DEG(etime) * 4.0; // in minutes of time
}

static double sunmodel_get_jd(struct tm *gmt) {
    int month, year;
    double a, b, jd;
    
    month = gmt->tm_mon + 1;
    year = gmt->tm_year + 1900;
  
    if (month <= 2) {
        year -= 1;
        month += 12;
    }
    
    a = floor(year * 0.01);
    b = 2.0 - a + floor(a * 0.25);
    jd = floor(365.25 * (year + 4716.0)) + floor(30.6001 * (month + 1.0)) 
            + gmt->tm_mday + b - 1524.5;
    return jd;
}

static double sunmodel_az_el(double T, double localtime, double latitude, double longitude, double zone) {
    double azimuth, eq_time, theta, solar_time_fix, earth_rad_vec;
    double true_solar_time, hour_angle, csz, zenith, az_denom;
    double exoatm_elevation, refraction_correction, te, solar_zen;
    
    eq_time = sunmodel_equation_of_time(T);
    theta  = sunmodel_sun_declination(T);
    
    printf("eq_time: %8.2f\n", floor(eq_time * 100.0 + 0.5) * 0.01);
    printf("theta  : %8.2f\n", floor(theta * 100.0 + 0.5) * 0.01);
    

    solar_time_fix = eq_time + 4.0 * longitude - 60.0 * zone;
    earth_rad_vec = sunmodel_sun_rad_vector(T);
    true_solar_time = localtime + solar_time_fix;
            
    while (true_solar_time > 1440.0) {
        true_solar_time -= 1440.0;
    }
    
    hour_angle = true_solar_time / 4.0 - 180.0;
    if (hour_angle < -180.0)  {
        hour_angle += 360.0;
    }
    
    hour_angle = DEG2RAD(hour_angle);
    csz = sin(DEG2RAD(latitude)) * sin(DEG2RAD(theta)) + 
            cos(DEG2RAD(latitude)) * cos(DEG2RAD(theta)) * cos(hour_angle);
    
    if (csz > 1.0) {
        csz = 1.0;
    }
    else {
        if (csz < -1.0) { 
            csz = -1.0;
        }
    }
    
    zenith = RAD2DEG(acos(csz));
    az_denom = cos(DEG2RAD(latitude)) * sin(DEG2RAD(zenith));
    
    if (fabs(az_denom) > DBL_EPSILON) {
        azimuth = ((sin(DEG2RAD(latitude)) * cos(DEG2RAD(zenith))) - sin(DEG2RAD(theta))) / az_denom;
    
        if (fabs(azimuth) > 1.0) {
            if (azimuth < 0.0) {
                azimuth = -1.0;
            } 
            else {
                azimuth = 1.0;
            }
        }
        
        azimuth = 180.0 - RAD2DEG(acos(azimuth));
        if (hour_angle > 0.0) {
            azimuth = -azimuth;
        }
    } 
    else {
        if (latitude > 0.0) {
            azimuth = 180.0;
        } else { 
            azimuth = 0.0;
        }
    }
    
    if (azimuth < 0.0) {
        azimuth += 360.0;
    }
    

// Atmospheric Refraction correction
    
    exoatm_elevation = 90.0 - zenith;

    if (exoatm_elevation > 85.0) {
          refraction_correction = 0.0;
    } 
    else {
        te = tan(DEG2RAD(exoatm_elevation));
        if (exoatm_elevation > 5.0) {
            refraction_correction = 
                    58.1 / te - 0.07 / (te * te * te) + 
                    0.000086 / (te * te * te * te * te);
        } 
        else {
            if (exoatm_elevation > -0.575) {
                refraction_correction = 
                        1735.0 + exoatm_elevation * (-518.2 + exoatm_elevation * 
                        (103.4 + exoatm_elevation * (-12.79 + exoatm_elevation * 0.711)));
            }
            else {
                refraction_correction = -20.774 / te;
            }
        }
        refraction_correction = refraction_correction / 3600.0;
    }

    solar_zen = zenith - refraction_correction;

    if (solar_zen > 108.0) {
        puts("A Night at the Roxbury");
    }
    
    printf("azm: %8.2f\n", floor(azimuth * 100.0 + 0.5) * 0.01);
    printf("elv: %8.2f\n", floor((90.0 - solar_zen) * 100.0 + 0.5) * 0.01);
    
    return azimuth;
}

void sunmodel_test(void) {
    time_t t;
    struct tm *gmt;
    double jday, tl, tz, total, T, lat, lon;

    t = time(NULL);
    gmt = gmtime(&t);
    puts(ctime(&t));
    if (gmt == NULL) {
        exit(-1);
    }
    
    jday = sunmodel_get_jd(gmt);
    tl = (gmt->tm_hour / 24.0 + gmt->tm_min / 1440.0 + gmt->tm_sec / 86400.0) * 1440.0;
    tz = 0.0;
    total = jday + tl / 1440.0 - tz / 24.0;
    T = sunmodel_time_julian_cent(total);
    lat = 58.0;
    lon = 63.0;
    sunmodel_az_el(T, tl, lat, lon, tz);
}


