/* 
 * File:   main.c
 * Author: taran
 *
 * Created on July 29, 2018, 4:00 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sunmodel.h"
 
#define RAD2DEG(X) (180.0 * (X) / M_PI)
#define DEG2RAD(X) (M_PI * (X) / 180.0)

/*
 * 
 */
int main(int argc, char** argv) {
    double lat, lon, azm, elv;
    time_t t;
    
    sunmodel_test();
    
    t = time(NULL);
    
    lat = DEG2RAD(64.133);
    lon = DEG2RAD(-21.9);
    
    sunmodel_ae(t, lat, lon, &azm, &elv, 'y');

    printf("azm: %8.2f\n", floor(RAD2DEG(azm) * 100.0 + 0.5) * 0.01);
    printf("elv: %8.2f\n", floor(RAD2DEG(elv) * 100.0 + 0.5) * 0.01);
    
    return (EXIT_SUCCESS);
}

