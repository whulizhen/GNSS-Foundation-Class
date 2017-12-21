//
//  sp3Tucl.cpp
//  GFC
//
//  Created by lizhen on 22/09/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include <stdio.h>

#include "GSpacecraftModel.hpp"
#include "GSpacecraft.hpp"

#include "GJPLEPH.h"

#include "GEarthOrientationParameter.hpp"

#include "GSensor.h"

#include "GVector.hpp"

using namespace gfc;


/* the parameters for the argv
 argv[1] sysName
 argv[2] prn
 argv[3] sp3file0
 argv[4] sp3file1
 argv[5] sp3file2
 argv[6] gps_year
 argv[7] gps_doy
 
 */
int main(int argc, char** argv)
{
    
    //argv[1] = "2014";
    //argv[2] = "90";
    
    int year = atoi(argv[1]);
    int doy = atoi(argv[2]);
    
    DOYTime dt(year, doy,0.0,"tsGPS");
    
    GTime t;
    t.SetFromDoyTime(dt);
    
    NavTime nt = GTime::GTime2NavTime(t);
    
    printf("%d %1d\n",nt.m_week, nt.getDOW());
    
    return 0;
}
