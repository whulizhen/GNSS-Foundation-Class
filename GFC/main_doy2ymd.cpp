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


int main(int argc, char** argv)
{
    
    //argv[1] = "2014";
    //argv[2] = "90";
    
    int year = atoi(argv[1]);
    int doy = atoi(argv[2]);
    
    DOYTime dt(year, doy,0.0,"tsGPS");
    
    GTime t;
    t.SetFromDoyTime(dt);
    CivilTime ct = GTime::GTime2CivilTime(t);
    
    //NavTime nt = GTime::GTime2NavTime(t);
    
    printf("%04d %02d %02d %02d %02d %.8f \n", ct.m_year, ct.m_month, ct.m_day, ct.m_hour, ct.m_minute, ct.m_second);
    
    return 0;
}
