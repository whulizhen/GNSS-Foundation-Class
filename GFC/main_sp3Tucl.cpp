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
    /*
    argv[1] = "ssGAL";
    argv[2] = "11";
    argv[3] = "/Users/lizhen/study/the_way_to_PhD/small_Program/sp3Tucl/codesp3/com17912.sp3";
    argv[4] = "/Users/lizhen/study/the_way_to_PhD/small_Program/sp3Tucl/codesp3/com17913.sp3";
    argv[5] = "/Users/lizhen/study/the_way_to_PhD/small_Program/sp3Tucl/codesp3/com17914.sp3";
    argv[6] = "2014";
    argv[7] = "127";
    argv[8] = "";
    */
    
    GString sp3file, outputdir, sysName=argv[1];
    int prn = atoi(argv[2]);
    outputdir = argv[8];
    
    int gps_year = atoi(argv[6]), gps_doy=atoi(argv[7]);
    
    DOYTime dt(gps_year,gps_doy,0.0,"tsGPS");
    
    
    //GSpacecraftModelMgr::initialiseModel( "spacecraft.model");
    
    GSpaceCraftMgr::Initializer();
    
    //GJPLEPH::loadEphFile_a("jpleph405.data");
    
    GEarthOrientationParameter::loadEOP("eopc04_IAU2000.62-now");
    
    GSensorID myid(sysName, prn); // C06 and C14 , G11, E11(IOV101)
    
    GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
    sp3file = argv[3];
    GSpaceCraftMgr::loadPreciseEphemeris(sp3file); // wum18000.sp3,97
    sp3file = argv[4];
    GSpaceCraftMgr::loadPreciseEphemeris(sp3file); // wum18000.sp3,97
    sp3file = argv[5];
    GSpaceCraftMgr::loadPreciseEphemeris(sp3file); // wum18000.sp3,97
    
    
    
    
    GVector pos_eci, vel_eci, pos_ecef, vel_ecef;
    
    GPreciseEphemeris pEph;
    
    GTime epoch_gps, epoch_utc;
    
    GEarthOrientationParameter eop;
    
    //epoch_gps = mysat.getEphValue(1).getEpoch();
    
    //epoch_utc = GTime::GPST2UTC(epoch_gps);
    CivilTime ct; //(utc_year,utc_month,utc_day,0,0,0.0,"tsUTC");
    epoch_gps.SetFromDoyTime(dt);
    
    NavTime nt = GTime::GTime2NavTime(epoch_gps);
    
    epoch_utc = GTime::GPST2UTC(epoch_gps);
    
    ct = GTime::GTime2CivilTime(epoch_gps);
    char tmp[1024];
    sprintf(tmp, "%s%02d-%04d%02d%02d",sysName.c_str(),prn,ct.m_year,ct.m_month,ct.m_day);
    
    GString outputFileName(tmp);
    outputFileName = outputdir + outputFileName;
    outputFileName += ".txt";
    ofstream dumpFile;
    dumpFile.open(outputFileName);
    double interval = 900.0;
    int num = int(86400.0/interval);
    GTime endtime = epoch_utc + 86400-interval;
    while(epoch_utc < endtime)
    {
        
        epoch_gps = GTime::UTC2GPST(epoch_utc);
        
        pEph = mysat.getEphValue(epoch_gps);
        
        ct = GTime::GTime2CivilTime( epoch_utc );
        
        if( pEph.isOK() == false )
        {
            printf("UTC: %s ,the current precise ephemeris is not available!\n",ct.TimeString().c_str());
            continue;
        }
        
        eop.setEpochTime(epoch_utc);
        
        pEph.getPV(pos_ecef, vel_ecef);
        
        eop.ECEF2ECI_pos(pos_ecef, pos_eci);
        
        eop.ECEF2ECI_vel(pos_ecef, vel_ecef, vel_eci);
        
        
        //dumpFile << 46<< ct.TimeString();
        dumpFile << std::setfill('0')<<prn
        << " " << std::setw(4) << (ct.m_year)
        << " " << std::setw(2) << (ct.m_month)
        << " " << std::setw(2) << ct.m_day
        << " " << std::setw(2) << ct.m_hour
        << " " << std::setw(2) << ct.m_minute
        
        << std::setprecision(2) << std::fixed
        << " " << std::setw(5) << ct.m_second
        
        << std::setprecision(6) << std::setfill(' ')
        << " " << std::setw(13) << pos_eci.x
        << " " << std::setw(13) << pos_eci.y
        << " " << std::setw(13) << pos_eci.z
        
        << std::setprecision(8)
        << " " << std::setw(11) << vel_eci.x
        << " " << std::setw(11) << vel_eci.y
        << " " << std::setw(11) << vel_eci.z
        << std::endl;
        
        epoch_utc = epoch_utc + interval;
        
    }
    
    dumpFile.close();
    
    return 0;
}
