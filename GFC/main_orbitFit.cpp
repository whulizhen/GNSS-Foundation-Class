//
//  testERP.cpp
//  GFC
//
//  Created by lizhen on 16/2/4.
//  Copyright © 2016年 lizhen. All rights reserved.
//

#include <stdio.h>
#include <algorithm>
#include <typeinfo>


//#include "../../src/GNetcdf.h"
#include "../../src/GEarthOrientationParameter.hpp"
#include "../../src/GRungeKutta.hpp"
#include "GEarthRadiationModel.h"
#include "GFMEarthRadiationPressure.hpp"
#include "GFMThermalRadiationForce.hpp"


#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Gnomonic.hpp>
#include <GeographicLib/PolygonArea.hpp>



#include "GOrbitPredictor.hpp"
#include "GSpaceEnv.hpp"
#include "GOrbitFitting.hpp"
#include "GGlobalPressureTemperature.hpp"
#include "GOBSData.h"


#include "GRinex.hpp"

#include "GSLRcrd.hpp"

#include "GTemperatureTest.hpp"

#include "GSLRValidation.hpp"


using namespace gfc;
using namespace GeographicLib;

//
//  testERP.cpp
//  GFC
//
//  Created by lizhen on 16/2/4.
//  Copyright © 2016年 lizhen. All rights reserved.
//

using namespace gfc;


/*
 
 sysName
 prn
 start_year start_doy
 end_year end_day
 output_filename
 number of sp3 ephemeris
 sp3file1
 sp3file2
 sp3file3
 sp3file4
 sp3file5
 sp3file6
 ...
 
 */

int main(int argc, char** argv)
{

    /*
     default setting
     */
    GString sysName,outputName;
    int satprn;
    GTime startTime, endTime;
    int num_sp3 ;
    
//    argv[1] = "ssGAL";  //system name
//    argv[2] = "11";      // satllite prn
//    argv[3] = "2015";    // start year
//    argv[4] = "186";      // start doy
//    argv[5] = "2015";    // end year
//    argv[6] = "187";      // end doy
//    argv[7] = "testOrbitFit.txt"; // output name
//    argv[8] = "3";  // num of sp3
//    argv[9] = "sp3/codsp3-2015/com18516.sp3";   // sp3 filename1
//    argv[10] = "sp3/codsp3-2015/com18520.sp3";  // sp3 filename2
//    argv[11] = "sp3/codsp3-2015/com18521.sp3";  // sp3 filename3

        argv[1] = "ssGAL";  //system name
        argv[2] = "11";      // satllite prn
        argv[3] = "2015";    // start year
        argv[4] = "2";      // start doy
        argv[5] = "2015";    // end year
        argv[6] = "3";      // end doy
        argv[7] = "testOrbitFit.txt"; // output name
        argv[8] = "3";  // num of sp3
        argv[9] = "/Users/lizhen/study/the_way_to_PhD/small_Program/data/whusp3-2015/wum18254.sp3";   // sp3 filename1
        argv[10] = "/Users/lizhen/study/the_way_to_PhD/small_Program/data/whusp3-2015/wum18255.sp3";  // sp3 filename2
        argv[11] = "/Users/lizhen/study/the_way_to_PhD/small_Program/data/whusp3-2015/wum18256.sp3";  // sp3 filename3
    
    
//    //during eclipse
//        argv[1] = "ssGAL";  //system name
//        argv[2] = "11";      // satllite prn
//        argv[3] = "2015";    // start year
//        argv[4] = "3";      // start doy
//        argv[5] = "2015";    // end year
//        argv[6] = "4";      // end doy
//        argv[7] = "testOrbitFit.txt"; // output name
//        argv[8] = "3";  // num of sp3
//        argv[9] = "sp3/codsp3-2015/com18255.sp3";   // sp3 filename1
//        argv[10] = "sp3/codsp3-2015/com18256.sp3";  // sp3 filename2
//        argv[11] = "sp3/codsp3-2015/com18260.sp3";  // sp3 filename3
    
    
    sysName = argv[1];
    satprn = atoi(argv[2]);
    DOYTime dt;
    dt.m_year = atoi(argv[3]); //year
    dt.m_doy = atoi(argv[4]); //day of year
    dt.m_ts = TimeSystem("tsGPS");
    startTime.SetFromDoyTime(dt);
    
    dt.m_year = atoi(argv[5]); // the end year
    dt.m_doy =  atoi(argv[6]); // the end day
    endTime.SetFromDoyTime(dt);
    
    outputName = argv[7]; // the output filename
    num_sp3 = atoi(argv[8]);
    
    //endTime = startTime + 3600*5;
    
    
    GSpacecraftModelMgr::initialiseModel( "spacecraft.model");
    GSpaceCraftMgr::Initializer();
    
    //load erp and srp grid files
    GEarthRadiationFlux::populateEMData();
    
//    CivilTime cct = GTime::GTime2CivilTime(startTime);
//    GString longwavefile, shortwavefile;
//    char tmpstr1[1024],tmpstr2[1024];
//    sprintf(tmpstr1, "flux/%02d/longwave%02d.grid",cct.m_month,cct.m_month);
//    sprintf(tmpstr2, "flux/%02d/shortwave%02d.grid",cct.m_month,cct.m_month);
//    longwavefile = tmpstr1;
//    shortwavefile = tmpstr2;
//    GEarthRadiationFlux::populateCERESGrid("../data/flux/04/longwave04.grid", "../data/flux/04/shortwave04.grid");
//    
    GJPLEPH::loadEphFile_a("jpleph405.data");
    
    GEarthOrientationParameter::loadEOP("eopc04_IAU2000.62-now");
    
    // choose a spacecraft, prn G11, svn46
    //these are for the gps test ()
    for(int i = 0 ; i< num_sp3; i++)
    {
        GSpaceCraftMgr::loadPreciseEphemeris(argv[9+i]);
    }
    
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19142.eph");
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19143.eph");
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19144.eph");
    
    //GSpaceCraftMgr::loadUCLEphemeris("sp3/ucl_svn46_01mar04.eci");
    
    GSensorID myid(sysName, satprn); // C08  C06 and C14 , G11, E11(IOV101)
    GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
    
    //startTime.SetFromCivilTime(ct0);
    //endTime.SetFromCivilTime(ct1);
    //GTime epoch = startTime;
    
    //std::vector<GString> forceList = { "GFMGravity","GFMGR","GFMERP","GFMSRP","GFMNbody","GFMANT" };
    
    std::vector<GString> forceList = { "GFMGravity","GFMGR","GFMNbody" ,"GFMSRP","GFMEMP"};
    
    // this predictor is just for satellite myid right now    
    
    GTime epoch_gps = startTime; // GPST
    
    GTime epoch_utc = GTime::GPST2UTC(epoch_gps);
    endTime   = GTime::GPST2UTC(endTime);
    
    GVector ps( -0.86956146241453E+04, 0.22187345110001E+05 , 0.11159092897201E+05 );
    GVector vs( -0.20432406470742E+01,-0.21481631147695E+01,0.25494488888535E+01 );
    mysat.getEphValue(epoch_gps).getPV(ps, vs);
    
    GSpaceEnv::updateSpaceEnvironment(epoch_utc);
    mysat.getStatePointer()->updateState_ecef(epoch_utc, ps, vs);
    
   // for initial orbit determination
    GOrbitFitting orbfit(&mysat);
    orbfit.setForceList(forceList);
    orbfit.setStepSize(300.0);
    orbfit.setInitialValue(epoch_utc, mysat.getStatePointer()->satpos_eci, mysat.getStatePointer()->satvel_eci);
    orbfit.setEndEpoch(endTime);
    orbfit.setOutputFile(outputName);
    
    orbfit.calculateInitialError();
    
    
    return 0;
}
