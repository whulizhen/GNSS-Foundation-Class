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
 outputfilename
 number of sp3 ephemeris
 number of normal point files
 
 normal point file1
 normal point file2
 normal point file3
 ...
 
 sp3file1
 sp3file2
 sp3file3
 sp3file4
 sp3file5
 sp3file6
 ...
 
 */

int main( int argc, char** argv )
{
    
//    std::map<GTime, double> testMap;
//    GTime t1,t2,t3,t4;
//    t1.SetData(TimeSystem("tsGPS"), 52000, 1000, 0.432);
//    
//    t2 = t1 + 1;
//    t3 = t2 + 1;
//    t4 = t1;
//    testMap[t1] = 0.0;
//    testMap[t2] = 1.0;
//    testMap[t3] = 2.0;
//    testMap[t4] = 3.0;
//    
//    bool bt1 = t4 < t1;
    
    GString sysName="ssGAL";
    int prn = 11;
    int sp3fileNum = 6;
    int npFileNum = 1;
    
    sysName = argv[1];
    prn = atoi(argv[2]);
    
    GString npFile = "galileo101_201501.npt";
    GString filename =  "slrtest.res";
    
    filename = argv[3];
    sp3fileNum = atoi(argv[4]);
    npFileNum  = atoi(argv[5]);
    
    // normal point file
    //argv[6] = "galileo101_201501.npt";
    //sp3 file test
    //argv[7] = "sp3/com18254.sp3";
    //argv[8] = "sp3/com18255.sp3";
    //argv[9] = "sp3/com18256.sp3";
    //argv[10] = "sp3/com18260.sp3";
    //argv[11] = "sp3/com18261.sp3";
    //argv[12] = "sp3/com18262.sp3";
    
    GSpacecraftModelMgr::initialiseModel( "spacecraft.model");
    GSpaceCraftMgr::Initializer();
    
    GJPLEPH::loadEphFile_a("jpleph405.data");
    
    GEarthOrientationParameter::loadEOP("eopc04_IAU2000.62-now");
    
    GSLRStation::loadStationCoordinate("SLRF2008_140210_2014.03.24.snx");
    GSLRStation::loadStationEccentricity("ecc_une.snx");
    
    for( int i = 0 ; i< sp3fileNum; i++ )
    {
        int index = 6 + npFileNum  + i;
        printf("ephemeris: %s\n",argv[index]);
        GSpaceCraftMgr::loadPreciseEphemeris(argv[ index]);
    }
    
    GSensorID myid(sysName, prn); // C08  C06 and C14 , G11, E11(IOV101)
    
    GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
    
    GSLRValidation slrvalidation;
    
    printf("outputfilename:%s\n",filename.c_str());
    
    slrvalidation.setOutputFile(filename);
    for(int i = 0 ; i< npFileNum; i++)
    {
        printf("NP file: %s\n",argv[6+i]);
        slrvalidation.loadDataFile( argv[6+i] ); // lageos2_201604.npt, slrTest.frd compassg1_201404.npt compassi3_201404.npt
    }
    
    slrvalidation.setSpacecraft(&mysat);
    
    slrvalidation.computeResidual();
    
       
}
