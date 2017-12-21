//
//  GSLRValidation.cpp
//  GFC
//
//  Created by lizhen on 05/09/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GSLRValidation.hpp"


namespace gfc
{
   void GSLRValidation::setSpacecraft(gfc::GSpaceCraft* spaceVehicle)
    {
        m_spacecraft = spaceVehicle;
        m_it  = m_slrdata.find(spaceVehicle->getSpaceCraftID());
        if( m_it == m_slrdata.end())
        {
            printf("ERROR: the current spacecraft does NOT exist in the slr obs file\n");
        }
    }
       
    void GSLRValidation::loadDataFile(gfc::GString filename)
    {
        m_datafile.setdatabase(filename);
        
        while ( m_datafile.nextDataBlock(m_slrdata)) { }
        
        m_datafile.closedatbase();
        
    }
    
    void GSLRValidation::computeResidual2()
    {
        
       
        std::map<GTime,slrResInfo> slrRes;
        
        FILE* testFile = fopen("slr.res", "w+");
        
        std::map<GSourceID,GSLRBlock>::iterator myit;
        double PI = GCONST("PI");
        GEllipsoid myellipsoid = GEllipsoidMgr::GetEllipsoid("WGS84");
        
        GVector satpos_ecef, satpos_eci, satvel_ecef, satvel_eci, stapos_ecef, stapos_eci;
        GVector sat_offset , satNEU;
        GPreciseEphemeris preEph;
        
        double xyz[3], blh[3];
        double latitude = 0.0;  // degree
        double height = 0.0; // meter
        double pressure = 0.0, temperature =0.0, wvp = 14.322 , lambda =0.0, ztd, zhd, zwd, ele;
        
        int count = 0;
        
        for( myit = m_it->second.begin(); myit != m_it->second.end(); myit++ )
        {
            //get the information about the station
            
            int station_code = myit->second.h2.cdp_pad_id;
            
            //according this station_code to get the coordinate and velocity information
            GSLRStation sta = GSLRStation::stationInfo[station_code];
            
            // how to get these eccentricity information for every station
            // ftp://cddis.gsfc.nasa.gov/slr/slrocc/ecc_une.snx
            // ftp://cddis.gsfc.nasa.gov/slr/slrocc/ecc_xyz.snx
            
            double dxyz[3] = {0.0};
            myellipsoid.NEU2XYZ(sta.neu, sta.staP, dxyz);
            
            stapos_ecef.x = sta.staP[0];stapos_ecef.y = sta.staP[1];stapos_ecef.z = sta.staP[2];
            
            myellipsoid.XYZ2BLH(sta.staP, blh);
            latitude = blh[0]*180.0/PI;
            height = blh[2];
            
            lambda = myit->second.c0.xmit_wavelength/1000.0;  // unit: micrometer
            myit->second.getTP(0.0, temperature, pressure);  // mba and kalvin
            
            // 1 hpa  = 100000 mpa, hectopascal
            
            for( int i = 0 ; i< myit->second.obsdata.size(); i++ )
            {
                // the epoch time at the laser sending time, need to be corrected
                GTime epoch_utc = myit->second.obsdata[i].obsepoch;
                
                // get satellite position and velocity in ECEF, the reference point is COM
                GTime gpsTime = GTime::UTC2GPST(epoch_utc);
                
                JDTime jdt_utc = GTime::GTime2JDTime(epoch_utc);
                JDTime jdt = GTime::GTime2JDTime(gpsTime);
                CivilTime ct_gps = GTime::JDTime2CivilTime(jdt);
                CivilTime ct_utc = GTime::JDTime2CivilTime(jdt_utc);
                
                double dxtide[3];
                
                if( myit->second.getSourceName() == "YARL")
                {
                    int testc = 0;
                }
                
                //GIERS::DEHANTTIDEINEL(sta.staP, jdt_utc.jdt(), GTime::getLeapSecond(ct.m_year, ct.m_month, ct.m_day),sunpos , moonpos, dxtide);
                
                //must update the attitde and state vector of the satellite
                GSpaceEnv::updateSpaceEnvironment(epoch_utc);
                
                GSpaceEnv::eop.ECEF2ECI_pos(stapos_ecef, stapos_eci);
                
                
                // need to make some corrections to the observation
                // THIS IS WRONG, the upside and downside are NOT euqal !!!!
                double obs = myit->second.obsdata[i].getDataValue() ;  //unit: meter
                preEph = m_spacecraft->getEphValue(gpsTime);
                if( preEph.isOK() == false)
                {
                    continue;  // next obs
                }
                // the sending time , used to calculate the satellite position when receiving the signal
                if( myit->second.obsdata[i].epoch_event == 2 )
                {
                    
                   
                    
                    double dt1 = obs/2.0; // initial value for uplink
                    double dt2 = obs/2.0; // initial value for downlink
                    GTime  t1, t2;
                   
                    while(1)
                    {
                        t1 = gpsTime + dt1;
                        preEph = m_spacecraft->getEphValue(t1);
                        preEph.getPV(satpos_ecef, satvel_ecef);
                        
                        //must update the attitde and state vector of the satellite
                        GSpaceEnv::updateSpaceEnvironment(GTime::GPST2UTC(t1));
                        satpos_ecef = satpos_ecef*1000.0;
                        satvel_ecef = satvel_ecef*1000.0;
                        GSpaceEnv::eop.ECEF2ECI_pos(satpos_ecef, satpos_eci);
                        
                        double dis = (satpos_eci - stapos_eci).norm();
                        
                       if( fabs(  dis - dt1*299792458.0) < 1.0E-4 )
                        {
                            break;
                        }
                        else
                        {
                            dt1 = dis/299792458.0;
                        }
                        
                    }  // the iteration for the uplink
                    
                    
                    //GSpaceEnv::updateSpaceEnvironment(epoch_utc + dt1);
                    
                    GSpaceEnv::eop.ECEF2ECI_vel(satpos_ecef, satvel_ecef, satvel_eci);
                    m_spacecraft->getStatePointer()->updateState_eci(epoch_utc+dt1, satpos_eci/1000.0, satvel_eci/1000.0);
                    
                    m_spacecraft->getOffsetCorrection(1, sat_offset);  // 1 for SLR, 0 for GNSS
                    ct_utc= GTime::GTime2CivilTime(epoch_utc + dt1);
                    
                    double sunpos[3] = {GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].x*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].y*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].z*1000.0};
                    double moonpos[3] = {GSpaceEnv::planetPos_ecef[GJPLEPH::MOON].x*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::MOON].y*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::MOON].z*1000.0};
                    double hours = ct_utc.m_hour+ct_utc.m_minute/60.0+ct_utc.m_second/3600.0;
                    GIERS::DEHANTTIDEINEL(sta.staP, ct_utc.m_year, ct_utc.m_month, ct_utc.m_day,hours , GTime::getLeapSecond(ct_utc.m_year, ct_utc.m_month, ct_utc.m_day), sunpos, moonpos, dxtide);
                    
                    
                    //the position of the satellite is determined, and the laser reflected time as well.
                    myellipsoid.XYZ2NEU(stapos_ecef, satpos_ecef, satNEU);
                    double len = satNEU.norm();
                    ele = asin(satNEU.z/len)*180.0/PI;
                    
                   
                    while(1)
                    {
                        GVector tmpstation;
                        t2 = t1 + dt2;
                        GSpaceEnv::updateSpaceEnvironment(GTime::GPST2UTC(t2) );
                        GSpaceEnv::eop.ECEF2ECI_pos(stapos_ecef, tmpstation);
                        
                        double dis = (satpos_eci - tmpstation).norm() ;
                        
                        if( fabs(dis- dt2*299792458.0 ) < 1.0E-4 )
                        {
                            break;
                        }
                        else
                        {
                            dt2 = dis/299792458.0;
                        }
                    }
                    
                    // error models: relativistic corrections, tropspheric delay, clock error , system calibration .et al
                    
                    GIERS::FCULZD_HPA(latitude, height, pressure, wvp, lambda, ztd, zhd, zwd);
                    double mapping =  GIERS::FCUL_A(latitude, height, temperature, ele);
                    
                    double tropdelay = ztd * mapping;
                    
                    obs = obs*299792458.0;
                    
                    GVector los = (satpos_eci - stapos_eci);
                    los.normalise();
                    
                    //convert ecef to eci
                    double tmp[3] = {0.0};
                    //GSpaceEnv::eop.ECEF2ECI_pos(dxtide, tmp);
                    GSpaceEnv::eop.ECEF2ECI(0,dxtide,tmp);
                    memcpy(dxtide,tmp,sizeof(double)*3);
                    
                    double dxyz1[3] = {0.0};
                    GSpaceEnv::eop.ECEF2ECI(0,dxyz,dxyz1);
                    
                    GVector tvector2;
                    GSpaceEnv::eop.ECEF2ECI_pos(sat_offset, tvector2);
                    
                    double tidecorrection =  -(dxtide[0]*los.x + dxtide[1]*los.y + dxtide[2]*los.z);
                    double sta_offset     =  -(dxyz1[0]*los.x + dxyz1[1]*los.y + dxyz1[2]*los.z);
                    double lra_offset     =  dotproduct(tvector2, los);
                    
                    // relativity is about 2-3 cm
                    double rel_error   = 0.0;
                    
                    double cal = (dt1 + dt2)*299792458.0 + 2.0*tropdelay + 2.0*lra_offset + 2.0*sta_offset + 2.0*tidecorrection ;
                    
                    // get rid of the gross error
                    //if( fabs(obs-cal) < 3.0  )
                    {
                        slrResInfo resinfo;
                        resinfo.res = obs - cal;
                        resinfo.obs = obs;
                        resinfo.geodis1 = dt1;
                        resinfo.geodis2 = dt2;
                        
                        resinfo.earthTide = tidecorrection;
                        resinfo.ele = ele;
                        resinfo.latitude = latitude;
                        resinfo.pressure = pressure;
                        
                        resinfo.satOffset = lra_offset;
                        resinfo.staOffset = sta_offset;
                        resinfo.temperature = temperature;
                        resinfo.stationName = myit->second.getSourceName();
                        resinfo.trop = tropdelay;
                        
                        resinfo.rel_error = rel_error;
                        
                        slrRes[t1] = resinfo;
                    }
                    
//                    fprintf(testFile, "%s %s %6.3f %6.3f %6.3f %8.3f %8.3f %8.3f %8.3f %16.12f %16.12f %16.12f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %8.4f \n", ct_utc.TimeString().c_str(),  myit->second.getSourceName().c_str(),
//                            temperature, pressure,ele,mapping*ztd,
//                            lra_offset, sta_offset, tidecorrection, obs, dt1, dt2,
//                            satpos_ecef.x, satpos_ecef.y, satpos_ecef.z,
//                            satpos_eci.x, satpos_eci.y, satpos_eci.z,
//                            obs -cal
//                            );
                    
                    int testc = 0;
                    
                    count++;
                }
            }
        }
        
        //output the slr residual
        std::map<GTime,slrResInfo>::iterator oit;
        for( oit = slrRes.begin(); oit != slrRes.end(); oit++)
        {
            // get satellite position and velocity in ECEF, the reference point is COM
            JDTime jdt = GTime::GTime2JDTime(oit->first);
            CivilTime ct = GTime::JDTime2CivilTime(jdt);
            
            // printf("%s , %.8f\n", ct.TimeString().c_str(), oit->second);
            
            fprintf(testFile, "%s,%s,%6.3f,%6.3f,%6.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f, %16.12f,%16.12f,%16.12f,%8.4f,\n", ct.TimeString().c_str(),
                    oit->second.stationName.c_str(),
                    oit->second.temperature,
                    oit->second.pressure,
                    oit->second.ele,
                    oit->second.trop,
                    oit->second.satOffset, oit->second.staOffset, oit->second.earthTide,oit->second.rel_error,
                    oit->second.obs,
                    oit->second.geodis1,
                    oit->second.geodis2,
                    oit->second.res
                    );
        }
        
        fclose(testFile);
        
        int testc =0;
        
    }
    
    
    
    void GSLRValidation::computeResidual()
    {
        
        std::map<GTime,slrResInfo> slrRes;
        
        FILE* testFile = fopen(m_outputFileName.c_str(), "w+");
        
        std::map<GSourceID,GSLRBlock>::iterator myit;
        double PI = GCONST("PI");
        GEllipsoid myellipsoid = GEllipsoidMgr::GetEllipsoid("WGS84");
        
        GVector satpos_ecef, satvel_ecef, stapos_ecef, stapos1_ecef, stapos2_ecef;
        GVector lrabias_ecef ;
        GPreciseEphemeris preEph;
        
        double  blh[3], stabias_ecef[3];
        double latitude = 0.0;  // degree
        double height = 0.0; // meter
        double pressure = 0.0, temperature =0.0, wvp = 14.322 , lambda =0.532, ztd, zhd, zwd, ele;
        int count = 0;
         double dxtide[3];
        DOYTime doy_epoch(2005,1,0.0,"tsUTC");
        GTime station_epoch;
        station_epoch.SetFromDoyTime(doy_epoch);
        
        for( myit = m_it->second.begin(); myit != m_it->second.end(); myit++ )
        {
            //get the information about the station
            
            int station_code = myit->second.h2.cdp_pad_id;
            
            //according this station_code to get the coordinate and velocity information
            GSLRStation sta = GSLRStation::stationInfo[station_code];
            
            if( fabs(sta.staP[0]) < 1.0 )  // current station is unavailable
            {
                continue;
            }
            
            //apply the corrections of reference frame, the velocity of the station coordinates
            
            printf("%s %.3f %.3f %.3f\n", myit->second.h2.stn_name, sta.staP[0], sta.staP[1], sta.staP[2]);
            
            // how to get these eccentricity information for every station
            // ftp://cddis.gsfc.nasa.gov/slr/slrocc/ecc_une.snx
            // ftp://cddis.gsfc.nasa.gov/slr/slrocc/ecc_xyz.snx
            
             // station ENU
            myellipsoid.NEU2XYZ( sta.neu, sta.staP,stabias_ecef);
            
            
            myellipsoid.XYZ2BLH(sta.staP, blh);
            latitude = blh[0]*180.0/PI;
            height = blh[2];
            
            lambda = myit->second.c0.xmit_wavelength/1000.0;  // unit: micrometer
            
            myit->second.getTP(0.0, temperature, pressure);  // mba and kalvin
            
            // 1 hpa  = 100000 mpa, hectopascal
            
            for( int i = 0 ; i< myit->second.obsdata.size(); i++ )
            {
                // the epoch time at the laser sending time, need to be corrected
                GTime epoch_utc = myit->second.obsdata[i].obsepoch;
                
                double yr_duration = (epoch_utc - station_epoch).toDays()/365.25;
                
                //yr_duration = 0.0;
                
                stapos_ecef.x =  sta.staP[0] + sta.staV[0]*yr_duration;
                stapos_ecef.y =  sta.staP[1] + sta.staV[1]*yr_duration;
                stapos_ecef.z =  sta.staP[2] + sta.staV[2]*yr_duration;
                
                // get satellite position and velocity in ECEF, the reference point is COM
                GTime gpsTime = GTime::UTC2GPST(epoch_utc);
                JDTime jdt_utc = GTime::GTime2JDTime(epoch_utc);
                JDTime jdt = GTime::GTime2JDTime(gpsTime);
                CivilTime ct_gps = GTime::JDTime2CivilTime(jdt);
                CivilTime ct_utc = GTime::JDTime2CivilTime(jdt_utc);
                
                //GIERS::DEHANTTIDEINEL(sta.staP, jdt_utc.jdt(), GTime::getLeapSecond(ct.m_year, ct.m_month, ct.m_day),sunpos , moonpos, dxtide);
                
                stapos1_ecef = stapos_ecef;
                stapos2_ecef = stapos_ecef;
                
                double obs = myit->second.obsdata[i].getDataValue() ;  //unit: meter
                
                preEph = m_spacecraft->getEphValue(gpsTime);
                
                if( preEph.isOK() == false )
                {
                    continue;  // next obs
                }
                // the sending time , used to calculate the satellite position when receiving the signal
                if( myit->second.obsdata[i].epoch_event == 2 )
                {
                    
                    double dt1 = obs/2.0; // initial value for uplink
                    double dt2 = obs/2.0; // initial value for downlink
                    GTime  t1, t2;
                    int test1 =0;
                    while(1)
                    {
                        t1 = gpsTime + dt1;
                        
                        //t1 = epoch_utc + dt1;
                        
                        preEph = m_spacecraft->getEphValue(t1);
                        preEph.getPV(satpos_ecef, satvel_ecef);
                        
                        //must update the attitde and state vector of the satellite
                        satpos_ecef *= 1000.0 ;
                        satvel_ecef *= 1000.0 ;
                        
                        // rotation for station coordinate
                        double rotationAngle = myellipsoid.AV() * dt1;
                        stapos1_ecef.x = stapos_ecef.x*cos(rotationAngle) + stapos_ecef.y*sin(rotationAngle);
                        stapos1_ecef.y = stapos_ecef.y*cos(rotationAngle) - stapos_ecef.x*sin(rotationAngle);
                        stapos1_ecef.z = stapos_ecef.z;
                        
//                        // rotation for the station offset
//                        dxyz1[0] = dxyz[0]*cos(rotationAngle) + dxyz[1]*sin(rotationAngle);
//                        dxyz1[1] = dxyz[1]*cos(rotationAngle) - dxyz[0]*sin(rotationAngle);
//                        dxyz1[2] = dxyz[2];
//
                        double dis = (satpos_ecef - stapos1_ecef).norm();
                        
                        if( fabs( dis - dt1*299792458.0) < 1.0E-4 )
                        {
                            break;
                        }
                        else
                        {
                            dt1 = dis/299792458.0;
                        }
                        
                    }  // the iteration for the uplink
                    
                    GTime  signalReceivedUTC = epoch_utc + dt1;
                    
                    ct_utc= GTime::GTime2CivilTime(signalReceivedUTC);
                    GSpaceEnv::updateSpaceEnvironment(signalReceivedUTC);
                    
                    // the updation must be in km
                    m_spacecraft->getStatePointer()->updateState_ecef(signalReceivedUTC, satpos_ecef/1000.0, satvel_ecef/1000.0);
                   
                    m_spacecraft->getOffsetCorrection(1, lrabias_ecef);  // 1 for SLR, 0 for GNSS
                    
                    GVector satpos_eci;
                    GSpaceEnv::eop.ECEF2ECI_pos(satpos_ecef, satpos_eci);
                    
                    double sunpos[3] = {GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].x*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].y*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].z*1000.0};
                    double moonpos[3] = {GSpaceEnv::planetPos_ecef[GJPLEPH::MOON].x*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::MOON].y*1000.0,
                        GSpaceEnv::planetPos_ecef[GJPLEPH::MOON].z*1000.0};
                    
                    //the position of the satellite is determined, and the laser reflected time as well.
                    GVector satNEU;
                    myellipsoid.XYZ2NEU(stapos1_ecef, satpos_ecef, satNEU);
                    double len = satNEU.norm();
                    ele = asin(satNEU.z/len)*180.0/PI;
                    
                    // the downlink
                    while(1)
                    {
                        // rotate from receiving time to reflected time
                        double rotationAngle = -myellipsoid.AV() * dt2;
                        
                        stapos2_ecef.x = stapos_ecef.x*cos(rotationAngle) + stapos_ecef.y*sin(rotationAngle);
                        stapos2_ecef.y = stapos_ecef.y*cos(rotationAngle) - stapos_ecef.x*sin(rotationAngle);
                        stapos2_ecef.z = stapos_ecef.z;
                        
                        double dis = (satpos_ecef - stapos2_ecef).norm();
                        
                        if( fabs( dis - dt2*299792458.0) < 1.0E-4 )
                        {
                            break;
                        }
                        else
                        {
                            dt2 = dis/299792458.0;
                        }
                    }
                    
                    // error models: relativistic corrections, tropspheric delay, clock error , system calibration .et al
                    
                    GIERS::FCULZD_HPA(latitude, height, pressure, wvp, lambda, ztd, zhd, zwd);
                    double mapping =  GIERS::FCUL_A(latitude, height, temperature, ele);
                    
                    double trop = mapping*ztd;
                    
                    obs = obs*299792458.0;
                    
                    double statmp[3] = {stapos1_ecef.x,stapos1_ecef.y,stapos1_ecef.z};
                    
                    double fhr = ct_utc.m_hour+ct_utc.m_minute/60.0+ct_utc.m_second/3600.0;
                    GIERS::DEHANTTIDEINEL(statmp, ct_utc.m_year, ct_utc.m_month, ct_utc.m_day, fhr, GTime::getLeapSecond(ct_utc.m_year, ct_utc.m_month, ct_utc.m_day), sunpos, moonpos, dxtide);
                    
                    GVector los = (satpos_ecef -  stapos1_ecef); // ecef vector
                    
                    double geodis = los.norm();
                    
                    los.normalise();
                    
                    // relativity makes distance shorter
                    double GM = 3.986004418E14;
                    double gamma = 1.0;
                    double rel_error1 = (1+gamma)*GM/299792458.0/299792458.0*log( (stapos1_ecef.norm() + satpos_ecef.norm() + dt1*299792458)/(stapos1_ecef.norm() + satpos_ecef.norm() - dt1*299792458) ) ;
                    
                    double tidecorrection =   -(dxtide[0]*los.x + dxtide[1]*los.y + dxtide[2]*los.z);
                    
                    //tidecorrection = sqrt(pow(stapos1_ecef.x + dxtide[0]- satpos_ecef.x ,2.0) +pow(stapos1_ecef.y + dxtide[1]- satpos_ecef.y ,2.0)+pow(stapos1_ecef.z + dxtide[2]- satpos_ecef.z ,2.0) ) - geodis;
                    
                    double sta_offset     =   -(stabias_ecef[0]*los.x + stabias_ecef[1]*los.y + stabias_ecef[2]*los.z);
                    
                   // sta_offset =sqrt(pow(stapos1_ecef.x + stabias_ecef[0]- satpos_ecef.x ,2.0) +pow(stapos1_ecef.y + stabias_ecef[1]- satpos_ecef.y ,2.0)+pow(stapos1_ecef.z + stabias_ecef[2]- satpos_ecef.z ,2.0) ) - geodis;
                    
                    double lra_offset     =  dotproduct(lrabias_ecef, los);
                    
                    //lra_offset = -0.251;
                    //lra_offset = (stapos1_ecef + lrabias_ecef - satpos_ecef).norm() - geodis;
                    
                    double cal = (dt1 + dt2 )*299792458.0 + 2.0*trop + 2.0*lra_offset + 2.0*sta_offset + 2.0*tidecorrection + 2.0*rel_error1;
                    
                    // get rid of the gross error
                    if( fabs(obs-cal) < 1.0  )
                    {
                        slrResInfo resinfo;
                        resinfo.res = obs - cal;
                        
                        resinfo.obs = obs;
                        resinfo.geodis1 = dt1*299792458.0;
                        resinfo.geodis2 = dt2*299792458.0;
                        
                        resinfo.earthTide = tidecorrection;
                        resinfo.ele = ele;
                        resinfo.beita = m_spacecraft->getStatePointer()->attitude_eci.beta*180.0/PI;
                        resinfo.eps = m_spacecraft->getStatePointer()->eps*180.0/PI;
                        resinfo.latitude = latitude;
                        resinfo.pressure = pressure;
                        
                        resinfo.satOffset = lra_offset;
                        resinfo.staOffset = sta_offset;
                        resinfo.temperature = temperature;
                        resinfo.stationName = myit->second.getSourceName();
                        resinfo.trop = trop;
                        
                        resinfo.los[0] = los.x;
                        resinfo.los[1] = los.y;
                        resinfo.los[2] = los.z;
                        
                        resinfo.dlra[0] = lrabias_ecef.x;
                        resinfo.dlra[1] = lrabias_ecef.y;
                        resinfo.dlra[2] = lrabias_ecef.z;
                        
                        resinfo.dxyz[0] = stabias_ecef[0];
                        resinfo.dxyz[1] = stabias_ecef[1];
                        resinfo.dxyz[2] = stabias_ecef[2];
                        
                        resinfo.dtide[0] = dxtide[0];
                        resinfo.dtide[1] = dxtide[1];
                        resinfo.dtide[2] = dxtide[2];
                        
                        resinfo.rel_error = rel_error1;
                        
                        GTime time_utc = GTime::GPST2UTC(t1);
                        
                        slrRes[time_utc] = resinfo;
                        
                    }
                    
//                    fprintf(testFile, "%s %s %6.3f %6.3f %6.3f %8.3f %8.3f %8.3f %8.3f %16.12f %16.12f %16.12f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %8.4f \n", ct_utc.TimeString().c_str(),  myit->second.getSourceName().c_str(),
//                            temperature, pressure,ele,mapping*ztd,
//                            lra_offset, sta_offset, tidecorrection, obs, dt1, dt2,
//                            satpos_ecef.x, satpos_ecef.y, satpos_ecef.z,
//                            satpos_eci.x,  satpos_eci.y,  satpos_eci.z,
//                            obs -cal
//                            );
//
                    
                    int testc = 0;
                    
                    count++;
                }
            }
        }
        
        
        
        fprintf(testFile, "timeUTC,JDT,staName,temperature,pressure,ele,beta,eps,trop,relativity,com_error,sta_eccentriciy,erathTide,obs,uplink,downlink,LOS_x,LOS_y,LOS_z,LRA_x,LRA_y,LRA_z,staOffset_x,staOffset_y,staOffset_z,tide_x,tide_y,tide_z,res,\n");
        
        //output the slr residual
        std::map<GTime,slrResInfo>::iterator oit;
        for( oit = slrRes.begin(); oit != slrRes.end(); oit++)
        {
            // get satellite position and velocity in ECEF, the reference point is COM
            JDTime jdt = GTime::GTime2JDTime(oit->first);
            CivilTime ct = GTime::JDTime2CivilTime(jdt);
            
            GString  timestr( jdt.jdt() );
            
            fprintf(testFile, "%s,%s,%s,%6.3f,%6.3f,%6.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%16.12f,%16.12f,%16.12f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.4f,\n",
                    ct.TimeString().c_str(),
                    timestr.c_str(),
                    oit->second.stationName.c_str(),
                    oit->second.temperature,
                    oit->second.pressure,
                    oit->second.ele,
                    oit->second.beita,
                    oit->second.eps,
                    oit->second.trop,
                    oit->second.rel_error,
                    oit->second.satOffset,oit->second.staOffset,oit->second.earthTide,
                    oit->second.obs,oit->second.geodis1,oit->second.geodis2,
                    
                    oit->second.los[0],oit->second.los[1],oit->second.los[2],
                    
                    oit->second.dlra[0],oit->second.dlra[1],oit->second.dlra[2],
                    oit->second.dxyz[0],oit->second.dxyz[1],oit->second.dxyz[2],
                    oit->second.dtide[0],oit->second.dtide[1],oit->second.dtide[2],
                    oit->second.res/2.0
                    
                    );
            
        }
        
        fclose(testFile);
        
        int testc =0;
        
    }
    
    
    
    
} // end of namespace
