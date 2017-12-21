//
//  GForceModelMgr.cpp
//  GFC
//
//  Created by lizhen on 22/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include <stdio.h>

#include "GForceModelMgr.hpp"

namespace gfc
{
    
    void GForceModelMgr::setForceList( std::vector<GString>& forceNames,std::vector<bool>& withPartials )
    {
        
        for( int i = 0 ; i< forceNames.size(); i++ )
        {
            
            bool withPartial = withPartials[i];
            
            // empirical force models
            if( forceNames[i] == "GFMEMP" )
            {
                
                m_forceModels["GFMEMP"] = new GFMSolarRadiationPressureEM();
                
                m_forceModels["GFMEMP"]->hasPartialDerivatives(withPartial);
                
                if( withPartial == true )
                {
                    totalNum_parameters += m_forceModels["GFMEMP"]->num_param;
                }
                
            }
            
            // the Y bias force
            else if (forceNames[i]=="GFMBFSbias")
            {
                 m_forceModels["GFMBFSbias"] = new GFMBFSbias();
            }
            else if( forceNames[i]=="GFMGravity" ) //gravity model
            {
                m_forceModels["GFMGravity"] = new GFMEarthGravity(gravity_degree,gravity_order);
                //m_forceModels["GFMGravity"] = new GFMGravity(gravity_degree,gravity_order);
                
                
                //((GFMGravity*)m_forceModels["GFMGravity"])->loadEGM2008(gravity_file);
                ((GFMEarthGravity*)m_forceModels["GFMGravity"])->loadEGM2008(gravity_file);
                
                //printf("end of loading gravity coefficient file\n");
                m_forceModels["GFMGravity"]->hasPartialDerivatives(withPartial);
                
            }
            else if( forceNames[i]=="GFMNbody" )
            {
                m_forceModels["GFMNbody"] = new GFMThirdBody();
                
                ((GFMThirdBody*)m_forceModels["GFMNbody"])->setBodies(GSpaceEnv::planetsUsed);
                
                m_forceModels["GFMNbody"]->hasPartialDerivatives(withPartial);
            }
            else if(forceNames[i] == "GFMSRP" ) // solar radiation pressure model
            {
                m_forceModels["GFMSRP"] = new GFMSolarRadiationPressure();
                m_forceModels["GFMSRP"]->hasPartialDerivatives(withPartial);
            }
            else if(forceNames[i] == "GFMERP" ) // earth radiation pressure model
            {
                m_forceModels["GFMERP"] = new GFMEarthRadiationPressure();
            }
            else if(forceNames[i] == "GFMANT") // antenna thrust recoil force
            {
                m_forceModels["GFMANT"] = new GFMAntennaThrust();
            }
            else if(forceNames[i] == "GFMTRR") // thermal reradiation recoil force, incluing SOLOR PANEL, the calculation of bus is inside grid file
            {
                m_forceModels["GFMTRR"] = new GFMThermalRadiationForce();
            }
            else if( forceNames[i] == "GFMGR")  // the general relativity
            {
                m_forceModels["GFMGR"] = new GFMGeneralRelativity();
            }
            else if(forceNames[i] == "GFMEarthTide") // the solid earth tide
            {
                m_forceModels["GFMEarthTide"] = new GFMEarthTide();
            }
            
        }
        
    } // end of function setForceList
    
    /*
     *  
        the structure of dydx matrix
     
     | u        v      w       ax      ay      az     |
     | dax/dx, dax/dy, dax/dz, dax/du, dax/dv, dax/dw |
     | day/dx, day/dy, day/dz, day/du, day/dv, day/dw |
     | daz/dx, daz/dy, daz/dz, daz/du, daz/dv, daz/dw |
     
     ct is in UTC
     */
    
    void GForceModelMgr::getDerivatives( GTime ct, int n, double x, double *y, double *dydx, GSpaceCraft* spacecraft )
    {
        
        GVector a_eci;
        GVector a_ecef;
        
        GMatrix dadr(3,3); // partial derivatives of accleration w.r.t position
        GMatrix dadv(3,3); // partial derivatives of accleration w.r.t velocity
        GMatrix dadp;//(3,num_parameters);      // partial derivatives of accelertion w.r.t parameters
        if( totalNum_parameters > 0 )
        {
            dadp.resize(3, totalNum_parameters);
        }
        //velocity
        dydx[0] = y[3];
        dydx[1] = y[4];
        dydx[2] = y[5];
        
        std::map< GString, GForceModel* >::iterator   it = m_forceModels.begin();
        for( ;it != m_forceModels.end(); ++it )
        {
            
            // empirical force models
            if( it->first == "GFMEMP" )
            {
                
                ((GFMSolarRadiationPressureEM *)it->second)->doCompute(spacecraft);
                
                GVector t_eci = it->second->getForce()/ spacecraft->getSpaceCraftGemotry()->m_mass;  //m/s^2
                
                a_eci += t_eci ;
                
                dadr += it->second->m_dadr;
                
                //need to consider the other force model parameters
                // may need to change in the future
                if( it->second->m_hasPartialDerivatives == true )
                {
                    dadp += it->second->m_dadp;
                }
                
            }
            else if( it->first == "GFMGravity" )
            {
                
               // ((GFMGravity*)it->second)->doCompute(ct, GSpaceEnv::eop.getUT1mUTC(), spacecraft->getStatePointer()->satpos_ecef);
                //GVector testP, testV, dif;
                //spacecraft->getEphValue(GTime::UTC2GPST(ct)).getPV(testP, testV);
                
                GVector p = spacecraft->getStatePointer()->satpos_ecef;
                
                ((GFMEarthGravity*)it->second)->doCompute(ct, GSpaceEnv::eop.getUT1mUTC(), p);
               // ((GFMGravity*)it->second)->doCompute(ct, GSpaceEnv::eop.getUT1mUTC(), p);
                
                
                //dif = testP - p;
                //dif = dif*1000.0;
                //printf("%s ,%.6f, %.6f, %.6f, \n",GTime::GTime2CivilTime(ct).TimeString().c_str(), dif.x, dif.y, dif.z);
                
                GVector f_eci = it->second->getForce();
                
                a_eci += f_eci;
                dadr += it->second->m_dadr;
            }
            else if(it->first == "GFMERP")
            {
                
                ((GFMEarthRadiationPressure*)it->second)->doCompute(spacecraft->getStatePointer(),
                                                                    spacecraft->getSpaceCraftGemotry()
                                                                    );
                
                GVector f = it->second->getForce() /spacecraft->getSpaceCraftGemotry()->m_mass;
                
                a_eci += f;
                
                dadr += it->second->m_dadr;
                
            }
            else if( it->first == "GFMNbody") // the thirdbody gravity
            {
                
                ((GFMThirdBody*)it->second)->doCompute(spacecraft->getStatePointer()->satpos_eci);
                
                GVector f_eci =  it->second->getForce();
                
                a_eci += (f_eci );
                dadr += it->second->m_dadr;
                
                //cout<<dadr;
                
            }
            else if(it->first == "GFMEarthTide")
            {
                ((GFMEarthTide*)it->second)->doCompute(spacecraft->getStatePointer()->satpos_eci,GSpaceEnv::planetPos_eci[GJPLEPH::SUN],GSpaceEnv::planetPos_eci[GJPLEPH::MOON]);
                
                GVector f_eci =  it->second->getForce();
                
                a_eci += (f_eci );
                dadr += it->second->m_dadr;
            }
            else if( it->first == "GFMSRP" ) // the solar radiation pressure acceleration
            {
                /*
                ((GFMSolarRadiationPressure *)it->second)->doCompute( spacecraft->getStatePointer()->solarFlux,
                                                                      spacecraft->getStatePointer()->getAttitude(),
                                                                      spacecraft->getSpaceCraftGemotry()
                                                                     );
                */
                
                ((GFMSolarRadiationPressure *)it->second)->doCompute( spacecraft );
                
                GVector t_eci = it->second->getForce()/spacecraft->getSpaceCraftGemotry()->m_mass;  //m/s^2
                
                a_eci += t_eci ;
                
                dadr += it->second->m_dadr;
                
            }
            
            // the Y bias force
            else if( it->first == "GFMBFSbias" )
            {
                GVector f_eci;
                ((GFMBFSbias *)it->second)->setParam(spacecraft->getSpaceCraftGemotry()->m_bias);
                GVector& nx = spacecraft->getStatePointer()->getAttitude()->xhat;
                GVector& ny = spacecraft->getStatePointer()->getAttitude()->yhat;
                GVector& nz = spacecraft->getStatePointer()->getAttitude()->zhat;
                ((GFMBFSbias *)it->second)->doCompute(nx,ny,nz);
                
                f_eci =  it->second->getForce();
                a_eci += (f_eci/spacecraft->getSpaceCraftGemotry()->m_mass);
                dadr += it->second->m_dadr;
            }
            
            else if( it->first == "GFMANT" ) // antenna thrust force
            {
                GVector f_eci;
                
                // the radiation of antenna should be in +z bfs direction
                GVector power_bfs(0,0,spacecraft->getSpaceCraftGemotry()->m_antennaPower);// (0,0,85) in bfs, should be transformed to ECI
                
                GVector power_eci;
                
                GVector& nx = spacecraft->getStatePointer()->getAttitude()->xhat;
                
                GVector& ny = spacecraft->getStatePointer()->getAttitude()->yhat;
                
                GVector& nz = spacecraft->getStatePointer()->getAttitude()->zhat;
                
                power_eci.x = nx.x * power_bfs.x + ny.x * power_bfs.y + nz.x * power_bfs.z;
                power_eci.y = nx.y * power_bfs.x + ny.y * power_bfs.y + nz.y * power_bfs.z;
                power_eci.z = nx.z * power_bfs.x + ny.z * power_bfs.y + nz.z * power_bfs.z;
                
                ((GFMAntennaThrust *)it->second)->setPower(power_eci);
                
                ((GFMAntennaThrust *)it->second)->doCompute();
                
                f_eci =  it->second->getForce();
                
                a_eci += (f_eci/spacecraft->getSpaceCraftGemotry()->m_mass);
                dadr += it->second->m_dadr;
                
            }
            else if( it->first == "GFMTRR" )  // the thermal reradiation pressure for solar panel
            {
                ((GFMThermalRadiationForce *)it->second)->doCompute(spacecraft);
                GVector f_eci = it->second->getForce();
                a_eci += (f_eci/spacecraft->getSpaceCraftGemotry()->m_mass);
                
            }
            else if(it->first == "GFMGR")  // general relativity
            {
                double gm = 398600.4415; // unit: in km
                ((GFMGeneralRelativity *)it->second)->doCompute(gm, spacecraft->getStatePointer()->satpos_eci,
                                                                spacecraft->getStatePointer()->satvel_eci);
                GVector f_eci = it->second->getForce();
                a_eci += f_eci;
                dadr += it->second->m_dadr;
            }
            
        }
        
        
        GVector t_eci ;
        GSpaceEnv::eop.ECEF2ECI_pos(a_ecef, t_eci);
        a_eci += t_eci;
        
        //acceleration must be in ECI, X Y and Z direction respectively.
        dydx[3] =  a_eci.x/1000.0;
        dydx[4] =  a_eci.y/1000.0;
        dydx[5] =  a_eci.z/1000.0;
        
        //set dadr,
        if( n > 6 ) // for orbit fitting problems
        {
            //dadr, because dadv =0
            dydx[6] = dadr[0];dydx[7] = dadr[1];dydx[8] = dadr[2];
            dydx[12] = dadr[3];dydx[13] = dadr[4];dydx[14] = dadr[5];
            dydx[18] = dadr[6];dydx[19] = dadr[7];dydx[20] = dadr[8];
        }
        if( totalNum_parameters > 0 ) // for the sensitivity matrix problem
        {
            
            
            //printf("dadp:\n");
            //cout<<dadp<<endl;
            
            //set dadp
            for(int i = 0 ; i< 3; i++ )
            {
                for( int j = 0 ; j<totalNum_parameters;j++)
                {
                    dydx[24+i*totalNum_parameters+j] = dadp(i,j);
                }
            }
        }
        
        
        //printf("dadr:\n");
        //dadr.print();
        
        int testc = 0;
        
        
    } // end of function getDerivatives
    
}
