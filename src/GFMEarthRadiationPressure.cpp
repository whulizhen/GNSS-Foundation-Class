//
//  GFMEarthRadiationPressure.cpp
//  GFC
//
//  Created by lizhen on 16/4/12.
//  Copyright © 2016年 lizhen. All rights reserved.
//


//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GFC is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GFC; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Copyright 2015, lizhen
//
//============================================================================

#include "GFMEarthRadiationPressure.hpp"
namespace gfc
{
    
    
    // fortran library in the code from Solano
    extern "C"
    {
        void myerptest_( int* ERM,int* ANT, int* GRD, int* REFF,double* YSAT,double* SUN,double* KAPPA,int* MONTH, int* BLKNUM,double* ACCEL);
    }
    
    
    GVector GFMEarthRadiationPressure::erp_solar_panel( GEarthRadiationFlux* earthRadiationflux, GSpaceCraftAttitude* svAttitude,GSpaceCraftModel* space_craft_geometry)
    {
        GVector force;
        
        //determine the specularity and reflectivity of solar panel according to the backside and foreside
        double specularity = 0.0, reflectivity = 0.0;
        double panel_area = space_craft_geometry->solarArray[0].area;
        
        for( int i = 0 ; i< earthRadiationflux->m_flux.size(); i++ )
        {
            double t = dotproduct(earthRadiationflux->m_flux[i].m_dir, svAttitude->phat);
            
            int test = 0;
            GVector f_short, f_long;
            GVector surfaceNomal;
            
            if( t < 0.0 )  // foreside
            {
                specularity = space_craft_geometry->solarArray[0].specularity;
                
                reflectivity = space_craft_geometry->solarArray[0].reflectivity;
                
                surfaceNomal = svAttitude->phat;
                
                //get force first
               // f = GRadiationFlux::radiationForce(surfaceNomal, earthRadiationflux->m_flux[i].m_dir,earthRadiationflux->m_flux[i].m_shortwave,
               //                                            panel_area, specularity, reflectivity);
                test = 1;
            }
            else //backside
            {
                specularity = space_craft_geometry->solarArray[1].specularity;
                
                reflectivity = space_craft_geometry->solarArray[1].reflectivity;
                surfaceNomal = -svAttitude->phat;
                //get force first
                //f = GRadiationFlux::radiationForce(surfaceNomal , earthRadiationflux->m_flux[i],
                //                                           panel_area, specularity, reflectivity);
                test = 0;
            }
            
            
            f_short = GRadiationFlux::radiationForce(surfaceNomal,
                                                     earthRadiationflux->m_flux[i].m_dir,
                                                     earthRadiationflux->m_flux[i].m_shortwave,
                                                     panel_area, specularity, reflectivity);
            
            f_long = GRadiationFlux::radiationForce(surfaceNomal,
                                                    earthRadiationflux->m_flux[i].m_dir,
                                                    earthRadiationflux->m_flux[i].m_longwave,
                                                    panel_area, specularity, reflectivity);
            
            //printf("%d\n",test);
            
        
            //get acceleration
            //a/= space_craft_geometry->m_mass;
            // the vector addition
            force += f_short;
            force += f_long;
            
        }
        
        return force;
        
    } // end of erp_solar_panel
    
    GVector GFMEarthRadiationPressure::erp_bus_simple( GEarthRadiationFlux& earthRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry )
    {
        
        GVector xhat = svAttitude->xhat, yhat=svAttitude->yhat, zhat=svAttitude->zhat;
        
        GVector acc;
        
        double areaX,areaY,areaZ;
        double speX,speY,speZ,refX,refY,refZ;
        
        areaX = space_craft_geometry->busX[0].area;
        areaY = space_craft_geometry->busY[0].area;
        areaZ = space_craft_geometry->busZ[0].area;
        
        speX = space_craft_geometry->busX[0].specularity;
        speY = space_craft_geometry->busY[0].specularity;
        speZ = space_craft_geometry->busZ[0].specularity;
        
        refX = space_craft_geometry->busX[0].reflectivity;
        refY = space_craft_geometry->busY[0].reflectivity;
        refZ = space_craft_geometry->busZ[0].reflectivity;
        
        
        for( int i = 0 ; i< earthRadiationflux.m_flux.size(); i++ )
        
        {
            
            GVector ax_s,ax_l,ay_s,ay_l,az_s,az_l;
            
            fluxdata& radiation = earthRadiationflux.m_flux[i];;
            
            
            //radiation.m_dir = normalise(earthRadiationflux.totalFlux_lw);
            //radiation.m_shortwave = earthRadiationflux.totalFlux_sw.norm();
            //radiation.m_longwave = earthRadiationflux.totalFlux_lw.norm();
            
            double tx = dotproduct(radiation.m_dir, svAttitude->xhat );
            double ty = dotproduct(radiation.m_dir, svAttitude->yhat );
            double tz = dotproduct(radiation.m_dir, svAttitude->zhat );
            
            if( tx > 0.0 )  // the angle between flux and surface normal is less than 90, we want greater than 90
            {
                xhat = -xhat;
                speX = space_craft_geometry->busX[1].specularity;
                refX = space_craft_geometry->busX[1].reflectivity;
            }
            if( ty > 0.0 )
            {
                yhat = -yhat;
                speY = space_craft_geometry->busY[1].specularity;
                refY = space_craft_geometry->busY[1].reflectivity;
            }
            if(tz > 0.0 )
            {
                zhat = -zhat;
                speZ = space_craft_geometry->busZ[1].specularity;
                refZ = space_craft_geometry->busZ[1].reflectivity;
            }
            
             ax_s = GRadiationFlux::radiationForce(xhat, radiation.m_dir, radiation.m_shortwave,
                                                       areaX,
                                                       speX,
                                                       refX);
             ax_l = GRadiationFlux::radiationForce(xhat, radiation.m_dir, radiation.m_longwave,
                                                  areaX,
                                                  speX,
                                                  refX);
            
            
             ay_s = GRadiationFlux::radiationForce(yhat, radiation.m_dir, radiation.m_shortwave,
                                                       areaY,
                                                       speY,
                                                       refY);
            
            ay_l = GRadiationFlux::radiationForce(yhat, radiation.m_dir, radiation.m_longwave,
                                                  areaY,
                                                  speY,
                                                  refY);
            
             az_s = GRadiationFlux::radiationForce(zhat, radiation.m_dir, radiation.m_shortwave,
                                                       areaZ,
                                                       speZ,
                                                       refZ);
            az_l = GRadiationFlux::radiationForce(zhat, radiation.m_dir, radiation.m_longwave,
                                                  areaZ,
                                                  speZ,
                                                  refZ);
            
            acc += (ax_s + ax_l + ay_s + ay_l + az_s + az_l );
           
        }
        
        return acc;
    }
    
    
    GVector GFMEarthRadiationPressure::erp_bus_grid(GEarthRadiationFlux& earthRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry)
    {
        GVector force_bfs, force;
        
        double R2D = GCONST("R2D");
        // make sure that the earthRadiationFlux is in ECI, NOT ECEF
        GVector body_frame;
        for( int i = 0 ; i< earthRadiationflux.m_flux.size(); i++ )
        {
            
            double W = earthRadiationflux.m_flux[i].m_longwave + earthRadiationflux.m_flux[i].m_shortwave;
            
            // the symmetric matrix, transfer from ECI to BFS
            body_frame.x = - dotproduct(svAttitude->xhat, earthRadiationflux.m_flux[i].m_dir);
            
            body_frame.y = - dotproduct(svAttitude->yhat, earthRadiationflux.m_flux[i].m_dir);
            
            body_frame.z = - dotproduct(svAttitude->zhat, earthRadiationflux.m_flux[i].m_dir);
            
            body_frame.normalise();
            
            double lambda = atan2(body_frame.y, body_frame.x )*R2D;
            //-pi to pi
            double phi = 0.0;  // -pi/2 to pi/2
            
            if( fabs(body_frame.z - 1.0)<1.0E-14 )
            {
                phi = 90.0;
            }
            else if(fabs(body_frame.z + 1.0)<1.0E-14)
            {
                phi = -90.0;
            }
            else
            {
                phi = R2D * std::asin(body_frame.z);
            }
            
            force_bfs += GRadiationGridMgr::getAcc_grid(space_craft_geometry->m_srpgriddata, lambda, phi)*W;  //uint: Newton
            
            int testc = 0;
           
        
        }
        
         force = force_bfs.x * svAttitude->xhat + force_bfs.y * svAttitude->yhat +force_bfs.z * svAttitude->zhat ;
        
        //need to transform from BFS to ECI
        // force is in BFS, transform from bfs to eci
//        force.x = force_bfs.x*svAttitude->xhat.x + force_bfs.y*svAttitude->yhat.x + force_bfs.z*svAttitude->zhat.x;
//        force.y = force_bfs.x*svAttitude->xhat.y + force_bfs.y*svAttitude->yhat.y + force_bfs.z*svAttitude->zhat.y;
//        force.z = force_bfs.x*svAttitude->xhat.z + force_bfs.y*svAttitude->yhat.z + force_bfs.z*svAttitude->zhat.z;
        
        return force;
        
    }
    
    
    GVector GFMEarthRadiationPressure::Solano_ERP( GVector& satpos_eci, GVector& satvel_eci)
    {
        //fortran test
        int ERM = 2;  //1 for analytical , 2 for CERES data
        int ANT = 0;
        int GRD = 1;
        int REFF = 0;
        int MONTH = 3;
        int BLKNUM = 4;
        double YSAT[6] = { satpos_eci.x*1000.0,satpos_eci.y*1000.0,satpos_eci.z*1000.0,
            satvel_eci.x*1000.0,satvel_eci.y*1000.0,satvel_eci.z*1000.0};
        
        double SUN[3] = {GSpaceEnv::planetPos_eci[GJPLEPH::SUN].x*1000.0,GSpaceEnv::planetPos_eci[GJPLEPH::SUN].y*1000.0,GSpaceEnv::planetPos_eci[GJPLEPH::SUN].z*1000.0};
        
        double tm[9] ={0.0};
        GSpaceEnv::eop.getECI2ECEFMatrix(tm);
        
//        double KAPPA[3][3] =
//        {
//            {tm[0],tm[1],tm[2]},
//            {tm[3],tm[4],tm[5]},
//            {tm[6],tm[7],tm[8]}
//        };
        
        double KAPPA[3][3] =
        {
            {tm[0],tm[3],tm[6]},
            {tm[1],tm[4],tm[7]},
            {tm[2],tm[5],tm[8]}
        };
        
        double FORCE[3] ={0.0};
        
        //this is the function in the fortran library
        //myerptest_(&ERM, &ANT, &GRD, &REFF, YSAT,SUN, &KAPPA[0][0],&MONTH,&BLKNUM,FORCE);
        
        GVector force;
        force.x = FORCE[0];
        force.y = FORCE[1];
        force.z = FORCE[2];
        
        return force;
        
    }
    
    
    /*
     * get the final earth radiation pressure acceleration
     *
     */
    void GFMEarthRadiationPressure::doCompute( GMotionState* statePointer ,  GSpaceCraftModel* space_craft_geometry )
    {
        //GEarthRadiationFlux& earthRadiationflux , GSpaceCraftAttitude* svAttitude,
        
        GVector force;
        
        force += erp_solar_panel(&(statePointer->earthFlux), statePointer->getAttitude(),space_craft_geometry );
        
        if(with_grid_on == false)
        {
            force += erp_bus_simple(statePointer->earthFlux ,statePointer->getAttitude(), space_craft_geometry);
        }
        else if(with_grid_on == true)
        {
            force += erp_bus_grid(statePointer->earthFlux ,statePointer->getAttitude(), space_craft_geometry);
        }
        
        // other choice
        //force =  Solano_ERP(statePointer.satpos_eci,statePointer.satvel_eci);
         
        setForce(force);
    }
    
    
    
} // end of namespace gfc
