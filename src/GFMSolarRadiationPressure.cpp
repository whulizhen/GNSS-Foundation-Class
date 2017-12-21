//
//  GFMSolarRadiationPressure.cpp
//  GFC
//
//  Created by lizhen on 18/05/2016.
//  Copyright © 2016 lizhen. All rights reserved.
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

#include "GFMSolarRadiationPressure.hpp"
#include "GSpaceCraftAttitude.hpp"

namespace   gfc
{
    
    
    /*
     the newly derived box-wing model used for analyzing empirical model
     only suitalbe for nominal attitude
     */
    GVector GFMSolarRadiationPressure::newBoxWing(GSolarRadiationFlux& solarRadiationflux, GSpaceCraftAttitude* svAttitude,GSpaceCraftModel* space_craft_geometry)
    {
        
        GVector force_bfs;
        
        double speX,speY,speZ,refX,refY,refZ,speP,refP;
        speX = space_craft_geometry->busX[1].specularity;
        speY = space_craft_geometry->busY[0].specularity;
        speZ = space_craft_geometry->busZ[0].specularity;
        speP = space_craft_geometry->solarArray[0].specularity;
        
        refX = space_craft_geometry->busX[1].reflectivity;
        refY = space_craft_geometry->busY[0].reflectivity;
        refZ = space_craft_geometry->busZ[0].reflectivity;
        refP = space_craft_geometry->solarArray[0].reflectivity;
        
        double areaX = space_craft_geometry->busX[1].area;
        double areaY = space_craft_geometry->busY[0].area;
        double areaZ = space_craft_geometry->busZ[0].area;
        double areaP = space_craft_geometry->solarArray[0].area;
        
        double W_C = solarRadiationflux.m_flux.m_shortwave/299792458.0;
        double phi = svAttitude->phi;
        
        
        //GVector f_solarpanel = W_C*areaP*(1.0+2.0/3.0*refP + 1.0/3.0*refP*speP)*solarRadiationflux.m_flux.m_dir;
        //GVector s(cos(phi),0,-sin(phi));
        //GVector nx(-1, 0,0);
        //GVector nz(0,0,1);
        //double nds = -sin(phi);
        //GVector fz = areaZ*( (refZ*speZ-1.0)*nds*s  + (2.0/3.0*refZ*(1.0-speZ)*nds -2.0*speZ*refZ*nds*nds)*nz  );
        
        
        force_bfs.x = areaP*(1.0+2.0/3.0*refP + 1.0/3.0*refP*speP)*cos(phi)
        
                       +2.0/3.0*areaX*refX*(1.0-speX)*cos(phi)
        
                       +0.5*areaX*(1.0+refX*speX)*(1.0+cos(2.0*phi));
        
        
        force_bfs.z =  -areaP*(1.0+2.0/3.0*refP + 1.0/3.0*refP*speP)*sin(phi)
        
        
                        - 0.5*areaX*(1.0-refX*speX)*sin(2.0*phi);
        
        
        
        if( phi >= 0.0 ) // +z
        {
            force_bfs.x += 0.5*areaZ*(1.0-speZ*refZ)*sin(2.0*phi);
            force_bfs.z += -(2.0/3.0*areaZ*refZ*(1.0-speZ))*sin(phi) + 0.5*areaZ*(1.0+refZ*speZ)*(cos(2.0*phi) - 1.0) ;
        }
        else  // -z
        {
            speZ = space_craft_geometry->busZ[1].specularity;
            refZ = space_craft_geometry->busZ[1].reflectivity;
            
            force_bfs.x -= 0.5*areaZ*(1.0-speZ*refZ)*sin(2.0*phi);
            force_bfs.z += -(2.0/3.0*areaZ*refZ*(1.0-speZ))*sin(phi) - 0.5*areaZ*(1.0+refZ*speZ)*(cos(2.0*phi) - 1.0) ;
        }
        
        //fz        = fz*W_C;
        force_bfs = force_bfs*W_C;
        
        GVector force_eci = svAttitude->xhat*force_bfs.x + svAttitude->yhat*force_bfs.y + svAttitude->zhat*force_bfs.z;
        GVector force_dyb;
        force_dyb.x = -cos(phi)*force_bfs.x + sin(phi)*force_bfs.z;
        force_dyb.y = -force_bfs.y;
        force_dyb.z = sin(phi)*force_bfs.x  + cos(phi)*force_bfs.z;
        
        GVector force_eci2 = svAttitude->eD*force_dyb.x + svAttitude->eY*force_dyb.y + svAttitude->eB*force_dyb.z;
        
        
        
        
        return  force_eci;
        
    }
    
    
    
    
    //return force
    GVector  GFMSolarRadiationPressure::srp_solar_panel( GSolarRadiationFlux& solarRadiationflux, GSpaceCraftAttitude* svAttitude,GSpaceCraftModel* space_craft_geometry )
    {
        GVector force;
        
        GVector n  = svAttitude->phat ;  // the normal of solar panel in ECI
        
        //determine the specularity and reflectivity of solar panel according to the backside of foreside
        double specularity = 0.0, reflectivity = 0.0;
        
        double t = dotproduct(solarRadiationflux.m_flux.m_dir, n);
        
        //correct the flux magnitude according to the incident angle (cos)
        //solarRadiationflux.m_flux.m_longwave = solarRadiationflux.m_flux.m_longwave * (-t);
        
        //solarRadiationflux.m_flux.m_shortwave = solarRadiationflux.m_flux.m_shortwave * (-t);
        double panel_area  = space_craft_geometry->solarArray[0].area;
        
        if( t <= 0.0 )  // foreside, sun light will always be on the frontsight of the solar panel
        {
            specularity = space_craft_geometry->solarArray[0].specularity;
            
            reflectivity = space_craft_geometry->solarArray[0].reflectivity;
        }
        else //backside
        {
            specularity = space_craft_geometry->solarArray[1].specularity;
            
            reflectivity = space_craft_geometry->solarArray[1].reflectivity;
        }
        
        //solar radiation is regarded as shortwave
        force += GRadiationFlux::radiationForce(n, solarRadiationflux.m_flux.m_dir, solarRadiationflux.m_flux.m_shortwave, panel_area, specularity, reflectivity);
        
        
        /*
        GVector testforce = solarRadiationflux.m_flux.m_shortwave*panel_area/299792458.0*(1.0+2.0/3.0*reflectivity+1.0/3.0*reflectivity*specularity)*solarRadiationflux.m_flux.m_dir;
        */
        
        return force;
    }
    
    
    GVector  GFMSolarRadiationPressure::srp_bus_simple(  GSolarRadiationFlux& solarRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry  )
    {
        
        GVector force;
        
        GVector xhat = svAttitude->xhat;
        GVector yhat = svAttitude->yhat;
        GVector zhat = svAttitude->zhat;
        double phi = svAttitude->phi;
        double tx = dotproduct(solarRadiationflux.m_flux.m_dir, xhat);
        double ty = dotproduct(solarRadiationflux.m_flux.m_dir, yhat);
        double tz = dotproduct(solarRadiationflux.m_flux.m_dir, zhat);
        
        double speX=0.0,speY=0.0,speZ=0.0,refX=0.0,refY=0.0,refZ=0.0,areaX=0.0,areaY=0.0,areaZ=0.0;
        
        // +x, -x, +y, -y, +z, -z
        int illumination[6]={0};
        
        /*
        speX = space_craft_geometry->busX[0].specularity;
        speY = space_craft_geometry->busY[0].specularity;
        speZ = space_craft_geometry->busZ[0].specularity;
        
        refX = space_craft_geometry->busX[0].reflectivity;
        refY = space_craft_geometry->busY[0].reflectivity;
        refZ = space_craft_geometry->busZ[0].reflectivity;
        
        double areaX = space_craft_geometry->busX[0].area;
        double areaY = space_craft_geometry->busY[0].area;
        double areaZ = space_craft_geometry->busZ[0].area;
        */
        
        // 卫星本体的受照情况
        if(fabs(tx) >1.0E-10)
        {
            if( tx > 0.0 )
            {
                xhat = -xhat;
                speX = space_craft_geometry->busX[1].specularity;
                refX = space_craft_geometry->busX[1].reflectivity;
                areaX =space_craft_geometry->busX[0].area;
                //printf("-1, %.3f\n", phi*180/3.1415926);
                illumination[0] = -1;
            }
            else
            {
                speX = space_craft_geometry->busX[0].specularity;
                refX = space_craft_geometry->busX[0].reflectivity;
                areaX =space_craft_geometry->busX[0].area;
                //printf("+1, %.3f \n",phi*180/3.1415926);
                illumination[1] = 1;
            }
            
        }
        
        if(fabs(ty) >1.0E-10)
        {
            if( ty > 0.0 )
            {
                yhat = -yhat;
                speY = space_craft_geometry->busY[1].specularity;
                refY = space_craft_geometry->busY[1].reflectivity;
                areaY =space_craft_geometry->busY[0].area;
                //printf("-2,%.3f\n",phi*180/3.1415926);
                illumination[2] = -3;
            }
            else
            {
                speY = space_craft_geometry->busY[0].specularity;
                refY = space_craft_geometry->busY[0].reflectivity;
                areaY =space_craft_geometry->busY[0].area;
                //printf("+2, %.3f\n",phi*180/3.1415926);
                illumination[3] = 3;
            }
            
        }
        
        if(fabs(tz) >1.0E-10)
        {
            if( tz > 0.0 )
            {
                zhat = -zhat;
                speZ = space_craft_geometry->busZ[1].specularity;
                refZ = space_craft_geometry->busZ[1].reflectivity;
                areaZ =space_craft_geometry->busZ[0].area;
                //printf("-3, %.3f\n",phi*180/3.1415926);
                illumination[4] = -2;
            }
            else
            {
                speZ = space_craft_geometry->busZ[0].specularity;
                refZ = space_craft_geometry->busZ[0].reflectivity;
                areaZ =space_craft_geometry->busZ[0].area;
                //printf("+3,%.3f\n",phi*180/3.1415926);
                illumination[5] = 2;
            }
            
        }
        
        
        printf("%d %d %d %d %d %d %.6f\n", illumination[0],illumination[1],
                                       illumination[2],illumination[3],
                                       illumination[4],illumination[5],
                                      svAttitude->eta*180/M_PI);
        
//        
//        if( tx > 0.0 )  // the angle between flux and surface normal is less than 90, we want greater than 90
//        {
//            xhat = -xhat;
//            speX = space_craft_geometry->busX[1].specularity;
//            refX = space_craft_geometry->busX[1].reflectivity;
//        }
//        if( ty > 0.0 )
//        {
//            yhat = -yhat;
//            speY = space_craft_geometry->busY[1].specularity;
//            refY = space_craft_geometry->busY[1].reflectivity;
//        }
//        if( tz > 0.0 )
//        {
//            zhat = -zhat;
//            speZ = space_craft_geometry->busZ[1].specularity;
//            refZ = space_craft_geometry->busZ[1].reflectivity;
//        }
//        
        
        
        
        
        
        
        // solar radiation is regarded as shortwave
        
        GVector fx = GRadiationFlux::radiationForce(xhat, solarRadiationflux.m_flux.m_dir,
                                                    solarRadiationflux.m_flux.m_shortwave,
                                                    areaX,
                                                    speX,
                                                    refX);
        
        
        GVector fy = GRadiationFlux::radiationForce(yhat, solarRadiationflux.m_flux.m_dir,
                                                    solarRadiationflux.m_flux.m_shortwave,
                                                    areaY,
                                                    speY,
                                                    refY);
        
        GVector fz = GRadiationFlux::radiationForce(zhat, solarRadiationflux.m_flux.m_dir,
                                                    solarRadiationflux.m_flux.m_shortwave,
                                                    areaZ,
                                                    speZ,
                                                    refZ);
        
        GVector fx_mli =  GFMThermalRadiationForce::ThermalRadiationForce_MLI(xhat,solarRadiationflux.m_flux.m_dir,
                                                                              solarRadiationflux.m_flux.m_shortwave,
                                                                              areaX,1.0-refX,0.84);
        
        GVector fy_mli =  GFMThermalRadiationForce::ThermalRadiationForce_MLI(yhat,solarRadiationflux.m_flux.m_dir,
                                                                              solarRadiationflux.m_flux.m_shortwave,
                                                                              areaY,1.0-refY,0.84);

        GVector fz_mli =  GFMThermalRadiationForce::ThermalRadiationForce_MLI(zhat,solarRadiationflux.m_flux.m_dir,
                                                                              solarRadiationflux.m_flux.m_shortwave,
                                                                              areaZ,1.0-refZ,0.84);

        
        
        
        force += (fx + fy + fz + fx_mli + fy_mli + fz_mli);
        
        return force;
        
    }
    
    
    
    GVector GFMSolarRadiationPressure::srp_bus_M(double e, GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry)
    {
        GVector force;
        double aC_a_d = 14.5;  // nm/s^2
        double aS_a_d = 5.0 ;  // nm/s^2
        
        double a_d = -aC_a_d*( fabs(cos(e)) + sin(e) + 2.0/3.0 ) - aS_a_d*( fabs(cos(e)) - sin(e) - 4.0/3.0*sin(e)*sin(e) + 2.0/3.0 ) ;
        
        double a_b = -4.0/3.0*aS_a_d*(cos(e)*sin(e));
        
        force = 1.0E-9*(svAttitude->eD * a_d + svAttitude->eB * a_b)*space_craft_geometry->m_mass;
        
        return force;
    }
    
    GVector  GFMSolarRadiationPressure::srp_bus_grid(GSolarRadiationFlux& solarRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry )
    {
        GVector force;
        
        double R2D = GCONST("R2D");
        
        //double radiationDir_bfs[3] = {0.0};
        GVector body_frame_sun;
        
        double W = solarRadiationflux.m_flux.m_longwave + solarRadiationflux.m_flux.m_shortwave;
        
        // the symmetric matrix eci to bfs
        //here minus sign, because the direction of solar radiation, the BFS is from satellite to sun
        // while the real solar radiation is from sun to satellite
        body_frame_sun.x = - dotproduct(svAttitude->xhat, solarRadiationflux.m_flux.m_dir);
        
        body_frame_sun.y = - dotproduct(svAttitude->yhat, solarRadiationflux.m_flux.m_dir);
        
        body_frame_sun.z = - dotproduct(svAttitude->zhat, solarRadiationflux.m_flux.m_dir);
        
        body_frame_sun.normalise();
        
        //longitude -pi to pi  , body_frame_sun.y 的正负决定-180和+180
        double lambda = atan2(body_frame_sun.y, body_frame_sun.x )*R2D;
        
        //latitude  -pi/2 to pi/2
        double phi = 0.0;
        
        if( fabs(body_frame_sun.z - 1.0)<1.0E-14 )
        {
            phi = 90.0;
        }
        else if(fabs(body_frame_sun.z + 1.0)<1.0E-14)
        {
            phi = -90.0;
        }
        else
        {
            phi = std::asin(body_frame_sun.z)*R2D ;
        }
        
        
        //printf("lat: %f lon: %f \n", phi, lambda);
        //lambda = 180.0;
        
        GVector force_bfs = GRadiationGridMgr::getAcc_grid( space_craft_geometry->m_srpgriddata, lambda, phi)*W ;
        
        // force is in BFS, transform from bfs to eci
        force =  force_bfs.x * svAttitude->xhat + force_bfs.y * svAttitude->yhat +force_bfs.z * svAttitude->zhat ;
        
        //need to transform from BFS to ECI
        // force is in BFS, transform from bfs to eci
        
        //force.x = force_bfs.x*svAttitude->xhat.x + force_bfs.y*svAttitude->yhat.x + force_bfs.z*svAttitude->zhat.x;
        //force.y = force_bfs.x*svAttitude->xhat.y + force_bfs.y*svAttitude->yhat.y + force_bfs.z*svAttitude->zhat.y;
        //force.z = force_bfs.x*svAttitude->xhat.z + force_bfs.y*svAttitude->yhat.z + force_bfs.z*svAttitude->zhat.z;
        
        return  force;
        
    } // end of srp_bus_grid
    
    
    //should account for the ECLIPSE situation
    void GFMSolarRadiationPressure::doCompute( GSpaceCraft* spacecraft)
    {
        
        GVector force;
        
        GVector tf;
        
        force += srp_solar_panel(spacecraft->getStatePointer()->solarFlux,
                                 &(spacecraft->getStatePointer()->attitude_eci),
                                 spacecraft->getSpaceCraftGemotry()
                                  );
        
        if(with_grid_on == false)
        {
            tf += srp_bus_simple( spacecraft->getStatePointer()->solarFlux,
                                  &(spacecraft->getStatePointer()->attitude_eci),
                                 spacecraft->getSpaceCraftGemotry());
            
            //cose = cosb * cosu, Yar bar-sever, this is just an approximation
            //double e = acos(cos(spacecraft->getStatePointer()->attitude_eci.beta)*cos(spacecraft->getStatePointer()->attitude_eci.eta));
            /*
            tf += srp_bus_M(
                            spacecraft->getStatePointer()->eps,
                            &(spacecraft->getStatePointer()->attitude_eci),
                            spacecraft->getSpaceCraftGemotry()
                           );
            */
        }
        else if(with_grid_on == true)
        {
           tf += srp_bus_grid( spacecraft->getStatePointer()->solarFlux ,
                              &(spacecraft->getStatePointer()->attitude_eci),
                              spacecraft->getSpaceCraftGemotry() );
        }
        
        force += tf;
        
        
        //GVector testforce = newBoxWing(solarRadiationflux,svAttitude,space_craft_geometry);
        
        //transfer force from ECI to BFS for test
        
        
        
        
        setForce(force);
        
    }
    
    
}
