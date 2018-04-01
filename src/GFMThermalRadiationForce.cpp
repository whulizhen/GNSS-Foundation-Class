//
//  GFMThermalRadiationForce.cpp
//  GFC
//
//  Created by lizhen on 10/07/2016.
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

#include "GFMThermalRadiationForce.hpp"


namespace gfc
{
    // Stefan-Boltzmann constant
   double GFMThermalRadiationForce::sigma =  5.670373E-8;
   
    // get the temperature of the solar panel
    
    GVector GFMThermalRadiationForce::ThermalRadiationForce_MLI(gfc::GVector &n, GVector& flux_dir, double flux,double area, double alfa, double emissivity)
    {
        GVector force;
        
        double eff_emssivity = 0.02;
        //double emssivity = 0.84;
        double Tsc = 298*298*298*298;  // the temperature inside satellite is 298K
        double cos_theta = -dotproduct(flux_dir, n);
        double t4 = (alfa*flux*cos_theta + eff_emssivity*sigma*Tsc)/(sigma*(eff_emssivity+emissivity));
        
        force = -2.0*sigma*emissivity*t4*area/3.0/299792458.0*n;
        
        return force;
    }
    
    
    
    void GFMThermalRadiationForce::getTemperatures( double radiation_front, double radiation_back, double power_draw,GSpaceCraftModel* space_craft_geometry)
    {
       
        // need to find a good way to get the approximate solution
        
        m_temperatures[0] = 315.0; // K, T0
        m_temperatures[1] = 315.0; // K, T1
        m_temperatures[2] = 308.0; // K, T2
        GMatrix X(3,1);
        
        GMatrix L(3,1); // the measurements
        GMatrix B(3,3); // the design matrix
        
        double emissivity0 = space_craft_geometry->solarArray[0].emmisivity; // emmisivity for the front side panel
        double emissivity1 = space_craft_geometry->solarArray[1].emmisivity; // emmisivity for the back side panel
        
        double conductivity0 = space_craft_geometry->solarArray[0].conductivity;
        double conductivity1 = space_craft_geometry->solarArray[1].conductivity;
        
        double absorptivity0 = space_craft_geometry->solarArray[0].absorbtivity;
        double absorptivity1 = space_craft_geometry->solarArray[1].absorbtivity;
        
        double thick0 = space_craft_geometry->solarArray[0].thickness;
        double thick1 = space_craft_geometry->solarArray[1].thickness;;
        
        double panel_area = space_craft_geometry->solarArray[0].area;
        
        
        int iteration = 0;
        
        while(1)
        {
            
            m_temperatures[0]+= X[0];
            m_temperatures[1]+= X[1];
            m_temperatures[2]+= X[2];
            
            L[0] = -sigma*emissivity0*pow(m_temperatures[0],4.0) - conductivity0*(m_temperatures[0]-m_temperatures[1])/thick0 + absorptivity0*radiation_front;
            
            L[1] = -conductivity0*(m_temperatures[1]-m_temperatures[0])/thick0 - conductivity1*(m_temperatures[1]-m_temperatures[2])/thick1 - power_draw/panel_area;
            
            L[2] = -sigma*emissivity1*pow(m_temperatures[2],4.0) - conductivity1*(m_temperatures[2]-m_temperatures[1])/thick1 + absorptivity1*radiation_back;
            
            //L = L*(-1.0);
            
            B(0,0) = 4.0*sigma*emissivity0*pow(m_temperatures[0],3.0) + conductivity0/thick0;
            
            B(0,1) = -conductivity0/thick0;
            
            B(0,2) = 0.0;
            
            B(1,0) = -conductivity0/thick0;
            
            B(1,1) = conductivity0/thick0 + conductivity1/thick1;
            
            B(1,2) = -conductivity1/thick1;
            
            B(2,0) = 0.0;
            
            B(2,1) = -conductivity1/thick1;
            
            B(2,2) = 4.0*sigma*emissivity1*pow(m_temperatures[2],3.0) + conductivity1/thick1;
            
            //cout<<"B:"<<endl;
            //cout<<B;
            
            //cout<<~X;
            //cout<< ~L;
            
            GMatrix XX = (!B)*L;
            
            iteration++;
            
            if(
               (fabs(L[0]) < 1.0E-10
               &&fabs(L[1]) < 1.0E-10
               &&fabs(L[2]) < 1.0E-10)
               || iteration >=5
               )
            {
                //cout<<"residual:";
                //cout<<~L;
                //cout<<(~XX);
                //printf("temperature: ");
                //cout<<"temperature: ";
                //cout<< radiation_front << " " <<m_temperatures[0]<<" "<<m_temperatures[1]<<" "<<m_temperatures[2]<<endl;
                //printf("%.6f %.6f %.6f %.6f \n", radiation_front, m_temperatures[0], m_temperatures[1], m_temperatures[2]);
                
                break;
            }
            else
            {
                //cout<<"X-XX"<<endl;
                //cout << ~(X-XX);
                X = XX;
            }
        }
        
    }
    
    
    void GFMThermalRadiationForce::doCompute( GSpaceCraft* spacecraft )
    {
         GEarthRadiationFlux& earthflux = spacecraft->getStatePointer()->earthFlux;
        
         GSolarRadiationFlux& solarflux = spacecraft->getStatePointer()->solarFlux;
        
         // the flux in normal, the solar panel normal is always in front side
         double front_flux =0.0, back_flux = 0.0;
        
        //get the flux on backside and front side, only the front side can generate electricity power
         double cos_theta_earth_flux_l = dotproduct(spacecraft->getStatePointer()->attitude_eci.phat, earthflux.totalFlux_lw ) ;
        if( cos_theta_earth_flux_l <= 0.0 ) // front side longwave
        {
             front_flux += -cos_theta_earth_flux_l;
        }
        else
        {
             back_flux += cos_theta_earth_flux_l;
        }
        
        double cos_theta_earth_flux_s = dotproduct(spacecraft->getStatePointer()->attitude_eci.phat, earthflux.totalFlux_sw ) ;
        if( cos_theta_earth_flux_l <= 0.0 ) // front side shortwave
        {
            front_flux += -cos_theta_earth_flux_s;
        }
        else
        {
            back_flux += cos_theta_earth_flux_s;
        }
        
        double cos_theta_solar_flux = dotproduct(spacecraft->getStatePointer()->attitude_eci.phat, solarflux.m_flux.m_dir ) ;
        if( cos_theta_solar_flux < 0.0 ) // front side
        {
            front_flux += -(solarflux.m_flux.m_longwave + solarflux.m_flux.m_shortwave)*cos_theta_solar_flux;
        }
        else
        {
            back_flux += (solarflux.m_flux.m_longwave + solarflux.m_flux.m_shortwave)*cos_theta_solar_flux;
        }
        
        // make sure that power_draw can not be larger than normal flux, especially in eclipse
        // still need to consider the efficientcy of the solar panel
        // only the front side of solar panel is used to collect energy
        double power_draw = spacecraft->getSpaceCraftGemotry()->m_solarPanelPowerDraw;
        
        // The efficiency of the solar cells should be considered
        double efficiency = 0.2;
        
        if( front_flux <= spacecraft->getSpaceCraftGemotry()->m_solarPanelPowerDraw )
        {
            power_draw = front_flux;
        }
        
        //front_flux = 1360.0;
        //back_flux = 0.0;
        double c = GCONST("CLIGHT");
        
        getTemperatures(front_flux, back_flux, power_draw, spacecraft->getSpaceCraftGemotry());
        
        
        //calculate the acceleration
        double magnitude =2.0/3.0*spacecraft->getSpaceCraftGemotry()->solarArray[0].area *sigma*
        (   spacecraft->getSpaceCraftGemotry()->solarArray[1].emmisivity*pow( m_temperatures[2],4.0)
         -  spacecraft->getSpaceCraftGemotry()->solarArray[0].emmisivity*pow( m_temperatures[0],4.0)
         )/c;
        
        GVector force = magnitude*spacecraft->getStatePointer()->attitude_eci.phat;
        
        setForce(force);
        
    }
    
    
} // end of namespace gfc

