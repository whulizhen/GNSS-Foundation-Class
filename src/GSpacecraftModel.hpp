//
//  GSpacecraftModel.hpp
//  GFC
//
//  Created by lizhen on 16/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GSpacecraftModel_hpp
#define GSpacecraftModel_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>

#include "GString.h"

#include "GSensor.h"

#include "GRadiationGrid.hpp"

#include <map>

namespace gfc
{
    
    /*
     * the geometry of space craft
     *
     */
    class GSpaceCraftModel
    {
        public:
        //this optical element is only used in SpaceCraft Model
        struct opticalElement
        {
            GString name;    // the name of this element
            GVector normal;  // the normal vector of this surface element in BFS
            double area;
            
            //these are for the solar radiation property
            double absorbtivity;
            double conductivity;
            double emmisivity;
            double reflectivity;
            double specularity;
            
            //these are for infrared optical property
            double absorbtivity_IR;
            double conductivity_IR;
            double emmisivity_IR;
            double reflectivity_IR;
            double specularity_IR;
            
            double thickness;
            
            opticalElement()
            {
                area = 0.0;
                absorbtivity =0.0;
                conductivity = 0.0;
                emmisivity = 0.0;
                reflectivity = 0.0;
                specularity = 0.0;
                
                absorbtivity_IR =0.0;
                conductivity_IR = 0.0;
                emmisivity_IR = 0.0;
                reflectivity_IR = 0.0;
                specularity_IR = 0.0;
                
                
                thickness = 0.0;
            }
            
        };
        
        double m_mass;       //  the real mass of spacecraft, kg
        
        double m_nomialMass; // the nomial mass for the srp grid file
        
        // center of mass in bfs
        double m_com[3];   // center of mass
        
        // offset of LRA in bfs
        double m_offsetLRA[3]; // the mount offset of Laser Retroreflectory Array
        
        // offset of GNSS Antenna in bfs
        double m_offsetGNSS[3];  // the antenna mount offset of the GNSS L band Antenna
        
        //the properties of the satellite goemetry
        std::vector<opticalElement> busX;
        std::vector<opticalElement> busY;
        std::vector<opticalElement> busZ;
        std::vector<opticalElement> solarArray;
        //this grid data is only for solar radiation
        GMatrix m_srpgriddata[3]; // for x , y and z seperately
        
        //this grid data is only for infrared radiation
        GMatrix m_erpgriddata[3]; // for x , y and z seperately
        
        double m_antennaPower;  // the power of antenna.
        
        double m_solarPanelPowerDraw; // the power draw of solar panel
        
        GVector m_bias;  // the value for y_bias in nanoNewton, pointing to the +y direction of the BFS
        
        void printSpacecraftModel()
        {
            
            printf("mass: %f\n", m_mass );
            printf("COM:  %f %f %f\n", m_com[0],m_com[1],m_com[2] );
            printf("LRAoffset:  %f %f %f\n", m_offsetLRA[0],m_offsetLRA[1],m_offsetLRA[2] );
            printf("GNSSoffset:  %f %f %f\n", m_offsetGNSS[0],m_offsetGNSS[1],m_offsetGNSS[2] );
            printf("nominalMass: %f\n", m_nomialMass);
            printf("AntennaPower: %f\n", m_antennaPower);
            for(int i = 0 ; i< 2; i++ )
            {
                printf("%f %f %f ", busX[i].normal.x,busX[i].normal.y,busX[i].normal.z);
                //printf();
            }
            
        }
        
        
    };
    
    //a class to store all the spacecraft model and geometry information
    class GSpacecraftModelMgr
    {
        
    public:
        
       static std::map<GString, GSpaceCraftModel> initialiseModel();
       static void initialiseModel(GString spacecraftModelFile);
        
       static GString sensorID2svType( GSensorID svID );
        
        // GString is the type of the spacecraft, for an example, GPSIIR, GPSIIF, GPSIII
       static std::map<GString, GSpaceCraftModel> spacecraftModelStore;
        
    };
    
    
    
}

#endif /* GSpacecraftModel_hpp */
