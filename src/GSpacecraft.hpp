//
//  GSpacecraft.hpp
//  GFC
//
//  Created by lizhen on 16/4/14.
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

#ifndef GSpacecraft_hpp
#define GSpacecraft_hpp

#include <stdio.h>
#include <fstream>

#include "GEphemeris.hpp"
#include "GMotionState.hpp"

#include "GSpacecraftModel.hpp"

namespace  gfc
{
    
    
    
    
    
    // the definition of spacecraft
    // in the future, this class should be inherited from class Sensor
    class GSpaceCraft
    {
        
    public:
        
        GSpaceCraft () {}
        
        GSpaceCraft (GSensorID id ) { m_scID = id; m_bodyModel = NULL;}
        
        GSpaceCraft& operator= (const GSpaceCraft& right);  //= override
        
        GSpaceCraft( const GSpaceCraft& right);   // copy constructor
        
        virtual ~GSpaceCraft() {}
        
        GString getSpaceCraftName()
        {
            return m_scID.getIDString();
        }
        
        //get the pointer of the m_state for the operation of the current state
        GMotionState* getStatePointer()
        {
            return &m_state;
        }
        
        GSpaceCraftModel* getSpaceCraftGemotry()
        {
            return m_bodyModel;
        }
        
        void pushPreciseEphemeris(GPreciseEphemeris pehp);
        
        // this function must be rewritten by the child class
        void setSpaceCraftID( GSensorID id ){ m_scID = id; };
        
        GSensorID getSpaceCraftID() {return m_scID;}
        
        
        int getEphLocation(GTime& ttag, const int& nhalf,
                                        std::map<GTime,GPreciseEphemeris>::const_iterator & it1,
                                        std::map<GTime,GPreciseEphemeris>::const_iterator& it2,
                                        bool& exactMatch,
                                        bool exactReturn
                                    );
        
        GPreciseEphemeris getEphValue1(GTime ct);
        
        GPreciseEphemeris getEphValue( int index );
        
        GPreciseEphemeris getEphValue( GTime ct );
        
        GPreciseEphemeris getEphValue_test(GTime ct);
        
        int getNumofPreEph();
        
        void bindSpacecraftModel(GSensorID);
        
        //get the corrections for LRA or GNSS Antenna PCO
        // the final result is in ECEF because all the measurements are in ECEF
        void getOffsetCorrection( int type, GVector& targetecef );
        
    private:
        
        // the broadcast ephemeris storage for this satellite
        std::vector< GBroadcastEphemeris > m_bephStore;
        
        std::map< GTime, GPreciseEphemeris > m_pephStore;
        
        GMotionState m_state;  // the state of the spacecraft in current time epoch, including pos vel and attitude
        
        GSensorID m_scID;  // spacecraft ID
        
        GSpaceCraftModel* m_bodyModel;  // the geometry of the current spacecraft
        
    };
    
    
    /*class for the management of all the spacecraft */
    class GSpaceCraftMgr
    {
        
        //make the constructor private forcely
    private:
        
        static int GPSMAXPRN ;
        
        static int BDSMAXPRN ;
        
        static int GLSMAXPRN ;
        
        static int GALMAXPRN ;
        
    public:
        
        GSpaceCraftMgr()  {}
        
        ~GSpaceCraftMgr() {}
        
        static std::map< GString,GSpaceCraft> gSpacecraft;
        
        static void loadPreciseEphemeris( GString ephemerisFile);
        
        // load the broadcast ephemeris
        static void loadBroadcastEphemeris( GString ephemerisFile );
        
        static void loadUCLEphemeris(GString ephemerisFile);
        
        static void loadGRACEEphemeris(GString ephemerisFile);
        
        static void Initializer();
        
    };
    
    
    
    
    
    
    
}  // end of namespace gfc






#endif /* GSpacecraft_hpp */
