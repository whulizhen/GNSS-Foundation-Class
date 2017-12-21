//
//  GEphemeris.cpp
//  GFC
//
//  Created by lizhen on 05/06/2016.
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

#include "GEphemeris.hpp"

namespace gfc
{
    //GPreciseEphemeris::GPreciseEphemeris(GTime epoch,GString info)
   // {
//        m_x = 0.0; m_dx =0.0;m_y=0.0; m_dy=0.0; m_z=0.0;m_dz=0.0;m_c=0.0;m_dc=0.0;
//        m_u =0.0; m_v =0.0; m_w=0.0; m_du=0.0; m_dv=0.0;m_dw=0.0;
//        
//        m_epoch = epoch;
//        if( info[0]=='P' )  //position info
//        {
//            int tag = 0;
//            
//            std::vector<GString> splitstr = info.split();
//            
//            if( info[1] == 'G' || info[1] == ' ' )  //gps
//            {
//                int prn  = -1;
//                //GString tmp =splitstr[0].substr(2,2);
//                if( splitstr[0] == "P" )
//                {
//                    prn = static_cast<int>(splitstr[1].asINT());
//                    tag = 1;
//                }
//                else
//                {
//                    prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                }
//                
//                m_scID = GSensorID("ssGPS",prn);
//            }
//            else if( info[1]=='R' ) //glonass
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssGLO",prn);
//            }
//            else if( info[1]=='C' ) //beidou
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssBDS",prn);
//            }
//            else if( info[1]=='E' ) //galileo
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssGAL",prn);
//            }
//            else if( info[1] == 'L')  // LAGEOS geodetic satellite
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssLAGEOS",prn);
//            }
//            
//            m_x = splitstr[1+tag].asDOUBLE();
//            m_dx = GCONST("MAXDEV");
//            m_y = splitstr[2+tag].asDOUBLE();
//            m_dy = GCONST("MAXDEV");
//            m_z = splitstr[3+tag].asDOUBLE();
//            m_dz = GCONST("MAXDEV");
//            if(splitstr.size() > 3 + tag)
//            {
//                m_c = splitstr[4+tag].asDOUBLE(); // unit:meter
//                m_dc = GCONST("MAXDEV");
//            }
//            
//        }
//        else if(info[0]=='V')  //velocity info
//        {
//            int tag = 0;
//            std::vector<GString> splitstr = info.split();
//            
//            if( info[1] == 'G' || info[1] == ' ' )  //gps
//            {
//                int prn  = -1;
//                //GString tmp =splitstr[0].substr(2,2);
//                if( splitstr[0] == "V" )
//                {
//                    prn = static_cast<int>(splitstr[1].asINT());
//                    tag = 1;
//                }
//                else
//                {
//                    prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                }
//                
//                m_scID = GSensorID("ssGPS",prn);
//            }
//            else if( info[1]=='R' ) //glonass
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssGLO",prn);
//            }
//            else if( info[1]=='C' ) //beidou
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssBDS",prn);
//            }
//            else if( info[1]=='E' ) //galileo
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssGAL",prn);
//            }
//            else if( info[1] == 'L')  // LAGEOS geodetic satellite
//            {
//                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
//                m_scID = GSensorID("ssLAGEOS",prn);
//            }
//            
//            m_u = splitstr[1+tag].asDOUBLE();
//            m_du = GCONST("MAXDEV");
//            m_v = splitstr[2+tag].asDOUBLE();
//            m_dv = GCONST("MAXDEV");
//            m_w = splitstr[3+tag].asDOUBLE();
//            m_dw = GCONST("MAXDEV");
//            
//            //m_c = splitstr[4+tag].asDOUBLE(); // unit:meter
//            //m_dc = GCONST("MAXDEV");
//            
//        }
//        
//        m_isOK = true;
        
 //   }
    
    void gfc::GPreciseEphemeris::setData( gfc::GString info)
    {
//        m_x = 0.0; m_dx =0.0;m_y=0.0; m_dy=0.0; m_z=0.0;m_dz=0.0;m_c=0.0;m_dc=0.0;
//        m_u =0.0; m_v =0.0; m_w=0.0; m_du=0.0; m_dv=0.0;m_dw=0.0;
        
        if( info[0]=='P' )  //position info
        {
            int tag = 0;
            
            std::vector<GString> splitstr = info.split();
            
            if( info[1] == 'G' || info[1] == ' ' )  //gps
            {
                int prn  = -1;
                //GString tmp =splitstr[0].substr(2,2);
                if( splitstr[0] == "P" )
                {
                    prn = static_cast<int>(splitstr[1].asINT());
                    tag = 1;
                }
                else
                {
                    prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                }
                
                m_scID = GSensorID("ssGPS",prn);
            }
            else if( info[1]=='R' ) //glonass
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssGLO",prn);
            }
            else if( info[1]=='C' ) //beidou
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssBDS",prn);
            }
            else if( info[1]=='E' ) //galileo
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssGAL",prn);
            }
            else if( info[1] == 'L')  // LAGEOS geodetic satellite
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssLAGEOS",prn);
            }
            else  // unrecognised satellite system
            {
                return ;
            }
            m_x = splitstr[1+tag].asDOUBLE();
            m_dx = GCONST("MAXDEV");
            m_y = splitstr[2+tag].asDOUBLE();
            m_dy = GCONST("MAXDEV");
            m_z = splitstr[3+tag].asDOUBLE();
            m_dz = GCONST("MAXDEV");
            if(splitstr.size() > 4 + tag)
            {
                m_c = splitstr[4+tag].asDOUBLE(); // unit:meter
                m_dc = GCONST("MAXDEV");
            }
            
            m_isOK = true;
            
        }
        else if(info[0]=='V')  //velocity info
        {
            int tag = 0;
            std::vector<GString> splitstr = info.split();
            double scale = 1.0;
            if( info[1] == 'G' || info[1] == ' ' )  //gps
            {
                int prn  = -1;
                //GString tmp =splitstr[0].substr(2,2);
                if( splitstr[0] == "V" )
                {
                    prn = static_cast<int>(splitstr[1].asINT());
                    tag = 1;
                }
                else
                {
                    prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                }
                
                m_scID = GSensorID("ssGPS",prn);
            }
            else if( info[1]=='R' ) //glonass
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssGLO",prn);
            }
            else if( info[1]=='C' ) //beidou
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssBDS",prn);
            }
            else if( info[1]=='E' ) //galileo
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssGAL",prn);
            }
            else if( info[1] == 'L')  // LAGEOS geodetic satellite
            {
                int prn = static_cast<int>(splitstr[0].substr(2,2).asINT()) ;
                m_scID = GSensorID("ssLAGEOS",prn);
                scale = 1.0/10000.0;
            }
            else // unrecognised satellite system
            {
                return ;
            }
            
            m_u = splitstr[1+tag].asDOUBLE()*scale;
            m_du = GCONST("MAXDEV");
            m_v = splitstr[2+tag].asDOUBLE()*scale;
            m_dv = GCONST("MAXDEV");
            m_w = splitstr[3+tag].asDOUBLE()*scale;
            m_dw = GCONST("MAXDEV");
            
            //m_c = splitstr[4+tag].asDOUBLE(); // unit:meter
            //m_dc = GCONST("MAXDEV");
            m_isOK = true;
        }
        
        
        
    }
    
    
    GPreciseEphemeris& GPreciseEphemeris::operator= (const GPreciseEphemeris& right)
    {
        this->m_x = right.m_x;this->m_y = right.m_y;this->m_z = right.m_z;this->m_c = right.m_c;
        this->m_dx = right.m_dx;this->m_dy = right.m_dy;this->m_dz = right.m_dz;this->m_dc = right.m_dc;
        this->m_u = right.m_u;this->m_v = right.m_v;this->m_w = right.m_w;
        this->m_du = right.m_du;this->m_dv = right.m_dv;this->m_dw = right.m_dw;
        this->m_epoch = right.m_epoch;
        this->m_scID = right.m_scID;
        this->m_isOK = right.m_isOK;
        
        return *this;
    }
    
    GPreciseEphemeris::GPreciseEphemeris( const GPreciseEphemeris& right)   // copy constructor
    {
        m_epoch = right.m_epoch;
        m_scID = right.m_scID;
        this->m_x = right.m_x;this->m_y = right.m_y;this->m_z = right.m_z;this->m_c = right.m_c;
        this->m_dx = right.m_dx;this->m_dy = right.m_dy;this->m_dz = right.m_dz;this->m_dc = right.m_dc;
        this->m_u = right.m_u;this->m_v = right.m_v;this->m_w = right.m_w;
        this->m_du = right.m_du;this->m_dv = right.m_dv;this->m_dw = right.m_dw;
        this->m_isOK = right.m_isOK;
    }
    
    
} // end of namespace
