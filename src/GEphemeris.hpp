//
//  GEphemeris.hpp
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

#ifndef GEphemeris_hpp
#define GEphemeris_hpp

#include <stdio.h>
#include "GTime.h"
#include "GSensor.h"
#include "GVector.hpp"

namespace gfc
{
    
    // the broad ephemeris for GNSS, except for GLONASS
    class GBroadcastEphemeris
    {
        
    public:
        
        GBroadcastEphemeris() {}
        
        virtual ~GBroadcastEphemeris() {}
        
    private:
        
        GSensorID m_scID;
        
        double m_sqrtA; //
        
    };
    
    
    class GPreciseEphemeris
    {
        
    public:
        
        GPreciseEphemeris()
        {
            m_x =0.0; m_y = 0.0; m_z =0.0; m_c =0.0;
            m_u =0.0; m_v =0.0; m_w =0.0;
            m_dx =0.0; m_dy =0.0; m_dz =0.0; m_dc =0.0;
            m_du =0.0; m_dv =0.0; m_dw =0.0;
            m_isOK = false;
        }
        
        //GPreciseEphemeris(GTime epoch,GString info);
        
        void setData(GString info);
        
        GTime getEpoch() {return m_epoch;}
        
        bool isOK()
        {
            return m_isOK;
        }
        
        void setOK(bool status)
        {
            m_isOK = status;
        }
        
        void getPV( double* pv )
        {
            pv[0] = m_x; pv[1] = m_y; pv[2] = m_z;
            pv[3] = m_u; pv[4] = m_v; pv[5] = m_w;
        }
        
        void getPV(GVector& p, GVector& v)
        {
            p.set(m_x, m_y, m_z);
            v.set(m_u, m_v, m_w);
        }
        
        //double getPX() {return m_x;}
        
        double getClock()
        {
            return m_c;
        }
        
        GSensorID getSensorID(){return m_scID;}
        
        void setEpoch(GTime t) {m_epoch = t;}
        void setPX(double value) {m_x = value;}
        void setPY(double value) {m_y = value;}
        void setPZ(double value) {m_z = value;}
        void setVX(double value) {m_u = value;}
        void setVY(double value) {m_v = value;}
        void setVZ(double value) {m_w = value;}
        void setC(double value) {m_c = value;}
        void setPos(double *v)
        {
            m_x = v[0];
            m_y = v[1];
            m_z = v[2];
        }
        
        void setVel(double *v)
        {
            m_u = v[0];
            m_v = v[1];
            m_w = v[2];
        }
        void setSCID(GSensorID scid)
        {
            m_scID = scid;
        }
        virtual ~ GPreciseEphemeris() {}
        
        GPreciseEphemeris& operator= (const GPreciseEphemeris& right);  //= override
        
        GPreciseEphemeris( const GPreciseEphemeris& right);   // copy constructor
        
        GTime m_epoch;
        GSensorID m_scID;
        double m_x; // x coordinate; km
        double m_dx; // sigma for x
        double m_y; // y coordinate
        double m_dy; // sigma for y
        double m_z; // z coordinate
        double m_dz; // sigma for z
        double m_u; // velocity x
        double m_du; // sigma for u
        double m_v;  // velocity y
        double m_dv; // sigma for v
        double m_w; // velocity z
        double m_dw; // sigma for w
        double m_c;  // clock difference at the moment , unit:meter
        double m_dc; // the sigma for the clock difference
        bool   m_isOK;  // whether this ephemeris is available ?
    };
    
    
    
} // end of namespace




#endif /* GEphemeris_hpp */
