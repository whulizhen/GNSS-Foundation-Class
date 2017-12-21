#ifndef GFC_ELLIPSOIDMGR_H
#define GFC_ELLIPSOIDMGR_H

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


#include "GEllipsoid.h"
#include <map>

namespace gfc
{
    class GEllipsoidMgr
    {
        
    public:
        
        //构造函数
        GEllipsoidMgr();
        //析构函数
        ~GEllipsoidMgr();
        //初始化函数
        static std::map< GString,  GEllipsoid > InitializeEllipsoidMgr();
        static GEllipsoid GetEllipsoid( GString ellipsoidName );
        double    A(GString ellipsoidName);
        double    A_km(GString ellipsoidName);
        double    GM(GString  ellipsoidName);
        double    F(GString ellipsoidName);
        double    E2(GString ellipsoidName);
        double    E(GString ellipsoidName);
        double    B(GString ellipsoidName);
        double    AV(GString ellipsoidName);
        double    J2(GString ellipsoidName);
        void dump(std::ostream& s)  ;
        GString getClassName();
        
    private:
        static std::map<GString ,GEllipsoid > mEllipsoidStorage;  //存储所有的椭球相关信息
        
    };
    
}

#endif