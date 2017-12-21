
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


#include "GEllipsoidMgr.h"

namespace gfc
{
    //调用常量的初始化函数；必须在类外调用
    std::map< GString, GEllipsoid > gfc::GEllipsoidMgr::mEllipsoidStorage = gfc::GEllipsoidMgr::InitializeEllipsoidMgr();
    
    GEllipsoidMgr::GEllipsoidMgr() {}
    
    GEllipsoidMgr::~GEllipsoidMgr() {}
    
    std::map< GString,  GEllipsoid > GEllipsoidMgr::InitializeEllipsoidMgr()
    {
        std::map< GString ,GEllipsoid > tmpEllipsoidStorage;
        //所有椭球的J2项(正常化二阶带球谐系数)均采用同一值；均为总地球椭球，全球平差计算椭球定向
        
        GEllipsoid ellipsoid1("WGS84",6378137.0,1.0/298.257223563,7.2921151467E-5,3.9860050E14,-4.8416685E-4);
        GEllipsoid ellipsoid2("CGCS2000",6378137.0,1.0/298.257222101,7.2921150e-5,3.986004418E14,1.0826257E-3);
        GEllipsoid ellipsoid3("PZ90",6378136.0,1.0/298.257222101,7.292115E-5,3.9860044E14,1.0826257E-3);
        
        GEllipsoid ellipsoid4("CERES",6408137.0,1.0/298.25680814310491,7.292115E-5,3.9860044E14,1.0826257E-3);
        
        GEllipsoid ellipsoid5("CERES_SPHERE",6408137.0,0.0,7.292115E-5,3.9860044E14,1.0826257E-3);
        
        GEllipsoid ellipsoid6("AVERAGE",6371000.0,0.0,7.2921151467E-5,3.9860050E14,-4.8416685E-4);
        
        
        tmpEllipsoidStorage["WGS84"] = ellipsoid1;
        tmpEllipsoidStorage["CGCS2000"] = ellipsoid2;
        tmpEllipsoidStorage["PZ90"] = ellipsoid3;
        tmpEllipsoidStorage["CERES"] = ellipsoid4;
        tmpEllipsoidStorage["CERES_SPHERE"] = ellipsoid5;
        tmpEllipsoidStorage["AVERAGE"] = ellipsoid6;
        
        return tmpEllipsoidStorage;
        
    }
    
    GEllipsoid GEllipsoidMgr::GetEllipsoid( GString ellipsoidName )
    {
        //需要对椭球的合法性进行检查
        return mEllipsoidStorage[ellipsoidName.upperCase()];
    }
    
    double GEllipsoidMgr::A(GString ellipsoidName)
    {
        return GetEllipsoid(ellipsoidName.upperCase()).A();
    }
    
    double GEllipsoidMgr::A_km(GString ellipsoidName)
    {
        return GetEllipsoid(ellipsoidName.upperCase()).A_km();
    }
    
    double GEllipsoidMgr::F(GString ellipsoidName)
    {
        return GetEllipsoid(ellipsoidName.upperCase()).F();
    }
    
    double GEllipsoidMgr::E(GString ellipsoidName)
    {
        return GetEllipsoid( ellipsoidName.upperCase()).E();
    }
    
    double GEllipsoidMgr::E2(GString ellipsoidName)
    {
        return GetEllipsoid( ellipsoidName.upperCase() ).E2();
    }
    
    double GEllipsoidMgr::B(GString ellipsoidName)
    {
        return GetEllipsoid( ellipsoidName.upperCase() ).B();
    }
    
    double GEllipsoidMgr::AV(GString ellipsoidName)
    {
        return GetEllipsoid( ellipsoidName.upperCase() ).AV();
    }
    
    double GEllipsoidMgr::GM(GString ellipsoidName)
    {
        return GetEllipsoid( ellipsoidName.upperCase() ).GM();
    }
    
    double GEllipsoidMgr::J2(GString ellipsoidName)
    {
        return GetEllipsoid( ellipsoidName.upperCase() ).J2();
    }
    
    GString GEllipsoidMgr::getClassName()
    {
        return "EllipsoidMgr";
    }
    
    void GEllipsoidMgr::dump(std::ostream& s)
    {
        s<<"ClassName:       "<<getClassName()<<std::endl;
        s<<"Total Ellipsoid Number: "<< mEllipsoidStorage.size()<<std::endl;
        //遍历所有的变量
        std::map<GString, GEllipsoid >::iterator myit = mEllipsoidStorage.begin();
        for(;myit!=mEllipsoidStorage.end();myit++)
        {
            myit->second.dump(s);
        }
    }
    
}
