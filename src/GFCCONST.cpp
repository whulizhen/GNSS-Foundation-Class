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

#include "GFCCONST.h"

namespace gfc
{
    
    //调用常量的初始化函数；必须在类外调用
    std::map< GString,  long double > gfc::GFCCONST::constantValue = gfc::GFCCONST::Initializer();
    
    /*
     静态函数的实现(初始化)
     所有的时间单位为：纳秒
     所有的距离单位为：米
     速度单位为：m/s
     */
    std::map< GString,  long double > GFCCONST::Initializer()
    {
        //printf("constant initializer\n");
        std::map< GString,long double>  tmpconsValue;  //map 内部用红黑树实现，查找效率是logN
        tmpconsValue["CLIGHT"]   = 299792458.0;
        tmpconsValue["GM"]       = 398600.4415;
        tmpconsValue["MAXDEV"]   = 999.9;   // the max standard deviation for all the data, you should square it if you want to get variance
        tmpconsValue["PI"]       = 3.14159265358979323846264338328;
        tmpconsValue["AU"]       = 149597870.7;  //天文单位, in Km
        tmpconsValue["MJD0"]     = 2400000.5;  //MJD0的儒略日
        tmpconsValue["UTC0"]     = 36934.0;      //UTC0零时刻的儒略日
        tmpconsValue["SECPMIN"]  = 60.0;  //sec per minute = 60
        tmpconsValue["SECPHOR"]  = 3600.0;  //sec per hour = 3600
        tmpconsValue["SECPDAY"]  = 86400.0;  //sec per day = 86400
        tmpconsValue["SECPWEK"]  = 604800.0;  //sec per week = 604800
        tmpconsValue["MINPHOR"]  = 60.0;  //sec per minute = 60
        tmpconsValue["MINPDAY"]  = 1440.0;  //sec per minute = 60
        tmpconsValue["MINPWEK"]  = 10080.0;  //sec per minute = 60
        tmpconsValue["HORPDAY"]  = 24.0;  //sec per minute = 60
        tmpconsValue["HORPWEK"]  = 168.0;  //sec per minute = 60
        tmpconsValue["DAYPWEK"]  = 7.0;  //sec per minute = 60
        tmpconsValue["SECPDEG"]  = 3600.0;
        tmpconsValue["NANO"]     = 1.0E9;                      //单位：纳 是10^9
        tmpconsValue["2PI"]      = 2.0*GetByName_internal("PI",tmpconsValue);
        tmpconsValue["IONCOEF"]  = 40.28E16/GetByName_internal("CLIGHT",tmpconsValue)/GetByName_internal("CLIGHT",tmpconsValue);
        tmpconsValue["D2R"]      = GetByName_internal("PI",tmpconsValue)/180.0;
        tmpconsValue["R2D"]      = 180.0/GetByName_internal("PI",tmpconsValue);
        tmpconsValue["SIN5"]     = sin(5.0*GetByName_internal("D2R",tmpconsValue));
        tmpconsValue["COS5"]     = cos(5.0*GetByName_internal("D2R",tmpconsValue));
        tmpconsValue["AS2R"]  = GetByName_internal("PI",tmpconsValue)/180.0/3600.0;  //角秒转换为弧度
        
        //constantValue.insert(std::map<std::string,long double >::value_type("3PI",9.0));
        
        return  tmpconsValue;
    }
    
    /*静态函数：取出变量的值*/
    long double GFCCONST::GetByName(GString variableName)
    {
        if( constantValue.count(variableName)>0 )
        {
            return constantValue[variableName];
        }
        else
        {
            constantUnexist e("constant unexist.",1001,gfc::GException::Severity::unrecoverable);
            GFC_THROW(e);
        }
    }
    
    /*静态函数：取出变量的值*/
    long double GFCCONST::GetByName_internal(GString variableName,std::map< GString,  long double >& myconstantValue)
    {
        if( myconstantValue.count(variableName)>0 )
        {
            return myconstantValue[variableName];
        }
        else
        {
            constantUnexist e("constant unexist.",1001,gfc::GException::Severity::unrecoverable);
            GFC_THROW(e);
        }
    }
    
    
    void GFCCONST::RegByName(GString variableName,long double variableValue)
    {
        if( constantValue.count(variableName) == 0 )
        {
            constantValue.insert(std::map<GString,long double >::value_type(variableName,variableValue));
        }
    }
    
    /*删除一个常量*/
    void GFCCONST::UnregByName(GString variableName)
    {
        if( constantValue.count(variableName)>0 )
        {
            constantValue.erase(variableName);
        }
        else
        {
            constantUnexist e("Delete ConstantValue: constant Unexist.",1002,gfc::GException::Severity::unrecoverable);
            GFC_THROW(e);
        }				
    }
    
    void GFCCONST::dump( std::ostream& s )
    {
        s<<"ClassName:    "<<"GFCCONST"<<std::endl;
        
        //遍历所有的变量
        std::map<GString,long double>::iterator myit = constantValue.begin();
        for(;myit!=constantValue.end(); myit++)
        {
            s<<"Name:  "<<myit->first<<"\t"<<"Value:  "<<GString(myit->second,6)<<std::endl;
        }		
    }
    
    GFCCONST::GFCCONST(void)
    {
        
    }
    
    GFCCONST::~GFCCONST(void)
    {
        
    }
    
    
    
}  //end of namespace



