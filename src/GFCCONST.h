#ifndef GFC_GNSSCONST_H
#define GFC_GNSSCONST_H

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


/*
 
 GNSS数据处理软件中的常量的存储
 
 其中的成员变量全部定义为static，方便调用
 
 该类中的常量全部在初始化时从文件中读取，
 未来替换掉代码中常量的定义文件constant.h
 
 创建：李桢 2015年2月3日
 
 */


/*
 每颗卫星的三维模型定义结构体
 目前初步想法是考虑用DTN来建立飞行器的三维模型
 包括飞行器的表面材料，飞行器发射日期，预计寿命等
 每颗飞行器的天线型号，以及星体本体坐标系的定义
 飞行器的实时姿态
 
 */

#include "GString.h"
#include <map>
#include <iomanip>
#include <math.h>

namespace gfc
{
    
    /*坐标参考系统的定义
     WGS-84坐标系是G1150周的第一天零时
     CGCS200坐标系的参考历元是？？？
     ITRS的Z轴指向为BIH1984.0
     但是ITRS采用何种参考椭球????
     
     */
    //	struct RefSys
    //	{
    //		GString  name;  //坐标系的名称
    //		Ellipsoid ellipsoid;  //每个坐标系都有一个参考椭球
    //		double T0;    // 每个坐标系都有一个参考起始历元时刻(MJD)，用于确定坐标轴指向
    //	};
    
    
    ///*GNSS卫星系统的定义*/
    //struct GNSSSYS
    //{
    //	int sys;    //系统标记，目前是GPS为0x01,BDS为0x02,GLONASS为0x04, GALLILUE为0x08
    //	std::string sysname;  //卫星系统的名称
    //	char signalType;  //信号类型，0为CDMA，1为FDMA，2为TDMA，目前只有GLONASS为FDMA，星间链路数据采用TDMA方式
    //	int satnum;  //该卫星系统的卫星颗数
    //	double T0;   //该卫星系统时间的起始时刻(MJD表示)
    //	RefSys cs;  //该卫星系统采用的坐标系统
    //	std::vector<SpaceCraftModel>  spaceCraft;  //该卫星系统所采用的所有的飞行器信息
    //};
    
    
    NEW_GEXCEPTION_CLASS( constantUnexist, gfc::GException );
    
    class GFCCONST
    {
        
    public:
        
        /*
         静态函数，用于对constantValue进行初始化
         */
        static std::map< GString,  long double > Initializer();
        static void RegByName(GString variableName,long double variableValue);
        static long double GetByName(GString variableName);
        
        static void UnregByName(GString variableName);
        static void dump( std::ostream& s ) ;
        
    private:
        
        // private construction function means this class can not be hesitated!
       	GFCCONST(void);
        virtual ~GFCCONST(void);
        
        static long double GetByName_internal( GString variableName,std::map< GString,  long double >& myconstantValue);
        
        static std::map< GString,  long double > constantValue;  //所有变量的值以及变量名均存储在这里
        
        //GString constFileName;  //常量配置文件(设计为文本文件)
        
        //std::map<Ellipsoid>   myEllipsoid;                        //所有的椭球信息
        
    };
    
    
    
    
    
    
}  // end of namespace


/*
 
 定义宏CONSTANT用于获取常量的值
 
 可以仿造LogStream类进行编写
 
 注意：常量名均为大写;
 
 调用格式如下：
 CONSTANT("CLIHT")
 CONSTANT("PI")
 
 */

#define GCONST(constantName) \
GFCCONST::GetByName(constantName)

#endif


