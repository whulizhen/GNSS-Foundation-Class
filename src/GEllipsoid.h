#ifndef GFC_ELLIPSOID_H
#define GFC_ELLIPSOID_H


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
 
 椭球基准类，可被用于继承（根据不同的需求）
 椭球的定义，用于定义大地坐标系统
 用椭球来模拟地球的实际情况
 
 只有地球坐标系才需要参考椭球参数
 天球坐标系不需要椭球参数（惯性坐标系）
 
 */

#include "GException.h"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <sstream>

#include "GString.h"
#include "GFCCONST.h"
#include "GVector.hpp"

using namespace std;

namespace gfc
{
    NEW_GEXCEPTION_CLASS(InvalidEllipsoid, gfc::GException);
    
    class GEllipsoid
    {
        
    public:
        
        /*椭球构造函数*/
        GEllipsoid()
        : name( "UNKNOWN" ), a(0.0), flattening(0.0), eccentricity(0.0), gm(0.0),\
        angVelocity(0.0), j2(0.0)
        {};
        
        /*参数构造函数*/
        GEllipsoid(  GString pmName, long double pmA,long double pmF,
                  long double pmAV, long double pmGM,  long double pmJ2
                  )
        throw ( InvalidEllipsoid )
        {
            
            if( pmName.empty() )
            {
                InvalidEllipsoid e("Invalid Ellipsoid:construction function",1101,gfc::GException::Severity::unrecoverable);
                e.addLocation(FILE_LOCATION);
                std::cout<<e<<endl;
                GFC_THROW(e);
                }
                else
                {
                    name = pmName.upperCase(); //椭球名均为大写字母
                    a = pmA;
                    flattening = pmF;
                    double b = a*(1.0 - flattening);
                    eccentricity = sqrt(a*a - b*b )/a;
                    angVelocity = pmAV;
                    gm = pmGM;
                    j2 = pmJ2;
                }
                };
                
                ///拷贝构造函数
                /// Copy constructor
                GEllipsoid(const GEllipsoid& s)
                : name(s.name), a(s.a) ,flattening(s.flattening),\
                eccentricity(s.eccentricity), gm(s.gm),j2(s.j2),\
                angVelocity(s.angVelocity)
                {};
                
                /// Assignment operator
                virtual GEllipsoid operator=(const GEllipsoid& right)
                {
                    if ( this == &right ) return (*this);
                    
                    (*this).name = right.name;(*this).a = right.a;
                    (*this).angVelocity = right.angVelocity;(*this).eccentricity = right.eccentricity;
                    (*this).flattening = right.flattening;(*this).gm = right.gm;(*this).j2 = right.j2;
                    return *this;
                };
                
                
                /// Equality requires all fields to be the same
                virtual bool operator==(const GEllipsoid& right) const
                { return name==right.name; };
                
                /// Inequality operator
                virtual bool operator!=(const GEllipsoid& right) const
                { return !(operator==(right)); };
                
                /// Greater than operator
                virtual bool operator>(const GEllipsoid& right) const
                {  return (!operator<(right) && !operator==(right)); };
                
                /// This ordering is somewhat arbitrary but is required to be able
                /// to use an TypeID as an index to a std::map. If an application
                /// needs some other ordering, inherit and override this function.
                virtual bool operator<(const GEllipsoid& right) const
                { return name < right.name; };
                
                /// Greater than or equal operator
                virtual bool operator>=(const GEllipsoid& right) const
                { return !(operator<(right)); };
                
                /// Less than or equal operator
                virtual bool operator<=(const GEllipsoid& right) const
                { return (operator<(right) || operator==(right)); };
                
                
                /*椭球析构函数*/
                ~GEllipsoid()  {};
                
                /// Return a string that will identify the derived class
                virtual GString getClassName(void) const
                { return GString("Ellpsoid"); }
                
                GString Name()  throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid: Invalid Ellipsoid Name!",1102,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return name;
                    }
                };
                
                //椭球的长半轴
                long double A()		throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid: A()",1103,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return a;
                    }
                };
                
                
                long double A_km()  throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid: A_km()",1104,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return a/1000.0;
                    }
                };
                
                //椭球扁率
                long double F()		throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid: F()",1105,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return flattening;
                    }
                };
                
                //椭球离心率
                long double E()		throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid: E()",1106,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return eccentricity;
                    }
                };
                
                long double E2()    throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid: E2()",1107,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return eccentricity*eccentricity;
                    }
                };
                
                long double B()     throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid:B()",1108,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return 0.0;
                    }
                };
                
                long double GM()    throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid:GM()",1109,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return gm;
                    }
                };
                
                //椭球转动角速度
                long double AV()    throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid:AV()",1110,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return angVelocity;
                    }
                };
                
                long double J2()    throw(InvalidEllipsoid)
                {
                    if( name == "UNKNOWN" || name.empty())
                    {
                        InvalidEllipsoid e("Ellipsoid Invalid:J2()",1111,gfc::GException::Severity::unrecoverable);
                        e.addLocation(FILE_LOCATION);
                        std::cout<<e<<endl;
                        GFC_THROW(e);
                    }
                    else
                    {
                        return j2;
                    }
                };
                
                // transfer from blh to xyz with the current elliposid parameters
                void BLH2XYZ(double blh[3], double* xyz )
                {
                    double e = eccentricity;   //第一偏心率
                    double sinb = sin(blh[0]);
                    double sinl = sin(blh[1]);
                    double cosb = cos(blh[0]);
                    double cosl = cos(blh[1]);
                    double N = a/sqrt(1-pow(e*sinb,2));
                    xyz[0] = (N + blh[2])*cosb*cosl;
                    xyz[1] = (N + blh[2])*cosb*sinl;
                    xyz[2] = (blh[2]+ N*(1-e*e))*sinb;
                    
                }
                // transfer from blh to xyz with the current elliposid parameters
                void BLH2XYZ(GVector blh, GVector& xyz )
                {
                    double e = eccentricity;   //第一偏心率
                    double sinb = sin(blh.x);
                    double sinl = sin(blh.y);
                    double cosb = cos(blh.x);
                    double cosl = cos(blh.y);
                    double N = a/sqrt(1-pow(e*sinb,2.0));
                    xyz.x = (N + blh.z)*cosb*cosl;
                    xyz.y = (N + blh.z)*cosb*sinl;
                    xyz.z = (blh.z+ N*(1-e*e))*sinb;
                    
                }

                void XYZ2BLH(double xyz[3], double* blh)
                {
                    double PI = GCONST("PI");
                    double e2= eccentricity*eccentricity,z,zk,v=a,sinp;
                    double r2 = xyz[0]*xyz[0] + xyz[1]*xyz[1];
                    for (z=xyz[2],zk=0.0;fabs(z-zk)>=1E-8;)
                    {
                        zk=z;
                        sinp=z/sqrt(r2+z*z);
                        v=a/sqrt(1.0-e2*sinp*sinp);
                        z=xyz[2]+v*e2*sinp;
                    }
                    
                    blh[0]=r2>1E-12?atan(z/sqrt(r2)):(xyz[2]>0.0?PI/2.0:-PI/2.0);
                    blh[1]=r2>1E-12?atan2(xyz[1],xyz[0]):0.0;
                    blh[2]=sqrt(r2+z*z)-v;
                    
                }
                
                void XYZ2BLH(GVector xyz, GVector& blh)
                {
                    double PI = GCONST("PI");
                    double e2= eccentricity*eccentricity,z,zk,v=a,sinp;
                    double r2 = xyz.x*xyz.x + xyz.y*xyz.y;
                    for (z=xyz.z,zk=0.0;fabs(z-zk)>=1E-6;)
                    {
                        zk=z;
                        sinp=z/sqrt(r2+z*z);
                        v=a/sqrt(1.0-e2*sinp*sinp);
                        z=xyz.z+v*e2*sinp;
                    }
                    
                    blh.x=r2>1E-12?atan(z/sqrt(r2)):(xyz.z>0.0?PI/2.0:-PI/2.0);
                    blh.y=r2>1E-12?atan2(xyz.y,xyz.x):0.0;
                    blh.z=sqrt(r2+z*z)-v;
                    
                }
                
                void XYZ2NEU( double* refXYZ, double* XYZ, double* NEU)
                {
                    double tmp[3] = {0.0}, blh[3] = {0.0};
                    for(int i = 0 ; i < 3; i++ )
                    {
                        tmp[i] = XYZ[i] - refXYZ[i];
                    }
                    
                    XYZ2BLH(refXYZ, blh);
                    
                    double sinl = sin(blh[1]);
                    double cosl = cos(blh[1]);
                    double cosb = cos(blh[0]);
                    double sinb = sin(blh[0]);
                    double m[9] = {0.0};//转换矩阵
                    m[0] = -sinb*cosl;m[1] = -sinb*sinl;m[2] = cosb;   //N
                    m[3] = -sinl;     m[4] = cosl;      m[5] = 0.0;    //E
                    m[6] = cosb*cosl; m[7] = cosb*sinl; m[8] = sinb;   //U
                    
                    NEU[0] = m[0]*tmp[0] + m[1]*tmp[1] + m[2]*tmp[2];
                    NEU[1] = m[3]*tmp[0] + m[4]*tmp[1] + m[5]*tmp[2];
                    NEU[2] = m[6]*tmp[0] + m[7]*tmp[1] + m[8]*tmp[2];
                    
                }
                
                void XYZ2NEU( GVector refXYZ, GVector XYZ, GVector& NEU)
                {
                    GVector tmp , blh;
                    tmp = XYZ - refXYZ;
                    
                    XYZ2BLH(refXYZ, blh);
                    
                    double sinl = sin(blh.y);
                    double cosl = cos(blh.y);
                    double cosb = cos(blh.x);
                    double sinb = sin(blh.x);
                    double m[9] = {0.0};//转换矩阵
                    m[0] = -sinb*cosl;m[1] = -sinb*sinl;m[2] = cosb;   //N
                    m[3] = -sinl;     m[4] = cosl;      m[5] = 0.0;    //E
                    m[6] = cosb*cosl; m[7] = cosb*sinl; m[8] = sinb;   //U
                    
                    NEU.x = m[0]*tmp.x + m[1]*tmp.y + m[2]*tmp.z;
                    NEU.y = m[3]*tmp.x + m[4]*tmp.y + m[5]*tmp.z;
                    NEU.z = m[6]*tmp.x + m[7]*tmp.y + m[8]*tmp.z;
                    
                }
                
                void NEU2XYZ(double* NEU, double* refXYZ, double* XYZ)
                {
                    double E[9] = {0.0};
                    double blh[3] ={0.0};
                    XYZ2BLH(refXYZ, blh);
                    double sinb = sin(blh[0]);
                    double cosb = cos(blh[0]);
                    double sinl = sin(blh[1]);
                    double cosl = cos(blh[1]);
                    E[0] =-sinb*cosl;	E[1] =-sinl ;		E[2]=cosb*cosl;
                    E[3] =-sinb*sinl;	E[4] =cosl ;		E[5]=cosb*sinl;
                    E[6] =cosb;			E[7] =0;			E[8]=sinb;
                    
                    XYZ[0] = E[0]*NEU[0] + E[1]*NEU[1] +E[2]*NEU[2];
                    XYZ[1] = E[3]*NEU[0] + E[4]*NEU[1] +E[5]*NEU[2];
                    XYZ[2] = E[6]*NEU[0] + E[7]*NEU[1] +E[8]*NEU[2];
                    
                    
                }
                
                void NEU2XYZ(GVector NEU, GVector refXYZ, GVector& XYZ)
                {
                    double E[9] = {0.0};
                    GVector blh ;
                    XYZ2BLH(refXYZ, blh);
                    double sinb = sin(blh.x);
                    double cosb = cos(blh.x);
                    double sinl = sin(blh.y);
                    double cosl = cos(blh.y);
                    E[0] =-sinb*cosl;	E[1] =-sinl ;		E[2]=cosb*cosl;
                    E[3] =-sinb*sinl;	E[4] =cosl ;		E[5]=cosb*sinl;
                    E[6] =cosb;			E[7] =0;			E[8]=sinb;
                    
                    XYZ.x = E[0]*NEU.x + E[1]*NEU.y +E[2]*NEU.z;
                    XYZ.y = E[3]*NEU.x + E[4]*NEU.y +E[5]*NEU.z;
                    XYZ.z = E[6]*NEU.x + E[7]*NEU.y +E[8]*NEU.z;
                    
                }
                
                void Geocentric(GVector pos,double& lat, double& lon)
                {
                    double xyz[3] = {pos.x,pos.y,pos.z};
                    Geocentric(xyz, lat, lon );
                }
                
                /* get the geocentric lat and lon*/
                void Geocentric(double xyz[3],double& lat, double& lon)
                {
                    double PI = GCONST("PI");
                    
                    double r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]+ xyz[2]*xyz[2]);
                    if( r< 1.0E-12 )
                    {
                        lat = 0.0;
                        lon = 0.0;
                        return ;
                    }
                    
                    lat = acos(xyz[2]/r);
                    double rr = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] );
                    if( rr < 10E-12)
                    {
                        if(xyz[2]> 0.0 )
                        {
                            lat = PI/2;
                        }
                        else
                        {
                            lat = -PI/2;
                        }
                        lon = 0.0;
                        return ;
                    }
                    
                    lon = atan2(xyz[1],xyz[0]);
                    if(lon < 0.0 )
                    {
                        lon += 2*PI;
                    }
                }
                
                /// Convenience output method.
                std::ostream& dump( std::ostream& s ) const
                {
                    
                    s<<"ClassName:       "<<"Ellipsoid"<<std::endl;
                    s<<"Ellipsoid Name:  "<<name<<std::endl;
                    s<<"a:               "<<GString(a,3)<<std::endl;
                    s<<"flattening:      "<<GString(flattening,6)<<std::endl;
                    s<<"eccentricity:    "<<GString(eccentricity,6)<<std::endl;
                    s<<"angVelocity:     "<<GString(angVelocity,6)<<std::endl;
                    s<<"gm:              "<<GString(gm,6)<<std::endl;
                    s<<"J2:              "<<GString(j2,6)<<std::endl;
                    
                    return s;
                };
                
            private:
                std::string name;  //椭球名称(均为大写字母)
                long double a;     //椭球长半轴
                long double flattening;  //椭球扁率
                long double eccentricity; //离心率
                long double angVelocity;  //自转角速度 unit:弧度/秒
                long double gm;   //地球重力常数
                long double j2;   //椭球的J2项
                
                };   //class Ellipsoid End
                
                
                
                inline std::string asString(const GEllipsoid& p)
                {
                    std::ostringstream oss;
                    p.dump(oss);
                    return oss.str();
                }
                
                
                /// stream output for Ellipsoid
                inline std::ostream& operator<<(std::ostream& s, const GEllipsoid& p)
                {
                    p.dump(s);
                    return s;
                }
                
                }
                
#endif
