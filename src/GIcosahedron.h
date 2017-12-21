//
//  GIcosahedron.hpp
//  GFC
//
//  Created by lizhen on 16/2/6.
//  Copyright © 2016年 lizhen. All rights reserved.
//



//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 2.1 of the License, or
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



#ifndef GIcosahedron_hpp
#define GIcosahedron_hpp

#include <stdio.h>
#include "GString.h"
#include <math.h>

#include "GFCCONST.h"
#include <map>
namespace gfc
{
    
    // class Vertex for the graph
    class GVertex
    {
        
    public:
        
        GVertex()
        {
            m_x = 0.0; m_y = 0.0; m_z = 0.0;m_lon =0.0; m_lat =0.0;
        }
        
        GVertex( GString name, double x, double y, double z)
        {
            double pi = GCONST("PI");
            m_name = name;
            m_x = x; m_y = y; m_z = z;
            
            //convert xyz to BLH, normaly, H = 0
            double xyz[3]={0}; xyz[0]=x;xyz[1]=y;xyz[2]=z;
            double blh[3]={0};
            
            //ellipsoid paraments for the flux data
            double m_radius = 6371000;
            //double a = 6408137;     // a = 6378137 for WGS84
            //double b = 6386651.7;   // b=6356752.3142 for WGS84
            //double alfa = (m_radius - m_radius)/m_radius;
            //xyz2blh(xyz, blh, alfa, m_radius);
            //double tmp = sqrt(m_x*m_x + m_y*m_y);
            //blh[1] = asin(m_z/m_radius);
            //blh[0] = atan2(m_y, m_x);
            
            //m_lat = blh[0]*180.0/pi;
            //m_lon = blh[1]*180.0/pi;
            
            
            double dot1 =  sqrt( pow(m_x,2.0) + pow(m_y,2.0));
            if(fabs(dot1)<1.0E-10)
            {
                if( m_z > 0 )
                {
                    m_lat = 90.0;
                    m_lon = 0.0;
                }
                else if(m_z < 0)
                {
                    m_lat = -90.0;
                    m_lon = 0.0;
                }
            }
            else
            {
                m_lat =  atan(m_z/dot1)*180.0/pi; // -pi/2 ~ pi/2
                m_lon =  atan2(m_y, m_x)*180.0/pi;
                if( m_lon < 0.0 )
                {
                    m_lon = m_lon + 360.0;
                }
            }
                
        }
        
        GVertex(GString name, double* v)
        {
            m_name = name ;
            m_x = v[0];m_y = v[1]; m_z = v[2];
        }
        
        //copy constructor
        GVertex( const GVertex& v)
        {
            m_name = v.m_name;
            m_x = v.m_x;
            m_y = v.m_y;
            m_z = v.m_z;
            m_lat = v.m_lat;
            m_lon = v.m_lon;
        }
        
        GVertex& operator=(const GVertex& v)
        {
            m_name = v.m_name;
            m_x = v.m_x;
            m_y = v.m_y;
            m_z = v.m_z;
            m_lon = v.m_lon;
            m_lat = v.m_lat;
            return *this;
        }
        
        virtual ~GVertex() {}
        
        double getX() {return m_x;}
        double getY() {return m_y;}
        double getZ() {return m_z;}
        GString getName() {return m_name;}
        
        void setBL(double lon, double lat)
        {
            m_lat = lat;
            m_lon = lon;
        }
        
        double getLat() {return m_lat;}
        double getLon() {return m_lon;}
        
        
        void blh2xyz( double* blh, double* xyz, double alfa,double RE)
        {
            double e = sqrt( 2 * alfa - alfa * alfa);   //第一偏心率
            double sinb = sin(blh[0]);
            double sinl = sin(blh[1]);
            double cosb = cos(blh[0]);
            double cosl = cos(blh[1]);
            double N = RE/sqrt(1-pow(e*sinb,2));
            xyz[0] = (N + blh[2])*cosb*cosl;
            xyz[1] = (N + blh[2])*cosb*sinl;
            xyz[2] = (blh[2]+ N*(1-e*e))*sinb;
        }

        
        void xyz2blh( double* xyz, double* blh, double alfa,double RE)
        {
            double PI = GCONST("PI");
            double e2=alfa*(2.0-alfa),z,zk,v=RE,sinp;
            double r2 = xyz[0]*xyz[0] + xyz[1]*xyz[1];
            for ( z=xyz[2],zk=0.0;fabs(z-zk)>=1E-4;)
            {
                zk=z;
                sinp=z/sqrt(r2+z*z);
                v=RE/sqrt(1.0-e2*sinp*sinp);
                z=xyz[2]+v*e2*sinp;
            }
            blh[0]=r2>1E-12?atan(z/sqrt(r2)):(xyz[2]>0.0?PI/2.0:-PI/2.0);
            blh[1]=r2>1E-12?atan2(xyz[1],xyz[0]):0.0;
            blh[2]=sqrt(r2+z*z)-v;
            
        }

        
        
        
    private:
        GString m_name;
        double m_x;
        double m_y;
        double m_z;
        double m_lon;
        double m_lat;
    };
    
    //class edge for the graph
    template<class POINT>
    class GEdge
    {
        
    public:
        
        GEdge() {}
        
        GEdge( POINT v_start, POINT v_end )
        {
            m_start = v_start; m_end = v_end;
            m_name = v_start.getName() + "-->" + v_end.getName();
        }
        
        GEdge( GString name_s, double *start, GString name_e,double* end)
        {
            m_start = POINT(name_s,start);
            m_end  =  POINT(name_e,end);
        }
        
        
        virtual ~GEdge() {}
        
    private:
        GString m_name;
        POINT m_start;
        POINT m_end;
    };
    
    
    //class Triangle
    //   A<------C
    //    \     /
    //     \   /
    //      \B/
    // the order: A->B->C
    
    template<class POINT>
    class GTriangle
    {
        
    public:
        
        GTriangle() { m_value[0] = 0.0;m_value[1] = 0.0;m_area = 0.0;m_visible = false;}
        
        // please be sure about the order of these points
        GTriangle( GString name, POINT v1, POINT v2 , POINT v3 )
        {
            m_name = name;
            m_pts[0] = v1;
            m_pts[1] = v2;
            m_pts[2] = v3;
            m_value[0] = 0.0;
            m_value[1] = 0.0;
            
            calculateArea(6371.0);
            
        }
        
        //copy constructor
        GTriangle( const GTriangle& t)
        {
            m_name = t.m_name;
            m_tag = t.m_tag;
            m_Ftag = t.m_Ftag;
            m_Ctag[0] = t.m_Ctag[0];m_Ctag[1] = t.m_Ctag[1];m_Ctag[2] = t.m_Ctag[2];m_Ctag[3] = t.m_Ctag[3];
            m_pts[0] = t.m_pts[0]; m_pts[1] = t.m_pts[1]; m_pts[2] = t.m_pts[2];
            m_value[0] = t.m_value[0];
            m_value[1] = t.m_value[1];
            m_area = t.m_area;
            m_visible = t.m_visible;
        }
        
        GTriangle& operator=(const GTriangle& t)
        {
            m_name = t.m_name;
            m_tag = t.m_tag;
            m_Ftag = t.m_Ftag;
            m_Ctag[0] = t.m_Ctag[0];m_Ctag[1] = t.m_Ctag[1];m_Ctag[2] = t.m_Ctag[2];m_Ctag[3] = t.m_Ctag[3];
            m_pts[0] = t.m_pts[0]; m_pts[1] = t.m_pts[1]; m_pts[2] = t.m_pts[2];
            m_value[0] = t.m_value[0];
            m_value[1] = t.m_value[1];
            m_area = t.m_area;
            
            m_visible =  t.m_visible;
            
            return *this;
        }

        
        
        POINT getPoint(int index)
        {
            return m_pts[index];
        }
        
        //get the centroid of this triangle
        POINT getCentroid()
        {
            double xyz[3]={0.0};
            for( int i = 0 ; i< 3; i++ )
            {
                xyz[0] += m_pts[i].getX();
                xyz[1] += m_pts[i].getY();
                xyz[2] += m_pts[i].getZ();
            }
            
            xyz[0] = xyz[0]/3.0;
            xyz[1] = xyz[1]/3.0;
            xyz[2] = xyz[2]/3.0;
            double len = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
            //extend this point onto the sphere
            xyz[0] = xyz[0]/len*6371000.0;
            xyz[1] = xyz[1]/len*6371000.0;
            xyz[2] = xyz[2]/len*6371000.0;
            
            POINT p("centroid",xyz[0],xyz[1],xyz[2]);
            
            return p;
        }
        
        
        double getArea()
        {
            return m_area;
        }
        
        
        void SetVisible(bool t)
        {
            m_visible = t;
        }
        
        void setArea(double area)
        {
            m_area = area;
        }
        /*return the area of spherical triangle
          * ref: http://www.doc88.com/p-9307125540394.html
          *
          *
          *
         */
        void calculateArea(double R)
        {
            double theta[3]={0.0};
            double phi[3] = {0.0};
            double h[3] = {0.0};
            double PI = GCONST("PI");
           
            for( int i = 0 ; i<3; i++ )
            {
                phi[i] = m_pts[i].getLat()*PI/180.0;
                if(m_pts[i].getLon()>=180.0)
                {
                    theta[i] = (m_pts[i].getLon()-360.0)*PI/180.0;
                }
                else
                {
                   theta[i] = (m_pts[i].getLon())*PI/180.0;
                }
            }
            
        h[0] = cos(phi[1])*cos(phi[2])*cos(theta[1]-theta[2]) + sin(phi[1])*sin(phi[2]);
        h[1] = cos(phi[2])*cos(phi[0])*cos(theta[2]-theta[0]) + sin(phi[2])*sin(phi[0]);
        h[2] = cos(phi[0])*cos(phi[1])*cos(theta[0]-theta[1]) + sin(phi[0])*sin(phi[1]);
            
            double A = acos( (h[0]-h[1]*h[2])/sqrt( (1.0-h[1]*h[1])*(1.0-h[2]*h[2])) );
            double B = acos( (h[1]-h[2]*h[0])/sqrt( (1.0-h[2]*h[2])*(1.0-h[0]*h[0])) );
            double C = acos( (h[2]-h[0]*h[1])/sqrt( (1.0-h[0]*h[0])*(1.0-h[1]*h[1])) );
            
            m_area = R*R*(A+B+C - PI);
            
        }
        
        double getX(int index)
        {
           return  m_pts[index].getX();
        }
        
        double getY(int index)
        {
            return  m_pts[index].getY();
        }
        
        double getZ(int index)
        {
            return  m_pts[index].getZ();
        }
        
        GString getName()
        {
            return m_name;
        }
        
        double getValue(int index)
        {
            return m_value[index];
        }
        
        void setValue( double value,int index )
        {
            m_value[index] = value;
        }
        
        void setTag(GString tag)
        {
            m_tag = tag;
        }
        
        GString getTag()
        {
            return m_tag;
        }
        
        //设置上一层的编号
        void setFtag( GString tag )
        {
            m_Ftag = tag;
        }
        
        GString getFtag()
        {
            return m_Ftag;
        }
        
        //设置下一层的子节点编号,因为有4个子节点
        void setCtag( GString tag[4])
        {
            m_Ctag[0] = tag[0];
            m_Ctag[1] = tag[1];
            m_Ctag[2] = tag[2];
            m_Ctag[3] = tag[3];
        }
        
        // get the children tag of this node
        void getCtag(std::vector<GString>& ptag )
        {
            ptag.push_back(m_Ctag[0]);
            ptag.push_back(m_Ctag[1]);
            ptag.push_back(m_Ctag[2]);
            ptag.push_back(m_Ctag[3]);
        }
        
        // get the point order of m_v1, m_v2 and m_v3
        // return 1 : in counterclockwise order
        // return -1 : in clockwise order
        int getPointOrder()
        {
            return 1;
        }
        
        void write2File(FILE* pf)
        {
            //std::vector<GString> childrenTag ;
            //m_TriGrid[i][j].getCtag(childrenTag);
            //double area = m_TriGrid[i][j].getArea();
            fprintf(pf, "%9s %9s %9s %9s %9s %9s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %16.8f %12.4f %12.4f\n",
                    m_tag.c_str(),
                    m_Ftag.c_str(),
                    m_Ctag[0].c_str(),
                    m_Ctag[1].c_str(),
                    m_Ctag[2].c_str(),
                    m_Ctag[3].c_str(),
                    m_pts[0].getX(),
                    m_pts[0].getY(),
                    m_pts[0].getZ(),
                    m_pts[1].getX(),
                    m_pts[1].getY(),
                    m_pts[1].getZ(),
                    m_pts[2].getX(),
                    m_pts[2].getY(),
                    m_pts[2].getZ(),
                    m_area,
                    m_value[0],
                    m_value[1]
                    );

        }
        
        
    private:
        // the vertex should be in counter-clockwise order looking from the normal direction
        // So the order of v1 , v2 and v3 is very important
        GString m_name;
        GString  m_Ftag;  // 当前节点的父节点
        GString  m_tag; // tag is the variable which record the father triangle and the number of it
        POINT  m_pts[3];  // the index of the points
        double m_area;  // the area of this triangle
        GString  m_Ctag[4]; // 4 children nodes
        
        double m_value[2];  // the flux data value in this triangle, m_value[0] longwave, m_value[1] shortwave
        
        bool  m_visible;  // it is a temp variable, whether visible for the satellite
        
        
        //POINT m_v1;
        //POINT m_v2;
        //POINT m_v3;
    };
    
    //
    template<class POINT>
    class GIcosahedron
    {
        public:
        GIcosahedron()
        {
            m_level = 0 ;
            originalShape();
            
            pf = fopen("triNet.data","w+");
        }
        
        void setRadius( double radius )
        {
            m_radius = radius;
        }
        
        
        void outputTriangleNetJS(GString filepath)
        {
            double pi = GCONST("PI");
            FILE* pf = fopen(filepath.c_str(),"w+");
            fprintf(pf, "var myfaces=[ \n");
            for( int k = 0 ; k< m_faceList.size() ; k++)
            {
                GVertex p1 = m_faceList[k].getPoint(0);
                GVertex p2 = m_faceList[k].getPoint(1);
                GVertex p3 = m_faceList[k].getPoint(2);
                int tmp[3]={0};
                tmp[0] = p1.getName().asINT()-1;
                tmp[1] = p2.getName().asINT()-1;
                tmp[2] = p3.getName().asINT()-1;
                
                fprintf(pf, "[%6d ,%6d, %6d ]",tmp[0], tmp[1],tmp[2]);
                if( k != m_faceList.size() -1 )
                {
                    fprintf(pf, ",");
                }
                fprintf(pf, "\n");
                
            }
            
            fprintf(pf, "];\n");
            fclose(pf);
            
            
        }
        
        void outputTriangleNet(GString filepath)
        {
            double pi = GCONST("PI");
            FILE* pf = fopen(filepath.c_str(),"w+");
            for( int k = 0 ; k < m_faceList.size() ; k++)
            {
                GVertex p1 = m_faceList[k].getPoint(0);
                GVertex p2 = m_faceList[k].getPoint(1);
                GVertex p3 = m_faceList[k].getPoint(2);
                
//                double dot1 =  sqrt(pow(p1.getX(),2.0) + pow(p1.getY(),2.0));
//                double lat1 =  atan(p1.getZ()/dot1)*180.0/pi; // -pi/2 ~ pi/2
//                double lon1 =  atan2(p1.getY(), p1.getX())*180.0/pi;
//                if( lon1 < 0.0 )
//                {
//                    lon1 = lon1 + 360.0;
//                }
//                
//                double dot2 = sqrt(pow(p2.getX(),2.0) + pow(p2.getY(),2.0));
//                double lat2 =  atan(p2.getZ()/dot2)*180.0/pi; // -pi/2 ~ pi/2
//                double lon2 =  atan2(p2.getY(), p2.getX())*180.0/pi;
//                if( lon2 < 0.0 )
//                {
//                    lon2 = lon2 + 360.0;
//                }
//                
//                double dot3 =  sqrt(pow(p3.getX(),2.0) + pow(p3.getY(),2.0));
//                double lat3 =  atan(p3.getZ()/dot3)*180.0/pi; // -pi/2 ~ pi/2
//                double lon3 =  atan2(p3.getY(), p3.getX())*180.0/pi;
//                if( lon3 <0.0 )
//                {
//                    lon3 = lon3 + 360.0;
//                }
                
                fprintf(pf, "%06d %24s %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", k,
                        m_faceList[k].getName().c_str(),
                        p1.getLat(),p1.getLon(),
                        p2.getLat(),p2.getLon(),
                        p3.getLat(),p3.getLon()
                        );
            }
            
            fclose(pf);
        }
        
        void outputVertexJS(GString filepath)
        {
            double pi = GCONST("PI");
            FILE* pf = fopen(filepath.c_str(),"w+");
            fprintf(pf, "var myvertex=[ \n");
            for( int k = 0 ; k< m_vertexList.size() ; k++ )
            {
                
//                double dot1 = sqrt(pow(m_vertexList[k].getX(),2.0) + pow(m_vertexList[k].getY(),2.0));
//                double lat1 =  atan(m_vertexList[k].getZ()/dot1)*180.0/pi; // -pi/2 ~ pi/2
//                double lon1 =  atan2(m_vertexList[k].getY(), m_vertexList[k].getX())*180.0/pi;
//                if( lon1 <0.0 )
//                {
//                    lon1 = lon1 + 360.0;
//                }
                fprintf(pf, "[%12.4f ,%12.4f, %12.4f, %12.4f, %12.4f ]",m_vertexList[k].getX(),m_vertexList[k].getY(),m_vertexList[k].getZ(),m_vertexList[k].getLon(),m_vertexList[k].getLat());
                if( k != m_vertexList.size() -1 )
                {
                    fprintf(pf, ",");
                }
                fprintf(pf, "\n");
            }
            
             fprintf(pf, "];\n");
            
            if(pf != NULL )
            {
                fclose(pf);
                pf = NULL;
            }
            
        }
        
        
        void outputVertex(GString filepath)
        {
            double pi = GCONST("PI");
            FILE* pf = fopen(filepath.c_str(),"w+");
            for( int k = 0 ; k< m_vertexList.size() ; k++ )
            {
             
                double dot1 = sqrt(pow(m_vertexList[k].getX(),2.0) + pow(m_vertexList[k].getY(),2.0));
                double lat1 =  atan(m_vertexList[k].getZ()/dot1)*180.0/pi; // -pi/2 ~ pi/2
                double lon1 =  atan2(m_vertexList[k].getY(), m_vertexList[k].getX())*180.0/pi;
                if( lon1 <0.0 )
                {
                    lon1 = lon1 + 360.0;
                }
                fprintf(pf, "%06d %s %12.6f %12.6f %12.6f %12.6f %12.6f\n", k,
                        m_vertexList[k].getName().c_str(),
                        m_vertexList[k].getX(),
                        m_vertexList[k].getY(),
                        m_vertexList[k].getZ(),
                        lat1,lon1
                        );
             
             }
             
            if(pf != NULL )
            {
                fclose(pf);
                pf = NULL;
            }
            
        }
        
        // create the original icosahedron with 12 points and 20 faces
        void originalShape()
        {
            
            m_index = 0;
            m_vertexList.clear();
            m_middlePointIndexCache.clear();
           
            int npoints = 6 ;
            int nfaces  = 8 ;
            /*
            POINT pts[6] = { POINT("1",0,0,m_b), POINT("2",m_a,0,0),POINT("3",0,m_a,0),
                             POINT("4",-m_a,0,0),POINT("5",0,-m_a,0),POINT("6",0,0,-m_b)
                           };
                */
            
            POINT pts[6] =
            {   POINT("1",0,0,m_radius), POINT("2",m_radius,0,0),POINT("3",0,m_radius,0),
                POINT("4",-m_radius,0,0),POINT("5",0,-m_radius,0),POINT("6",0,0,-m_radius)
            };
            
            
            GTriangle<POINT> faces[8] =
            {
                GTriangle<POINT>("1-2-3",pts[0],pts[1],pts[2]),
                GTriangle<POINT>("1-3-4",pts[0],pts[2],pts[3]),
                GTriangle<POINT>("1-4-5",pts[0],pts[3],pts[4]),
                GTriangle<POINT>("1-5-2",pts[0],pts[4],pts[1]),
                GTriangle<POINT>("2-6-3",pts[1],pts[5],pts[2]),
                GTriangle<POINT>("3-6-4",pts[2],pts[5],pts[3]),
                GTriangle<POINT>("4-6-5",pts[3],pts[5],pts[4]),
                GTriangle<POINT>("5-6-2",pts[4],pts[5],pts[1])
            };
            
            //设定triangle 的m_tag 属性，tag 的格式为Nddnnnn-Fddnnn－Addnnnn-Bddnnnn-Cddnnnn-Dddnnnn
            // C表示当前，dd表示level，nnnn表示这个三角形在这个level中的序号
            // L表示上一层级，格式一致
            
            
//            int npoints = 3 ;
//            int nfaces  = 1 ;
//            // for the xyz coordinates
//            POINT pts[3] = { POINT("1",0,0,m_b), POINT("2",m_a,0,0),POINT("3",0,m_a,0) };
//            //POINT pts[3] = { POINT("1",90,0,0), POINT("2", 0,0,0),POINT("3",0,90,0)};
//            GTriangle<POINT> faces[1] =
//            {
//                GTriangle<POINT>("1-2-3",pts[0],pts[1],pts[2])
//            };
//
            
            for( int i = 0 ; i< npoints ; i++ )
            {
                AddVertex(pts[i]);
            }
            
            // for the origin shape, there is no last level, so just set the last tag as F000000
            for( int i = 0 ; i< nfaces; i++ )
            {
                char tmp[9] = {0}; // 必须是8位
                sprintf(tmp, "%08d",i+1);
                GString str1 = "0";  //level 0
                GString str2 = tmp;
                //set father
                faces[i].setFtag("N00000000"); // the very begining node
                faces[i].setTag(str1+str2);
                
                m_faceList.push_back(faces[i]);
            }
            
        }
        
        
        std::vector< std::vector<GTriangle<POINT> > > createGeometry( int level, std::vector< std::vector<POINT> >& myvertices )
        {
            std::vector< std::vector<GTriangle<POINT> > > retval;
            
            //retval.push_back(m_faceList);
            std::vector<GTriangle<POINT> > faces2;
            for( int i = 0; i < level; ++i )
            {
                int count = 0;  // count for this level
                char tmp1[2] ={0}; // level str
                char tmp2[9] ={0}; // number str
                sprintf(tmp1, "%01d",i+1);
                GString strlevel = tmp1;
                // every level should have its faces
                //std::vector<GTriangle<POINT> > faces2;
                faces2.clear();
                
                if( i == 0 )
                {
                    myvertices.push_back(m_vertexList);
                }
                
                for( int j = 0 ; j< m_faceList.size(); j++ )
                {
                    POINT p1 = m_faceList[j].getPoint(0);
                    POINT p2 = m_faceList[j].getPoint(1);
                    POINT p3 = m_faceList[j].getPoint(2);
                    
                    GString strF = m_faceList[j].getTag();
                    
                    
                    int ip1 = m_faceList[j].getPoint(0).getName().asINT();
                    int ip2 = m_faceList[j].getPoint(1).getName().asINT();
                    int ip3 = m_faceList[j].getPoint(2).getName().asINT();
                    
                    int a = GetMiddlePoint( ip1-1, ip2-1);
                    int b = GetMiddlePoint( ip2-1, ip3-1);
                    int c = GetMiddlePoint( ip3-1, ip1-1);
                    
                    GTriangle<POINT> t1( p1.getName()+"-"+m_vertexList[a].getName()+"-"+m_vertexList[c].getName(),p1,m_vertexList[a],m_vertexList[c]);
                    GTriangle<POINT> t2( p2.getName()+"-"+m_vertexList[b].getName()+"-"+m_vertexList[a].getName(),p2,m_vertexList[b],m_vertexList[a]);
                    GTriangle<POINT> t3( p3.getName()+"-"+m_vertexList[c].getName()+"-"+m_vertexList[b].getName(),p3,m_vertexList[c],m_vertexList[b]);
                    GTriangle<POINT> t4( m_vertexList[a].getName()+"-"+m_vertexList[b].getName()+"-"+m_vertexList[c].getName(),m_vertexList[a],m_vertexList[b],m_vertexList[c]);
                    
                    //set father
                    t1.setFtag(strF); t2.setFtag(strF);t3.setFtag(strF);t4.setFtag(strF);
                    GString ctag[4];
                    sprintf(tmp2, "%08d",count+1); ctag[0] = tmp2; ctag[0]= strlevel + ctag[0];
                    sprintf(tmp2, "%08d",count+2); ctag[1] = tmp2; ctag[1]= strlevel + ctag[1];
                    sprintf(tmp2, "%08d",count+3); ctag[2] = tmp2; ctag[2]= strlevel + ctag[2];
                    sprintf(tmp2, "%08d",count+4); ctag[3] = tmp2; ctag[3]= strlevel + ctag[3];
                    
                    m_faceList[j].setCtag(ctag);  //为父节点设置子节点
                    
                    t1.setTag(ctag[0]);t2.setTag(ctag[1]);t3.setTag(ctag[2]);t4.setTag(ctag[3]);
                    faces2.push_back(t1);faces2.push_back(t2);faces2.push_back(t3);faces2.push_back(t4);
                    
                
                   // faces2.push_back( GTriangle<POINT>(GString(p1.getName()+"-"+m_vertexList[a].getName()+"-"+m_vertexList[c].getName()),p1,m_vertexList[a],m_vertexList[c] ));
                    
                   // faces2.push_back( GTriangle<POINT>(p2.getName()+"-"+m_vertexList[b].getName()+"-"+m_vertexList[a].getName(),p2,m_vertexList[b],m_vertexList[a]) );
                    
                   // faces2.push_back( GTriangle<POINT>(p3.getName()+"-"+m_vertexList[c].getName()+"-"+m_vertexList[b].getName(),p3,m_vertexList[c],m_vertexList[b]));
                    
                  //  faces2.push_back( GTriangle<POINT>(m_vertexList[a].getName()+"-"+m_vertexList[b].getName()+"-"+m_vertexList[c].getName(),m_vertexList[a],m_vertexList[b],m_vertexList[c]));
                   
                    count = count +4 ;
                }
               
                 retval.push_back(m_faceList);
                 m_faceList = faces2;
                 myvertices.push_back(m_vertexList);
            }
            
            if( level == 0 )
            {
                retval.push_back(m_faceList);
            }
            else
            {
               retval.push_back(faces2);
            }
            
            
            
            /*
            for( int j = 0 ; j< m_faceList.size(); j++ )
            {
                POINT p1 = m_faceList[j].getPoint(0);
                POINT p2 = m_faceList[j].getPoint(1);
                POINT p3 = m_faceList[j].getPoint(2);
                
                
                int ip1 = m_faceList[j].getPoint(0).getName().asINT();
                int ip2 = m_faceList[j].getPoint(1).getName().asINT();
                int ip3 = m_faceList[j].getPoint(2).getName().asINT();
                m_indices.push_back(ip1-1);
                m_indices.push_back(ip2-1);
                m_indices.push_back(ip3-3);
                
            }
            */
        
            return retval;
        }
        
        virtual ~GIcosahedron() {}
        
        private:
        // private function
        int AddVertex( POINT position)
        {
            double xyz[3]={0.0};
            double pi = GCONST("PI");
            /* for xyz coordinates */
            
            //double length = sqrt(position.x * position.x + position.y * position.y + position.z * position.z);
            double length = sqrt( pow(position.getX(),2.0)+pow(position.getY(),2.0)+pow(position.getZ(),2.0));
            if( fabs(length)<0.00001 )
            {
                length = 0.00000001;
            }
            
            // this is for the spherical surface
            double r = m_radius;
            
            //double sin_alfa_2 = pow(position.getZ(),2.0)/( position.getX()*position.getX() + position.getY()*position.getY() + position.getZ()*position.getZ());
            //double r = m_a*m_b/sqrt((m_a*m_a-m_b*m_b)*sin_alfa_2+m_b*m_b);
            
            // this is the error way
            // double tan_alfa_2 =0.0, r =0.0;
            // if((position.getX()*position.getX() + position.getY()*position.getY())<0.0001)
            // {
            // r = m_b;
            // }
            // else
            // {
            // tan_alfa_2 = pow(position.getZ(),2.0)/( position.getX()*position.getX() + position.getY()*position.getY() );
            // r = sqrt((m_a*m_a +m_b*m_b*tan_alfa_2)/(1.0+tan_alfa_2));
            // }
            
            // make the point on the spherical surface, BUT, how abou the ellipsoidal surface ???
            xyz[0] = position.getX()/length*r;
            xyz[1] = position.getY()/length*r;
            xyz[2] = position.getZ()/length*r;
            
//            /* for spherical surface */
//            xyz[0] = position.getX() ;
//            xyz[1] = position.getY() ;
            
            // form new point
            POINT pt(GString(m_vertexList.size()+1),xyz[0],xyz[1],xyz[2]);
//            double dot1 =  sqrt(pow(pt.getX(),2.0) + pow(pt.getY(),2.0));
//            double lat1 =  atan(pt.getZ()/dot1)*180.0/pi; // -pi/2 ~ pi/2
//            double lon1 =  atan2(pt.getY(), pt.getX())*180.0/pi;
//            if( lon1 < 0.0 )
//            {
//                lon1 = lon1 + 360.0;
//            }
//            pt.setBL(lon1,lat1);
            
            m_vertexList.push_back(pt);
            return m_index++;
        }
        
        // p1 and p2 are the index of points in m_vertexList
        int GetMiddlePoint(int p1, int p2 )
        {
            bool firstPointIsSmaller = p1 < p2;
            double xyz[3] = {0.0};
            int64_t smallerIndex = firstPointIsSmaller ? p1 : p2;
            int64_t greaterIndex = firstPointIsSmaller ? p2 : p1;
            int64_t key = (smallerIndex << 32) + greaterIndex;
            
            auto foundValueIterator = m_middlePointIndexCache.find(key);
            if( foundValueIterator !=  m_middlePointIndexCache.end())
            {
                return foundValueIterator->second;
            }
            
            POINT point1 = m_vertexList[p1]; // m_vertices[p1];
            POINT point2 = m_vertexList[p2]; // m_vertices[p2];
            /*for xyz coordinates*/
            xyz[0] = (point1.getX() + point2.getX()) / 2.0 ;
            xyz[1] = (point1.getY() + point2.getY()) / 2.0 ;
            xyz[2] = (point1.getZ() + point2.getZ()) / 2.0 ;
            
            POINT middle = POINT("",  xyz[0],xyz[1],xyz[2]);
            
            int i = this->AddVertex(middle);
            
            this->m_middlePointIndexCache.insert(std::make_pair(key, i));
            
            return i;
        }
        
        double m_radius = 6371000;
        //a = 6408.1370 km and b = 6386.6517
        
        double m_a = 6408137;     // a = 6378137 for WGS84
        
        double m_b = 6386651.7;   //b=6356752.3142 for WGS84
        
        FILE* pf;
        int m_index;
        //std::vector<unsigned int > m_indices;
        std::map< int64_t, int > m_middlePointIndexCache;
        std::vector<POINT>  m_vertexList;
        std::vector< GTriangle<POINT> >  m_faceList;
        int                           m_level;
        
    };
    
    
    
    
    namespace PE
    {
        
        struct TriangleIndices
        {
            int v1, v2, v3;
            TriangleIndices(int v1, int v2, int v3)
            {
                this->v1 = v1; this->v2 = v2; this->v3 = v3;
            }
           
        };
        
        class IcoSphere
        {
        public:
            
            IcoSphere()
            {
                
            }
            
            ~IcoSphere()
            {
                
            }
            
            void Create(int recursionLevel)
            {
                middlePointIndexCache.clear();
                vertices.clear();
                indices.clear();
                index = 0;
                std::vector<TriangleIndices> faces;
                
                AddVertex(gfc::GVertex("0", 0, 0, 1));
                AddVertex(gfc::GVertex("1" ,1, 0, 0));
                AddVertex(gfc::GVertex("2", 0, 1, 0));
                AddVertex(gfc::GVertex("3", -1, 0, 0));
                AddVertex(gfc::GVertex("4", 0, -1, 0));
                AddVertex(gfc::GVertex("5", 0, 0, -1));
                
                faces.push_back(TriangleIndices(0, 1, 2));
                faces.push_back(TriangleIndices(0, 2, 3));
                faces.push_back(TriangleIndices(0, 3, 4));
                faces.push_back(TriangleIndices(0, 4, 1));
                faces.push_back(TriangleIndices(1, 5, 2));
                faces.push_back(TriangleIndices(2, 5, 3));
                faces.push_back(TriangleIndices(3, 5, 4));
                faces.push_back(TriangleIndices(4, 5, 1));
               
            /*
                auto t = (1.0 + sqrt(5.0)) / 2.0;
                
                AddVertex(gfc::GVertex("", -1, t, 0));
                AddVertex(gfc::GVertex("" ,1, t, 0));
                AddVertex(gfc::GVertex("", -1, -t, 0));
                AddVertex(gfc::GVertex("", 1, -t, 0));
                
                AddVertex(gfc::GVertex("", 0, -1, t));
                AddVertex(gfc::GVertex("", 0, 1, t));
                AddVertex(gfc::GVertex("", 0, -1, -t));
                AddVertex(gfc::GVertex("", 0, 1, -t));
                
                AddVertex(gfc::GVertex("", t, 0, -1));
                AddVertex(gfc::GVertex("", t, 0, 1));
                AddVertex(gfc::GVertex("", -t, 0, -1));
                AddVertex(gfc::GVertex("", -t, 0, 1));
               
               
               
                
                faces.push_back(TriangleIndices(0, 11, 5));
                faces.push_back(TriangleIndices(0, 5, 1));
                faces.push_back(TriangleIndices(0, 1, 7));
                faces.push_back(TriangleIndices(0, 7, 10));
                faces.push_back(TriangleIndices(0, 10, 11));
                
                faces.push_back(TriangleIndices(1, 5, 9));
                faces.push_back(TriangleIndices(5, 11, 4));
                faces.push_back(TriangleIndices(11, 10, 2));
                faces.push_back(TriangleIndices(10, 7, 6));
                faces.push_back(TriangleIndices(7, 1, 8));
                
                faces.push_back(TriangleIndices(3, 9, 4));
                faces.push_back(TriangleIndices(3, 4, 2));
                faces.push_back(TriangleIndices(3, 2, 6));
                faces.push_back(TriangleIndices(3, 6, 8));
                faces.push_back(TriangleIndices(3, 8, 9));
                
                faces.push_back(TriangleIndices(4, 9, 5));
                faces.push_back(TriangleIndices(2, 4, 11));
                faces.push_back(TriangleIndices(6, 2, 10));
                faces.push_back(TriangleIndices(8, 6, 7));
                faces.push_back(TriangleIndices(9, 8, 1));
                */
                for(int i = 0; i < recursionLevel; ++i)
                {
                    //auto faces2 = std::vector<std::shared_ptr<triangleindices>>();
                    std::vector<TriangleIndices> faces2;
                    std::vector<TriangleIndices>::iterator tri = faces.begin();
                    for( tri = faces.begin(); tri!= faces.end(); tri++)
                    {
                        int a = GetMiddlePoint(tri->v1, tri->v2);
                        int b = GetMiddlePoint(tri->v2, tri->v3);
                        int c = GetMiddlePoint(tri->v3, tri->v1);
                        
                        faces2.push_back(TriangleIndices(tri->v1, a, c));
                        faces2.push_back(TriangleIndices(tri->v2, b, a));
                        faces2.push_back(TriangleIndices(tri->v3, c, b));
                        faces2.push_back(TriangleIndices(a, b, c));
                    }
                    faces = faces2;
                }
                
                std::vector<TriangleIndices>::iterator tri = faces.begin();
                for( tri = faces.begin(); tri!= faces.end(); tri++)
                //for(auto tri : faces)
                {
                    this->indices.push_back(tri->v1);
                    this->indices.push_back(tri->v2);
                    this->indices.push_back(tri->v3);
                }
                
                int testc = 0;
                FILE* pf = fopen("testICO.data","w+");
                double pi = gfc::GCONST("PI");
                /*
                for( int k = 0 ; k< vertices.size() ; k++ )
                {
                    
                    double dot1 = sqrt(pow(vertices[k].getX(),2.0) + pow(vertices[k].getY(),2.0));
                    double lat1 =  atan(vertices[k].getZ()/dot1)*180.0/pi; // -pi/2 ~ pi/2
                    double lon1 =  atan2(vertices[k].getY(), vertices[k].getX())*180.0/pi;
                    if( lon1 <0.0 )
                    {
                        lon1 = lon1 + 360.0;
                    }
                fprintf(pf, "%06d %12.6f %12.6f %12.6f %12.6f %12.6f\n", k,
                                            vertices[k].getX(),
                                            vertices[k].getY(),
                                            vertices[k].getZ(),
                                            lat1,lon1
                                            );
                    
                }
                */
                //just for output
                for( int k = 0 ; k< faces.size() ; k++)
                {
                    GVertex p1 = vertices[faces[k].v1];
                    GVertex p2 = vertices[faces[k].v2];
                    GVertex p3 = vertices[faces[k].v3];
                    
                    double dot1 = sqrt(pow(p1.getX(),2.0) + pow(p1.getY(),2.0));
                    double lat1 =  atan(p1.getZ()/dot1)*180.0/pi; // -pi/2 ~ pi/2
                    double lon1 =  atan2(p1.getY(), p1.getX())*180.0/pi;
                    if( lon1 < 0.0 )
                    {
                        lon1 = lon1 + 360.0;
                    }
                    
                    double dot2 = sqrt(pow(p2.getX(),2.0) + pow(p2.getY(),2.0));
                    double lat2 =  atan(p2.getZ()/dot2)*180.0/pi; // -pi/2 ~ pi/2
                    double lon2 =  atan2(p2.getY(), p2.getX())*180.0/pi;
                    if( lon2 < 0.0 )
                    {
                        lon2 = lon2 + 360.0;
                    }

                    double dot3 =  sqrt(pow(p3.getX(),2.0) + pow(p3.getY(),2.0));
                    double lat3 =  atan(p3.getZ()/dot3)*180.0/pi; // -pi/2 ~ pi/2
                    double lon3 =  atan2(p3.getY(), p3.getX())*180.0/pi;
                    if( lon3 <0.0 )
                    {
                        lon3 = lon3 + 360.0;
                    }

                    fprintf(pf, "%06d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", k,
                            lat1,lon1,
                            lat2,lon2,
                            lat3,lon3
                            );
                }
                
                fclose(pf);
                
            }
            
        private:
            
            int AddVertex(gfc::GVertex position)
            {
                
                //double length = sqrt(position.x * position.x + position.y * position.y + position.z * position.z);
                double length = sqrt(pow(position.getX(),2.0)+pow(position.getY(),2.0)+pow(position.getZ(),2.0));
                if(fabs(length)<0.00001)
                {
                    length = 0.00000001;
                }
                GVertex pt("",position.getX()/length,position.getY()/length,position.getZ()/length);
                //vertices.push_back(glm::vec3(position.x/length, position.y/length, position.z/length));
                vertices.push_back(pt);
                return index++;
            }
            
            int GetMiddlePoint(int p1, int p2 )
            {
                bool firstPointIsSmaller = p1 < p2;
                int64_t smallerIndex = firstPointIsSmaller ? p1 : p2;
                int64_t greaterIndex = firstPointIsSmaller ? p2 : p1;
                int64_t key = (smallerIndex << 32) + greaterIndex;
                
                auto foundValueIterator = middlePointIndexCache.find(key);
                if( foundValueIterator != middlePointIndexCache.end())
                {
                    return foundValueIterator->second;
                }
                
                GVertex point1 = vertices[p1];
                GVertex point2 = vertices[p2];
                GVertex middle = GVertex("", (point1.getX() + point2.getX()) / 2.0,
                                             (point1.getY() + point2.getY()) / 2.0,
                                             (point1.getZ() + point2.getZ()) / 2.0);
                
                int i = this->AddVertex(middle);
                
                this->middlePointIndexCache.insert(std::make_pair(key, i));
                return i;
            }
            
        public:
            
            std::vector<gfc::GVertex> vertices;
            std::vector<unsigned int > indices;
            
        private:
            
            int index;
            std::map< int64_t, int > middlePointIndexCache;
        };
        
    };
    
    
} // end of namespace gfc









#endif /* GIcosahedron_hpp */
