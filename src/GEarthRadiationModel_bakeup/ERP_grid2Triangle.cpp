//
//  ERP_grid2Triangle.cpp
//  GFC
//
//  Created by lizhen on 16/3/11.
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

#include "ERP_grid2Triangle.h"

namespace gfc
{
    
    GEarthRadiationModel::GEarthRadiationModel()
    {
        
    }
    
    GEarthRadiationModel::~GEarthRadiationModel()
    {
        
    }
    
    //dump original monthly averaged ceres flux data
    void GEarthRadiationModel::dumpMonToaData()
    {
        for( int k =0 ; k< m_monToaData.size() ; k++ )
        {
            char tmpstr[4]={0};
            sprintf(tmpstr, "%02d",k+1);
            GString dataFileName(tmpstr);
            dataFileName.append(".data");
            GString dirpath = "../data/";
            
            FILE* meanFile = fopen((dirpath+dataFileName).c_str(),"w+");
            //fprintf(meanFile, "Month: %02d\n",k+1);
            for( int i = 0 ; i< 360; i++  ) //360
            {
                for( int j = 0 ; j< 180; j++ ) //180
                {
                    //m_monToaData[k].data[j][i] = m_monToaData[k].data[j][i] / dataTag[k];
                    fprintf(meanFile, "%8.3f ",m_monToaData[k].data[j][i]);
                }
                fprintf(meanFile, "\n");
            }
            fclose(meanFile);
        }
    }
    
    /*dump the triangle grids and the flux data in each triangle grid*/
    void GEarthRadiationModel::dumpTriGridFlux()
    {
        GString faceFilePath = "../data/";
        GString levelstr = "level";
        std::vector<GString> childrenTag ;
        for( int i = 0; i < m_level+1 ; i++ )
        {
            char charlevel[2] ={0};
            sprintf(charlevel, "%1d",i);
            GString tmpstr(charlevel);
            levelstr = levelstr + tmpstr;
            faceFilePath = faceFilePath +levelstr;
            FILE* faceFile = fopen(faceFilePath.c_str(),"w+");
            
            for( int j = 0 ; j < m_TriGrid[i].size(); j++ )
            {
                m_TriGrid[i][j].getCtag(childrenTag);
                fprintf(faceFile, "%9s %9s %9s %9s %9s %9s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                        m_TriGrid[i][j].getTag().c_str(),
                        m_TriGrid[i][j].getFtag().c_str(),
                        childrenTag[0].c_str(),
                        childrenTag[1].c_str(),
                        childrenTag[2].c_str(),
                        childrenTag[3].c_str(),
                        m_TriGrid[i][j].getX(0),
                        m_TriGrid[i][j].getY(0),
                        m_TriGrid[i][j].getZ(0),
                        m_TriGrid[i][j].getX(1),
                        m_TriGrid[i][j].getY(1),
                        m_TriGrid[i][j].getZ(1),
                        m_TriGrid[i][j].getX(2),
                        m_TriGrid[i][j].getY(2),
                        m_TriGrid[i][j].getZ(2),
                        m_TriGrid[i][j].getValue()
                        );
            }
            
            if( faceFile != NULL )
            {
                fclose(faceFile);
            }
        }
    }
    
    //output the triangle grid flux into python script
    void GEarthRadiationModel::outputTriFlux_py(int mylevel)
    {
        //for( int i = 0; i < m_level+1 ; i++ )
        //{
            int i = mylevel;
            FILE* fluxData_py = fopen("../data/myflux.py", "w+");
            fprintf(fluxData_py, "myflux=[ \n");
            for( int j = 0 ; j < m_TriGrid[i].size(); j++ )
            {
                double flux = m_TriGrid[i][j].getValue();
                fprintf(fluxData_py, "%.4f ",flux);
                if( j != m_TriGrid[i].size() -1 )
                {
                    fprintf(fluxData_py, ",");
                }
                fprintf(fluxData_py, "\n");
                
            }
            
            fprintf(fluxData_py, "];\n");
            if(fluxData_py != NULL)
            {
                fclose(fluxData_py);
                fluxData_py = NULL;
            }
        //}
        
    }
    
    
    void GEarthRadiationModel::createGrid()
    {
        //generating the triangle grid data
        GIcosahedron<GVertex> myico;
        std::vector< std::vector<GVertex> > allVer; // all the vertex
        m_TriGrid =  myico.createGeometry( m_level,allVer);
    }
    
    
    
    
    void GEarthRadiationModel::attachFlux()
    {
        int imonth = 2; // 0 -11, Jan to Dec
        FILE*  polygonFile = fopen("../data/mypolygon.py","w+");
        fprintf( polygonFile, "Points=[\n");
        // I need to judege which grid belonds to a triangle, and add up all the flux.
        for( int i = 0; i < m_level+1 ; i++ )
        {
            for( int j = 0 ; j < m_TriGrid[i].size(); j++ )
            {
                GVertex p[3];
                for( int k = 0 ; k< 3; k++ )
                {
                    p[k] = m_TriGrid[i][j].getPoint(k);
                }
                
                double Ra = 6371000;//6408137;     // a = 6378137 for WGS84
                double Rb = 6371000;//6386651.7;
                double perimeter = 0.0;
                double myarea = 0.0;
                GeodesicExact geod(Ra, (Ra-Rb)/Ra);
                //Alternatively: const Geodesic& geod = Geodesic::WGS84();
                PolygonAreaExact poly(geod);
                
                double flux =  countingStars(poly ,p, m_monToaData[imonth],polygonFile);
                m_TriGrid[i][j].setValue(flux);
        }
    }
        
    fprintf(polygonFile, "];\n");
    if( polygonFile != NULL )
    {
        fclose(polygonFile);
        polygonFile = NULL;
    }

}
    
    /* lat1A lon1A , lat1B , lon1B should be the scanning line
     from the start point to the end point, the arc should be in counterclockwise order
     So the intersection point should be between the start point and end point whatever
     in latitude and longitude!
     
     A is for start point , B is for end point
     */
bool  GEarthRadiationModel::SphereIntersection(double lat1a, double lon1a, double lat1b,
                            double lon1b, double lat2a, double lon2a,
                            double lat2b, double lon2b, double& lat3A,
                            double& lon3A )
    {
        
        double lat1A= lat1a;
        double lat1B= lat1b;
        double lat2A= lat2a;
        double lat2B= lat2b;
        
        double lon1A= lon1a;
        double lon1B= lon1b;
        double lon2A= lon2a;
        double lon2B= lon2b;
        
        double pi = GCONST("PI");
        double de2ra = pi/180.0;
        double eps = 1E-8;
        bool isInside = false;
        double radius = sphere_radius; //km
        double v1[3], v2[3], v3a[3], v3b[3],v3[3], n1[3], n2[3];
        double m=0.0,tmp =0.0;
        double lat3B;
        double lon3B;
        double dd[2]={0};
        bool test1 = false, test2 = false;
        double d1 = GCDistance(lat1A, lon1A, lat1B, lon1B,radius);
        double d2 = GCDistance(lat2A, lon2A, lat2B, lon2B,radius);
        
        //for the scanning line 1
        if( fabs(fabs(lat1A)-90)<eps_angle &&fabs(fabs(lat1B)-90)<eps_angle && fabs(lon1A-lon1B)<eps_angle  )
        {
            if( lat1A <= lat1B )  // start point < end point
            {
                lat1B = 0;
                lon1B = lon1B+180; // change the end point
                
                //             if(lon1B >= 180)
                //             {
                //               lon1B = lon1B-180;  // change the end point
                //             }
                //             else
                //             {
                //                 lon1B = lon1B+180; // change the end point
                //             }
                if(lon1B>=360)
                {
                    lon1B = lon1B -360;
                }
            }
            else if(lat1B<=lat1A) // end point <= start point
            {
                lat1B = 0;
                //lon1B = 0;  // always change the end point
            }
            test1 = true;
        }
        else if(fabs(fabs(lat2A)-90)<eps_angle &&fabs(fabs(lat2B)-90)<eps_angle && fabs(lon2A-lon2B)<eps_angle)
        {
            if( lat2A <= lat2B )  // start point < end point
            {
                lat2B = 0;  // change the end point
                lon2B = lon2B+180;
                if( lon2B>=360 )
                {
                    lon2B = lon2B -360;
                }
            }
            else if(lat2B<=lat2A) // end point <= start point
            {
                lat2B = 0;  // always change the end point
            }
            test2 = true;
        }
        
        //
        // for path 1, setting up my 2 vectors, v1 is vector
        // from center of the Earth to point A, v2 is vector
        // from center of the Earth to point B.
        //
        v1[0] = cos(lat1A*de2ra)*cos(lon1A*de2ra)*radius;
        v1[1] = cos(lat1A*de2ra)*sin(lon1A*de2ra)*radius;
        v1[2] = sin(lat1A*de2ra)*radius;
        
        v2[0] = cos(lat1B*de2ra)*cos(lon1B*de2ra)*radius;
        v2[1] = cos(lat1B*de2ra)*sin(lon1B*de2ra)*radius;
        v2[2] = sin(lat1B*de2ra)*radius;
        
        //
        // n1 is the normal to the plane formed by the three points:
        // center of the Earth, point 1A, and point 1B
        //
        bool testLat1 = false, testLat2= false;
        n1[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
        n1[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
        n1[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
        tmp = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
        n1[0] = n1[0]/tmp;n1[1] = n1[1]/tmp;n1[2] = n1[2]/tmp;
        // if the latitudes are equal, there are tow situations
        // One is that the latitude lines , Second is the big circle goes through the polar point
        // Luckily, Only situation one exist in this program.
        //    if( fabs(lat1a-lat1b)<eps)  // for the latitudes
        //    {
        //        n1[0] = 0.0; n1[1]=0.0;
        //        if(lon1a <= lon1b)
        //        {
        //            n1[2]=1.0;
        //        }
        //        else
        //        {
        //            n1[2]= -1.0;
        //        }
        //
        //        testLat1 = true;
        //    }
        
        dd[0] = -(v1[0]*n1[0] + v1[1]*n1[1] + v1[2]*n1[2]);
        
        //
        // for path 2, setting up my 2 vectors, v1 is vector
        // from center of the Earth to point A, v2 is vector
        // from center of the Earth to point B.
        //
        v1[0] = cos(lat2A * de2ra) * cos(lon2A * de2ra)*radius;
        v1[1] = cos(lat2A * de2ra) * sin(lon2A * de2ra)*radius;
        v1[2] = sin(lat2A * de2ra)*radius;
        
        v2[0] = cos(lat2B * de2ra) * cos(lon2B * de2ra)*radius;
        v2[1] = cos(lat2B * de2ra) * sin(lon2B * de2ra)*radius;
        v2[2] = sin(lat2B * de2ra)*radius;
        
        //
        // n2 is the normal to the plane formed by the three points:
        // center of the Earth, point 2A, and point 2B
        //
        n2[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
        n2[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
        n2[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
        tmp = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
        n2[0] = n2[0]/tmp;n2[1] = n2[1]/tmp;n2[2] = n2[2]/tmp;
        //    if( fabs(lat2a-lat2b)<eps)  // for the latitudes
        //    {
        //        n2[0] = 0.0; n2[1]=0.0;
        //        if(lon2a <= lon2b)
        //        {
        //            n2[2]=1.0;
        //        }
        //        else
        //        {
        //            n2[2]= -1.0;
        //        }
        //        testLat2 = true;
        //    }
        
        dd[1] = -(v1[0]*n2[0] + v1[1]*n2[1] + v1[2]*n2[2]);
        
        //check whether n1 is paralled with n2
        
        double A=0.0,B=0.0,C=0.0;
        //t1 = n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2];
        //t2 = n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2];
        //t3 = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];  // t3 can determined whether n1 is normal to n2
        double pos0[3]={0.0}; // the special point for the line，it is very troublesome to determine pos0
        v3[0] = (n2[1]*n1[2]) - (n1[1]*n2[2]);
        v3[1] = (n2[2]*n1[0]) - (n1[2]*n2[0]);
        v3[2] = (n1[1]*n2[0]) - (n2[1]*n1[0]);
        tmp = sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2]);
        v3[0] = v3[0]/tmp;v3[1] = v3[1]/tmp;v3[2] = v3[2]/tmp;
        
        if(fabs(v3[2])<eps)
        {
            if( fabs(fabs(n1[2])-1.0)<eps)
            {
                pos0[2] = -dd[0];
            }
            else if(fabs(fabs(n2[2])-1.0)<eps)
            {
                pos0[2] = -dd[1];
            }
        }
        else
        {
            pos0[2] = 0.0; // just assuming z component is zero
        }
        
        
        if(fabs(n1[0])<eps)
        {
            pos0[0] = 0.0;
            if(fabs(n1[1])<eps)
            {
                pos0[1] = 0.0;
            }
            else
            {
                pos0[1] = -(dd[0]+n1[2]*pos0[2])/n1[1];
            }
        }
        else if(fabs(n2[0])<eps)
        {
            pos0[0] = 0.0;
            if(fabs(n2[1])<eps)
            {
                pos0[1] = 0.0;
            }
            else
            {
                pos0[1] = -(dd[1]+n2[2]*pos0[2])/n2[1];
            }
        }
        else
        {
            if(fabs(v3[2])<eps)
            {
                pos0[1] = 0.0;
                pos0[0] = n2[0]*(-dd[1]-n2[2]*pos0[2])/n1[0]/n2[0];
            }
            else
            {
                pos0[1] = -(n2[0]*(-dd[0]-n1[2]*pos0[2]) - n1[0]*(-dd[1]-n2[2]*pos0[2]))/v3[2];
                pos0[0] = (n1[0]*(-dd[1]-n2[2]*pos0[2])-n1[0]*n2[1]*pos0[1])/n1[0]/n2[0];
            }
        }
        
        //if treat the earth as an sphere, then the pos0 can be the origin point
        // that is:
        pos0[0] = 0.0;pos0[1] = 0.0;pos0[2] = 0.0;
        
        
        double mya =1.0, myb=1.0;
        A = (v3[0]*v3[0]+v3[1]*v3[1])/mya/mya + v3[2]*v3[2]/myb/myb;
        B = 2*( (v3[0]*pos0[0]+v3[1]*pos0[1])/mya/mya + v3[2]*pos0[2]/myb/myb);
        C = (pos0[0]*pos0[0]+pos0[1]*pos0[1])/mya/mya + pos0[2]*pos0[2]/myb/myb-1;
        double test = (B*B-4*A*C)*10000.0;
        double x1 =0.0, x2=0.0;
        if( test >= 0 )
        {
            x1 = (-B+sqrt(B*B-4*A*C))/2.0/A;
            x2 = (-B-sqrt(B*B-4*A*C))/2.0/A;
            
            // choose the right solution
            // the rule is the intersection point should be between the 2 sections
            
            v3a[0] = v3[0]*x1 + pos0[0];
            v3a[1] = v3[1]*x1 + pos0[1];
            v3a[2] = v3[2]*x1 + pos0[2];
            
            v3b[0] = v3[0]*x2 + pos0[0];
            v3b[1] = v3[1]*x2 + pos0[1];
            v3b[2] = v3[2]*x2 + pos0[2];
            m = sqrt(v3a[0]*v3a[0]+v3a[1]*v3a[1]+v3a[2]*v3a[2]  );
            lat3A = asin(v3a[2]/m);
            lon3A = atan2(v3a[1],v3a[0]);
            
            m = sqrt(v3b[0]*v3b[0]+v3b[1]*v3b[1]+v3b[2]*v3b[2]  );
            lat3B = asin(v3b[2]/m);
            lon3B = atan2(v3b[1],v3b[0]);
            
            lat3A *= (1.0/de2ra);
            lon3A *= (1.0/de2ra);
            lat3B *= (1.0/de2ra);
            lon3B *= (1.0/de2ra);
            
            if(lon3B < 0.0 ) { lon3B = lon3B+360;}
            if(lon3A < 0.0 ) { lon3A = lon3A+360;}
            
            // However, for latitude, the start point should be larger than end poind
            lat1A= lat1a; lon1A= lon1a;
            lat1B= lat1b; lon1B= lon1b;
            lat2A= lat2a; lon2A= lon2a;
            lat2B= lat2b; lon2B= lon2b;
            
            double d1a3a = GCDistance(lat1A, lon1A, lat3A, lon3A,radius);
            double d1b3a = GCDistance(lat1B, lon1B, lat3A, lon3A,radius);
            double d2a3a = GCDistance(lat2A, lon2A, lat3A, lon3A,radius);
            double d2b3a = GCDistance(lat2B, lon2B, lat3A, lon3A,radius);
            
            double d1a3b = GCDistance(lat1A, lon1A, lat3B, lon3B,radius);
            double d1b3b = GCDistance(lat1B, lon1B, lat3B, lon3B,radius);
            double d2a3b = GCDistance(lat2A, lon2A, lat3B, lon3B,radius);
            double d2b3b = GCDistance(lat2B, lon2B, lat3B, lon3B,radius);
            
            double test1 = d1a3a + d1b3a-d1;
            double test2 = d2a3a + d2b3a-d2;
            double test3 = d1a3b + d1b3b-d1;
            double test4 = d2a3b + d2b3b-d2;
            if( fabs(test1)<eps_distance && fabs(test2)<eps_distance )
            {
                isInside = true;
            }
            else if( fabs(test3)<eps_distance && fabs(test4)<eps_distance)
            {
                m = lat3A;
                lat3A = lat3B;
                lat3B = m;
                m = lon3A;
                lon3A = lon3B;
                lon3B = m;
                isInside = true;
            }
        }
        
        return isInside;
    }

    
    
    
    /*treat it as sphere*/
bool GEarthRadiationModel::IsPointInTri_sph( GVertex p[3] , double lat, double lon)
    {
        double pi = GCONST("PI");
        bool test = false;
        //transfer lat and lon to xyz
        double xyz[3] ={0.0};
        //double a = 6408137;     // a = 6378137 for WGS84
        //double b = 6386651.7;   // b=6356752.3142 for WGS84
        //double m_radius = 6371000;
        double eps = 1.0E-10;
        
        double blh[3]={lat/180.0*pi,lon/180.0*pi,0.0};
        //p[0].blh2xyz(blh, xyz, 0.0, m_radius);
        xyz[2] = sin(blh[0])*sphere_radius;
        xyz[1] = cos(blh[0])*sin(blh[1])*sphere_radius;
        xyz[0] = cos(blh[0])*cos(blh[1])*sphere_radius;
        
        double len = sqrt( pow(xyz[0],2.0) + pow(xyz[1],2.0)+pow(xyz[2],2.0) );
        xyz[0] = xyz[0]/len;
        xyz[1] = xyz[1]/len;
        xyz[2] = xyz[2]/len;
        
        double v[3][3]={{0.0}};
        for( int i = 0 ; i< 3; i++ )
        {
            double len = sqrt( pow(p[i].getX(),2.0)+pow(p[i].getY(),2.0)+pow(p[i].getZ(),2.0) );
            v[i][0] = p[i].getX()/len;
            v[i][1] = p[i].getY()/len;
            v[i][2] = p[i].getZ()/len;
        }
        
        double n[3][3] ={{0.0}};
        for( int i = 0 ; i< 3 ; i++ )
        {
            int j =i+1;
            if( j==3 ) {j = 0;}
            n[i][0] = v[i][1]*v[j][2] - v[i][2]*v[j][1];
            n[i][1] = v[i][2]*v[j][0] - v[j][2]*v[i][0];
            n[i][2] = v[i][0]*v[j][1] - v[i][1]*v[j][0];
            
            double sum = sqrt(n[i][0]*n[i][0] + n[i][1]*n[i][1]+ n[i][2]*n[i][2]);
            n[i][0] = n[i][0]/sum;
            n[i][1] = n[i][1]/sum;
            n[i][2] = n[i][2]/sum;
            
        }
        
        double test1 = xyz[0]*n[0][0]+xyz[1]*n[0][1]+xyz[2]*n[0][2];
        double test2 = xyz[0]*n[1][0]+xyz[1]*n[1][1]+xyz[2]*n[1][2];
        double test3 = xyz[0]*n[2][0]+xyz[1]*n[2][1]+xyz[2]*n[2][2];
        
        bool t1 = ( test1 >0 || fabs(test1)<eps );
        bool t2 = ( test2 >0 || fabs(test2)<eps );
        bool t3 = ( test3 >0 || fabs(test3)<eps );
        
        test = t1&&t2&&t3;
        
        
        return test;
    }
    
    
    
double GEarthRadiationModel::countingStars(PolygonAreaExact poly, gfc::GVertex *p, gfc::GEarthRadiationModel::FLUXDATA fluxdata, FILE *polygonFile)
    {
        //process 270 and 0 problem
        if( fabs(p[1].getLon()-270)<eps_angle && fabs(p[2].getLon())<eps_angle )
        {
            if(fabs(fabs(p[2].getLat())-90)>eps_angle)  // not the polar
            {
                p[2].setBL(360, p[2].getLat());
            }
        }
        if( fabs(p[2].getLon()-270)<eps_angle && fabs(p[1].getLon())<eps_angle )
        {
            if(fabs(fabs(p[1].getLat())-90)>eps_angle)  // not the polar
            {
                p[1].setBL(360, p[1].getLat());
            }
        }
        if( fabs(p[0].getLon()-270)<eps_angle && fabs(p[1].getLon())<eps_angle )
        {
            if(fabs(fabs(p[1].getLat())-90)>eps_angle)  // not the polar
            {
                p[1].setBL(360, p[1].getLat());
            }
        }
        if( fabs(p[1].getLon()-270)<eps_angle && fabs(p[0].getLon())<eps_angle )
        {
            if(fabs(fabs(p[0].getLat())-90)>eps_angle)  // not the polar
            {
                p[0].setBL(360, p[0].getLat());
            }
        }
        if( fabs(p[0].getLon()-270)<eps_angle && fabs(p[2].getLon())<eps_angle)
        {
            if(fabs(fabs(p[2].getLat())-90)>eps_angle)  // not the polar
            {
                p[2].setBL(360, p[2].getLat());
            }
        }
        if( fabs(p[2].getLon()-270)<eps_angle && fabs(p[0].getLon())<eps_angle )
        {
            if(fabs(fabs(p[0].getLat())-90)>eps_angle)  // not the polar
            {
                p[0].setBL(360, p[0].getLat());
            }
        }
        
        //process the 360 and 0 problems
        if(  p[0].getLon() - 270 > eps_angle  )
        {
            if( fabs(p[1].getLon())<eps_angle)
            {
                p[1].setBL(360, p[1].getLat());
            }
            if( fabs(p[2].getLon())<eps_angle)
            {
                p[2].setBL(360, p[2].getLat());
            }
        }
        
        if(  p[1].getLon() - 270 > eps_angle  )
        {
            if( fabs(p[2].getLon())<eps_angle)
            {
                p[2].setBL(360, p[2].getLat());
            }
            if( fabs(p[0].getLon())<eps_angle)
            {
                p[0].setBL(360, p[0].getLat());
            }
        }
        
        if(  p[2].getLon() - 270 > eps_angle  )
        {
            if( fabs(p[1].getLon())<eps_angle)
            {
                p[1].setBL(360, p[1].getLat());
            }
            if( fabs(p[0].getLon())<eps_angle)
            {
                p[0].setBL(360, p[0].getLat());
            }
        }
        
        
        // double Ra = 6371000;//6408137;     // a = 6378137 for WGS84
        // double Rb = 6371000;//6386651.7;
        double perimeter = 0.0;
        // Geodesic geod(Ra, (Ra-Rb)/Ra);
        
        //Alternatively: const Geodesic& geod = Geodesic::WGS84();
        // PolygonArea poly(geod);
        
        // just choose the minimum
        int miniIndex = 0;    // the minimum one
        int middleIndex = -1; // the mididle one
        int maxIndex = -1;
        //param coef[3]={0};
        double tmpLon[3] ={0.0};
        int   tmpindex[2]={-1};
        for( int i = 0 ; i< 3; i++ )
        {
            tmpLon[i] = p[i].getLon();
            if( p[i].getLon() < p[miniIndex].getLon() )
            {
                miniIndex = i;
            }
        }
        
        //miniIndex should not be the polar
        if( (fabs(fabs(p[miniIndex].getLat())-90))<eps_angle )
        {
            if( miniIndex == 0 )
            {
                miniIndex = p[1].getLon()<=p[2].getLon()?1:2;
            }
            else if(miniIndex == 1)
            {
                miniIndex = p[0].getLon()<=p[2].getLon()?0:2;
            }
            else if(miniIndex == 2)
            {
                miniIndex = p[0].getLon()<=p[1].getLon()?0:1;
            }
        }
        
        
        for( int i = 0,j =0  ; i< 3; i++ )
        {
            if( i != miniIndex )
            {
                tmpindex[j] = i;
                j++;
            }
        }
        
        middleIndex = p[tmpindex[0]].getLon()< p[tmpindex[1]].getLon() ? tmpindex[0] : tmpindex[1];
        
        for( int k = 0 ; k < 3; k++ )
        {
            if( k != miniIndex && k!= middleIndex)
            {
                maxIndex = k;
                break;
            }
        }
        
        // starting the counting process
        double min[2] = {-89.5,0.5}; // lat lon
        double max[2] = {89.5, 359.5}; // lat lon
        double interval[2] = {1.0,1.0}; // lat lon
        int pos[2] ={-1};
        double xy[2] ={0.0};
        double iter = 0.0;
        double starter[2]={0.0};
        //double eps = 1.0E-10;
        determinGrid(p[miniIndex].getLat(), p[miniIndex].getLon(),interval, min, max, starter);
        
        double starterLon = starter[1];
        double starterLat = starter[0];
        
        iter = starterLon;
        double area = 0.0;
        while(1)
        {
            if(iter >360.0 )
            {
                break;
            }
            
            if( iter-interval[1]  > p[maxIndex].getLon() )
            {
                break;
            }
            
            double up = 999.0, down = -999.0;  // for latitude
            
            // according to the righthand law, the start point scanlineP1 should be
            // arc one
            double scanlineP1[2] ={ -90, iter};
            double scanlineP2[2] ={  90, iter};
            double scanlineP3[2] ={-90,iter-interval[1]};
            double scanlineP4[2] ={90,iter-interval[1]};
            // arc two
            //double scanlinePP1[2] ={ 90, iter};
            //double scanlinePP2[2] ={ -90, iter};
            
            // from the start point to the end point , should be in counterclockwise order
            double pointC[2] ={0}, pointD[2] = {0};
            double intersectionPoint[2]={0.0};
            std::vector<double> tmpvalue;
            for( int g = 0 ; g< 3; g++ )
            {
                int k = g + 1;
                if( k == 3 ) { k = 0; }
                pointC[0] = p[g].getLat(); pointC[1] = p[g].getLon();
                pointD[0] = p[k].getLat(); pointD[1] = p[k].getLon();
                //bool test = intersection_2D(scanlineP1, scanlineP2, pointC, pointD, intersectionPoint);
                bool test1 = SphereIntersection(    scanlineP1[0], scanlineP1[1], scanlineP2[0],
                                                scanlineP2[1], pointC[0], pointC[1],
                                                pointD[0], pointD[1], intersectionPoint[0],
                                                intersectionPoint[1]);
                if( test1 == true )
                {
                    tmpvalue.push_back(intersectionPoint[0]);
                }
                
                bool test2 = SphereIntersection(     scanlineP2[0], scanlineP2[1], scanlineP1[0],
                                                scanlineP1[1], pointC[0], pointC[1],
                                                pointD[0], pointD[1], intersectionPoint[0],
                                                intersectionPoint[1]);
                if( test2 == true )
                {
                    tmpvalue.push_back(intersectionPoint[0]);
                }
                
                bool test3 = SphereIntersection(     scanlineP3[0], scanlineP3[1], scanlineP4[0],
                                                scanlineP4[1], pointC[0], pointC[1],
                                                pointD[0], pointD[1], intersectionPoint[0],
                                                intersectionPoint[1]);
                if( test3 == true )
                {
                    tmpvalue.push_back(intersectionPoint[0]);
                }
                
            }
            
            // that means the triangle is smaller than the grid
            if( tmpvalue.size() == 0 )
            {
                poly.Clear();
                double centerPoint[2]={0.0};
                for( int i = 0 ; i< 3; i++ )
                {
                    centerPoint[0] += p[i].getLat();
                    centerPoint[1] += p[i].getLon();
                    poly.AddPoint(p[i].getLat(), p[i].getLon());
                }
                
                centerPoint[0] = centerPoint[0]/3;
                centerPoint[1] = centerPoint[1]/3;
                
                double perimeter = 0.0,area1,area2;
                poly.Compute(false, true, perimeter, area1);
                poly.Clear();
                double starter[2] ={0.0};
                determinGrid(centerPoint[0], centerPoint[1],interval, min, max, starter);
                
                std::vector<std::vector<double> >  mypolypoints;
                std::vector<double>  tmp;
                tmp.push_back(starter[0]);tmp.push_back(starter[1]);
                mypolypoints.push_back(tmp);
                tmp.clear();
                tmp.push_back(starter[0]+interval[0]);tmp.push_back(starter[1]);
                mypolypoints.push_back(tmp);
                tmp.clear();
                tmp.push_back(starter[0]+interval[0]);tmp.push_back(starter[1]-interval[1]);
                mypolypoints.push_back(tmp);
                tmp.clear();
                tmp.push_back(starter[0]);tmp.push_back(starter[1]-interval[1]);
                mypolypoints.push_back(tmp);
                
                mypolypoints =  polygonPointsProcess(mypolypoints);
                
                ClockwiseSortPoints(mypolypoints);
                for(int i = 0 ; i< mypolypoints.size();i++)
                {
                    poly.AddPoint(mypolypoints[i][0], mypolypoints[i][1]);
                }
                
                poly.Compute(false, true, perimeter, area2);
                poly.Clear();
                
                int gridpos[2]={0};
                
                gridIndex( min, max,interval, centerPoint, gridpos  );
                
                double ratio = area1/area2;
                
                //area = area1;
                //area = ratio* fluxdata.data[gridpos[0]][gridpos[1]];
                area = 0.0;
                
                break;
            }
            
            for( int i = 0 ; i< 3; i++ )
            {
                if( p[i].getLon()>= iter-interval[1] && p[i].getLon()<iter)
                {
                    tmpvalue.push_back(p[i].getLat());
                }
            }
            
            up =  *max_element(tmpvalue.begin(),tmpvalue.end());
            down = *min_element(tmpvalue.begin(),tmpvalue.end());
            
            double upStarter = ( up/interval[0] )*interval[0];
            double downStarter = ( floor( down/interval[0] ))*interval[0];
            if(fabs(fabs(downStarter-down)-interval[0])<eps_angle)
            {
                downStarter = downStarter + interval[0];
            }
            
            double area1 = 0.0;
            double latlooper = downStarter;
            
            printf("latitude domain: downstarter: %f  upstarter: %f\n", downStarter,upStarter);
            for( ; latlooper < upStarter;  latlooper = latlooper + interval[0] )
            {
                double point1[2] ={latlooper,iter};
                double point2[2] ={latlooper+interval[0],iter};
                double point3[2] ={latlooper+interval[0],iter-interval[1]};
                double point4[2] ={latlooper,iter-interval[1] };
                std::vector<std::vector<double> >  allpoints;
                //check the 3 vertices of the triangle
                for(int r= 0; r<3; r++)
                {
                    if( p[r].getLat()>=point4[0] && p[r].getLat()<=point3[0]
                       &&p[r].getLon()>=point4[1] && p[r].getLon()<=point1[1]
                       )
                    {
                        std::vector<double> tmp;
                        tmp.push_back(p[r].getLat());
                        tmp.push_back(p[r].getLon());
                        allpoints.push_back(tmp);
                    }
                }
                //std::vector<std::vector<double> >  tmppoints;
                // point 1 is at the right down side
                
                // the start point situation
                if(  starterLat >  latlooper && starterLat <latlooper + interval[0] )
                {
                    
                }
                
                //bool test = IsPointInTri(p, point1[0], point1[1]);
                bool test = IsPointInTri_sph( p , point1[0], point1[1] );
                if( test == true )
                {
                    std::vector<double> myp;
                    myp.push_back(point1[0]);
                    myp.push_back(point1[1]);
                    allpoints.push_back(myp);
                }
                for(int k = 0 ; k< 3 ; k++ )
                {
                    int t = k+1 ;
                    if( k == 2 ) {t = 0;}
                    double pointC[2] ={0}, pointD[2] = {0};
                    pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                    pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                    double section[2]={0.0};
                    //bool IsIntersection = intersection_2D(point1, point2, pointC, pointD, section);
                    bool IsIntersection = SphereIntersection( point2[0], point2[1], point1[0],
                                                             point1[1], pointC[0], pointC[1],
                                                             pointD[0], pointD[1], section[0],
                                                             section[1]);
                    if( IsIntersection == true )
                    {
                        std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                        allpoints.push_back(tmp);
                        //tmppoints.push_back(tmp);
                    }
                }
                
                
                //test = IsPointInTri(p, point2[0], point2[1]);
                test = IsPointInTri_sph( p , point2[0], point2[1] );
                if( test == true )
                {
                    std::vector<double> myp;
                    myp.push_back(point2[0]);
                    myp.push_back(point2[1]);
                    allpoints.push_back(myp);
                }
                
                for(int k = 0 ; k< 3 ; k++ )
                {
                    int t = k+1 ;
                    if( k == 2 ) {t = 0;}
                    double pointC[2] ={0}, pointD[2] = {0};
                    pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                    pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                    double section[2]={0.0};
                    //bool IsIntersection = intersection_2D(point2, point3, pointC, pointD, section);
                    bool IsIntersection = SphereIntersection( point3[0], point3[1], point2[0],
                                                             point2[1], pointC[0], pointC[1],
                                                             pointD[0], pointD[1], section[0],
                                                             section[1]);
                    if( IsIntersection == true )
                    {
                        std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                        allpoints.push_back(tmp);
                    }
                }
                
                // test = IsPointInTri(p, point3[0], point3[1]);
                test = IsPointInTri_sph( p , point3[0], point3[1] );
                if( test == true )
                {
                    std::vector<double> myp;
                    myp.push_back(point3[0]);
                    myp.push_back(point3[1]);
                    allpoints.push_back(myp);
                }
                
                for(int k = 0 ; k< 3 ; k++ )
                {
                    int t = k+1 ;
                    if( k == 2 ) {t = 0;}
                    double pointC[2] ={0}, pointD[2] = {0};
                    pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                    pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                    double section[2]={0.0};
                    //bool IsIntersection = intersection_2D(point3, point4, pointC, pointD, section);
                    bool IsIntersection = SphereIntersection( point3[0], point3[1], point4[0],
                                                             point4[1], pointC[0], pointC[1],
                                                             pointD[0], pointD[1], section[0],
                                                             section[1]);
                    if( IsIntersection == true )
                    {
                        std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                        allpoints.push_back(tmp);
                    }
                }
                
                //test = IsPointInTri(p, point4[0], point4[1]);
                test = IsPointInTri_sph( p , point4[0], point4[1] );
                if( test == true )
                {
                    std::vector<double> myp;
                    myp.push_back(point4[0]);
                    myp.push_back(point4[1]);
                    allpoints.push_back(myp);
                }
                
                for(int k = 0 ; k< 3 ; k++ )
                {
                    int t = k+1 ;
                    if( k == 2 ) {t = 0;}
                    double pointC[2] ={0}, pointD[2] = {0};
                    pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                    pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                    double section[2]={0.0};
                    bool IsIntersection = false ;
                    
                    //bool IsIntersection = intersection_2D(point4, point1, pointC, pointD, section);
                    IsIntersection = SphereIntersection( point4[0], point4[1], point1[0],
                                                        point1[1], pointC[0], pointC[1],
                                                        pointD[0], pointD[1], section[0],
                                                        section[1]);
                    if( IsIntersection == true )
                    {
                        std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                        allpoints.push_back(tmp);
                    }
                }
                
                std::vector<std::vector<double> > mypolypoints = polygonPointsProcess(allpoints);
                
                if( mypolypoints.size() <= 2 )
                {
                    continue;
                }
                
                ClockwiseSortPoints(mypolypoints);
                
                fprintf(polygonFile, "[");
                //output this polygon
                for( int t = 0 ; t< mypolypoints.size(); t++)
                {
                    // 输出 经度,纬度
                    fprintf(polygonFile, "[%8.5f,%8.5f] ",mypolypoints[t][1],mypolypoints[t][0]);
                    
                    if(t == mypolypoints.size() -1 )
                    {
                        fprintf(polygonFile, "],");
                    }
                    else
                    {
                        fprintf(polygonFile, ",");
                    }
                }
                //in order to form a closed polygon
                //fprintf(polygonFile, "[%8.5f,%8.5f], ",polypoints[0][1],polypoints[0][0]);
                fprintf(polygonFile, "\n");
                
                poly.Clear();
                // then calculate the area of this polygon
                double myarea1 = 0.0;
                for( int t = 0  ;t < mypolypoints.size() ;t++ )
                {
                    poly.AddPoint(mypolypoints[t][0], mypolypoints[t][1]);
                }
                poly.Compute(false, true, perimeter, myarea1);
                if( myarea1 < 0 )
                {
                    printf("warning: area1 is negative!!\n");
                }
                poly.Clear();
                poly.AddPoint(point1[0], point1[1]);poly.AddPoint(point2[0], point2[1]);poly.AddPoint(point3[0], point3[1]);poly.AddPoint(point4[0], point4[1]);
                double myarea2 =0.0;
                poly.Compute(false, true, perimeter, myarea2);
                if( myarea2 < 0 )
                {
                    printf("warning: area1 is negative!!\n");
                }
                poly.Clear();
                
                //myarea =  polygon_area(allpoints);
                // then get the flux data of this grid
                int gridpos[2]={0};
                gridIndex( min, max,interval, point1, gridpos  );
                double ratio = myarea1/myarea2;
                double certain_flux =fluxdata.data[gridpos[0]][gridpos[1]];
                double testflux = fluxdata.data[168][0];
                double myflux = ratio * certain_flux;
                //double myflux = myarea1;
                
                area1 = area1+ myflux;
                printf("Lat:%8.4f Lon%8.4f myflux:%.4f ratio: %.4f\n",point1[0],point1[1],myflux,ratio);
                int testc = 0;
                
                
            }
            
            printf("area1: %f\n",area1);
            area = area + area1;
            iter = iter + interval[1];// the upper boundary
        }
        
        return  area;
        
        
    }
    
    
    void GEarthRadiationModel::loadData()
    {
        m_ceresDataFile="../data/CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Subset_200503-201503.nc";
        m_varableName = "toa_sw_all_mon";
        
        GNetcdf mync(m_ceresDataFile,GNetcdf::NCMODE::MODE_NOWRITE);
        
        int ndim = mync.getNdim();
        std::vector<double> timeData = mync.getDimensionData("time");
        std::vector<double> lonData =  mync.getDimensionData("lon");
        std::vector<double> latData =  mync.getDimensionData("lat");
        int dim_time_len = timeData.size();
        int dim_lon_len =  lonData.size();
        int dim_lat_len =  latData.size();
        // first step get the average data monthly;
        int dataTag[12] ={0};
        
        size_t* start = new size_t[ndim];
        memset(start,0,sizeof(size_t)*ndim);
        size_t* count = new size_t[ndim];
        memset( count,0,sizeof(size_t)*ndim);
        count[0] = 1; count[1] = dim_lat_len; count[2] = dim_lon_len;
        
        double toaData[180][360]={{0}};
        //std::vector<FLUXDATA> totalData;
        m_monToaData.resize(12);
        //totalData.resize(12);
        CivilTime ct( 2000,3,1,0,0,0,"tsUTC" ); //starter time
        m_startTime.SetFromCivilTime(ct);
        
        for( int rec = 0 ; rec< dim_time_len; rec++ )
        {
            start[0] = rec;
            mync.getData(m_varableName, start, count,  &toaData[0][0] );
            GTime curtime ;
            curtime.SetData(TimeSystem::GetByName("tsUTC"), timeData[rec], 0, 0);
            curtime = m_startTime + curtime;
            JDTime jt = GTime::GTime2JDTime(curtime);
            CivilTime  myct = GTime::JDTime2CivilTime(jt);
            dataTag[myct.m_month-1]++;
            for( int i = 0 ; i< dim_lat_len; i++ )
            {
                for( int j =0; j< dim_lon_len; j++ )
                {
                    m_monToaData[myct.m_month-1].data[i][j] += toaData[i][j];
                }
            }
        }
        
        // get the monthly average data
        
        for( int k =0 ; k< m_monToaData.size() ; k++ )
        {
            for( int i = 0 ; i<dim_lon_len; i++  ) //360
            {
                for( int j = 0 ; j< dim_lat_len; j++ ) //180
                {
                    m_monToaData[k].data[j][i] = m_monToaData[k].data[j][i] / dataTag[k];
                }
            }
        }
    }
    
}









/*
 1 is for start point and 2 is for end point
 
 the range of longitude is -PI to PI
 */
double GCDistance(double lat1, double lon1, double lat2, double lon2,double radius)
{
    
   
    
    double pi = gfc::GCONST("PI");
    double delLon = fabs(lon2-lon1);
    if( fabs(delLon) >=180.0 )
    {
        delLon = 360 - delLon;
    }
    
    
    // the same point
    if( fabs(lat1-lat2)<eps_angle && fabs(lon1-lon2)<eps_angle)
    {
        return 0.0;
    }
    //lon1 = lon1 - 180.0;
    //lon2 = lon2 - 180.0;
    //总是计算优弧长
   
    double de2ra = pi/180.0;
    lat1 *= de2ra;
    lon1 *= de2ra;
    lat2 *= de2ra;
    lon2 *= de2ra;
    delLon*=de2ra;
    
    double d = 0.0;
    
    long double cosD = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(delLon);
    d =  ( radius*acos(cosD));
    
    return d;
}



/* lat1A lon1A , lat1B , lon1B should be the scanning line
 from the start point to the end point, the arc should be in counterclockwise order
 So the intersection point should be between the start point and end point whatever
 in latitude and longitude!
 
 A is for start point , B is for end point
 */
bool SphereIntersection(double lat1a, double lon1a, double lat1b,
                        double lon1b, double lat2a, double lon2a,
                        double lat2b, double lon2b, double& lat3A,
                        double& lon3A )
{
    
    double lat1A= lat1a;
    double lat1B= lat1b;
    double lat2A= lat2a;
    double lat2B= lat2b;
    
    double lon1A= lon1a;
    double lon1B= lon1b;
    double lon2A= lon2a;
    double lon2B= lon2b;
    
    double pi = GCONST("PI");
    double de2ra = pi/180.0;
    double eps = 1E-8;
    bool isInside = false;
    double radius = sphere_radius; //km
    double v1[3], v2[3], v3a[3], v3b[3],v3[3], n1[3], n2[3];
    double m=0.0,tmp =0.0;
    double lat3B;
    double lon3B;
    double dd[2]={0};
    bool test1 = false, test2 = false;
    double d1 = GCDistance(lat1A, lon1A, lat1B, lon1B,radius);
    double d2 = GCDistance(lat2A, lon2A, lat2B, lon2B,radius);
    
    //for the scanning line 1
    if( fabs(fabs(lat1A)-90)<eps_angle &&fabs(fabs(lat1B)-90)<eps_angle && fabs(lon1A-lon1B)<eps_angle  )
    {
        if( lat1A <= lat1B )  // start point < end point
        {
            lat1B = 0;
            lon1B = lon1B+180; // change the end point
            
            //             if(lon1B >= 180)
            //             {
            //               lon1B = lon1B-180;  // change the end point
            //             }
            //             else
            //             {
            //                 lon1B = lon1B+180; // change the end point
            //             }
            if(lon1B>=360)
            {
                lon1B = lon1B -360;
            }
        }
        else if(lat1B<=lat1A) // end point <= start point
        {
            lat1B = 0;
            //lon1B = 0;  // always change the end point
        }
        test1 = true;
    }
    else if(fabs(fabs(lat2A)-90)<eps_angle &&fabs(fabs(lat2B)-90)<eps_angle && fabs(lon2A-lon2B)<eps_angle)
    {
        if( lat2A <= lat2B )  // start point < end point
        {
            lat2B = 0;  // change the end point
            lon2B = lon2B+180;
            if( lon2B>=360 )
            {
                lon2B = lon2B -360;
            }
        }
        else if(lat2B<=lat2A) // end point <= start point
        {
            lat2B = 0;  // always change the end point
        }
        test2 = true;
    }
    
    //
    // for path 1, setting up my 2 vectors, v1 is vector
    // from center of the Earth to point A, v2 is vector
    // from center of the Earth to point B.
    //
    v1[0] = cos(lat1A*de2ra)*cos(lon1A*de2ra)*radius;
    v1[1] = cos(lat1A*de2ra)*sin(lon1A*de2ra)*radius;
    v1[2] = sin(lat1A*de2ra)*radius;
    
    v2[0] = cos(lat1B*de2ra)*cos(lon1B*de2ra)*radius;
    v2[1] = cos(lat1B*de2ra)*sin(lon1B*de2ra)*radius;
    v2[2] = sin(lat1B*de2ra)*radius;
    
    //
    // n1 is the normal to the plane formed by the three points:
    // center of the Earth, point 1A, and point 1B
    //
    bool testLat1 = false, testLat2= false;
    n1[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
    n1[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
    n1[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
    tmp = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
    n1[0] = n1[0]/tmp;n1[1] = n1[1]/tmp;n1[2] = n1[2]/tmp;
    // if the latitudes are equal, there are tow situations
    // One is that the latitude lines , Second is the big circle goes through the polar point
    // Luckily, Only situation one exist in this program.
    //    if( fabs(lat1a-lat1b)<eps)  // for the latitudes
    //    {
    //        n1[0] = 0.0; n1[1]=0.0;
    //        if(lon1a <= lon1b)
    //        {
    //            n1[2]=1.0;
    //        }
    //        else
    //        {
    //            n1[2]= -1.0;
    //        }
    //
    //        testLat1 = true;
    //    }
    
    dd[0] = -(v1[0]*n1[0] + v1[1]*n1[1] + v1[2]*n1[2]);
    
    //
    // for path 2, setting up my 2 vectors, v1 is vector
    // from center of the Earth to point A, v2 is vector
    // from center of the Earth to point B.
    //
    v1[0] = cos(lat2A * de2ra) * cos(lon2A * de2ra)*radius;
    v1[1] = cos(lat2A * de2ra) * sin(lon2A * de2ra)*radius;
    v1[2] = sin(lat2A * de2ra)*radius;
    
    v2[0] = cos(lat2B * de2ra) * cos(lon2B * de2ra)*radius;
    v2[1] = cos(lat2B * de2ra) * sin(lon2B * de2ra)*radius;
    v2[2] = sin(lat2B * de2ra)*radius;
    
    //
    // n2 is the normal to the plane formed by the three points:
    // center of the Earth, point 2A, and point 2B
    //
    n2[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
    n2[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
    n2[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
    tmp = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
    n2[0] = n2[0]/tmp;n2[1] = n2[1]/tmp;n2[2] = n2[2]/tmp;
    //    if( fabs(lat2a-lat2b)<eps)  // for the latitudes
    //    {
    //        n2[0] = 0.0; n2[1]=0.0;
    //        if(lon2a <= lon2b)
    //        {
    //            n2[2]=1.0;
    //        }
    //        else
    //        {
    //            n2[2]= -1.0;
    //        }
    //        testLat2 = true;
    //    }
    
    dd[1] = -(v1[0]*n2[0] + v1[1]*n2[1] + v1[2]*n2[2]);
    
    //check whether n1 is paralled with n2
    
    double A=0.0,B=0.0,C=0.0;
    //t1 = n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2];
    //t2 = n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2];
    //t3 = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];  // t3 can determined whether n1 is normal to n2
    double pos0[3]={0.0}; // the special point for the line，it is very troublesome to determine pos0
    v3[0] = (n2[1]*n1[2]) - (n1[1]*n2[2]);
    v3[1] = (n2[2]*n1[0]) - (n1[2]*n2[0]);
    v3[2] = (n1[1]*n2[0]) - (n2[1]*n1[0]);
    tmp = sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2]);
    v3[0] = v3[0]/tmp;v3[1] = v3[1]/tmp;v3[2] = v3[2]/tmp;
    
    if(fabs(v3[2])<eps)
    {
        if( fabs(fabs(n1[2])-1.0)<eps)
        {
            pos0[2] = -dd[0];
        }
        else if(fabs(fabs(n2[2])-1.0)<eps)
        {
            pos0[2] = -dd[1];
        }
    }
    else
    {
        pos0[2] = 0.0; // just assuming z component is zero
    }
    
    
    if(fabs(n1[0])<eps)
    {
        pos0[0] = 0.0;
        if(fabs(n1[1])<eps)
        {
            pos0[1] = 0.0;
        }
        else
        {
            pos0[1] = -(dd[0]+n1[2]*pos0[2])/n1[1];
        }
    }
    else if(fabs(n2[0])<eps)
    {
        pos0[0] = 0.0;
        if(fabs(n2[1])<eps)
        {
            pos0[1] = 0.0;
        }
        else
        {
            pos0[1] = -(dd[1]+n2[2]*pos0[2])/n2[1];
        }
    }
    else
    {
        if(fabs(v3[2])<eps)
        {
            pos0[1] = 0.0;
            pos0[0] = n2[0]*(-dd[1]-n2[2]*pos0[2])/n1[0]/n2[0];
        }
        else
        {
            pos0[1] = -(n2[0]*(-dd[0]-n1[2]*pos0[2]) - n1[0]*(-dd[1]-n2[2]*pos0[2]))/v3[2];
            pos0[0] = (n1[0]*(-dd[1]-n2[2]*pos0[2])-n1[0]*n2[1]*pos0[1])/n1[0]/n2[0];
        }
    }
    
    //if treat the earth as an sphere, then the pos0 can be the origin point
    // that is:
    pos0[0] = 0.0;pos0[1] = 0.0;pos0[2] = 0.0;
    
    
    double mya =1.0, myb=1.0;
    A = (v3[0]*v3[0]+v3[1]*v3[1])/mya/mya + v3[2]*v3[2]/myb/myb;
    B = 2*( (v3[0]*pos0[0]+v3[1]*pos0[1])/mya/mya + v3[2]*pos0[2]/myb/myb);
    C = (pos0[0]*pos0[0]+pos0[1]*pos0[1])/mya/mya + pos0[2]*pos0[2]/myb/myb-1;
    double test = (B*B-4*A*C)*10000.0;
    double x1 =0.0, x2=0.0;
    if( test >= 0 )
    {
        x1 = (-B+sqrt(B*B-4*A*C))/2.0/A;
        x2 = (-B-sqrt(B*B-4*A*C))/2.0/A;
        
        // choose the right solution
        // the rule is the intersection point should be between the 2 sections
        
        v3a[0] = v3[0]*x1 + pos0[0];
        v3a[1] = v3[1]*x1 + pos0[1];
        v3a[2] = v3[2]*x1 + pos0[2];
        
        v3b[0] = v3[0]*x2 + pos0[0];
        v3b[1] = v3[1]*x2 + pos0[1];
        v3b[2] = v3[2]*x2 + pos0[2];
        m = sqrt(v3a[0]*v3a[0]+v3a[1]*v3a[1]+v3a[2]*v3a[2]  );
        lat3A = asin(v3a[2]/m);
        lon3A = atan2(v3a[1],v3a[0]);
        
        m = sqrt(v3b[0]*v3b[0]+v3b[1]*v3b[1]+v3b[2]*v3b[2]  );
        lat3B = asin(v3b[2]/m);
        lon3B = atan2(v3b[1],v3b[0]);
        
        lat3A *= (1.0/de2ra);
        lon3A *= (1.0/de2ra);
        lat3B *= (1.0/de2ra);
        lon3B *= (1.0/de2ra);
        
        if(lon3B < 0.0 ) { lon3B = lon3B+360;}
        if(lon3A < 0.0 ) { lon3A = lon3A+360;}
        
        // However, for latitude, the start point should be larger than end poind
        lat1A= lat1a; lon1A= lon1a;
        lat1B= lat1b; lon1B= lon1b;
        lat2A= lat2a; lon2A= lon2a;
        lat2B= lat2b; lon2B= lon2b;
        
        double d1a3a = GCDistance(lat1A, lon1A, lat3A, lon3A,radius);
        double d1b3a = GCDistance(lat1B, lon1B, lat3A, lon3A,radius);
        double d2a3a = GCDistance(lat2A, lon2A, lat3A, lon3A,radius);
        double d2b3a = GCDistance(lat2B, lon2B, lat3A, lon3A,radius);
        
        double d1a3b = GCDistance(lat1A, lon1A, lat3B, lon3B,radius);
        double d1b3b = GCDistance(lat1B, lon1B, lat3B, lon3B,radius);
        double d2a3b = GCDistance(lat2A, lon2A, lat3B, lon3B,radius);
        double d2b3b = GCDistance(lat2B, lon2B, lat3B, lon3B,radius);
        
        double test1 = d1a3a + d1b3a-d1;
        double test2 = d2a3a + d2b3a-d2;
        double test3 = d1a3b + d1b3b-d1;
        double test4 = d2a3b + d2b3b-d2;
        if( fabs(test1)<eps_distance && fabs(test2)<eps_distance )
        {
            isInside = true;
        }
        else if( fabs(test3)<eps_distance && fabs(test4)<eps_distance)
        {
            m = lat3A;
            lat3A = lat3B;
            lat3B = m;
            m = lon3A;
            lon3A = lon3B;
            lon3B = m;
            isInside = true;
        }
    }
    
    return isInside;
}



/* compute the intersection point of two geodesic lines */
void GeodesicIntersection(double* lat, double* lon, GeographicLib::Geodesic geod, double* intersectionPoint)
{
    class vector3
    {
    public:
        double _x, _y, _z;
        vector3(double x, double y, double z = 1) throw()
        : _x(x)
        , _y(y)
        , _z(z) {}
        vector3 cross(const vector3& b) const throw()
        {
            return vector3(_y * b._z - _z * b._y,
                           _z * b._x - _x * b._z,
                           _x * b._y - _y * b._x);
        }
        void norm() throw()
        {
            _x /= _z;
            _y /= _z;
            _z = 1;
        }
    };
    
    
    const GeographicLib::Gnomonic gn(geod);//an equidistant azimuthal projection
    double lat0 =0.0, lon0 = 0.0; // these the initial results
    for( int i = 0; i< 4; i++ )
    {
        lat0+= lat[i];
        lon0+= lon[i];
    }
    
    lat0 = lat0/4;
    lon0 = lat0/4;
    
    for (int i = 0; i < 10; ++i) // the number of iteration is 10
    {
        double xa1, ya1, xa2, ya2;
        double xb1, yb1, xb2, yb2;
        gn.Forward(lat0, lon0, lat[0], lon[0], xa1, ya1);
        gn.Forward(lat0, lon0, lat[1], lon[1], xa2, ya2);
        gn.Forward(lat0, lon0, lat[2], lon[2], xb1, yb1);
        gn.Forward(lat0, lon0, lat[3], lon[3], xb2, yb2);
        // See Hartley and Zisserman, Multiple View Geometry, Sec. 2.2.1
        vector3 va1(xa1, ya1); vector3 va2(xa2, ya2);
        vector3 vb1(xb1, yb1); vector3 vb2(xb2, yb2);
        
        // la is homogeneous representation of line A1,A2
        // lb is homogeneous representation of line B1,B2
        vector3 la = va1.cross(va2);
        vector3 lb = vb1.cross(vb2);
        
        // p0 is homogeneous representation of intersection of la and lb
        vector3 p0 = la.cross(lb);
        p0.norm();
        
        double lat1, lon1;
        gn.Reverse(lat0, lon0, p0._x, p0._y, lat1, lon1);
        //std::cout << "Increment " << lat1-lat0 << " " << lon1-lon0 << "\n";
        lat0 = lat1;
        lon0 = lon1;
    }
    
    intersectionPoint[0] = lat0;
    intersectionPoint[1] = lon0;
    //std::cout  << "Final result " << lat0 << " " << lon0 << "\n";
}


/* get the index of current xy in the grid
 *  pay attention to the min and max; they should be the middle of the grid
 *
 * input:  min[0]-> xmin ;  min[1]->ymin
 *         max[0]-> xmax;   max[1]->ymax
 interval[0] -> xinterval; interval[1] -> yinterval
 xy[0]: latitude; xy[1] longtitude
 * output: pos[2]: the index of x and y;pos[0] for x. pos[1] for y
 */
void gridIndex( double min[2], double max[2],double interval[2], double xy[2], int* pos  )
{
    int tmp1 = ceil(xy[1]/interval[1])-1;
    // if( fabs(tmp1*interval[1]-xy[1])<eps_angle )
    //  {
    //      tmp1 = tmp1 -1 ;
    //  }
    //for lon
    //pos[1] = floor(xy[1]) -1 ;
    pos[1] = tmp1;
    
    //for lat
    double tmp = xy[0] - ( min[0] - interval[0]/2);
    int tmp2 = floor(tmp);
    pos[0] = tmp2;
    
//    if( xy[0] <= 0.0 )
//    {
//       if( fabs(tmp2-tmp)<eps_angle && fabs(tmp)>eps_angle )
//       {
//           pos[0] = tmp2 -1 ;
//       }
//    }
//    else if( xy[1] > 0.0 )
//    {
//        if( fabs(fabs(tmp)-fabs(tmp2)*interval[1] -interval[1])<eps_angle  )
//        {
//           pos[0] = tmp2+1;
//        }
//        else
//        {
//            pos[0] = tmp2;
//        }
//        
//    }
    
    if( pos[0] <0 || pos[0] >= 180 )
    {
        int testc = 0;
    }
    if( pos[1]<0 || pos[1]>=360 )
    {
        int testc = 0;
    }
    
    //pos[0] = floor(xy[0] -( min[0]-interval[0]/2));
    
    
}


/*treat it as sphere*/
bool IsPointInTri_sph( GVertex p[3] , double lat, double lon)
{
    double pi = GCONST("PI");
    bool test = false;
    //transfer lat and lon to xyz
    double xyz[3] ={0.0};
    //double a = 6408137;     // a = 6378137 for WGS84
    //double b = 6386651.7;   // b=6356752.3142 for WGS84
    //double m_radius = 6371000;
    double eps = 1.0E-10;
    
    double blh[3]={lat/180.0*pi,lon/180.0*pi,0.0};
    //p[0].blh2xyz(blh, xyz, 0.0, m_radius);
    xyz[2] = sin(blh[0])*sphere_radius;
    xyz[1] = cos(blh[0])*sin(blh[1])*sphere_radius;
    xyz[0] = cos(blh[0])*cos(blh[1])*sphere_radius;
    
    double len = sqrt( pow(xyz[0],2.0) + pow(xyz[1],2.0)+pow(xyz[2],2.0) );
    xyz[0] = xyz[0]/len;
    xyz[1] = xyz[1]/len;
    xyz[2] = xyz[2]/len;
    
    double v[3][3]={{0.0}};
    for( int i = 0 ; i< 3; i++ )
    {
        double len = sqrt( pow(p[i].getX(),2.0)+pow(p[i].getY(),2.0)+pow(p[i].getZ(),2.0) );
        v[i][0] = p[i].getX()/len;
        v[i][1] = p[i].getY()/len;
        v[i][2] = p[i].getZ()/len;
    }
    
    double n[3][3] ={{0.0}};
    for( int i = 0 ; i< 3 ; i++ )
    {
        int j =i+1;
        if( j==3 ) {j = 0;}
        n[i][0] = v[i][1]*v[j][2] - v[i][2]*v[j][1];
        n[i][1] = v[i][2]*v[j][0] - v[j][2]*v[i][0];
        n[i][2] = v[i][0]*v[j][1] - v[i][1]*v[j][0];
        //        if( fabs(p[i].getLat()-p[j].getLat())<0.00001 ) // the same latitude
        //        {
        //            n[i][0] = 0; n[i][1] = 0;
        //            n[i][2] = 1;
        //        }
        
        double sum = sqrt(n[i][0]*n[i][0] + n[i][1]*n[i][1]+ n[i][2]*n[i][2]);
        n[i][0] = n[i][0]/sum;
        n[i][1] = n[i][1]/sum;
        n[i][2] = n[i][2]/sum;
        
    }
    
    double test1 = xyz[0]*n[0][0]+xyz[1]*n[0][1]+xyz[2]*n[0][2];
    double test2 = xyz[0]*n[1][0]+xyz[1]*n[1][1]+xyz[2]*n[1][2];
    double test3 = xyz[0]*n[2][0]+xyz[1]*n[2][1]+xyz[2]*n[2][2];
    
    bool t1 = ( test1 >0 || fabs(test1)<eps );
    bool t2 = ( test2 >0 || fabs(test2)<eps );
    bool t3 = ( test3 >0 || fabs(test3)<eps );
    
    test = t1&&t2&&t3;
    
    
    return test;
}


// test whether a point is in the triangle
// treat it as plane
bool IsPointInTri( GVertex p[3] , double lat, double lon)
{
    double v0[2]= {0},v1[2]={0},v2[2]={0};
    //v0[0] = p3[1] -p1[1];
    //v0[1] = p3[0] -p1[0];
    
    v0[0] = p[2].getLon() -p[0].getLon();
    v0[1] = p[2].getLat() -p[0].getLat();
    
    v1[0] = p[1].getLon() -p[0].getLon();
    v1[1] = p[1].getLat() -p[0].getLat();
    
    //v1[0] = p2[1] -p1[1];
    //v1[1] = p2[0] -p1[0];
    
    
    v2[0] = lon - p[0].getLon();
    v2[1] = lat - p[0].getLat();
    
    double dot00 = v0[0]*v0[0] + v0[1]*v0[1];
    double dot01 = v0[0]*v1[0] + v0[1]*v1[1];
    double dot02 = v0[0]*v2[0] + v0[1]*v2[1];
    double dot11 = v1[0]*v1[0] + v1[1]*v1[1];
    double dot12 = v1[0]*v2[0] + v1[1]*v2[1];
    
    double inverDeno = 1.0/( dot00*dot11 - dot01*dot01);
    
    double u = (dot11*dot02 - dot01*dot12)*inverDeno;
    
    if( u<0 || u > 1 )
    {
        return false;
    }
    
    double v = (dot00*dot12 -dot01*dot02)*inverDeno;
    
    if( v< 0 || v > 1)
    {
        return false;
    }
    
    return u+v<=1;
}


/* compute the area of a ploygon
 * reference: http://wenku.baidu.com/link?url=Yf36bTpcX0-Swd8NL67aj-zUsghEb0QBwS9iemjkL68yDEldnpuQLykGnopP0U_BC0X6orgWZS4hjqzeJq5wSwKTAJZBskHiTYs6kkLN7Sq
 * points[k][0] : the x coordinate
 * points[k][1] : the y coordinate
 * just for convex polygon
 */
double polygon_area( std::vector< std::vector<double> >& points )
{
    double area = 0.0;
    int npoints = (int)points.size();
    int j = 0;
    for( int i = 0 ; i< npoints ; i++ )
    {
        j = i+ 1;
        if( i == npoints -1 )
        {
            j = 0;
        }
        
        area += points[i][0]*points[j][1] - points[j][0]*points[i][1];
    }
    
    area = fabs(area * 0.5);
    
    return area;
}


double getAngle(double* pA, double* pC)
{
    double test[2]={1,0},test1[2]={0}, test2[2]={0};
    double pi = GCONST("PI");
    double a= 0 ;
    for( int i = 0 ; i< 2; i++ )
    {
        test1[i] = pA[i] - pC[i];
    }
    
    //double b1 = sqrt( test1[0]*test1[0] + test1[1]*test1[1] );
    a = atan2( test1[1],test1[0] );
    if( a < 0.0 )
    {
        a = a + 2.0*pi;
    }
    
    return a;
}


void ClockwiseSortPoints(std::vector<std::vector<double> > &Points)
{
    if( Points.size() <= 0 )
    {
        return;
    }
    
    double pi = GCONST("PI");
    
    int count1=0; //>=0
    int count2=0; //<=0
    for(int i = 0 ; i< Points.size(); i++ )
    {
        if(Points[i][0]>=0.0) // the maxmum latitude
        {
            count1++;
        }
        if(Points[i][0]<=0.0)
        {
            count2++;
        }
    }
    
    bool testDirection = false;
    if(count1>count2)
    {
        testDirection = true;
    }
    
    
    std::vector<std::vector<double> > vPoints ;
    //project to xoy plane
    for( int i = 0 ; i< Points.size(); i++ )
    {
        //double z = sin(Points[i][0]*pi/180.0);
        double x = std::cos(Points[i][0]*pi/180.0)*std::cos(Points[i][1]*pi/180.0);
        double y = std::cos(Points[i][0]*pi/180.0)*std::sin(Points[i][1]*pi/180.0);
        std::vector<double> tmp;
        tmp.push_back(x);
        tmp.push_back(y);
        vPoints.push_back(tmp);
    }
    
    //计算重心
    double center[2]={0.0};
    double x = 0,y = 0;
    for (int i = 0;i < vPoints.size();i++)
    {
        x += vPoints[i][0];
        y += vPoints[i][1];
    }
    if( vPoints.size() > 0 )
    {
        center[0] = (double)(x/vPoints.size()); //x
        center[1] = (double)(y/vPoints.size()); //y
    }
    
    //double testData[6]={3,2,4,6,1,5};
    
    //buble sort
    for( int j = 0;j < vPoints.size();j++)
    {
        for ( int i = j;i < vPoints.size();i++)
        {
            //bool test = VertexCmp(&vPoints[j][0],&vPoints[j-1][0],center);
            double angel1 = getAngle(&vPoints[j][0],center);
            double angel2 = getAngle(&vPoints[i][0],center);
            bool test = false;
            if( testDirection==true)
            {
                test = angel1 > angel2 ;
            }
            else if(testDirection == false)
            {
                test = angel1 < angel2;
            }
            if ( test )
            {
                std::vector<double> tmp = vPoints[j];
                vPoints[j] = vPoints[i];
                vPoints[i] = tmp;
                // the same operation to Points
                std::vector<double> mytmp = Points[j];
                Points[j] = Points[i];
                Points[i] = mytmp;
            }
        }
    }
    
    //int testc = 0;
}


/**/
void determinGrid(double lat, double lon, double* interval, double* min, double* max,double* starter)
{
    int tag[2] ={0};
    if(lon <0 ) {lon = lon + 360;}
    
    if( fabs( lon - min[1]+interval[1]/2.0 )<eps_angle )  // 0
    {
        starter[1] = lon + interval[1];
        tag[1] = 1;
    }
    
    if( fabs( lon - max[1] - interval[1]/2.0 ) < eps_angle) // 360
    {
        starter[1] = lon;
        tag[1] = 1;
    }
    
    if( fabs( lat - min[0]+interval[0]/2.0 )<eps_angle )  //-90
    {
        starter[0] = lat;
        tag[0] = 1;
    }
    
    if( fabs( lat - max[0] - interval[0]/2.0 ) <eps_angle) //90
    {
        starter[0] = lat-interval[0];
        tag[0] = 1;
    }
    
    if( tag[0] == 0 )
    {
        int int_tmp = std::floor( lat/interval[0] );
        double tmp = int_tmp*interval[0]; // lat
        if( fabs(fabs(tmp-lat)-interval[0])<eps_angle )
        {
            starter[0] = tmp+interval[0];
        }
        else
        {
            starter[0] = tmp;
        }
        
        //starter[0] = tmp;
        
    }
    if( tag[1] == 0 )
    {
        double tmp = ceil( lon/interval[1] )*interval[1];  // lon
        if( fabs(tmp - lon) <eps_angle )
        {
            tmp = lon + interval[1];
        }
        //        if( lon > 180 )
        //        {
        //            lon = lon -360;
        //        }
        starter[1] = tmp;
    }
}




std::vector<std::vector<double> > polygonPointsProcess(std::vector<std::vector<double> >& allpoints)
{
    
    std::vector<std::vector<double> > polypoints;
    for(int i = 0 ; i< allpoints.size() ; i++ )
    {
        bool test = false;
        for( int j = 0 ; j<polypoints.size();j++ )
        {
            double t1 = allpoints[i][0]- polypoints[j][0]; //latitude
            double t2 = allpoints[i][1] - polypoints[j][1]; //longtitude
            if( (fabs(t1)<eps_angle && fabs(t2)<eps_angle)
               || (fabs(t1)<eps_angle && fabs(allpoints[i][1])<eps_angle&&fabs(allpoints[j][1]-360)<eps_angle)
               || (fabs(t1)<eps_angle && fabs(allpoints[j][1])<eps_angle&&fabs(allpoints[i][1]-360)<eps_angle)
               )
            {
                test = true;
                break;
            }
        }
        if(test == false)
        {
            polypoints.push_back(allpoints[i]);
        }
    }
    
    bool testP = false;
    bool testQ = false;
    std::vector<std::vector<double> > mypolypoints;
    for( int i = 0 ; i< polypoints.size(); i++ )
    {
        if( fabs(polypoints[i][0] - 90)<eps_angle )
        {
            if( testP == false)
            {
                mypolypoints.push_back(polypoints[i]);
                testP=true;
            }
        }
        else if( fabs(polypoints[i][0]+90)<eps_angle)
        {
            if( testQ == false)
            {
                mypolypoints.push_back(polypoints[i]);
                testQ=true;
            }
        }
        else
        {
            mypolypoints.push_back(polypoints[i]);
        }
    }
    
    //allpoints = mypolypoints;
    
    return mypolypoints;
}




/*choose the starter point
 *
 * the input triangle p[3] should be always in counterclockwise order
 *
 */
double countingStars(PolygonAreaExact poly, GVertex p[3] ,FLUXDATA fluxdata,FILE* polygonFile)
{
    //process 270 and 0 problem
    if( fabs(p[1].getLon()-270)<eps_angle && fabs(p[2].getLon())<eps_angle )
    {
        if(fabs(fabs(p[2].getLat())-90)>eps_angle)  // not the polar
        {
            p[2].setBL(360, p[2].getLat());
        }
    }
    if( fabs(p[2].getLon()-270)<eps_angle && fabs(p[1].getLon())<eps_angle )
    {
        if(fabs(fabs(p[1].getLat())-90)>eps_angle)  // not the polar
        {
            p[1].setBL(360, p[1].getLat());
        }
    }
    if( fabs(p[0].getLon()-270)<eps_angle && fabs(p[1].getLon())<eps_angle )
    {
        if(fabs(fabs(p[1].getLat())-90)>eps_angle)  // not the polar
        {
            p[1].setBL(360, p[1].getLat());
        }
    }
    if( fabs(p[1].getLon()-270)<eps_angle && fabs(p[0].getLon())<eps_angle )
    {
        if(fabs(fabs(p[0].getLat())-90)>eps_angle)  // not the polar
        {
            p[0].setBL(360, p[0].getLat());
        }
    }
    if( fabs(p[0].getLon()-270)<eps_angle && fabs(p[2].getLon())<eps_angle)
    {
        if(fabs(fabs(p[2].getLat())-90)>eps_angle)  // not the polar
        {
            p[2].setBL(360, p[2].getLat());
        }
    }
    if( fabs(p[2].getLon()-270)<eps_angle && fabs(p[0].getLon())<eps_angle )
    {
        if(fabs(fabs(p[0].getLat())-90)>eps_angle)  // not the polar
        {
            p[0].setBL(360, p[0].getLat());
        }
    }
    
    //process the 360 and 0 problems
    if(  p[0].getLon() - 270 > eps_angle  )
    {
        if( fabs(p[1].getLon())<eps_angle)
        {
            p[1].setBL(360, p[1].getLat());
        }
        if( fabs(p[2].getLon())<eps_angle)
        {
            p[2].setBL(360, p[2].getLat());
        }
    }
    
    if(  p[1].getLon() - 270 > eps_angle  )
    {
        if( fabs(p[2].getLon())<eps_angle)
        {
            p[2].setBL(360, p[2].getLat());
        }
        if( fabs(p[0].getLon())<eps_angle)
        {
            p[0].setBL(360, p[0].getLat());
        }
    }
    
    if(  p[2].getLon() - 270 > eps_angle  )
    {
        if( fabs(p[1].getLon())<eps_angle)
        {
            p[1].setBL(360, p[1].getLat());
        }
        if( fabs(p[0].getLon())<eps_angle)
        {
            p[0].setBL(360, p[0].getLat());
        }
    }
    
    
   // double Ra = 6371000;//6408137;     // a = 6378137 for WGS84
   // double Rb = 6371000;//6386651.7;
    double perimeter = 0.0;
   // Geodesic geod(Ra, (Ra-Rb)/Ra);
    
    //Alternatively: const Geodesic& geod = Geodesic::WGS84();
   // PolygonArea poly(geod);
    
    // just choose the minimum
    int miniIndex = 0;    // the minimum one
    int middleIndex = -1; // the mididle one
    int maxIndex = -1;
    //param coef[3]={0};
    double tmpLon[3] ={0.0};
    int   tmpindex[2]={-1};
    for( int i = 0 ; i< 3; i++ )
    {
        tmpLon[i] = p[i].getLon();
        if( p[i].getLon() < p[miniIndex].getLon() )
        {
            miniIndex = i;
        }
    }
    
    //miniIndex should not be the polar
    if( (fabs(fabs(p[miniIndex].getLat())-90))<eps_angle )
    {
        if( miniIndex == 0 )
        {
            miniIndex = p[1].getLon()<=p[2].getLon()?1:2;
        }
        else if(miniIndex == 1)
        {
            miniIndex = p[0].getLon()<=p[2].getLon()?0:2;
        }
        else if(miniIndex == 2)
        {
            miniIndex = p[0].getLon()<=p[1].getLon()?0:1;
        }
    }
    
    
    for( int i = 0,j =0  ; i< 3; i++ )
    {
        if( i != miniIndex )
        {
            tmpindex[j] = i;
            j++;
        }
    }
    
    middleIndex = p[tmpindex[0]].getLon()< p[tmpindex[1]].getLon() ? tmpindex[0] : tmpindex[1];
    
    for( int k = 0 ; k < 3; k++ )
    {
        if( k != miniIndex && k!= middleIndex)
        {
            maxIndex = k;
            break;
        }
    }
    
    // starting the counting process
    double min[2] = {-89.5,0.5}; // lat lon
    double max[2] = {89.5, 359.5}; // lat lon
    double interval[2] = {1.0,1.0}; // lat lon
    int pos[2] ={-1};
    double xy[2] ={0.0};
    double iter = 0.0;
    double starter[2]={0.0};
    //double eps = 1.0E-10;
    determinGrid(p[miniIndex].getLat(), p[miniIndex].getLon(),interval, min, max, starter);
    
    double starterLon = starter[1];
    double starterLat = starter[0];
    
    iter = starterLon;
    double area = 0.0;
    while(1)
    {
        if(iter >360.0 )
        {
            break;
        }
        
        if( iter-interval[1]  > p[maxIndex].getLon() )
        {
            break;
        }
        
        double up = 999.0, down = -999.0;  // for latitude
        
        // according to the righthand law, the start point scanlineP1 should be
        // arc one
        double scanlineP1[2] ={ -90, iter};
        double scanlineP2[2] ={  90, iter};
        double scanlineP3[2] ={-90,iter-interval[1]};
        double scanlineP4[2] ={90,iter-interval[1]};
        // arc two
        //double scanlinePP1[2] ={ 90, iter};
        //double scanlinePP2[2] ={ -90, iter};
        
        // from the start point to the end point , should be in counterclockwise order
        double pointC[2] ={0}, pointD[2] = {0};
        double intersectionPoint[2]={0.0};
        std::vector<double> tmpvalue;
        for( int g = 0 ; g< 3; g++ )
        {
            int k = g + 1;
            if( k == 3 ) { k = 0; }
            pointC[0] = p[g].getLat(); pointC[1] = p[g].getLon();
            pointD[0] = p[k].getLat(); pointD[1] = p[k].getLon();
            //bool test = intersection_2D(scanlineP1, scanlineP2, pointC, pointD, intersectionPoint);
            bool test1 = SphereIntersection(    scanlineP1[0], scanlineP1[1], scanlineP2[0],
                                            scanlineP2[1], pointC[0], pointC[1],
                                            pointD[0], pointD[1], intersectionPoint[0],
                                            intersectionPoint[1]);
            if( test1 == true )
            {
                tmpvalue.push_back(intersectionPoint[0]);
            }
            
            bool test2 = SphereIntersection(     scanlineP2[0], scanlineP2[1], scanlineP1[0],
                                            scanlineP1[1], pointC[0], pointC[1],
                                            pointD[0], pointD[1], intersectionPoint[0],
                                            intersectionPoint[1]);
            if( test2 == true )
            {
                tmpvalue.push_back(intersectionPoint[0]);
            }
            
            bool test3 = SphereIntersection(     scanlineP3[0], scanlineP3[1], scanlineP4[0],
                                            scanlineP4[1], pointC[0], pointC[1],
                                            pointD[0], pointD[1], intersectionPoint[0],
                                            intersectionPoint[1]);
            if( test3 == true )
            {
                tmpvalue.push_back(intersectionPoint[0]);
            }
            
        }
        
        // that means the triangle is smaller than the grid
        if( tmpvalue.size() == 0 )
        {
            poly.Clear();
            double centerPoint[2]={0.0};
            for( int i = 0 ; i< 3; i++ )
            {
                centerPoint[0] += p[i].getLat();
                centerPoint[1] += p[i].getLon();
                poly.AddPoint(p[i].getLat(), p[i].getLon());
            }
            
            centerPoint[0] = centerPoint[0]/3;
            centerPoint[1] = centerPoint[1]/3;
            
            double perimeter = 0.0,area1,area2;
            poly.Compute(false, true, perimeter, area1);
            poly.Clear();
            double starter[2] ={0.0};
            determinGrid(centerPoint[0], centerPoint[1],interval, min, max, starter);
            
            std::vector<std::vector<double> >  mypolypoints;
            std::vector<double>  tmp;
            tmp.push_back(starter[0]);tmp.push_back(starter[1]);
            mypolypoints.push_back(tmp);
            tmp.clear();
            tmp.push_back(starter[0]+interval[0]);tmp.push_back(starter[1]);
            mypolypoints.push_back(tmp);
            tmp.clear();
            tmp.push_back(starter[0]+interval[0]);tmp.push_back(starter[1]-interval[1]);
            mypolypoints.push_back(tmp);
            tmp.clear();
            tmp.push_back(starter[0]);tmp.push_back(starter[1]-interval[1]);
            mypolypoints.push_back(tmp);
            
            mypolypoints =  polygonPointsProcess(mypolypoints);
            
            ClockwiseSortPoints(mypolypoints);
            for(int i = 0 ; i< mypolypoints.size();i++)
            {
                poly.AddPoint(mypolypoints[i][0], mypolypoints[i][1]);
            }
            
            poly.Compute(false, true, perimeter, area2);
            poly.Clear();
            
            int gridpos[2]={0};
            
            gridIndex( min, max,interval, centerPoint, gridpos  );
            
            double ratio = area1/area2;
            
            //area = area1;
            //area = ratio* fluxdata.data[gridpos[0]][gridpos[1]];
            area = 0.0;
            
            break;
        }
        
        for( int i = 0 ; i< 3; i++ )
        {
            if( p[i].getLon()>= iter-interval[1] && p[i].getLon()<iter)
            {
                tmpvalue.push_back(p[i].getLat());
            }
        }
        
        up =  *max_element(tmpvalue.begin(),tmpvalue.end());
        down = *min_element(tmpvalue.begin(),tmpvalue.end());
        
        double upStarter = ( up/interval[0] )*interval[0];
        double downStarter = ( floor( down/interval[0] ))*interval[0];
        if(fabs(fabs(downStarter-down)-interval[0])<eps_angle)
        {
            downStarter = downStarter + interval[0];
        }
        
        double area1 = 0.0;
        double latlooper = downStarter;
        
        printf("latitude domain: downstarter: %f  upstarter: %f\n", downStarter,upStarter);
        for( ; latlooper < upStarter;  latlooper = latlooper + interval[0] )
        {
            double point1[2] ={latlooper,iter};
            double point2[2] ={latlooper+interval[0],iter};
            double point3[2] ={latlooper+interval[0],iter-interval[1]};
            double point4[2] ={latlooper,iter-interval[1] };
            std::vector<std::vector<double> >  allpoints;
            //check the 3 vertices of the triangle
            for(int r= 0; r<3; r++)
            {
                if( p[r].getLat()>=point4[0] && p[r].getLat()<=point3[0]
                   &&p[r].getLon()>=point4[1] && p[r].getLon()<=point1[1]
                   )
                {
                    std::vector<double> tmp;
                    tmp.push_back(p[r].getLat());
                    tmp.push_back(p[r].getLon());
                    allpoints.push_back(tmp);
                }
            }
            //std::vector<std::vector<double> >  tmppoints;
            // point 1 is at the right down side
            
            // the start point situation
            if(  starterLat >  latlooper && starterLat <latlooper + interval[0] )
            {
                
            }
            
            //bool test = IsPointInTri(p, point1[0], point1[1]);
            bool test = IsPointInTri_sph( p , point1[0], point1[1] );
            if( test == true )
            {
                std::vector<double> myp;
                myp.push_back(point1[0]);
                myp.push_back(point1[1]);
                allpoints.push_back(myp);
            }
            for(int k = 0 ; k< 3 ; k++ )
            {
                int t = k+1 ;
                if( k == 2 ) {t = 0;}
                double pointC[2] ={0}, pointD[2] = {0};
                pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                double section[2]={0.0};
                //bool IsIntersection = intersection_2D(point1, point2, pointC, pointD, section);
                bool IsIntersection = SphereIntersection( point2[0], point2[1], point1[0],
                                                         point1[1], pointC[0], pointC[1],
                                                         pointD[0], pointD[1], section[0],
                                                         section[1]);
                if( IsIntersection == true )
                {
                    std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                    allpoints.push_back(tmp);
                    //tmppoints.push_back(tmp);
                }
            }
            
            
            //test = IsPointInTri(p, point2[0], point2[1]);
            test = IsPointInTri_sph( p , point2[0], point2[1] );
            if( test == true )
            {
                std::vector<double> myp;
                myp.push_back(point2[0]);
                myp.push_back(point2[1]);
                allpoints.push_back(myp);
            }
            
            for(int k = 0 ; k< 3 ; k++ )
            {
                int t = k+1 ;
                if( k == 2 ) {t = 0;}
                double pointC[2] ={0}, pointD[2] = {0};
                pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                double section[2]={0.0};
                //bool IsIntersection = intersection_2D(point2, point3, pointC, pointD, section);
                bool IsIntersection = SphereIntersection( point3[0], point3[1], point2[0],
                                                         point2[1], pointC[0], pointC[1],
                                                         pointD[0], pointD[1], section[0],
                                                         section[1]);
                if( IsIntersection == true )
                {
                    std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                    allpoints.push_back(tmp);
                }
            }
            
            // test = IsPointInTri(p, point3[0], point3[1]);
            test = IsPointInTri_sph( p , point3[0], point3[1] );
            if( test == true )
            {
                std::vector<double> myp;
                myp.push_back(point3[0]);
                myp.push_back(point3[1]);
                allpoints.push_back(myp);
            }
            
            for(int k = 0 ; k< 3 ; k++ )
            {
                int t = k+1 ;
                if( k == 2 ) {t = 0;}
                double pointC[2] ={0}, pointD[2] = {0};
                pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                double section[2]={0.0};
                //bool IsIntersection = intersection_2D(point3, point4, pointC, pointD, section);
                bool IsIntersection = SphereIntersection( point3[0], point3[1], point4[0],
                                                         point4[1], pointC[0], pointC[1],
                                                         pointD[0], pointD[1], section[0],
                                                         section[1]);
                if( IsIntersection == true )
                {
                    std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                    allpoints.push_back(tmp);
                }
            }
            
            //test = IsPointInTri(p, point4[0], point4[1]);
            test = IsPointInTri_sph( p , point4[0], point4[1] );
            if( test == true )
            {
                std::vector<double> myp;
                myp.push_back(point4[0]);
                myp.push_back(point4[1]);
                allpoints.push_back(myp);
            }
            
            for(int k = 0 ; k< 3 ; k++ )
            {
                int t = k+1 ;
                if( k == 2 ) {t = 0;}
                double pointC[2] ={0}, pointD[2] = {0};
                pointC[0] = p[k].getLat(); pointC[1] = p[k].getLon();
                pointD[0] = p[t].getLat(); pointD[1] = p[t].getLon();
                double section[2]={0.0};
                bool IsIntersection = false ;
                
                //bool IsIntersection = intersection_2D(point4, point1, pointC, pointD, section);
                IsIntersection = SphereIntersection( point4[0], point4[1], point1[0],
                                                    point1[1], pointC[0], pointC[1],
                                                    pointD[0], pointD[1], section[0],
                                                    section[1]);
                if( IsIntersection == true )
                {
                    std::vector<double> tmp;  tmp.push_back(section[0]); tmp.push_back(section[1]);
                    allpoints.push_back(tmp);
                }
            }
            
            std::vector<std::vector<double> > mypolypoints = polygonPointsProcess(allpoints);
            
            if( mypolypoints.size() <= 2 )
            {
                continue;
            }
            
            ClockwiseSortPoints(mypolypoints);
            
            fprintf(polygonFile, "[");
            //output this polygon
            for( int t = 0 ; t< mypolypoints.size(); t++)
            {
                // 输出 经度,纬度
                fprintf(polygonFile, "[%8.5f,%8.5f] ",mypolypoints[t][1],mypolypoints[t][0]);
                
                if(t == mypolypoints.size() -1 )
                {
                    fprintf(polygonFile, "],");
                }
                else
                {
                    fprintf(polygonFile, ",");
                }
            }
            //in order to form a closed polygon
            //fprintf(polygonFile, "[%8.5f,%8.5f], ",polypoints[0][1],polypoints[0][0]);
            fprintf(polygonFile, "\n");
            
            poly.Clear();
            // then calculate the area of this polygon
            double myarea1 = 0.0;
            for( int t = 0  ;t < mypolypoints.size() ;t++ )
            {
                poly.AddPoint(mypolypoints[t][0], mypolypoints[t][1]);
            }
            poly.Compute(false, true, perimeter, myarea1);
            if( myarea1 < 0 )
            {
                printf("warning: area1 is negative!!\n");
            }
            poly.Clear();
            poly.AddPoint(point1[0], point1[1]);poly.AddPoint(point2[0], point2[1]);poly.AddPoint(point3[0], point3[1]);poly.AddPoint(point4[0], point4[1]);
            double myarea2 =0.0;
            poly.Compute(false, true, perimeter, myarea2);
            if( myarea2 < 0 )
            {
                printf("warning: area1 is negative!!\n");
            }
            poly.Clear();
            
            //myarea =  polygon_area(allpoints);
            // then get the flux data of this grid
            int gridpos[2]={0};
            gridIndex( min, max,interval, point1, gridpos  );
            double ratio = myarea1/myarea2;
            double certain_flux =fluxdata.data[gridpos[0]][gridpos[1]];
            double testflux = fluxdata.data[168][0];
            double myflux = ratio * certain_flux;
            //double myflux = myarea1;
            
            area1 = area1+ myflux;
            printf("Lat:%8.4f Lon%8.4f myflux:%.4f ratio: %.4f\n",point1[0],point1[1],myflux,ratio);
            int testc = 0;
            
            
        }
        
        printf("area1: %f\n",area1);
        area = area + area1;
        iter = iter + interval[1];// the upper boundary
    }
    
    return  area;
    
}






/*
 *
 * visiability test for one triangle
 *
 * test the 3 points of the triangle, if any one of the point is visialbe, then return true
 */
bool visiabilityTest( GTriangle<GVertex>& tri ,double* pos)
{
    double PI = GCONST("PI");
    bool test = false;
    double m[9]={0.0};
    double neu_sat[3]={0.0};// the neu coordinate of satellite relative to the station
    double tmp[3]={0.0};// the difference in coordinate between satellite and station
    double azel[2] = {0.0}; // azimuth and elevation of the satellite
    
    GVertex tpoints[3];
    tpoints[0] = tri.getPoint(0);tpoints[1] = tri.getPoint(1);tpoints[2] = tri.getPoint(2);
    
    // get the sub-satellite point first
    double len = sqrt( pos[0]*pos[0]+pos[1]*pos[1]);
    double subsatpoint[2]={0.0};  //星下点的经纬
    subsatpoint[0] = atan(pos[2]/sqrt(len))*180.0/PI;  // latitude
    subsatpoint[1] = atan2(pos[1],pos[0])*180.0/PI;  // longitude
    if( subsatpoint[1] < 0.0 )  { subsatpoint[1] = subsatpoint[1]+360.0; }
    
    bool test1 = IsPointInTri_sph(tpoints, subsatpoint[0], subsatpoint[1]);
    if( test1 == true )  // if the subsatpoint is inside the triangle, that means this triangle must be visiable
    {
        test =  true;
    }
    else if( test1 == false ) // 如果不在，则距离星下点最近的点的高度角最大，因为星下点高度角是90度,还要考虑到三角形边上的点
    {
        int minIndex1 = 0;
        double minDistance1 = 6371000; //至少40000km 最下的距离意味着最大的高度角
        double minDistance2 = 6371000; //至少40000km 最下的距离意味着最大的高度角
        double nearestPointBL[2]={0.0} ;  // the longtitude and latitude of the nearest point
        double neareatPointXYZ[3]={0.0} ; // the xyz coordinates of the nearest point
        
        for( int i = 0 ; i< 3; i++ )
        {
           double dis = GCDistance( tpoints[i].getLat() , tpoints[i].getLon(), subsatpoint[0], subsatpoint[1], 6371);
           if( dis < minDistance1 )
           {
               minIndex1 = i;
               minDistance1 = dis;
           }
        }
        
        //找到了三角形3个顶点中离星下点最近的点
        //然后寻找三条边上离星下点最近的点
        // first, find out the direction vector of the intersection line
        // all these planes should include the origin of the sphere
        for( int j = 0 ; j< 3; j++ )
        {
            int k = j+1;
            if( k== 3 ) { k = 0;}
            
            double a[3] = {0.0}; // the vector of substellar point to the origin
            double b[3] = {0.0}; // the normal vector of the known edge plane
            double c[3] = {0.0}; // the normal vecor of the new plane
            double d[3] = {0.0};// the director vector of the intersection line
            
            double t1[3] = {0.0}; // the vector of the point A and origin
            double t2[3] = {0.0}; // the vector of the point B and origin
            t1[0] = tpoints[j].getX()/1000.0;t1[1] = tpoints[j].getY()/1000.0;t1[2] = tpoints[j].getZ()/1000.0;
            t2[0] = tpoints[k].getX()/1000.0;t2[1] = tpoints[k].getY()/1000.0;t2[2] = tpoints[k].getZ()/1000.0;
            memcpy(a,pos,sizeof(double)*3);
            a[0] = a[0]/1000.0;a[1] = a[1]/1000.0;a[2] = a[2]/1000.0;
            double t1_len = sqrt(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2] );
            double t2_len = sqrt(t2[0]*t2[0]+t2[1]*t2[1]+t2[2]*t2[2] );
            double a_len = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2] );
            a[0] = a[0]/a_len;a[1] = a[1]/a_len;a[2] = a[2]/a_len;
            t1[0] = t1[0]/t1_len;t1[1] = t1[1]/t1_len;t1[2] = t1[2]/t1_len;
            t2[0] = t2[0]/t2_len;t2[1] = t2[1]/t2_len;t2[2] = t2[2]/t2_len;
            
            b[0] = t1[1]*t2[2]-t1[2]*t2[1];
            b[1] = t1[2]*t2[0]-t1[0]*t2[2];
            b[2] = t1[0]*t2[1]-t1[1]*t2[0];
            
            c[0] = a[1]*b[2]-a[2]*b[1];
            c[1] = a[2]*b[0]-a[0]*b[2];
            c[2] = a[0]*b[1]-a[1]*b[0];
            
            d[0] = c[1]*b[2]-c[2]*b[1];
            d[1] = c[2]*b[0]-c[0]*b[2];
            d[2] = c[0]*b[1]-c[1]*b[0];
            
            //the vector d should be between vector t1 and t2
            // this can help to determine the direction of the intersection line
            double  len1 = sqrt(t1[0]*t1[0]+t1[1]*t1[1] + t1[2]*t1[2]);
            double  len2 = sqrt(t2[0]*t2[0]+t2[1]*t2[1] + t2[2]*t2[2]);
            double theta1 = acos((t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2])/len1/len2);
            double len3 = sqrt(d[0]*d[0]+d[1]*d[1] + d[2]*d[2]);
            double theta2 = acos((t1[0]*d[0]+t1[1]*d[1]+t1[2]*d[2])/len1/len3);
            double theta3 = acos((t2[0]*d[0]+t2[1]*d[1]+t2[2]*d[2])/len2/len3);
            if( theta2 > theta1 || theta3 > theta1)
            {
                d[0] = -d[0];
                d[1] = -d[1];
                d[2] = -d[2];
            }
            
            
            
            //将d坐标转换为经纬度
            len1 = sqrt(a[0]*a[0]+a[1]*a[1] + a[2]*a[2]);
            double theta = acos((a[0]*d[0]+a[1]*d[1]+a[2]*d[2])/len1/len3);
            double dis  = theta*6371.0;
            
            //check whether this point is inside the triangle
            neareatPointXYZ[0] = d[0]/len3*6371000;
            neareatPointXYZ[1] = d[1]/len3*6371000;
            neareatPointXYZ[2] = d[2]/len3*6371000;
            len = sqrt(neareatPointXYZ[0]*neareatPointXYZ[0] + neareatPointXYZ[1]*neareatPointXYZ[1]);
            nearestPointBL[0] = atan(neareatPointXYZ[2]/sqrt(len))*180.0/PI;  // latitude
            nearestPointBL[1] = atan2(neareatPointXYZ[1],neareatPointXYZ[0])*180.0/PI;  // longitude
            if( nearestPointBL[1] < 0.0 )  { nearestPointBL[1] = nearestPointBL[1]+360.0; }
            bool t = IsPointInTri_sph(tpoints, nearestPointBL[0], nearestPointBL[1]);
            if(t == false )
            {
                dis = 6371000;
            }
            
            if( dis < minDistance2 )
            {
                minDistance2 = dis;
                neareatPointXYZ[0] = d[0]/len2*6371000;
                neareatPointXYZ[1] = d[1]/len2*6371000;
                neareatPointXYZ[2] = d[2]/len2*6371000;
            }
        }
        
        if( minDistance2 < minDistance1 )
        {
            len = sqrt(neareatPointXYZ[0]*neareatPointXYZ[0] + neareatPointXYZ[1]*neareatPointXYZ[1]);
            nearestPointBL[0] = atan(neareatPointXYZ[2]/sqrt(len))*180.0/PI;  // latitude
            nearestPointBL[1] = atan2(neareatPointXYZ[1],neareatPointXYZ[0])*180.0/PI;  // longitude
            if( nearestPointBL[1] < 0.0 )  { nearestPointBL[1] = nearestPointBL[1]+360.0; }
        }
        else if(minDistance1< minDistance2)
        {
            neareatPointXYZ[0] = tpoints[minIndex1].getX();
            neareatPointXYZ[1] = tpoints[minIndex1].getY();
            neareatPointXYZ[2] = tpoints[minIndex1].getZ();
            nearestPointBL[0] =  tpoints[minIndex1].getLat();
            nearestPointBL[1] =  tpoints[minIndex1].getLon();
        }
        
        
        //start calculate the elevation
        double tmp[3] ={0.0};
        double L = nearestPointBL[1]/180.0*PI;
        double B = nearestPointBL[0]/180.0*PI;
        double sinl = sin(L);
        double cosl = cos(L);
        double sinb = sin(B);
        double cosb = cos(B);
        tmp[0] = pos[0] - /*tri.getPoint(minIndex1).getX();*/  neareatPointXYZ[0];
        tmp[1] = pos[1] -       /*tri.getPoint(minIndex1).getY();*/ neareatPointXYZ[1];
        tmp[2] = pos[2] -       /*tri.getPoint(minIndex1).getZ();*/ neareatPointXYZ[2];
        
        m[0] = -sinb*cosl;m[1] = -sinb*sinl;m[2] = cosb;   //N
        m[3] = -sinl;     m[4] = cosl;      m[5] = 0.0;    //E
        m[6] = cosb*cosl; m[7] = cosb*sinl; m[8] = sinb;   //U
        
        neu_sat[0] = m[0]*tmp[0] + m[1]*tmp[1] + m[2]*tmp[2];
        neu_sat[1] = m[3]*tmp[0] + m[4]*tmp[1] + m[5]*tmp[2];
        neu_sat[2] = m[6]*tmp[0] + m[7]*tmp[1] + m[8]*tmp[2];
        
        double len = sqrt( pow(neu_sat[0] ,2) + pow(neu_sat[1] ,2) + pow(neu_sat[2],2)  );
        
        azel[0] = atan2(neu_sat[1],neu_sat[0])*180.0/PI;
        
        if( azel[0] < 0.0 )
        {
            azel[0]+= 360.0;
        }
        
        azel[1] = asin(neu_sat[2]/len)*180.0/PI;  //取值范围 (-PI/2 , PI/2) 高度角
        
        // visiable if elevation is greater than 0.0
        if(  azel[1] < -1.0E-6 )  // 最大的高度角小于0, 说明完全不可见
        {
            test = false;
        }
        else
        {
            test = true;
        }
        
    }

    return test;
}





/*
 * ********the searching strategy which is regarded faster******
 * According to the position of satellite, find out which triangles are visiable.
 * inputs: pos, alltri
 * output: triCode, the code string of the visiable triangles
 *
 */
void visibleArea( double* pos, std::vector< std::vector<GTriangle<GVertex> > >& alltri,std::vector<GTriangle<GVertex> >& myres )
{
    
    int maxlevel = alltri.size(); // the maxlevel of triangle division
    int levelID  = 0;  // 总是从第1层
    
    std::vector<GTriangle<GVertex> >  mytri;
    
    while (1)
    {
        if( levelID == 0 )
        {
            // the 8 triangles in level 0
            for( int k = 0 ; k< 8 ; k++ )
            {
                bool t = visiabilityTest(alltri[levelID][k], pos );
                if( t == true )
                {
                    myres.push_back(alltri[levelID][k]);
                }
            }
        }
        else
        {
            //every triangle has 4 children
            for( int k = 0 ; k< mytri.size(); k++ )
            {
                std::vector<GString> strcode;
                mytri[k].getCtag(strcode); // get str code of the children
                
                for( int j = 0 ; j< 4; j++ ) // check its 4 children
                {
                    int level =  strcode[j].substr(0,1).asINT();
                    int num   =  strcode[j].substr(1,8).asINT();
                    bool t = visiabilityTest( alltri[level][num-1], pos );
                    if( t == true )
                    {
                        myres.push_back(alltri[level][num-1]);
                    }
                }
            }
        }
        
        levelID++;
        if( levelID >= maxlevel )
        {
            break;
        }
        mytri = myres;
        myres.clear();
        
    }
    
    
}



/*
 
 get the certain netcdf flux data
 
*/
void getFluxDataTest(gfc::GString netcdfname)
{
    
    GNetcdf mync("../data/CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Subset_200503-201503.nc",GNetcdf::NCMODE::MODE_NOWRITE);
    int ndim = mync.getNdim();
    
    std::vector<double> timeData = mync.getDimensionData("time");
    std::vector<double> lonData =  mync.getDimensionData("lon");
    std::vector<double> latData =  mync.getDimensionData("lat");
    
    int dim_time_len = timeData.size();
    int dim_lon_len =  lonData.size();
    int dim_lat_len =  latData.size();
    
    // first step get the average data monthly;
    int dataTag[12] ={0};
    
    size_t* start = new size_t[ndim];
    memset(start,0,sizeof(size_t)*ndim);
    size_t* count = new size_t[ndim];
    memset( count,0,sizeof(size_t)*ndim);
    
    GString  varName = "toa_sw_all_mon";
    count[0] = 1; count[1] = dim_lat_len; count[2] = dim_lon_len;
    
    double toaData[180][360]={{0}};
    std::vector<FLUXDATA> totalData;
    totalData.resize(12);
    CivilTime ct( 2000,3,1,0,0,0,"tsUTC" ); //starter time
    GTime  time0;
    time0.SetFromCivilTime(ct);
    for( int rec = 0 ; rec< dim_time_len; rec++ )
    {
        start[0] = rec;
        mync.getData(varName, start, count, &toaData[0][0]);
        GTime curtime ;
        curtime.SetData(TimeSystem::GetByName("tsUKN"), timeData[rec], 0, 0);
        curtime = time0 + curtime;
        JDTime jt = GTime::GTime2JDTime(curtime);
        CivilTime  myct = GTime::JDTime2CivilTime(jt);
        dataTag[myct.m_month-1]++;
        for( int i = 0 ; i< dim_lat_len; i++)
        {
            for( int j =0; j< dim_lon_len; j++ )
            {
                totalData[myct.m_month-1].data[i][j] += toaData[i][j];
            }
        }
    }
    
    // get the monthly average data
    
    for( int k =0 ; k< totalData.size() ; k++ )
    {
        char tmpstr[4]={0};
        sprintf(tmpstr, "%02d",k+1);
        GString dataFileName(tmpstr);
        dataFileName.append(".data");
        GString dirpath = "../data/";
        
        FILE* meanFile = fopen((dirpath+dataFileName).c_str(),"w+");
        //fprintf(meanFile, "Month: %02d\n",k+1);
        for( int i = 0 ; i<dim_lon_len; i++  ) //360
        {
            for( int j = 0 ; j< dim_lat_len; j++ ) //180
            {
                totalData[k].data[j][i] = totalData[k].data[j][i] / dataTag[k];
                fprintf(meanFile, "%8.3f ",totalData[k].data[j][i]);
            }
            fprintf(meanFile, "\n");
        }
        fclose(meanFile);
    }
    
    //the data is stored in totalData variable monthly
    
    //generating the triangle grid data
    GIcosahedron<GVertex> myico;
    
    int level = 6;   //the maximum may be 6
    
    std::vector< std::vector<GTriangle<GVertex> > > allTri; // 第一层是level,第二层是三角形的vector
    
    std::vector< std::vector<GVertex> > allVer;
    allTri =  myico.createGeometry( level,allVer);
    
    //visibility test
    //given the position of the satellite,find out which triangles are visible
    double satpos[3] = { -7247484.740, 14497341.270, -22684596.470};
    std::vector<GTriangle<GVertex> > myres;
    clock_t startTime,endTime;
    startTime = clock();
    visibleArea( satpos, allTri,myres );
    //…calculating…
    endTime = clock();
    printf("time=%f\n",( (double)(endTime-startTime))/CLK_TCK);
    
    //check the visiable area
    //output the visiable triangles
    FILE* visiableTriFile = fopen("../data/vt.py","w+");
    fprintf(visiableTriFile, "Points=[\n");
    for( int i = 0 ; i< myres.size(); i++ )
    {
        fprintf(visiableTriFile, "[ [%.5f, %.5f],[%.5f, %.5f],[%.5f, %.5f] ],\n",
                myres[i].getPoint(0).getLon(),
                myres[i].getPoint(0).getLat(),
                myres[i].getPoint(1).getLon(),
                myres[i].getPoint(1).getLat(),
                myres[i].getPoint(2).getLon(),
                myres[i].getPoint(2).getLat()
                );
    }
    fprintf(visiableTriFile, "];\n");
    fclose(visiableTriFile);
    
    GString triNetfileName("../data/triNet.data");
    myico.outputTriangleNet(triNetfileName);
    myico.outputTriangleNetJS("../data/myfaces.js");
    myico.outputVertexJS("../data/myvertex.js");
    FILE*  polygonFile = fopen("../data/mypolygon.py","w+");
    fprintf( polygonFile, "Points=[\n");
    
    // I need to judege which grid belonds to a triangle, and add up all the flux.
    for( int i = 0; i < level+1 ; i++ )
    {
        char charlevel[2] ={0};
        GString levelstr = "level";
        sprintf(charlevel, "%1d",i);
        GString tmpstr(charlevel);
        levelstr = levelstr + tmpstr;
        GString faceFilePath = "../data/";
        faceFilePath = faceFilePath +levelstr;
        FILE* faceFile = fopen(faceFilePath.c_str(),"w+");
        
        FILE* fluxFile_js = fopen("../data/myflux.js","w+");
        FILE* fluxData_py = fopen("../data/myflux.py", "w+");
        
        fprintf(fluxFile_js, "var myflux=[ \n");
        fprintf(fluxData_py, "myflux=[ \n");
        double totalFluxData =0.0, totalFluxDataTest = 0.0;
        for( int g = 0 ; g< 180; g++ )
        {
            for( int t = 0 ;t< 360; t++ )
            {
                double value = totalData[2].data[g][t];
                totalFluxData+= value;
                //totalFluxData+= 1.0;
            }
        }
        //area
        //totalFluxData = 510064471909788.255;
        int imax = 0;
        double maxratio =0.0;
        printf("correct flux data value: %.3f\n",totalFluxData);
        for( int j = 0 ; j < allTri[i].size() ; j++ )
            //for( int j = 7 ; j < 8 ; j++ )
        {
            //j = 17409 ;
            GVertex p[3];
            for( int k = 0 ; k< 3; k++ )
            {
                p[k] = allTri[i][j].getPoint(k);
            }
            
            long double Ra = 6371000;//6408137;     // a = 6378137 for WGS84
            double Rb = 6371000;//6386651.7;
            double perimeter = 0.0;
            double myarea = 0.0;
            GeodesicExact geod(Ra, (Ra-Rb)/Ra);
            //Alternatively: const Geodesic& geod = Geodesic::WGS84();
            PolygonAreaExact poly(geod);
            
            //poly.AddPoint(90, 0);
            //poly.AddPoint(0, 0);
            //poly.AddPoint(0, 90);
            //poly.AddPoint(0, 180);
            //poly.Compute(false, true, perimeter, myarea);
            //long double areaA = 4*GeographicLib::Math::pi()*Ra*Ra;
            
            for( int t = 0  ;t < 3 ;t++ )
            {
                poly.AddPoint(p[t].getLat(), p[t].getLon());
            }
            poly.Compute(false, true, perimeter, myarea);
            poly.Clear();
            
            double flux =  countingStars(poly ,p,totalData[2],polygonFile);
            
            //double ratio = (flux-myarea)/myarea;
            //if(fabs(ratio) > fabs(maxratio))
            //{
            //    imax = j;
            //    maxratio = ratio;
            //}
            //if(fabs(ratio*100) > 1E-6)
            //{
            //    int c = 0;
            //}
            
//            fprintf(polygonFile, "];\n");
//            if( polygonFile != NULL )
//            {
//                fclose(polygonFile);
//                polygonFile = NULL;
//            }
            
            allTri[i][j].setValue(flux);
            totalFluxDataTest += flux;
            
            std::vector<GString> childrenTag ;
            allTri[i][j].getCtag(childrenTag);
            fprintf(faceFile, "%9s %9s %9s %9s %9s %9s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                    allTri[i][j].getTag().c_str(),
                    allTri[i][j].getFtag().c_str(),
                    childrenTag[0].c_str(),
                    childrenTag[1].c_str(),
                    childrenTag[2].c_str(),
                    childrenTag[3].c_str(),
                    allTri[i][j].getX(0),
                    allTri[i][j].getY(0),
                    allTri[i][j].getZ(0),
                    allTri[i][j].getX(1),
                    allTri[i][j].getY(1),
                    allTri[i][j].getZ(1),
                    allTri[i][j].getX(2),
                    allTri[i][j].getY(2),
                    allTri[i][j].getZ(2),
                    allTri[i][j].getValue()
                    );
            
            
            //printf("J_Number:%d flux value: %f  ratio: %f\n",j,flux,ratio);
            
            fprintf(fluxFile_js, "[ %.4f ]",flux);
            fprintf(fluxData_py, "%.4f ",flux);
            if( j != allTri[i].size() -1 )
            {
                fprintf(fluxFile_js, ",");
                fprintf(fluxData_py, ",");
            }
            fprintf(fluxFile_js, "\n");
            fprintf(fluxData_py, "\n");
        }
        
        if(faceFile != NULL )
        {
            fclose(faceFile);
        }
        
        
        fprintf( fluxFile_js, "];\n");
        fprintf(fluxData_py, "];\n");
        if( fluxFile_js != NULL )
        {
            fclose(fluxFile_js);
            fluxFile_js = NULL;
        }
        
        if( fluxData_py != NULL )
        {
            fclose(fluxData_py);
            fluxData_py = NULL;
        }
        
        
        
        printf("difference: %e\n",(totalFluxDataTest- totalFluxData)/totalFluxData);
        
    }
    
    fprintf(polygonFile, "];\n");
    if( polygonFile != NULL )
    {
        fclose(polygonFile);
        polygonFile = NULL;
    }


    
}




