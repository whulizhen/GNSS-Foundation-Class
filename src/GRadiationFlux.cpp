
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

//
//  GRadiationFlux.cpp
//  GFC
//
//  Created by lizhen on 07/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GRadiationFlux.hpp"
namespace gfc
{
    
    
    double GRadiationFlux::sun_rad  = 695700.0;
    double GRadiationFlux::sun_rad2 = 695700.0*695700.0;
    
    
    // the static variables in GEarthRadiationFlux
    double GEarthRadiationFlux::Re = 6371.0; //km
    
    int GEarthRadiationFlux::maxLevel = 4; // tlevelhe max level
    
    
    bool GEarthRadiationFlux::simple = true;
    bool GEarthRadiationFlux::ceres_tri = false;
    bool GEarthRadiationFlux::ceres_original = false;
    
    GString GEarthRadiationFlux::triFluxPath = "../data/flux";
    
    std::vector< std::vector< std::vector<GEarthRadiationFlux::triGridFlux > > > GEarthRadiationFlux::ERPgrid_t;
    GEarthRadiationFlux::ceresGridFlux GEarthRadiationFlux::ERPgrid_g[12][360][180];
    
    
    
    
    /*
     inputs: n is the normal vector of the plane, unit vector
     flux: is the flux that acted on that plane
     area: is the area of that plane, in m^2
     v: reflectivity
     u: specularity
     */
    GVector GRadiationFlux::radiationForce( gfc::GVector &n, GVector& flux_dir, double flux,double area, double u, double v )
    {
        GVector force(0.0,0.0,0.0);
        double c = GCONST("CLIGHT");
        double cos_theta = -dotproduct(flux_dir, n);
        GVector reflection_dir = flux_dir + 2.0*cos_theta*n;
        double WA_C = flux*area*cos_theta/c;
        
        force = flux_dir - u*v*reflection_dir - 2.0/3.0*v*(1.0-u)*n;
        force = force*WA_C;
        
        /*
        double dp = -dotproduct(flux_dir, n);
        //only when the angle between flux and n is larger than 90, it will be illuminated
        if( dp < 0.0 )
        {
            return force;
        }
        
        double cos_theta = dp/n.norm()/flux_dir.norm();
        
        
        double W = flux*cos_theta; // w/m^2
        
        // This is one way to calculate the radiation force
        
        GVector h = dp / n.norm() * n;
        
        GVector reflect_direction = 2.0 * h + flux_dir;
        
        reflect_direction.normalise();
        
        // force due to directly incident ray
        GVector f1 = area * W / c * flux_dir;
        
        // force due to directly incident ray
        //GVector f1 = area * W*(1.0-v) / c * flux.m_dir;
        
        // force due to specular reflection
        GVector f2 = -area * W / c * u * v * reflect_direction;
        // force due to diffuse reflection
        GVector f3 = -2.0 / 3.0 * area * W * cos_theta / c * v * (1.0 - u) * n;
        
        force = f1 + f2 + f3;
        */
        
        
        
        //        // Another way to calculate the radiation force
        //          double dp = -dotproduct(flux.m_dir, n);
        //          double cos_theta = dp/n.norm()/flux.m_dir.norm();
        //          double coef = W*area/c;
        //          GVector f1 = coef*(1.0-u*v)*flux.m_dir;
        //          GVector f2 = coef*(2.0/3.0*v*(u-1.0) - 2*u*v*cos_theta)*n;
        //          force = f1 + f2;
        
        
        // the third way of making calculation
        
        //double dp = -dotproduct(flux.m_dir, n);
        
//        double area_v = area*v;
//        double spec_f = area - area_v*u;
//        double spec_n = area_v*u;
//        double diff = area_v*(u - 1.0);
//        
//        spec_f = spec_f / c;
//        spec_n = spec_n / (-0.5*c);
//        diff = diff / (1.5*c);
//        
//        GVector f1 = (W * spec_f) * flux.m_dir;
//        GVector f2 = (W * (spec_n * dp + diff)) * n; ///1086.4512
//        
//        force = f1 + f2;
//        
        
        
        return force;
    }
    
    
    /*
     * get this function by Wolfram Alpha
     */
    double GRadiationFlux::Integral(double theta, double R, double r)
    {
        double coef = 0.0;
        
        coef = 2.0*(R*cos(theta) - r)/sqrt(r*r - 2.0*r*R*cos(theta) + R*R);
        
        return coef;
    }
    
    
   
    
    
    
    //http://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-in-space
    // http://education.gsfc.nasa.gov/experimental/all98invproject.site/pages/science-briefs/ed-stickler/ed-irradiance.html
    // inverse square law
    
    void GSolarRadiationFlux::makeFlux( double tsi, GVector& sat_sunHat_eci, double dis_sat_sun)
    {
        //double au = GCONST("AU"); // au should be in km
        
        m_flux.m_dir = -sat_sunHat_eci;
        
        //m_flux.m_dir.setData(sat_sun_hat, 3, 1);
        
        m_flux.m_longwave = 0.0;
        
        //double dis_sat_sun2 = dis_sat_sun * dis_sat_sun;
        
        //        double flux_scale = au * au /((sun_rad2 + au * au) * std::log((au + sun_rad) /(au - sun_rad)) +2.0 * sun_rad * (sun_rad - au) );
        
        
        
        //        m_flux.m_shortwave = tsi * flux_scale *((sun_rad2 + dis_sat_sun2) *std::log((dis_sat_sun + sun_rad) /(dis_sat_sun - sun_rad)) +
        //         2.0 * sun_rad * (sun_rad - dis_sat_sun)) /dis_sat_sun2;
        
        
        //need to account for the eclipse
        
       //double shadowfactor= GRadiationFlux::shadowFactor(GVector& sunpos_eci, GVector& satpos_eci);
        
        
        m_flux.m_shortwave = tsi;
        
        
    }
    
    
    GEarthRadiationFlux::GEarthRadiationFlux()
    {
        m_level = maxLevel;
    }
    
    //compute the polygon area on spheres
    // vertices are in CLOCKWISE order, NOT counter-clockwise
    //[lat,lon]
    double GEarthRadiationFlux::SphericalPolygonArea(std::vector<double* > vertices, double radius )
    {
        double area = 0.0;
        unsigned long nv = vertices.size();
        double flat,flon,blat,blon,sum=0.0;
        double pi = 3.14159265357;
        double dtr = pi/180.0;
        for(int iv = 0 ; iv< nv ; iv++)
        {
            // the first vertex
            if(iv == 0)
            {
                flat = *vertices[1];
                flon = *(vertices[1]+1);
                blat = *(vertices[nv-1]);
                blon = *(vertices[nv-1]+1);
            }
            // the last vertex
            else if (iv == nv -1)
            {
                flat = *vertices[0];
                flon = *(vertices[0]+1);
                blat = *(vertices[nv-2]);
                blon = *(vertices[nv-2]+1);
            }
            
            // in the middle
            else
            {
                flat = *vertices[iv+1];
                flon = *(vertices[iv+1]+1);
                blat = *(vertices[iv-1]);
                blon = *(vertices[iv-1]+1);
            }
            
            double lat = *(vertices[iv]);
            double lon = *(vertices[iv]+1);
            
            double t_f = sin( (flon - lon)*dtr ) * cos( flat*dtr );
            double c_f = sin( flat*dtr )*cos(lat*dtr) - cos(flat*dtr)*sin(lat*dtr)*cos( (lon-flon)*dtr );
            
            double t_b = sin( (blon - lon)*dtr ) * cos( blat*dtr );
            double c_b = sin( blat*dtr )*cos(lat*dtr) - cos(blat*dtr)*sin(lat*dtr)*cos( (lon-blon)*dtr );
            
            double tranlon_f = atan2(t_f, c_f);
            double tranlon_b = atan2(t_b, c_b);
            
            double fvb = tranlon_b - tranlon_f;
            if(fvb < 0.0 ) fvb += (2*pi);
            sum += fvb;
        }
        
        area = ( sum - (nv-2)*pi )*radius*radius;
        
        return area;
    }
    
    void GEarthRadiationFlux::readingTriGridFile( gfc::GString filename, std::vector< GEarthRadiationFlux::triGridFlux >&  mydata)
    {
        //std::vector< triFGridFlux > test;
        
        fstream fs(filename.c_str(),ios::in);
        if( !fs )
        {
            printf("open file: %s failed!\n",filename.c_str());
            exit(0);
        }
        
        while( !fs.eof() )
        {
            //char tmpstr[1024]={0};
            GString tmpstr;
            getline(fs,tmpstr);
            if(tmpstr=="")
            {
                break;
            }
            
            std::vector<GString> split_str = tmpstr.split();
            
            if( split_str.size() == 14 )  // the last level, without children
            {
                
                triGridFlux mytri;
                mytri.name = split_str[0];
                mytri.fname = split_str[1];
                //without children
                
                mytri.vertice[0].xyz[0] = split_str[2].asDOUBLE()/1000.0;
                mytri.vertice[0].xyz[1] = split_str[3].asDOUBLE()/1000.0;
                mytri.vertice[0].xyz[2] = split_str[4].asDOUBLE()/1000.0;
                
                mytri.vertice[1].xyz[0] = split_str[5].asDOUBLE()/1000.0;
                mytri.vertice[1].xyz[1] = split_str[6].asDOUBLE()/1000.0;
                mytri.vertice[1].xyz[2] = split_str[7].asDOUBLE()/1000.0;
                
                mytri.vertice[2].xyz[0] = split_str[8].asDOUBLE()/1000.0;
                mytri.vertice[2].xyz[1] = split_str[9].asDOUBLE()/1000.0;
                mytri.vertice[2].xyz[2] = split_str[10].asDOUBLE()/1000.0;
                
                mytri.area = split_str[11].asDOUBLE(); //km^2
                
                mytri.longwave = split_str[12].asDOUBLE();
                mytri.shortwave = split_str[13].asDOUBLE();
                
                mytri.centroid = getCentroid(mytri);
                
                mydata.push_back(mytri);
                
                int testc = 0;
                
            }
            else if( split_str.size() == 18 ) // the normal level
            {
                
                triGridFlux mytri;
                mytri.name = split_str[0];
                mytri.fname = split_str[1];
                
                mytri.cname[0] = split_str[2];
                mytri.cname[1] = split_str[3];
                mytri.cname[2] = split_str[4];
                mytri.cname[3] = split_str[5];
                
                
                mytri.vertice[0].xyz[0] = split_str[6].asDOUBLE()/1000.0;
                mytri.vertice[0].xyz[1] = split_str[7].asDOUBLE()/1000.0;
                mytri.vertice[0].xyz[2] = split_str[8].asDOUBLE()/1000.0;
                
                mytri.vertice[1].xyz[0] = split_str[9].asDOUBLE()/1000.0;
                mytri.vertice[1].xyz[1] = split_str[10].asDOUBLE()/1000.0;
                mytri.vertice[1].xyz[2] = split_str[11].asDOUBLE()/1000.0;
                
                mytri.vertice[2].xyz[0] = split_str[12].asDOUBLE()/1000.0;
                mytri.vertice[2].xyz[1] = split_str[13].asDOUBLE()/1000.0;
                mytri.vertice[2].xyz[2] = split_str[14].asDOUBLE()/1000.0;
                
                mytri.area = split_str[15].asDOUBLE(); //km^2
                
                mytri.longwave = split_str[16].asDOUBLE();
                mytri.shortwave = split_str[17].asDOUBLE();
                
                mytri.centroid = getCentroid(mytri);
                
                mydata.push_back(mytri);
                
                
            }
            else
            {
                printf("format error of erp data files\n");
            }
            
        }
    } // end of function readingTrigridFile
    
    
    void GEarthRadiationFlux::populateEMData()
    {
        GString datadir= triFluxPath;
        
        int level = maxLevel;
        
        ERPgrid_t.resize(12);
        
        //m_lwdata.resize(12);
        //m_swdata.resize(12);
       
        for( int i = 0 ; i< 12; i++ )
        {
            ERPgrid_t[i].resize(level+1);
            
            for( int j = 0 ; j< level+1; j++ )
            {
                char filepathstr[100] = {0};
                
                sprintf(filepathstr, "%s/%02d/level%1d.txt",datadir.c_str(),i+1,j);
                
                GString filepath(filepathstr);
                //start reading files
                readingTriGridFile(filepath, ERPgrid_t[i][j]);
                
                //double total_longwave  =0.0;
                //for(int k = 0 ; k< ERPgrid[i][j].size(); k++)
               // {
               //     total_longwave += ERPgrid[i][j][k].longwave *ERPgrid[i][j][k].area;
               // }
                
                int testc = 0;
                
            }
            
        }
    
    } // end of populateEMdata
    
    
    void GEarthRadiationFlux::readingCERESGridFile( int type, gfc::GString filename, ceresGridFlux ceresgrid[360][180] )
    {
        fstream fs(filename.c_str(),ios::in);
        if( !fs )
        {
            printf("open file: %s failed!\n",filename.c_str());
            exit(0);
        }
        
        int lineNum = 0;
        while( !fs.eof() )
        {
            //char tmpstr[1024]={0};
            GString tmpstr;
            getline(fs,tmpstr);
            if( tmpstr=="" )
            {
                break;
            }
            
            std::vector<GString> split_str = tmpstr.split();
            if( split_str.size() != 180 )
            {
                printf("ERROR: CERES GRID FLUX FILE ERROR!\n");
            }
            
            for( int i = 0 ; i< 180 ; i++)
            {
                
                if( type == 0 ) // longwave part
                {
                    //ceresFlux[lineNum][i].longwave = split_str[i].asDOUBLE();
                    ceresgrid[lineNum][i].longwave = split_str[i].asDOUBLE();
                }
                else if(type == 1 )
                {
                     //ceresFlux[lineNum][i].shortwave = split_str[i].asDOUBLE();
                     ceresgrid[lineNum][i].shortwave = split_str[i].asDOUBLE();
                }
            }
            
            lineNum++;
        }

        
    } // end of function readingCERESGridFile
    
    
    void GEarthRadiationFlux::populateCERESGrid()
    {
        
        
        //calculate the area of every grid, the area is the same for all months
        
        //GString ellipoidName = "AVERAGE";//"CERES_SPHERE";//"CERES_SPHERE";
        //GEllipsoid myellipsoid = GEllipsoidMgr::GetEllipsoid(ellipoidName);
        
        double perimeter =0.0;
        double area = 0.0;
        double D2R = GCONST("D2R");
        double hgt = 0.0; // 30 km
        double xyz[3]={0.0};
        
        double radius = 6408.137; // in km
        //double radius = 6371.0;
        
        //double totalArea = 0.0;
        
        //double total_longwave = 0.0;
        
        
        //first figure out the area for every grid
        //GeodesicExact geod(6371000.0,0.0 );
        //Alternatively: const Geodesic& geod = Geodesic::WGS84();
        //PolygonAreaExact poly(geod);
        
        
        for( int i = 0 ; i< 360 ; i++ )
        {
            for( int j = 0 ; j< 180 ; j++ )
            {
                double lon  = 0.5 + i;
                
                double lat =  -89.5 + j;
                
                double blh[3] = {lat*D2R,lon*D2R,hgt};
                
                
                double sinb = sin(blh[0]);
                double sinl = sin(blh[1]);
                double cosb = cos(blh[0]);
                double cosl = cos(blh[1]);
                xyz[0] = (radius + blh[2])*cosb*cosl;
                xyz[1] = (radius + blh[2])*cosb*sinl;
                xyz[2] = (blh[2]+ radius)*sinb;
                
                
                ERPgrid_g[0][i][j].centroid.xyz[0]= xyz[0];
                ERPgrid_g[0][i][j].centroid.xyz[1]= xyz[1];
                ERPgrid_g[0][i][j].centroid.xyz[2]= xyz[2];
                
                std::vector<double*> cell;
                double p1[2] = {lat-0.5, lon-0.5};
                double p2[2] = {lat-0.5, lon+0.5};
                double p3[2] = {lat+0.5, lon+0.5};
                double p4[2] = {lat+0.5, lon-0.5};
                
                cell.push_back(p4);cell.push_back(p3);cell.push_back(p2);cell.push_back(p1);
                
                //poly.AddPoint(lat-0.5, lon-0.5);
                //poly.AddPoint(lat-0.5, lon+0.5);
                //poly.AddPoint(lat+0.5, lon+0.5);
                //poly.AddPoint(lat+0.5, lon-0.5);
                //poly.Compute(false, true, perimeter, area);  // unit: m^2
                
                double area = SphericalPolygonArea(cell, radius);
                
                //ceresFlux[i][j].area = area*1.0E-6; // unit: km^2
                ERPgrid_g[0][i][j].area = area;
                
                // totalArea += ceresFlux[i][j].area;
                //ceresFlux[i][j].area = ceresFlux[i][j].area;
                
                //printf("%.3f ",ceresFlux[i][j].area);
                
                if(ERPgrid_g[0][i][j].area< 0 )
                {
                    printf("WARNING: area < 0 \n");
                }
                //poly.Clear();
                
                //total_longwave += ceresFlux[i][j].longwave*ceresFlux[i][j].area;
                
            }  // end of for( int j = 0 ; j< 180 ; j++ )
            
        } // end of for( int i = 0 ; i< 360 ; i++ )
        
        
        //copy area and centroid information to the other 11 months
        for(int i = 1 ; i< 12; i++)
        {
            memcpy(ERPgrid_g[i], ERPgrid_g[0], sizeof(GEarthRadiationFlux::ceresGridFlux)*360*180);
        }
        
        
        GString datadir= triFluxPath;
        for( int imonth = 0 ; imonth< 12; imonth++ )
        {
            //"../data/flux/06/longwave06.grid"
            char filepathstr_l[100] = {0};
            char filepathstr_s[100] = {0};
            
            sprintf(filepathstr_l, "%s/%02d/longwave%02d.grid",datadir.c_str(),imonth+1,imonth+1);
            sprintf(filepathstr_s, "%s/%02d/shortwave%02d.grid",datadir.c_str(),imonth+1,imonth+1);
            GString longwaveFile(filepathstr_l);
            GString shortwaveFile(filepathstr_s);
            
            readingCERESGridFile(0,longwaveFile,ERPgrid_g[imonth]);
            readingCERESGridFile(1,shortwaveFile,ERPgrid_g[imonth]);
        }
        
        int testc =0;
        
    } // end of function populateCERESGrid
    
    
    void GEarthRadiationFlux::visibleArea_origin( GVector& satpos_ecef, std::vector<std::vector<triGridFlux>> &alltri,
                                                       std::vector<triGridFlux> &myres)
    {
        int levelID = m_level;
        
        std::vector<triGridFlux> mytri;
        
        for ( auto tri : alltri[levelID])
        {
            point centre = tri.centroid;
            
            GVector satmtri (satpos_ecef.x - centre.xyz[0],satpos_ecef.y - centre.xyz[1],satpos_ecef.z - centre.xyz[2]);
            
            //double rsrp = satpos_ecef.x*centre.xyz[0] + satpos_ecef.y*centre.xyz[1] + satpos_ecef.z*centre.xyz[2];
            //double rprp = centre.xyz[0]*centre.xyz[0] + centre.xyz[1]*centre.xyz[1] + centre.xyz[2]*centre.xyz[2];
            
            double dot = satmtri.x * centre.xyz[0] + satmtri.y * centre.xyz[1] + satmtri.z * centre.xyz[2];
            
            if (dot >= 0) // visible
            {
                myres.push_back(tri);
            }
        }
    } // end of function visibleArea_origin

    
    /*
     * a function to check whether a plane triangle and a circle
     * return : 0 : no intersection, 1 partly intersected or triangle is inside the circle
     * ref:http://www.phatcode.net/articles.php?id=459
     */
    bool GEarthRadiationFlux::tri_circle_intersection( double r, double x[3], double y[3])
    {
        double r2 = r*r;
        bool state = false;
        double csqr[3] = {0.0};
        //test1: vertices within circle
        for( int i = 0 ; i< 3; i++ )
        {
            csqr[i] = x[i]*x[i] + y[i]*y[i] - r2;
            if( csqr[i] <= 0.0 )
            {
                state = true;
                return true;
            }
        }
        
        //test2: circle centre within triangle
        double ex[3] ={0.0}, ey[3]={0.0};
        
        ex[0] = x[1] - x[0];
        ey[0] = y[1] - y[0];
        
        ex[1] = x[2] - x[1];
        ey[1] = y[2] - y[1];
        
        ex[2] = x[0] - x[2];
        ey[2] = y[0] - y[2];
        
        if(   ey[0]*x[0] >= ex[0]*y[0]
           && ey[1]*x[1] >= ex[1]*y[1]
           && ey[2]*x[2] >= ex[2]*y[2]
           )
        {
            state = true;
            return true;
        }
        
        //test3: circle intersects edge
        
        for( int i = 0 ; i< 3; i++ )
        {
            double k = x[i]*ex[i] + y[i]*ey[i];
            
            if( k >= 0.0 )
            {
                double len = ex[i]*ex[i] + ey[i]*ey[i];
                if( k <= len )
                {
                    if( csqr[i]*len <= k*k )
                    {
                        state = true;
                        return true;
                    }
                }
            }
        }
        
        return state;
        
    } // end of bool tri_circle_intersection
    
    
    
    bool GEarthRadiationFlux::visibilityTest_new( triGridFlux& tri ,GVector& satpos_ecef, double rv)
    {
        double x[3] = {0.0}, y[3] = {0.0}, z[3] = {0.0};
        bool t[3] = {false};
        //double points[2][3] ={{0.0}};
        GVector points[2];
        GVector subpos_ecef;
        
        //the equation of the plane is: subpos[0]*x+subpos[1]*y+subpos[2]*z = R^2
        double  xx =0.0,yy =0.0,zz =0.0;
        
        GVector normal = normalise(satpos_ecef) ;
        
        subpos_ecef = normal*Re;
        
        double c = GCONST("CLIGHT");
        double PI = GCONST("PI");
        double dt = (satpos_ecef - subpos_ecef).norm()/c;
        
        //double rotationAngle = - dt*1.002737909350795*2*PI/86400.0;  // in radians, negative means the opposite direction
        
        //exclude the triangle that is not at the same semi-sphere as subpos
        for( int i = 0 ; i< 3; i++ )
        {
            //should rotate x and y coordinate , z keep the same
//            xx = tri.vertice[i].xyz[0]*cos(rotationAngle) - tri.vertice[i].xyz[1]*sin(rotationAngle);
//            yy = tri.vertice[i].xyz[0]*sin(rotationAngle) + tri.vertice[i].xyz[1]*cos(rotationAngle);
//            zz = tri.vertice[i].xyz[2];
//            
            xx = tri.vertice[i].xyz[0] - tri.vertice[i].xyz[1];
            yy = tri.vertice[i].xyz[0] + tri.vertice[i].xyz[1];
            zz = tri.vertice[i].xyz[2];
            
            if( (xx*normal.x + yy*normal.y + zz*normal.z ) >=0.0 )
            {
                t[i] = true;
            }
            else
            {
                t[i] = false;
            }
        }
        
        if( t[0]==false && t[1] == false && t[2] == false )
        {
            return false;
        }
        
        //find out two points which are perpendicular on the projection plane
        if(normal.x!= 0.0)
        {
            points[0].x = Re /subpos_ecef.x;
            points[0].y = 0.0;
            points[0].z = 0.0;
        }
        else if(normal.y!=0.0 )
        {
            points[0].x = 0.0;
            points[0].y = Re /subpos_ecef.y;
            points[0].z = 0.0;
        }
        else if(normal.z!=0.0 )
        {
            points[0].x = 0.0;
            points[0].y = 0.0;
            points[0].z = Re /subpos_ecef.z;
        }
        
        points[0].x -=  subpos_ecef.x;
        points[0].y -=  subpos_ecef.y;
        points[0].z -=  subpos_ecef.z;
        
        points[0].normalise();
        points[1] = crossproduct(normal, points[0]);
        
        //here, choose a sphecial projection, just ignore the z cooridinate
        // !!! pay attention to the projection transformation !!!
        for( int i = 0 ; i< 3; i++ )
        {
            //should rotate x and y coordinate , z keep the same
            //xx = tri.vertice[i].xyz[0]*cos(rotationAngle) - tri.vertice[i].xyz[1]*sin(rotationAngle);
            //yy = tri.vertice[i].xyz[0]*sin(rotationAngle) + tri.vertice[i].xyz[1]*cos(rotationAngle);
            //zz = tri.vertice[i].xyz[2];
            
            //get the distance from the points to the projection plane
            //double d = fabs(xx*normal[0] + yy*normal[1] + zz*normal[2]);
            //figure out the coordiante of projection points in ECEF
            //xx = xx + normal[0]*d - subpos[0];
            //yy = yy + normal[1]*d - subpos[1];
            //zz = zz + normal[2]*d - subpos[2];
            
            double t = -((xx-subpos_ecef.x)*normal.x + (yy-subpos_ecef.y)*normal.y + (zz-subpos_ecef.z)*normal.z);
            xx= t*normal.x + xx;
            yy= t*normal.y + yy;
            zz= t*normal.z + zz;
            
            //get the vector from the projection point to subpos
            xx = xx - subpos_ecef.x;
            yy = yy - subpos_ecef.y;
            zz = zz - subpos_ecef.z;
            
            x[i]  = xx*points[0].x + yy*points[0].y + zz*points[0].z;
            y[i]  = xx*points[1].x + yy*points[1].y + zz*points[1].z;
            z[i]  = xx*normal.x  + yy*normal.y + zz*normal.z;
            
        }
        
        bool res = tri_circle_intersection(rv, x, y);
        
        return res;
        
    }  // end of function visibilityTest_new
    
    
    /*
     * ********the searching strategy which is regarded faster******
     * According to the position of satellite, find out which triangles are visiable.
     * inputs: pos, alltri
     * output: triCode, the code string of the visiable triangles
     *
     */
    void  GEarthRadiationFlux::visibleArea( GVector& satpos_ecef, std::vector< std::vector<triGridFlux > >& alltri,std::vector<triGridFlux>& myres )
    {
        
        int maxlevel = m_level + 1; // the maxlevel of triangle division
        
        int levelID  = 0;  // 总是从第1层
        
        double r = satpos_ecef.norm();
        
        //GVector subpos;
        
        //subpos = normalise(satpos_ecef)*Re;
        
        // the radius of visible circle area
        //double rv = m_Ra*sqrt(r*r - m_Ra*m_Ra)/(r + m_Ra);
        double rv = Re/r*sqrt(r*r - Re*Re); // km
        
        std::vector<triGridFlux>  mytri;
        
        while (1)
        {
            if( levelID == 0 )
            {
                // the 8 triangles in level 0
                for( int k = 0 ; k< 8 ; k++ )
                {
                    //bool t = visibilityTest(alltri[levelID][k], pos );
                    bool t = visibilityTest_new(alltri[levelID][k],satpos_ecef,rv);
                    
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
                    //GString strcode[4]= mytri[k].cname;
                    //mytri[k].getCtag(strcode); // get str code of the children
                    
                    for( int j = 0 ; j< 4; j++ ) // check its 4 children
                    {
                        int level = static_cast<int>(mytri[k].cname[j].substr(0,1).asINT()) ;
                        
                        int num   = static_cast<int>(mytri[k].cname[j].substr(1,8).asINT()) ;
                        
                        //bool t = visibilityTest( alltri[level][num-1], pos );
                        bool t = visibilityTest_new(alltri[level][num-1],satpos_ecef,rv);
                        
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
        
    } // end of function visibleArea
    
    
    
    //get the centroid of this triangle
    GEarthRadiationFlux::point GEarthRadiationFlux::getCentroid(triGridFlux& tri)
    {
        //double xyz[3]={0.0};
        
        GVector xyz;
        
        for( int i = 0 ; i< 3; i++ )
        {
            xyz.x += tri.vertice[i].xyz[0];
            xyz.y += tri.vertice[i].xyz[1];
            xyz.z += tri.vertice[i].xyz[2];
        }
        
        //xyz[0] = xyz[0]/3.0;
        //xyz[1] = xyz[1]/3.0;
        //xyz[2] = xyz[2]/3.0;
        double len = xyz.norm();
        
        xyz = xyz/len*Re;
        //extend this point onto the sphere
        //xyz[0] = xyz[0]/len*Re;
        //xyz[1] = xyz[1]/len*Re;
        //xyz[2] = xyz[2]/len*Re;
        
        //POINT p("centroid",xyz[0],xyz[1],xyz[2]);
        
        point p(xyz.x,xyz.y,xyz.z);
        
        return p;
    }
    
    
    
    void GEarthRadiationFlux::outputVisibleTri_py(std::vector<triGridFlux>& vis_triangles)
    {
        FILE* faceFile = fopen("visible_triangle.txt", "w+");
        for( int j = 0 ; j < vis_triangles.size(); j++ )
        {
            double area = vis_triangles[j].area;
            // the output of area , flux must be very high precision
            fprintf(faceFile, "%9s %9s %9s %9s %9s %9s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %20.16f %20.12f %20.12f\n",
                    vis_triangles[j].name.c_str(),
                    vis_triangles[j].fname.c_str(),
                    vis_triangles[j].cname[0].c_str(),
                    vis_triangles[j].cname[1].c_str(),
                    vis_triangles[j].cname[2].c_str(),
                    vis_triangles[j].cname[3].c_str(),
                    
                    vis_triangles[j].vertice[0].xyz[0],
                    vis_triangles[j].vertice[0].xyz[1],
                    vis_triangles[j].vertice[0].xyz[2],
                    
                    vis_triangles[j].vertice[1].xyz[0],
                    vis_triangles[j].vertice[1].xyz[1],
                    vis_triangles[j].vertice[1].xyz[2],
                    
                    vis_triangles[j].vertice[2].xyz[0],
                    vis_triangles[j].vertice[2].xyz[1],
                    vis_triangles[j].vertice[2].xyz[2],
                    area,
                    vis_triangles[j].longwave,
                    vis_triangles[j].shortwave
                    );
        }
        
        fclose(faceFile);
        
//
//        FILE* pf = fopen("vt.py", "w+");
//        double pi = GCONST("PI");
//
//        fprintf(pf, "Points=[\n");
//
//        for( int i = 0 ; i< vis_triangles.size() ; i++ )
//        {
//            fprintf(pf, "[ ");
//            for(int j = 0 ; j<3 ; j++ )
//            {
//                double lat =0.0, lon =0.0;
//
//                double dot1 =  sqrt( pow(vis_triangles[i].vertice[j].xyz[0],2.0) + pow(vis_triangles[i].vertice[j].xyz[1],2.0));
//                if(fabs(dot1)<1.0E-10)
//                {
//                    if( vis_triangles[i].vertice[j].xyz[2] >= 0 )
//                    {
//                        lat = 90.0;
//                        lon = 0.0;
//                    }
//                    else if(vis_triangles[i].vertice[j].xyz[2] < 0)
//                    {
//                        lat = -90.0;
//                        lon = 0.0;
//                    }
//                }
//                else
//                {
//                    lat =  atan(vis_triangles[i].vertice[j].xyz[2]/dot1)*180.0/pi; // -pi/2 ~ pi/2
//                    lon =  atan2(vis_triangles[i].vertice[j].xyz[1], vis_triangles[i].vertice[j].xyz[0])*180.0/pi;
//                    if( lon < 0.0 )
//                    {
//                        lon = lon + 360.0;
//                    }
//                }
//
//                if(j == 2 )
//                {
//                   fprintf(pf, "[%9.6f, %9.6f] ",lon, lat);
//                }
//                else
//                {
//                    fprintf(pf, "[%9.6f, %9.6f], ",lon, lat);
//                }
//
//
//            }
//
//           if( i == vis_triangles.size() - 1)
//           {
//               fprintf(pf, "]\n");
//           }
//           else
//           {
//              fprintf(pf, "],\n");
//           }
//
//        }
//
//        fprintf(pf, "];");
//
//        fclose(pf);
        
    }
    
    
    
    double GEarthRadiationFlux::fluxCoef( double theta, double R, double r)
    {
        double coef = 0.0;
        coef = (r*r*(R*cos(theta)-r))/(R*R*sqrt(r*r-2.0*r*R*cos(theta) + R*R) );
        //coef = Integral( theta_up, R, r) - Integral(theta_down, R, r);
        
        return coef;
    }
    
    
    void GEarthRadiationFlux::Flux_simple(GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos)
    {
        double pi = GCONST("PI");
        
        // Average long and short wavelength fluxes from CERES data quality summary
        double avg_earth_flux_lw = 239.6;
        double avg_earth_flux_sw = 99.6;
        
        double a_ceres = 6408.1370;
        double b_ceres = 6386.6517;
        double ceres_earth_rad = (a_ceres + b_ceres) / 2.0;
        double ceres_earth_rad2 = ceres_earth_rad * ceres_earth_rad;
        
        
        // the half angle of the satellite's visible field
        double cos_theta = ceres_earth_rad / satpos_ecef.norm();
        
        //圆锥的顶角和立体角的关系
        // OM is the solid angle https://en.wikipedia.org/wiki/Solid_angle
        double OM = 2.0*pi*(1.0 - cos_theta);
        
        //the visible area
        double A = OM*ceres_earth_rad2;
        
        double coef = A*2.0/(3.0 *pi*satpos_ecef.norm2() );
        //double coef = ceres_earth_rad2/satpos_ecef.norm2();
        
        
        
        
        //double erp_flux_lw_scale = 0.5 * avg_earth_flux_lw;
        //double erp_flux_sw_scale = 0.5 * avg_earth_flux_sw;
        
        //double erp_flux_lw_scale = -2.0 * avg_earth_flux_lw;
        //double erp_flux_sw_scale = -2.0 * avg_earth_flux_sw;
        
        //double r2 = dis_satpos* dis_satpos;
        
//        double scale_factor =
//        ( (ceres_earth_rad2 + r2) *
//         std::log((dis_satpos + ceres_earth_rad) / (dis_satpos - ceres_earth_rad))
//         + 2.0*ceres_earth_rad * (ceres_earth_rad - dis_satpos) )/r2;
//
        // get this formular according to the integral
        //double scale_factor_lw = 2.0*( 1.0-std::sqrt(satpos_ecef.norm2() - ceres_earth_rad2) / satpos_ecef.norm() );
        
        fluxdata ecef_flux;
        
        //should use another scale_factor for shortwave,because the position of sun should be considered
        //GVector psun = normalise(sunpos_ecef);
        //GVector psat = normalise(satpos_ecef);
        
        //double gamma = acos(dotproduct(psun, psat)/psun.norm()/psat.norm()); // (0,pi/2)
        //if( gamma > pi/2.0 ) gamma = gamma - pi/2.0;
        
        //double e = pi/2 - gamma;
        
        //double ee = acos(Re/dis_satpos);
        
//        if( e < ee )
//        {
//            double t1 = fluxCoef(ee, Re, dis_satpos);
//            double t2 = fluxCoef(e, Re, dis_satpos);
//            scale_factor_sw = t1 - t2 ;
//        }
//        else
//        {
//            scale_factor_sw = scale_factor_lw;
//        }
        
        
        

        
        //double az = atan2(psun.y,psun.x);
        //if( az < 0) az += 2*pi;
       // double el = asin(psun.z/psun.norm());
        
       // double p1[2]={0.0}, p2[2]={0.0};
        
        ecef_flux.m_dir = normalise(satpos_ecef);
        
        ecef_flux.m_longwave = coef * avg_earth_flux_lw;
        ecef_flux.m_shortwave = coef* avg_earth_flux_sw;
        
        // Apply additional scaling to the shortwave flux based on sun position:
        double test = dotproduct(normalise(sunpos_ecef), normalise(satpos_ecef));
        if( test <= 0.0 )
        {
            ecef_flux.m_shortwave = 0.0;
        }
        
        
        /*
        ecef_flux.m_shortwave *=
        (0.5 * dotproduct(normalise(sunpos_ecef), normalise(satpos_ecef)) + 0.5);
        */
        
        
        
        // here m_flux is the ECEF flux
        m_flux.push_back(ecef_flux);
        
    }
    
    
    /*
     * this function is used for the original ceres data
     *  satpos_ecef is in WGS84 ellipsoid
     *  for sure that, CERES FLUX is just a scaling of WGS84, thus, do NOT need a converstion between them,
     *   that means they have same longitude and latitude,with hight: 30km
     */
    void GEarthRadiationFlux::Flux_ceresOriginal( int imonth, GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos)
    {
        /*
        struct polygon
        {
            std::vector<double> lat;
            std::vector<double> lon;
        } ;
        std::vector<polygon> visibleArea;
        */
        
        double pi = 3.14159265357;
        
        fluxdata myflux;
        
        double coef = 0.0 ;
        
        //the position of the satellite
        //GVector r = satpos_ecef;
        //GVector s = sunpos_ecef;
        
       // GVector totalFLux;
        int count = 0;
        
        //FILE* pf = fopen("visible.txt", "w+");
        
        for( int i = 0 ; i< 360 ; i++ )
        {
            for( int j = 0; j< 180; j++ )
            {
                // ECEF pos of the cells (the centroid of the cell)
                GVector gridpos(ERPgrid_g[imonth][i][j].centroid.xyz[0],ERPgrid_g[imonth][i][j].centroid.xyz[1],ERPgrid_g[imonth][i][j].centroid.xyz[2]);
                
                double r_grid = gridpos.norm();
                
                GVector satmgrid = satpos_ecef - gridpos;
                
                GVector sunmgrid = sunpos_ecef - gridpos;
                
                double r_satmgrid = satmgrid.norm();
                
                double cos_theta = dotproduct(gridpos, satmgrid)/(r_grid*satmgrid.norm());
                
                double cos_theta_s = dotproduct(gridpos, sunmgrid)/ (r_grid*sunmgrid.norm());
                
                // the direction of the earth radiation flux
                myflux.m_dir = normalise(gridpos);
                
                if( cos_theta > 0 ) // visible for the satellite
                {
                    coef = ERPgrid_g[imonth][i][j].area*cos_theta/(pi*r_satmgrid*r_satmgrid);
                    
                    myflux.m_longwave = ERPgrid_g[imonth][i][j].longwave*coef;
                    
                    myflux.m_shortwave = ERPgrid_g[imonth][i][j].shortwave*coef;
                    
                    if( cos_theta_s <= 0) // invisible for the sun
                    {
                        myflux.m_shortwave = 0.0;
                    }
                    
                    count++;
                    
                    //m_flux.push_back(myflux);
                    
                    totalFlux_lw += myflux.m_dir*myflux.m_longwave;
                    totalFlux_sw += myflux.m_dir*myflux.m_shortwave;
                    
                    //totalFLux += myflux.m_dir*(myflux.m_longwave + myflux.m_shortwave);
                }
            }
        }
        
    //    fprintf(pf, "Points=[\n");
    //
    //    for( int i = 0 ; i< visibleArea.size() ; i++ )
    //    {
    //        fprintf(pf, "[ ");
    //        for(int j = 0 ; j<visibleArea[i].lat.size() ; j++ )
    //        {
    //            double lat =0.0, lon =0.0;
    //
    //            lon = visibleArea[i].lon[j];
    //            lat = visibleArea[i].lat[j];
    //
    //            if(j == visibleArea[i].lat.size() - 1 )
    //            {
    //                fprintf(pf, "[%9.6f, %9.6f] ",lon, lat);
    //            }
    //            else
    //            {
    //                fprintf(pf, "[%9.6f, %9.6f], ",lon, lat);
    //            }
    //
    //
    //        }
    //
    //        if( i == visibleArea.size() - 1)
    //        {
    //            fprintf(pf, "]\n");
    //        }
    //        else
    //        {
    //            fprintf(pf, "],\n");
    //        }
   //
   //     }
   //
   //     fprintf(pf, "];");
   //
   //     fclose(pf);
   //
   //     int testc = 0;
    

        
        
        
//        double t1 = totalFLux.norm();
//        double total = 0.0;
//        for( int i = 0 ; i< m_flux.size() ; i++ )
//        {
//            double test = m_flux[i].m_dir.norm();
//            
//            if(fabs(test-1)>0.01 )
//            {
//                
//                printf("WARNING\n");
//            }
//            
//            double t = (m_flux[i].m_longwave + m_flux[i].m_shortwave);
//            total += t;
//        }
//        
//        double norm = totalFLux.norm();
//        
//        printf("Flux Norm: %f\n",total);
//        
//        fluxdata tflux;
//        tflux.m_dir = normalise(r);
//        tflux.m_longwave = total;
//        
//        m_flux.clear();
//        
//        m_flux.push_back(tflux);
//        
//        int testc = 0;
        
    }
    
    
    
    void gfc::GEarthRadiationFlux::Flux_ceres(int imonth, GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos)
    {
        double pi = GCONST("PI");
        
        double c = GCONST("CLIGHT");
        
        fluxdata ecef_flux;
        
        
        std::vector< triGridFlux > vis_triangles;
        vis_triangles.reserve(4*pow(4,m_level));
        m_flux.reserve(4*pow(4,m_level));
        //m_level = 0;
        //visibleArea( satpos_ecef, ERPgrid[imonth],vis_triangles);
        
        visibleArea_origin(satpos_ecef, ERPgrid_t[imonth], vis_triangles);
        
        outputVisibleTri_py( vis_triangles );
        
        double sin_sun_elev = 0.0, r_minus_p_mag=0.0, r_minus_p_mag2 = 0.0, coef =0.0;
        
        GVector p(0.0,0.0,0.0);
        
        //GVector totalFlux;
        
        for( int i =0; i< vis_triangles.size(); i++ )
        {
            p.set( vis_triangles[i].centroid.xyz[0], vis_triangles[i].centroid.xyz[1], vis_triangles[i].centroid.xyz[2]);
            
            GVector r_p =  satpos_ecef - p;
            
            GVector s_p =  sunpos_ecef - p;
            
            ecef_flux.m_dir = normalise(p);
            
            double cos_theta = dotproduct(r_p, p)/(r_p.norm()*p.norm());
            
            coef = vis_triangles[i].area * cos_theta /( pi*r_p.norm2() );
            
            ecef_flux.m_longwave =  vis_triangles[i].longwave*coef;
            ecef_flux.m_shortwave = vis_triangles[i].shortwave*coef;
            
            //if invisible for the Sun, the shortwave flux is set to be 0
            if(dotproduct(s_p, p) < 0.0 )
            {
                ecef_flux.m_shortwave = 0.0;
            }
            
            //m_flux.push_back(ecef_flux);
            
            totalFlux_lw += ecef_flux.m_dir*ecef_flux.m_longwave;
            totalFlux_sw += ecef_flux.m_dir*ecef_flux.m_shortwave;
        }
        
       // double test = 0.0;
       // for( int i = 0 ; i< m_flux.size(); i++ )
       // {
       //     test += (m_flux[i].m_longwave+ m_flux[i].m_shortwave);
       // }
        
        //totalFlux_lw = totalFlux_lw.norm()*normalise(satpos_ecef);
        //totalFlux_sw = totalFlux_sw.norm()*normalise(satpos_ecef);
        
        int test_int = 0;
    }
    
    
    //use the simple flux model
    void gfc::GEarthRadiationFlux::makeFlux( int imonth, GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos  )
    {
        
        totalFlux_lw.x =0.0;
        totalFlux_lw.y =0.0;
        totalFlux_lw.z =0.0;
        totalFlux_sw.x =0.0;
        totalFlux_sw.y =0.0;
        totalFlux_sw.z =0.0;
        
        if( ceres_tri == true)
        {
            Flux_ceres( imonth, sunpos_ecef, satpos_ecef, dis_satpos);
        }
        else if(ceres_original == true)
        {
             Flux_ceresOriginal( imonth, sunpos_ecef, satpos_ecef, dis_satpos);
        }
        else if(simple==true)
        {
            Flux_simple( sunpos_ecef,satpos_ecef, dis_satpos);
        }
        
        
        
//        for( int i =0; i< m_flux.size() ; i++ )
//        {
//            totalFlux_lw += m_flux[i].m_dir*m_flux[i].m_longwave;
//            totalFlux_sw += m_flux[i].m_dir*m_flux[i].m_shortwave;
//        }
        
        //convert the Earth flux from ECEF to ECI
        GVector t;
        GSpaceEnv::eop.ECEF2ECI_pos(totalFlux_lw, t);
        totalFlux_lw = t;
        GSpaceEnv::eop.ECEF2ECI_pos(totalFlux_sw, t);
        totalFlux_sw = t;
        
        int testc = 0;
        
    } // end of function makeFlux
    
} // end of namespace
