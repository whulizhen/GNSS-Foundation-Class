
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
//  GRadiationFlux.hpp
//  GFC
//
//  Created by lizhen on 07/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GRadiationFlux_hpp
#define GRadiationFlux_hpp

#include <stdio.h>
#include <fstream>
#include "GMatrix.h"
#include "GMath.hpp"
#include "GVector.hpp"
#include "GFCCONST.h"
#include "GEllipsoidMgr.h"
#include "GSpaceEnv.hpp"


//#include <GeographicLib/PolygonArea.hpp>
//#include <GeographicLib/Geodesic.hpp>
//#include <GeographicLib/Constants.hpp>
//#include <GeographicLib/Gnomonic.hpp>

//using namespace GeographicLib;

namespace gfc
{
    
    struct fluxdata
    {
        GVector m_dir; //unit vector of the flux direction
        double  m_longwave;
        double  m_shortwave;
        fluxdata()
        {
            m_longwave =  0.0;
            m_shortwave = 0.0;
        }
    };
    
    
    // This is the class for the modelling of all the flux in the space
    class GRadiationFlux
    {
        
    public:
        
        static double sun_rad ;
        
        static double sun_rad2 ;
        
        GRadiationFlux() {}
        
        virtual ~ GRadiationFlux() {}
        
        static GVector radiationForce( gfc::GVector &n, GVector& flux_dir, double flux,double area, double u, double v );
        //a test function
        static double Integral( double theta ,double R, double r );
        
    };
    
    
    /*
     *
     * this is the class for providing the radiation flux from sun
     * that acted on the spacecraft
     *
    */
    class GSolarRadiationFlux : public GRadiationFlux
    {
        
    public:
        
        fluxdata m_flux;
        
        // sat_sun must be unit vector from sat to sun, here should consider the ECLIPSE factor for the radiation flux
        void makeFlux( double tsi, GVector& sat_sunHat_eci, double dis_sat_sun );
        
    };
    
    
    /*
     * the radiation flux from the earth
     *
     */
    class GEarthRadiationFlux: public GRadiationFlux
    {
        
    private:
        struct point
        {
            double xyz[3];
            int visited;
            point()
            {
                visited = 0;
                memset(xyz,0,sizeof(double)*3);
            }
            point(double x, double y, double z, int isvisited=0)
            {
                xyz[0] = x; xyz[1] =y; xyz[2] = z; visited = isvisited;
            }
            
        };
        
        struct triGridFlux
        {
            GString name;
            GString fname;
            GString cname[4];
            point   vertice[3];
            point   centroid;
            double  area;
            double  longwave;
            double  shortwave;
        };
        
       
        // this is for the original ceres grid data
        /*the struct for storing the original ceres grid data*/
        struct ceresGridFlux
        {
            point centroid;
            
            double longwave;
            
            double shortwave;
            
            double area;
            
            ceresGridFlux()
            {
                longwave = 0.0;
                shortwave =0.0;
                area = 0.0;
            }
        };
        
    public:
        
        // the choices for the earth radiation flux
        static bool ceres_original; // ceres original grid
        static bool ceres_tri; // ceres trianle grid
        static bool simple; // simple earth flux
        
        //member variable section
        std::vector<fluxdata> m_flux;
        
        GVector totalFlux_lw;
        GVector totalFlux_sw;
        
        int m_level;   // the level for the current choice
        
        GEarthRadiationFlux();
        
        void makeFlux( int imonth, GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos);
        
        void outputVisibleTri_py( std::vector<GEarthRadiationFlux::triGridFlux>& vis_triangles);
        
        static void populateEMData();
        
        static void populateCERESGrid();
        
        
        
        static GString  triFluxPath;  // this should be set in the configure file
        
        static int maxLevel;  // the max level for the subdivision of the grid file
        
    private:
        
        //this is for the original ceres grid data
        
       //static ceresGridFlux ceresFlux[360][180];
       static ceresGridFlux ERPgrid_g[12][360][180];
       //static double        ceresGridArea[360][180];
        
        
        static void readingCERESGridFile( int type, gfc::GString filename, ceresGridFlux ceresgrid[360][180]);
        
        // these are functions and variables for the triangular grid
        
        static double Re;  // the radius of the TOA sphere for the CERES flux
        
        static std::vector< std::vector< std::vector<triGridFlux> > >   ERPgrid_t;  // 12 months data for all the level
        
        static void readingTriGridFile( gfc::GString filename, std::vector< GEarthRadiationFlux::triGridFlux >&  mydata);
        
        static point getCentroid(triGridFlux& tri);
        
        void visibleArea_origin( GVector& satpos_ecef, std::vector<std::vector<triGridFlux> > &alltri,
                                std::vector<triGridFlux> &myres);
        
        bool tri_circle_intersection(double r, double x[3], double y[3]);
        
        bool visibilityTest_new( triGridFlux& tri ,GVector& satpos_ecef, double rv );
        
        void  visibleArea( GVector& satpos_ecef, std::vector< std::vector<triGridFlux > >& alltri,std::vector<triGridFlux>& myres );
        
        // the function for making the ceres flux
        void Flux_ceres(int imonth, GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos);
        
        void Flux_simple(GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos);
        
        // for the orininal ceres grid data
        void Flux_ceresOriginal(int imonth, GVector& sunpos_ecef, GVector& satpos_ecef, double dis_satpos);
        
        
        //this is the integral function for theta
        double fluxCoef(double theta, double R, double r);
        
    };
    
    
    
    
} // end of namespace




#endif /* GRadiationFlux_hpp */
