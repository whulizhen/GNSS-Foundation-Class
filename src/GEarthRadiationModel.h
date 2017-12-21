//
//  ERP_grid2Triangle.h
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



#ifndef ERP_grid2Triangle_h
#define ERP_grid2Triangle_h

#include <stdio.h>

#include "GNetcdf.h"

#include "GEarthRadiationModel.h"

#include "GIcosahedron.h"
#include <algorithm>
#include <fstream>
#include "GTime.h"
#include "GLegendre.hpp"
#include "GMatrix.h"

#include <GeographicLib/PolygonArea.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Gnomonic.hpp>


using namespace GeographicLib;

using namespace std;
using  std::sqrt;
using  std::sin;
using  std::cos;
using  std::atan2;
using  std::acos;
using  std::asin;

//const long double eps_angle = 1.0E-6;   //degree
//const long double eps_distance = eps_angle*110000;  //meter
//const long double eps_general = 1.0E-10;  // for general use
//const double sphere_radius = 6371000;  // the radius of sphere

//struct FLUXDATA
//{
//    double data[180][360];
//    FLUXDATA()
//    {
//        memset(data,0,sizeof(double)*180*360);
//    }
//};
//
//double GCDistance(double lat1, double lon1, double lat2, double lon2,double radius);
//double getAngle(double* pA, double* pC );
//
//bool SphereIntersection(double lat1a, double lon1a, double lat1b,
//                        double lon1b, double lat2a, double lon2a,
//                        double lat2b, double lon2b, double& lat3A,
//                        double& lon3A );
//
//void gridIndex( double min[2], double max[2],double interval[2], double xy[2], int* pos  );
//
///*treat it as sphere*/
//bool IsPointInTri_sph( GVertex p[3] , double lat, double lon);
//
//void ClockwiseSortPoints(std::vector<std::vector<double> > &Points);
//
//void determinGrid(double lat, double lon, double* interval, double* min, double* max,double* starter);
//
//std::vector<std::vector<double> > polygonPointsProcess(std::vector<std::vector<double> >& allpoints);
//
//
//double countingStars( GVertex p[3] ,FLUXDATA fluxdata,FILE* polygonFile);
//
//
//
////the main test function
//void getFluxDataTest(gfc::GString netcdfname);




//not used
//void GeodesicIntersection(double* lat, double* lon, GeographicLib::Geodesic geod, double* intersectionPoint);
//bool IsPointInTri( GVertex p[3] , double lat, double lon);
//double polygon_area( std::vector< std::vector<double> >& points );



namespace gfc
{
    
    
    

    
//the definition of the class
// this class is used to modelling the earth radiation using CERES data
class GEarthRadiationModel
{
    
public:
    
    constexpr  static   double eps_angle = 1.0E-6;   //degree
    constexpr  static   double eps_distance = eps_angle*110000;  //meter
    
    struct FLUXDATA
    {
        double data[180][360];
        FLUXDATA()
        {
            memset(data,0,sizeof(double)*180*360);
        }
    };
    
    GEarthRadiationModel();
    
    virtual ~GEarthRadiationModel();
    
    void loadData();
    
    void setVariableName(GString varName);
    
    void createGrid();  //creat the triangle grids
    void attachFlux(int imonth, int type);  // attach the original CERES flux data onto the triangle girds
    
    void dumpMonToaData(GString dir);
    
    void dumpTriGridFlux(int imonth,GString dir);
    void outputTriFlux_py(int mylevel);
    void outputTriFlux_js();
    
    void outputVertice_js( std::vector<GVertex>& vertice, FILE* pf);
    void outputTri_js( std::vector<GTriangle<GVertex> >& tri, FILE* pf);
    
    
    // tool functions
    double countingStars(PolygonAreaExact poly, GVertex p[3] ,FLUXDATA fluxdata,FILE* polygonFile);
    bool  IsPointInTri_sph( GVertex p[3] , double lat, double lon);
    bool  SphereIntersection(double lat1a, double lon1a, double lat1b,
                                                   double lon1b, double lat2a, double lon2a,
                                                   double lat2b, double lon2b, double& lat3A,
                                                   double& lon3A );
    
    double GCDistance( double lat1, double lon1, double lat2, double lon2,double radius);
    void  gridIndex( double min[2], double max[2],double interval[2], double xy[2], int* pos  );
    void  ClockwiseSortPoints(std::vector<std::vector<double> > &Points);
    void  determinGrid(double lat, double lon, double* interval, double* min, double* max,double* starter);
    std::vector<std::vector<double> > polygonPointsProcess(std::vector<std::vector<double> >& allpoints);
    double getAngle( double* pA, double* pC);
    
/*visibility test*/
bool visiabilityTest( GTriangle<GVertex>& tri ,double* pos);
void visibleArea( double* pos, std::vector< std::vector<GTriangle<GVertex> > >& alltri,std::vector<GTriangle<GVertex> >& myres );
    
void createFluxModel();
    
    
/* for the on-line calculation */
    void readingTriGridFile(gfc::GString filename, std::vector< GTriangle<GVertex> >&  mydata);
    
    void SphericalHarmonicsTest();
    
    
    //store for all the 12 months
   // std::vector< std::vector< std::vector<GTriangle<GVertex> > > > m_allstorage;
    
private:
    
    double m_Ra;
    
    GString m_ceresDataFile;
    GString m_varableName;  // the varable extracted
    GString m_dataDirectory;  // the directory for storing the data files
    
    //double m_toaData[180][360]; // the monthly averaged data
    
    std::vector<FLUXDATA> m_monToaData[2]; //// the monthly averaged data
    
    //std::vector<FLUXDATA> m_monToaData; //// the monthly averaged data
    
    GTime m_startTime;  //utc
    
    //generating the triangle grids
    //GIcosahedron<GVertex> m_myico;
    int m_level ;
    //for the single month
    std::vector< std::vector<GTriangle<GVertex> > >  m_TriGrid; //第一层是level,第二层是三角形的vector
    
    
    
};

    
void RadiationModelTest();
    
    
    
} // end of namespace

#endif /* ERP_grid2Triangle_h */
