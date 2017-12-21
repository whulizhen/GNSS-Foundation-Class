//
//  GEarthTide.hpp
//  GFC
//
//  Created by lizhen on 22/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GEarthTide_hpp
#define GEarthTide_hpp
#include <math.h>

#include <stdio.h>
#include <complex.h>
#include "GVector.hpp"
#include "GEllipsoid.h"
#include "GEllipsoidMgr.h"
#include "GIERS.hpp"

using namespace std;

namespace gfc
{
    
    //need to include ocean tide
    
    class GTidalCorrection
    {
        
    public:
        
       void getSolidEarthTideCorrection( double JD_UT1_I,double JD_UT1_F,double JD_TT, GVector& sunpos_ecef, GVector& moonpos_ecef,double dC[], double dS[]);
       
      void getSolidEarthTideCorrection_1(double JD_UT1_I,double JD_UT1_F,double JD_TT, GVector &sunpos_ecef, GVector &moonpos_ecef, double dC[], double dS[]);
        
       static double normFactor ( int n, int m ) ;
       
        //http://iers-conventions.obspm.fr/2010/2010_update/chapter6/additional_info/tidemodels/fes2004_Cnm-Snm.dat
       static void loadFES2004(GString ocentidefile);
        
    private:
        
        ///  Legendre polynomial
        //  relevant formula can be found in "satellite orbits"(3.23,3.24,3.25 in chapter 3.2.4)
        static double legendrePoly( int n,int m,double u);
        
        /// Objects to hold parameters, coefficients for solid earth tide
        static const double Argu_C20[21][13];
        static const double Argu_C21[48][13];
        static const double Argu_C22[2][13];
        
        //the degree and order of ocean tide
        int m_n;
        int m_m;
        //in the ocean tide file: fes2004_Cnm-Snm.dat
        static const std::vector<double> dCp;  // delC+
        static const std::vector<double> dCm;  // delC-
        static const std::vector<double> dSp;  // delS+
        static const std::vector<double> dSm;  // delS-
        
    };
    
    
    
}  // end of namespace



#endif /* GEarthTide_hpp */
