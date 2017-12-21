//
//  GGlobalPressureTemperature.hpp
//  GFC
//
//  Created by lizhen on 14/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GGlobalPressureTemperature_hpp
#define GGlobalPressureTemperature_hpp

#include <stdio.h>
#include "GString.h"
#include <fstream>
#include <math.h>
namespace gfc
{
    // the class to calculate the Global Pressure and Temperature
    /*
     ref: IERS 2010, GPT2.F
     */
    class GGlobalPressureTemperature
    {
        
        
    public:
        
        static void loadGridFile( GString gridfile);
        
        void gpt2( double MJD, double DLAT, double DLON, double HELL,
                  double &P, double &T, double &DT, double &E,
                  double &AH, double &AW, double &UNDU, int IT = 0);
        
        
        void gpt( double MJD, double DLAT, double DLON, double HGT,
                 double& PRES, double& TEMP, double& UNDU);
        
        
        void getPressureTemperature(double MJD,double LAT,double LON, double HGT,
                                    double& P,double& T,double& UNDU);
        
    private:
        
        //these are variables to store the grid file, gpt2_5.grid
        static double PGRID[2592][5];
        static double TGRID[2592][5];
        static double QGRID[2592][5];
        static double DTGRID[2592][5];
        static double U[2592];
        static double HS[2592];
        static double AHGRID[2592][5];
        static double AWGRID[2592][5];
        static bool gridFileLoaded;
        
    };
    
    
    
} // end of namespace gfc



#endif /* GGlobalPressureTemperature_hpp */
