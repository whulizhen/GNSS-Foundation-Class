//
//  GGlobalPressureTemperature.cpp
//  GFC
//
//  Created by lizhen on 14/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GGlobalPressureTemperature.hpp"

namespace gfc
{
    //these are variables to store the grid file, gpt2_5.grid
    double GGlobalPressureTemperature::PGRID[2592][5];
    double GGlobalPressureTemperature::TGRID[2592][5];
    double GGlobalPressureTemperature::QGRID[2592][5];
    double GGlobalPressureTemperature::DTGRID[2592][5];
    double GGlobalPressureTemperature::U[2592];
    double GGlobalPressureTemperature::HS[2592];
    double GGlobalPressureTemperature::AHGRID[2592][5];
    double GGlobalPressureTemperature::AWGRID[2592][5];
    bool GGlobalPressureTemperature::gridFileLoaded = false;
    void GGlobalPressureTemperature::loadGridFile(gfc::GString gridfile)
    {
        char LINE[1024] = {0};
        std::fstream gridFile(gridfile);
        if( gridFile.is_open() == false )
        {
            printf("GPT: global grid file open error!\n");
            return;
        }
        gridFile.getline(LINE, 1024);
        int lineNum = 0;
        while(!gridFile.eof())
        {
            memset(LINE,0,sizeof(char)*1024);
            gridFile.getline(LINE, 1024);
            GString myline(LINE);
            std::vector<GString> splitstr = myline.split();
            
            for( int i = 0 ; i< 5; i++ )
            {
                PGRID[lineNum][i] = splitstr[2+i].asDOUBLE(); // pressure in Pascal
                TGRID[lineNum][i] = splitstr[7+i].asDOUBLE(); // temperature in Kelvin
                QGRID[lineNum][i] = splitstr[12+i].asDOUBLE()/1000.0; // specific humidity in kg/kg
                DTGRID[lineNum][i] = splitstr[17+i].asDOUBLE()/1000.0; // temperature lapse rate in Kelvin/m
                
                AHGRID[lineNum][i] = splitstr[24+i].asDOUBLE()/1000.0; // hydrostatic mapping function coefficient, dimensionless
                AWGRID[lineNum][i] = splitstr[29+i].asDOUBLE()/1000.0; // wet mapping function coefficient, dimensionless
            }
            
            U[lineNum] = splitstr[22].asDOUBLE(); //  geoid undulation in m
            HS[lineNum] = splitstr[23].asDOUBLE(); // orthometric grid height in m
            
            lineNum++;
        }
        
        gridFile.close();
        
        gridFileLoaded = true;
        
    }  // end of function loadGridFile
    
    
    // get all the parameters at certain time and position
    void GGlobalPressureTemperature::gpt2(double MJD, double DLAT, double DLON, double HELL,
                   double &P, double &T, double &DT, double &E,
                   double &AH, double &AW, double &UNDU, int IT)
    {
        
        static double PI = 3.1415926535;
        //Define the mean gravity in m/s**2
        static double GM = 9.80665;
        //  Define the molar mass of dry air in kg/mol
        static double DMTR = 28.965E-3;
        //  Universal gas constant in J/K/mol
        static double RG = 8.3143;
        
        double UNDUL[4] = {0.0};
        double QL[4] = {0.0};
        double DTL[4] = {0.0};
        double TL[4] = {0.0};
        double PL[4] = {0.0};
        double AHL[4] = {0.0};
        double AWL[4] = {0.0};
        
        int INDX[4] = {-1};
        
        //Change the reference epoch to January 1 2000
        double DMJD1 = MJD-51544.5;
        
        double COSFY =0.0, COSHY =0.0, SINFY =0.0,SINHY =0.0;
        
        // Define factors for amplitudes
        if ( IT == 1)   // constant parameters
        {
            COSFY = 0.0;
            COSHY = 0.0;
            SINFY = 0.0;
            SINHY = 0.0;
        }
        else
        {
            COSFY = cos(DMJD1/365.25*2*PI);
            COSHY = cos(DMJD1/365.25*4*PI);
            SINFY = sin(DMJD1/365.25*2*PI);
            SINHY = sin(DMJD1/365.25*4*PI);
        }
        
        
        double PLON =0.0, PLAT =0.0;
        double PPOD =0.0;
        //! only positive longitude in degrees
        if( DLON < 0.0)
        {
            PLON = ( DLON + 2.0*PI)*180.0/PI;
        }
        else
        {
            PLON = DLON*180.0/PI;
        }
        
        //! transform to polar distance in degrees
        PPOD = (- DLAT + PI/2.0)*180.0/PI;
        //! find the index (line in the grid file) of the nearest point
        int IPOD = floor((PPOD+5.0)/5.0);
        int ILON = floor((PLON+5.0)/5.0);
        
        // ! normalized (to one) differences, can be positive or negative
        double DIFFPOD = (PPOD - (IPOD*5.0 - 2.50))/5.0 ;
        double DIFFLON = (PLON - (ILON*5.0 - 2.50))/5.0 ;
        
        // ! added by HCY
        if ( IPOD == 37)
        {
            IPOD = 36;
        }
        
        // ! get the number of the corresponding line, !!minus one
        INDX[0] = (IPOD - 1)*72 + ILON - 1;
        
        //! near the poles: nearest neighbour interpolation, otherwise: bilinear
        int IBILINEAR = 0;
        if ((PPOD > 2.50 ) && (PPOD < 177.50))
        {
            IBILINEAR = 1;
        }
        
        
        //! case of nearest neighbourhood
        if( IBILINEAR == 0 )
        {
            int IX = INDX[0];
            //! transforming ellipsoidial height to orthometric height
            UNDU = U[IX];
            
            double HGT = HELL - UNDU;
            
            // ! pressure, temperature at the height of the grid
            double T0 =   TGRID[IX][0]
            +  TGRID[IX][1]*COSFY + TGRID[IX][2]*SINFY
            +  TGRID[IX][3]*COSHY + TGRID[IX][4]*SINHY ;
            
            double P0 = PGRID[IX][0]
            + PGRID[IX][1]*COSFY + PGRID[IX][2]*SINFY
            + PGRID[IX][3]*COSHY + PGRID[IX][4]*SINHY;
            
            // ! specific humidity
            double Q =  QGRID[IX][0]
            + QGRID[IX][1]*COSFY + QGRID[IX][2]*SINFY
            + QGRID[IX][3]*COSHY + QGRID[IX][4]*SINHY ;
            
            // ! lapse rate of the temperature
            DT =  DTGRID[IX][0]
            + DTGRID[IX][1]*COSFY + DTGRID[IX][2]*SINFY
            + DTGRID[IX][3]*COSHY + DTGRID[IX][4]*SINHY;
            
            //! station height - grid height
            double REDH = HGT - HS[IX];
            //! temperature at station height in Celsius
            T = T0 + DT*REDH - 273.150;
            
            //! temperature lapse rate in degrees / km
            DT = DT*1000.0;
            //! virtual temperature in Kelvin
            double TV = T0*(1.0 + 0.60770*Q);
            
            double C = GM*DMTR/(RG*TV);
            
            // ! pressure in hPa
            P = (P0*exp(-C*REDH))/100.0;
            
            // ! water vapour pressure in hPa
            E = (Q*P)/(0.6220 + 0.3780*Q);
            
            // ! hydrostatic coefficient ah
            AH =  AHGRID[IX][0]
            + AHGRID[IX][1]*COSFY + AHGRID[IX][2]*SINFY
            + AHGRID[IX][3]*COSHY + AHGRID[IX][4]*SINHY ;
            
            //! wet coefficient aw
            AW = AWGRID[IX][0]
            + AWGRID[IX][1]*COSFY + AWGRID[IX][2]*SINFY
            + AWGRID[IX][3]*COSHY + AWGRID[IX][4]*SINHY;
            
        }
        else   // ! bilinear interpolation
        {
            //int t1 =  DIFFPOD>0?1:-1;
            //int t2 = DIFFLON>0?1:-1;
            int IPOD1 = IPOD +  int(DIFFPOD>0?1:-1);
            int ILON1 = ILON +  int(DIFFLON>0?1:-1);// int( SIGN(1.0,DIFFLON));
            
            if (ILON1 == 73)
            {
                ILON1 = 1 ;
            }
            if (ILON1 == 0 )
            {
                ILON1 = 72 ;
            }
            
            //! get the number of the line
            INDX[1] = (IPOD1 - 1)*72 + ILON -1; // !% along same longitude
            INDX[2] = (IPOD  - 1)*72 + ILON1 -1  ; // !% along same polar distance
            INDX[3] = (IPOD1 - 1)*72 + ILON1 -1 ; // !% diagonal
            
            for( int L = 0 ; L< 4; L++)
            {
                // ! transforming ellipsoidial height to orthometric height:
                // ! Hortho = -N + Hell
                UNDUL[L] = U[INDX[L]];
                double HGT = HELL-UNDUL[L];
                // ! pressure, temperature at the heigtht of the grid
                double T0 = TGRID[ INDX[L] ][0]
                + TGRID[ INDX[L] ][1]*COSFY + TGRID[INDX[L]][2]*SINFY
                + TGRID[INDX[L]][3]*COSHY + TGRID[INDX[L]][4]*SINHY;
                
                double P0 = PGRID[INDX[L]][0] +
                PGRID[INDX[L]][1]*COSFY + PGRID[INDX[L]][2]*SINFY +
                PGRID[INDX[L]][3]*COSHY + PGRID[INDX[L]][4]*SINHY;
                
                // ! humidity
                QL[L] =    QGRID[INDX[L]][0]
                + QGRID[INDX[L]][1]*COSFY + QGRID[INDX[L]][2]*SINFY
                + QGRID[INDX[L]][3]*COSHY + QGRID[INDX[L]][4]*SINHY;
                
                // ! reduction = stationheight - gridheight
                double HS1 = HS[INDX[L]];
                double REDH = HGT - HS1 ;
                
                // ! lapse rate of the temperature in degree / m
                DTL[L] = DTGRID[INDX[L]][0] +
                + DTGRID[INDX[L]][1]*COSFY + DTGRID[INDX[L]][2]*SINFY
                + DTGRID[INDX[L]][3]*COSHY + DTGRID[INDX[L]][4]*SINHY;
                
                // ! temperature reduction to station height
                TL[L] = T0 + DTL[L]*REDH - 273.150;
                
                //! virtual temperature
                double TV = T0*(1.0+0.60770*QL[L]);
                double C = GM*DMTR/(RG*TV);
                
                //! pressure in hPa
                PL[L] = (P0*exp(-C*REDH))/100.0;
                
                // ! hydrostatic coefficient ah
                AHL[L] = AHGRID[INDX[L]][0]
                +AHGRID[INDX[L]][1]*COSFY + AHGRID[INDX[L]][2]*SINFY
                +AHGRID[INDX[L]][3]*COSHY + AHGRID[INDX[L]][4]*SINHY;
                
                //! wet coefficient aw
                AWL[L] = AWGRID[INDX[L]][0]
                + AWGRID[INDX[L]][1]*COSFY + AWGRID[INDX[L]][2]*SINFY
                + AWGRID[INDX[L]][3]*COSHY + AWGRID[INDX[L]][4]*SINHY;
            }
            
            
            double DNPOD1 = fabs(DIFFPOD); // !% distance nearer point
            double DNPOD2 = 1.0 - DNPOD1 ;//!% distance to distant point
            double DNLON1 = fabs(DIFFLON);
            double DNLON2 = 1.0 - DNLON1;
            
            //! pressure
            double R1 = DNPOD2*PL[0]+DNPOD1*PL[1];
            double R2 = DNPOD2*PL[2]+DNPOD1*PL[3];
            P = DNLON2*R1+DNLON1*R2;
            
            // ! temperature
            R1 = DNPOD2*TL[0]+DNPOD1*TL[1];
            R2 = DNPOD2*TL[2]+DNPOD1*TL[3];
            T = DNLON2*R1+DNLON1*R2;
            
            //! temperature in degree per km
            R1 = DNPOD2*DTL[0]+DNPOD1*DTL[1];
            R2 = DNPOD2*DTL[2]+DNPOD1*DTL[3];
            DT = (DNLON2*R1+DNLON1*R2)*1000.0;
            
            //! humidity
            R1 = DNPOD2*QL[0]+DNPOD1*QL[1];
            R2 = DNPOD2*QL[2]+DNPOD1*QL[3];
            double Q = DNLON2*R1+DNLON1*R2;
            E = (Q*P)/(0.6220+0.3780*Q);
            
            //! hydrostatic
            R1 = DNPOD2*AHL[0]+DNPOD1*AHL[1];
            R2 = DNPOD2*AHL[2]+DNPOD1*AHL[3];
            AH = DNLON2*R1+DNLON1*R2;
            
            //! wet
            R1 = DNPOD2*AWL[0]+DNPOD1*AWL[1];
            R2 = DNPOD2*AWL[2]+DNPOD1*AWL[3];
            AW = DNLON2*R1+DNLON1*R2;
            
            // ! undulation
            R1 = DNPOD2*UNDUL[0]+DNPOD1*UNDUL[1];
            R2 = DNPOD2*UNDUL[2]+DNPOD1*UNDUL[3];
            UNDU = DNLON2*R1+DNLON1*R2;
            
        }
        
    } // end of function getParameters
    
    /*
     
     ref: IERS2000 gpt.f
     
     */
    
    void GGlobalPressureTemperature::gpt( double MJD, double DLAT, double DLON, double HGT,
                  double& PRES, double& TEMP, double& UNDU)
    {
        double V[10][10] = {{0.0}};
        double W[10][10] = {{0.0}};
        
        //double AP_MEAN[55], BP_MEAN[55],AP_AMP[55],BP_AMP[55];
        //double AT_MEAN[55], BT_MEAN[55],AT_AMP[55],BT_AMP[55];
        //double A_GEOID[55], B_GEOID[55];
        
        int I, N, M, NMAX, MMAX;
        
        double TWOPI = 6.283185307179586476925287;
        double DOY,TEMP0,PRES0,APM,APA,ATM,ATA,HORT,X,Y,Z;
        
        double A_GEOID[55] =
        {-5.6195E-001,-6.0794E-002,-2.0125E-001,-6.4180E-002,-3.6997E-002,
            +1.0098E+001,+1.6436E+001,+1.4065E+001,+1.9881E+000,+6.4414E-001,
            -4.7482E+000,-3.2290E+000,+5.0652E-001,+3.8279E-001,-2.6646E-002,
            +1.7224E+000,-2.7970E-001,+6.8177E-001,-9.6658E-002,-1.5113E-002,
            +2.9206E-003,-3.4621E+000,-3.8198E-001,+3.2306E-002,+6.9915E-003,
            -2.3068E-003,-1.3548E-003,+4.7324E-006,+2.3527E+000,+1.2985E+000,
            +2.1232E-001,+2.2571E-002,-3.7855E-003,+2.9449E-005,-1.6265E-004,
            +1.1711E-007,+1.6732E+000,+1.9858E-001,+2.3975E-002,-9.0013E-004,
            -2.2475E-003,-3.3095E-005,-1.2040E-005,+2.2010E-006,-1.0083E-006,
            +8.6297E-001,+5.8231E-001,+2.0545E-002,-7.8110E-003,-1.4085E-004,
            -8.8459E-006,+5.7256E-006,-1.5068E-006,+4.0095E-007,-2.4185E-008};
        
        double B_GEOID[55] =
        {
            +0.0000E+000,+0.0000E+000,-6.5993E-002,+0.0000E+000,+6.5364E-002,
            -5.8320E+000,+0.0000E+000,+1.6961E+000,-1.3557E+000,+1.2694E+000,
            +0.0000E+000,-2.9310E+000,+9.4805E-001,-7.6243E-002,+4.1076E-002,
            +0.0000E+000,-5.1808E-001,-3.4583E-001,-4.3632E-002,+2.2101E-003,
            -1.0663E-002,+0.0000E+000,+1.0927E-001,-2.9463E-001,+1.4371E-003,
            -1.1452E-002,-2.8156E-003,-3.5330E-004,+0.0000E+000,+4.4049E-001,
            +5.5653E-002,-2.0396E-002,-1.7312E-003,+3.5805E-005,+7.2682E-005,
            +2.2535E-006,+0.0000E+000,+1.9502E-002,+2.7919E-002,-8.1812E-003,
            +4.4540E-004,+8.8663E-005,+5.5596E-005,+2.4826E-006,+1.0279E-006,
            +0.0000E+000,+6.0529E-002,-3.5824E-002,-5.1367E-003,+3.0119E-005,
            -2.9911E-005,+1.9844E-005,-1.2349E-006,-7.6756E-009,+5.0100E-008
        };
        
        double AP_MEAN[55] =
        {
            +1.0108E+003,+8.4886E+000,+1.4799E+000,-1.3897E+001,+3.7516E-003,
            -1.4936E-001,+1.2232E+001,-7.6615E-001,-6.7699E-002,+8.1002E-003,
            -1.5874E+001,+3.6614E-001,-6.7807E-002,-3.6309E-003,+5.9966E-004,
            +4.8163E+000,-3.7363E-001,-7.2071E-002,+1.9998E-003,-6.2385E-004,
            -3.7916E-004,+4.7609E+000,-3.9534E-001,+8.6667E-003,+1.1569E-002,
            +1.1441E-003,-1.4193E-004,-8.5723E-005,+6.5008E-001,-5.0889E-001,
            -1.5754E-002,-2.8305E-003,+5.7458E-004,+3.2577E-005,-9.6052E-006,
            -2.7974E-006,+1.3530E+000,-2.7271E-001,-3.0276E-004,+3.6286E-003,
            -2.0398E-004,+1.5846E-005,-7.7787E-006,+1.1210E-006,+9.9020E-008,
            +5.5046E-001,-2.7312E-001,+3.2532E-003,-2.4277E-003,+1.1596E-004,
            +2.6421E-007,-1.3263E-006,+2.7322E-007,+1.4058E-007,+4.9414E-009
        };
        
        double BP_MEAN[55] =
        {
            +0.0000E+000,+0.0000E+000,-1.2878E+000,+0.0000E+000,+7.0444E-001,
            +3.3222E-001,+0.0000E+000,-2.9636E-001,+7.2248E-003,+7.9655E-003,
            +0.0000E+000,+1.0854E+000,+1.1145E-002,-3.6513E-002,+3.1527E-003,
            +0.0000E+000,-4.8434E-001,+5.2023E-002,-1.3091E-002,+1.8515E-003,
            +1.5422E-004,+0.0000E+000,+6.8298E-001,+2.5261E-003,-9.9703E-004,
            -1.0829E-003,+1.7688E-004,-3.1418E-005,+0.0000E+000,-3.7018E-001,
            +4.3234E-002,+7.2559E-003,+3.1516E-004,+2.0024E-005,-8.0581E-006,
            -2.3653E-006,+0.0000E+000,+1.0298E-001,-1.5086E-002,+5.6186E-003,
            +3.2613E-005,+4.0567E-005,-1.3925E-006,-3.6219E-007,-2.0176E-008,
            +0.0000E+000,-1.8364E-001,+1.8508E-002,+7.5016E-004,-9.6139E-005,
            -3.1995E-006,+1.3868E-007,-1.9486E-007,+3.0165E-010,-6.4376E-010
        };
        
        double AP_AMP[55] =
        {
            -1.0444E-001,+1.6618E-001,-6.3974E-002,+1.0922E+000,+5.7472E-001,
            -3.0277E-001,-3.5087E+000,+7.1264E-003,-1.4030E-001,+3.7050E-002,
            +4.0208E-001,-3.0431E-001,-1.3292E-001,+4.6746E-003,-1.5902E-004,
            +2.8624E+000,-3.9315E-001,-6.4371E-002,+1.6444E-002,-2.3403E-003,
            +4.2127E-005,+1.9945E+000,-6.0907E-001,-3.5386E-002,-1.0910E-003,
            -1.2799E-004,+4.0970E-005,+2.2131E-005,-5.3292E-001,-2.9765E-001,
            -3.2877E-002,+1.7691E-003,+5.9692E-005,+3.1725E-005,+2.0741E-005,
            -3.7622E-007,+2.6372E+000,-3.1165E-001,+1.6439E-002,+2.1633E-004,
            +1.7485E-004,+2.1587E-005,+6.1064E-006,-1.3755E-008,-7.8748E-008,
            -5.9152E-001,-1.7676E-001,+8.1807E-003,+1.0445E-003,+2.3432E-004,
            +9.3421E-006,+2.8104E-006,-1.5788E-007,-3.0648E-008,+2.6421E-010
        };
        
        
        double BP_AMP[55] =
        {
            +0.0000E+000,+0.0000E+000,+9.3340E-001,+0.0000E+000,+8.2346E-001,
            +2.2082E-001,+0.0000E+000,+9.6177E-001,-1.5650E-002,+1.2708E-003,
            +0.0000E+000,-3.9913E-001,+2.8020E-002,+2.8334E-002,+8.5980E-004,
            +0.0000E+000,+3.0545E-001,-2.1691E-002,+6.4067E-004,-3.6528E-005,
            -1.1166E-004,+0.0000E+000,-7.6974E-002,-1.8986E-002,+5.6896E-003,
            -2.4159E-004,-2.3033E-004,-9.6783E-006,+0.0000E+000,-1.0218E-001,
            -1.3916E-002,-4.1025E-003,-5.1340E-005,-7.0114E-005,-3.3152E-007,
            +1.6901E-006,+0.0000E+000,-1.2422E-002,+2.5072E-003,+1.1205E-003,
            -1.3034E-004,-2.3971E-005,-2.6622E-006,+5.7852E-007,+4.5847E-008,
            +0.0000E+000,+4.4777E-002,-3.0421E-003,+2.6062E-005,-7.2421E-005,
            +1.9119E-006,+3.9236E-007,+2.2390E-007,+2.9765E-009,-4.6452E-009
        };
        
        double AT_MEAN[55] =
        {
            +1.6257E+001,+2.1224E+000,+9.2569E-001,-2.5974E+001,+1.4510E+000,
            +9.2468E-002,-5.3192E-001,+2.1094E-001,-6.9210E-002,-3.4060E-002,
            -4.6569E+000,+2.6385E-001,-3.6093E-002,+1.0198E-002,-1.8783E-003,
            +7.4983E-001,+1.1741E-001,+3.9940E-002,+5.1348E-003,+5.9111E-003,
            +8.6133E-006,+6.3057E-001,+1.5203E-001,+3.9702E-002,+4.6334E-003,
            +2.4406E-004,+1.5189E-004,+1.9581E-007,+5.4414E-001,+3.5722E-001,
            +5.2763E-002,+4.1147E-003,-2.7239E-004,-5.9957E-005,+1.6394E-006,
            -7.3045E-007,-2.9394E+000,+5.5579E-002,+1.8852E-002,+3.4272E-003,
            -2.3193E-005,-2.9349E-005,+3.6397E-007,+2.0490E-006,-6.4719E-008,
            -5.2225E-001,+2.0799E-001,+1.3477E-003,+3.1613E-004,-2.2285E-004,
            -1.8137E-005,-1.5177E-007,+6.1343E-007,+7.8566E-008,+1.0749E-009
            
        };
        
        double BT_MEAN[55] =
        {
          	 +0.0000E+000,+0.0000E+000,+1.0210E+000,+0.0000E+000,+6.0194E-001,
            +1.2292E-001,+0.0000E+000,-4.2184E-001,+1.8230E-001,+4.2329E-002,
            +0.0000E+000,+9.3312E-002,+9.5346E-002,-1.9724E-003,+5.8776E-003,
            +0.0000E+000,-2.0940E-001,+3.4199E-002,-5.7672E-003,-2.1590E-003,
            +5.6815E-004,+0.0000E+000,+2.2858E-001,+1.2283E-002,-9.3679E-003,
            -1.4233E-003,-1.5962E-004,+4.0160E-005,+0.0000E+000,+3.6353E-002,
            -9.4263E-004,-3.6762E-003,+5.8608E-005,-2.6391E-005,+3.2095E-006,
            -1.1605E-006,+0.0000E+000,+1.6306E-001,+1.3293E-002,-1.1395E-003,
            +5.1097E-005,+3.3977E-005,+7.6449E-006,-1.7602E-007,-7.6558E-008,
            +0.0000E+000,-4.5415E-002,-1.8027E-002,+3.6561E-004,-1.1274E-004,
            +1.3047E-005,+2.0001E-006,-1.5152E-007,-2.7807E-008,+7.7491E-009
            
        };
        
        double AT_AMP[55] =
        {
            -1.8654E+000,-9.0041E+000,-1.2974E-001,-3.6053E+000,+2.0284E-002,
            +2.1872E-001,-1.3015E+000,+4.0355E-001,+2.2216E-001,-4.0605E-003,
            +1.9623E+000,+4.2887E-001,+2.1437E-001,-1.0061E-002,-1.1368E-003,
            -6.9235E-002,+5.6758E-001,+1.1917E-001,-7.0765E-003,+3.0017E-004,
            +3.0601E-004,+1.6559E+000,+2.0722E-001,+6.0013E-002,+1.7023E-004,
            -9.2424E-004,+1.1269E-005,-6.9911E-006,-2.0886E+000,-6.7879E-002,
            -8.5922E-004,-1.6087E-003,-4.5549E-005,+3.3178E-005,-6.1715E-006,
            -1.4446E-006,-3.7210E-001,+1.5775E-001,-1.7827E-003,-4.4396E-004,
            +2.2844E-004,-1.1215E-005,-2.1120E-006,-9.6421E-007,-1.4170E-008,
            +7.8720E-001,-4.4238E-002,-1.5120E-003,-9.4119E-004,+4.0645E-006,
            -4.9253E-006,-1.8656E-006,-4.0736E-007,-4.9594E-008,+1.6134E-009
        };
        
        double BT_AMP[55] =
        {
            +0.0000E+000,+0.0000E+000,-8.9895E-001,+0.0000E+000,-1.0790E+000,
            -1.2699E-001,+0.0000E+000,-5.9033E-001,+3.4865E-002,-3.2614E-002,
            +0.0000E+000,-2.4310E-002,+1.5607E-002,-2.9833E-002,-5.9048E-003,
            +0.0000E+000,+2.8383E-001,+4.0509E-002,-1.8834E-002,-1.2654E-003,
            -1.3794E-004,+0.0000E+000,+1.3306E-001,+3.4960E-002,-3.6799E-003,
            -3.5626E-004,+1.4814E-004,+3.7932E-006,+0.0000E+000,+2.0801E-001,
            +6.5640E-003,-3.4893E-003,-2.7395E-004,+7.4296E-005,-7.9927E-006,
            -1.0277E-006,+0.0000E+000,+3.6515E-002,-7.4319E-003,-6.2873E-004,
            -8.2461E-005,+3.1095E-005,-5.3860E-007,-1.2055E-007,-1.1517E-007,
            +0.0000E+000,+3.1404E-002,+1.5580E-002,-1.1428E-003,+3.3529E-005,
            +1.0387E-005,-1.9378E-006,-2.7327E-007,+7.5833E-009,-9.2323E-009
        };
        
        /*
         *     Reference day is 28 January 1980
         *     This is taken from Niell (1996) to be consistent (See References)
         *     For constant values use: doy = 91.3125
         */
        DOY = MJD  - 44239 + 1 - 28 ;
        
        /*     Define degree n and order m EGM */
        NMAX = 9;
        MMAX = 9;
        
        /**     Define unit vector */
        X = cos(DLAT)*cos(DLON);
        Y = cos(DLAT)*sin(DLON);
        Z = sin(DLAT);
        
        /*     Legendre polynomials */
        V[0][0] = 1.0;
        W[0][0] = 0.0;
        V[0][1] = Z * V[0][0];
        W[0][1] = 0.0;
        
        for( N = 2; N<= NMAX; N++ )
        {
            V[0][N] = ((2*N-1) * Z * V[0][N-1] - (N-1) * V[0][N-2]) / N ;
            W[0][N] = 0.0;
        }
        
        for( M =1; M <= NMAX; M++ )
        {
            V[M][M] = (2*M-1) * (X*V[M-1][M-1] - Y*W[M-1][M-1]);
            W[M][M] = (2*M-1) * (X*W[M-1][M-1] + Y*V[M-1][M-1]);
            
            if( M < NMAX )
            {
                V[M][M+1] = (2*M+1) * Z * V[M][M];
                W[M][M+1] = (2*M+1) * Z * W[M][M];
            }
            
            for( N = M + 2; N <= NMAX; N++)
            {
                V[M][N] = ((2*N-1)*Z*V[M][N-1] - (N+M-1)*V[M][N-2]) / (N-M);
                W[M][N] = ((2*N-1)*Z*W[M][N-1] - (N+M-1)*W[M][N-2]) / (N-M);
            }
            
        }
        
        /*     Geoidal height  */
        UNDU = 0.0;
        I = 0;
        for ( N = 0 ; N<=NMAX; N++ )
        {
            for( M =0; M <=N; M++ )
            {
                UNDU = UNDU  + (A_GEOID[I]*V[M][N] + B_GEOID[I]*W[M][N] );
                I = I +1;
            }
        }
        
        
        /*     orthometric height */
        HORT = HGT - UNDU;
        
        /*     Surface pressure on the geoid */
        APM = 0.0;
        APA = 0.0;
        I = 0;
        for(N = 0 ; N<=NMAX; N++ )
        {
            for( M =0; M<=N; M++ )
            {
                APM = APM + (AP_MEAN[I]*V[M][N] + BP_MEAN[I]*W[M][N]);
                APA = APA + (AP_AMP[I] *V[M][N] + BP_AMP[I] *W[M][N]);
                I = I+1;
            }
        }
        
        PRES0  = APM + APA*cos(DOY/365.250*TWOPI);
        
        /*     height correction for pressure */
        PRES = PRES0*pow((1.0-0.00002260*HORT), 5.2250);
        
        /*     Surface temperature on the geoid */
        ATM = 0;
        ATA = 0;
        I = 0;
        for( N =0; N<=NMAX; N++)
        {
            for( M =0;M<=N; M++)
            {
                ATM = ATM + (AT_MEAN[I]*V[M][N] + BT_MEAN[I]*W[M][N]);
                ATA = ATA + (AT_AMP[I] *V[M][N] + BT_AMP[I] *W[M][N]);
                I = I+1;
            }
        }
        
        TEMP0 =  ATM + ATA*cos(DOY/365.250*TWOPI);
        
        /*     height correction for temperature */
        TEMP = TEMP0 - 0.0065*HORT;
        
    } // end of function gpt
    
    
    void GGlobalPressureTemperature::getPressureTemperature(double MJD, double LAT, double LON, double HGT, double &P, double &T, double &UNDU)
    {
        if(gridFileLoaded == true) // call gpt2
        {
            double dt, e,ah,aw ;
            gpt2(MJD, LAT, LON, HGT, P, T, dt, e, ah, aw, UNDU);
        }
        else if(gridFileLoaded == false) // call gpt
        {
            gpt(MJD, LAT, LON, HGT, P, T, UNDU);
        }
        
    }
    
    
} // end of namespace gfc