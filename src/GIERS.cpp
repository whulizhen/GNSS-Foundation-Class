//
//  GIERS.cpp
//  GFC
//
//  Created by lizhen on 12/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GIERS.hpp"

namespace gfc
{
    
    
    
    
    
    /*
     *+
     *  - - - - - - - - - - -
     *   F U N D A R G
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine computes the lunisolar fundamental arguments.
     *  The model used is from Simon et al. (1994) as recommended by the IERS
     *  Conventions (2010).  Refer to IERS Conventions (2010) Chapter 5
     *  Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Canonical model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     T           d      TT, Julian centuries since J2000 (Note 1)
     *
     *  Returned:
     *     L           d      Mean anomaly of the Moon (Note 2)
     *     LP          d      Mean anomaly of the Sun (Note 2)
     *     F           d      L - OM (Notes 2 and 3)
     *     D           d      Mean elongation of the Moon from the Sun
     *                                                         (Note 2)
     *     OM          d      Mean longitude of the ascending node of
     *                                                the Moon (Note 2)
     *
     *  Notes:
     *
     *  1) Though T is strictly TDB, it is usually more convenient to use
     *     TT, which makes no significant difference.  Julian centuries since
     *     J2000 is (JD - 2451545.0)/36525.
     *
     *  2) The expression used is as adopted in IERS Conventions (2010) and
     *     is from Simon et al. (1994).  Arguments are in radians.
     *
     *  3) L in this instance is the Mean Longitude of the Moon. OM is the
     *     Mean longitude of the ascending node of the Moon.
     *
     *  Test case:
     *     given input: T = 0.07995893223819302 Julian centuries since J2000
     *                  (MJD = 54465)
     *     expected output:  L = 2.291187512612069099 radians
     *                       LP = 6.212931111003726414 radians
     *                       F = 3.658025792050572989 radians
     *                       D = 4.554139562402433228 radians
     *                       OM = -0.5167379217231804489 radians
     *
     *  References:
     *
     *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     *     Francou, G., Laskar, J., 1994, Astron.Astrophys. 282, 663-683
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  2008 January 18 B.E.Stetzler  Initial changes to header
     *               and used 2PI instead of PI as parameter
     *  2008 January 25 B.E. Stetzler Additional changes to header
     *               and defined fundamental arguments
     *  2008 January 28 B.E. Stetzler Additional changes to header
     *  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
     *  2008 April   03 B.E. Stetzler Provided example test case.
     *  2009 February 11 B.E. Stetzler Corrected term in OM from 6962890.2665
     *                                 to 6962890.5431 and updated test case
     *  2009 May     07 B.E. Stetzler Code formatting changes based on
     *                                client recommendations
     *  2009 May     07 B.E. Stetzler Updated test case due to above changes
     *  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
     *-----------------------------------------------------------------------
     ref : http://iers-conventions.obspm.fr/2010/2010_official/chapter5/software/FUNDARG.F
     */
    void GIERS::FUNDARG(double T, double &L, double &LP, double &F, double &D, double &OM)
    {
        
        //*  Arcseconds to radians
        double  DAS2R = 4.848136811095359935899141E-6 ;
        
        //*  Arcseconds in a full circle
        double TURNAS = 1296000E0;
        
        //*  2Pi
        double D2PI = 6.283185307179586476925287E0;
        
        //*  Compute the fundamental argument L.
        L = fmod ( 485868.249036E0 + T*( 1717915923.2178E0 + T*( 31.8792E0 + T*( 0.051635E0 + T*( - 0.00024470E0 )))), TURNAS ) * DAS2R;
        
        //*  Compute the fundamental argument LP.
        LP = fmod ( 1287104.79305E0 + T*( 129596581.0481E0 + T*( - 0.5532E0 + T*( 0.000136E0 + T*( - 0.00001149E0 )))), TURNAS ) * DAS2R;
        
        //*  Compute the fundamental argument F.
        F  = fmod ( 335779.526232E0 + T*( 1739527262.8478E0 + T*( - 12.7512E0 + T*( - 0.001037E0 + T*( 0.00000417E0 )))), TURNAS ) * DAS2R;
        
        //*  Compute the fundamental argument D.
        D = fmod ( 1072260.70369E0 + T*( 1602961601.2090E0 + T*( - 6.3706E0 + T*( 0.006593E0 + T*( - 0.00003169E0 )))), TURNAS ) * DAS2R;
        
        //*  Compute the fundamental argument OM.
        OM = fmod (450160.398036E0 + T*( - 6962890.5431E0 + T*(  7.4722E0 + T*( 0.007702E0 + T*(  - 0.00005939E0 )))), TURNAS ) * DAS2R ;
        
        // normalisze the Angle to -pi to pi
        /*
        L = normalizeAngle(L);
        LP = normalizeAngle(LP);
        F = normalizeAngle(F);
        D = normalizeAngle(D);
        OM = normalizeAngle(OM);
        */
    } // end of function FUNDARG
    
    
    
    //BETA are : (τ, s, h, p, N′, ps)
    void GIERS::DOODSONARG(double JD_UT1_I, double JD_UT1_F, double JD_TT,double FNUT[5], double BETA[6])
    {
        const double PI = GCONST("PI");
        
        double T = ( JD_TT - 2451545.0)/36525.0;
        double THETA = getGMST(JD_UT1_I, JD_UT1_F, JD_TT);
        //double FNUT[5] = {0.0};
        FUNDARG(T, FNUT[0], FNUT[1], FNUT[2], FNUT[3], FNUT[4]);
        
        double S = FNUT[2] + FNUT[4];
        BETA[0] = THETA+PI-S;
        BETA[1] = FNUT[2]+FNUT[4];
        BETA[2] = S-FNUT[3];
        BETA[3] = S-FNUT[0];
        BETA[4] = -FNUT[4];
        BETA[5] = S-FNUT[3]-FNUT[1];
        
    }
    
    
    
    
    
    
    
    // this function is used to calculate the corrections for x, y and s
    //  input: MJD in TT
    //  output: array H with 12 elements, partials of the tidal variation
    //  with respect to the orthoweights
    /*
     Test case:
     *     given input: dmjd = 54964.0D0
     *
     *     expected output: h(1) = 15.35873641938967360D0
     *                      h(2) = 9.784941251812741214D0
     *                      h(3) = -5.520740128266865554D0
     *                      h(4) = 3.575314211234633888D0
     *                      h(5) = -13.93717453496387648D0
     *                      h(6) = -9.167400321705855504D0
     *                      h(7) = 5.532815475865292321D0
     *                      h(8) = 9.558741883500834646D0
     *                      h(9) = -10.22541212627272600D0
     *                      h(10)= 0.8367570529461261231D0
     *                      h(11)= 1.946355176475630611D0
     *                      h(12)= -13.55702062247304696D0
     
     References:
     *
     *     Ray,R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
     *     "Diurnal and Semidiurnal Variations in the Earth's Rotation
     *     Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     * original fortran code: ftp://tai.bipm.org/iers/convupdt/chapter8/CNMTX.F
     */
    
    void GIERS::CNMTX(double MJD_TT, double *H)
    {
        
        int I = 0, J = 0, K = 0, M = 0, N = 0, NLINES = 71, NMAX = 0;
        //char NUMARG[71];
        memset(H, 0, sizeof(double)*12);
        double DMJD = MJD_TT;
        double DT60 =0.0, D1960 = 0.0 , DT =0.0 ;
        //double TWOPI = 2.0*GCONST("PI");
        static double PI = 3.141592653589793238462643E0 ;
        static double TWOPI = 6.283185307179586476925287E0;
        
        double  HS[71], PHASE[71], FREQ[71];
        double AP, AM, BP, BM, PINM, ALPHA;
        double ANM[3][4][2]={{{0.0}}};  /* 2:3 , 0:3, -1:1 */
        double BNM[3][4][2]={{{0.0}}}; /* 2:3 , 0:3, -1:1 */
        double P[2][3]={{0.0}}, Q[2][3]={{0.0}} /* 0:2, 1:2*/;
        int NJ[71]={0}, MJ[71]={0};
        /*  Define the orthotide weight factors */
        double SP[2][6] = {
            { 0.0298E0,0.1408E0,+0.0805E0, 0.6002E0,+0.3025E0, 0.1517E0 },
            { 0.0200E0,0.0905E0,+0.0638E0, 0.3476E0,+0.1645E0, 0.0923E0 }
        };
        DT = 2.0;
        NMAX = 2;
        
        /* tidal potential model for 71 diurnal and semidiurnal lines*/
        D1960 =  37076.5;
        
        static double data[71][6] = {
            2E0,1E0,  -1.94E0,  9.0899831E0,  5.18688050E0, 117.655 ,
            2E0,1E0,  -1.25E0,  8.8234208E0,  5.38346657E0, 125.745 ,
            2E0,1E0,  -6.64E0, 12.1189598E0,  5.38439079E0, 125.755 ,
            2E0,1E0,  -1.51E0,  1.4425700E0,  5.41398343E0, 127.545 ,
            2E0,1E0,  -8.02E0,  4.7381090E0,  5.41490765E0, 127.555 ,
            2E0,1E0,  -9.47E0,  4.4715466E0,  5.61149372E0, 135.645 ,
            2E0,1E0, -50.20E0,  7.7670857E0,  5.61241794E0, 135.655 ,
            2E0,1E0,  -1.80E0, -2.9093042E0,  5.64201057E0, 137.445 ,
            2E0,1E0,  -9.54E0,  0.3862349E0,  5.64293479E0, 137.455 ,
            2E0,1E0,   1.52E0, -3.1758666E0,  5.83859664E0, 145.535 ,
            2E0,1E0, -49.45E0,  0.1196725E0,  5.83952086E0, 145.545 ,
            2E0,1E0,-262.21E0,  3.4152116E0,  5.84044508E0, 145.555 ,
            2E0,1E0,   1.70E0, 12.8946194E0,  5.84433381E0, 145.755 ,
            2E0,1E0,   3.43E0,  5.5137686E0,  5.87485066E0, 147.555 ,
            2E0,1E0,   1.94E0,  6.4441883E0,  6.03795537E0, 153.655 ,
            2E0,1E0,   1.37E0, -4.2322016E0,  6.06754801E0, 155.445 ,
            2E0,1E0,   7.41E0, -0.9366625E0,  6.06847223E0, 155.455 ,
            2E0,1E0,  20.62E0,  8.5427453E0,  6.07236095E0, 155.655 ,
            2E0,1E0,   4.14E0, 11.8382843E0,  6.07328517E0, 155.665 ,
            2E0,1E0,   3.94E0,  1.1618945E0,  6.10287781E0, 157.455 ,
            2E0,1E0,  -7.14E0,  5.9693878E0,  6.24878055E0, 162.556 ,
            2E0,1E0,   1.37E0, -1.2032249E0,  6.26505830E0, 163.545 ,
            2E0,1E0,-122.03E0,  2.0923141E0,  6.26598252E0, 163.555 ,
            2E0,1E0,   1.02E0, -1.7847596E0,  6.28318449E0, 164.554 ,
            2E0,1E0,   2.89E0,  8.0679449E0,  6.28318613E0, 164.556 ,
            2E0,1E0,  -7.30E0,  0.8953321E0,  6.29946388E0, 165.545 ,
            2E0,1E0, 368.78E0,  4.1908712E0,  6.30038810E0, 165.555 ,
            2E0,1E0,  50.01E0,  7.4864102E0,  6.30131232E0, 165.565 ,
            2E0,1E0,  -1.08E0, 10.7819493E0,  6.30223654E0, 165.575 ,
            2E0,1E0,   2.93E0,  0.3137975E0,  6.31759007E0, 166.554 ,
            2E0,1E0,   5.25E0,  6.2894282E0,  6.33479368E0, 167.555 ,
            2E0,1E0,   3.95E0,  7.2198478E0,  6.49789839E0, 173.655 ,
            2E0,1E0,  20.62E0, -0.1610030E0,  6.52841524E0, 175.455 ,
            2E0,1E0,   4.09E0,  3.1345361E0,  6.52933946E0, 175.465 ,
            2E0,1E0,   3.42E0,  2.8679737E0,  6.72592553E0, 183.555 ,
            2E0,1E0,   1.69E0, -4.5128771E0,  6.75644239E0, 185.355 ,
            2E0,1E0,  11.29E0,  4.9665307E0,  6.76033111E0, 185.555 ,
            2E0,1E0,   7.23E0,  8.2620698E0,  6.76125533E0, 185.565 ,
            2E0,1E0,   1.51E0, 11.5576089E0,  6.76217955E0, 185.575 ,
            2E0,1E0,   2.16E0,  0.6146566E0,  6.98835826E0, 195.455 ,
            2E0,1E0,   1.38E0,  3.9101957E0,  6.98928248E0, 195.465 ,
            2E0,2E0,   1.80E0, 20.6617051E0, 11.45675174E0, 225.855 ,
            2E0,2E0,   4.67E0, 13.2808543E0, 11.48726860E0, 227.655 ,
            2E0,2E0,  16.01E0, 16.3098310E0, 11.68477889E0, 235.755 ,
            2E0,2E0,  19.32E0,  8.9289802E0, 11.71529575E0, 237.555 ,
            2E0,2E0,   1.30E0,  5.0519065E0, 11.73249771E0, 238.554 ,
            2E0,2E0,  -1.02E0, 15.8350306E0, 11.89560406E0, 244.656 ,
            2E0,2E0,  -4.51E0,  8.6624178E0, 11.91188181E0, 245.645 ,
            2E0,2E0, 120.99E0, 11.9579569E0, 11.91280603E0, 245.655 ,
            2E0,2E0,   1.13E0,  8.0808832E0, 11.93000800E0, 246.654 ,
            2E0,2E0,  22.98E0,  4.5771061E0, 11.94332289E0, 247.455 ,
            2E0,2E0,   1.06E0,  0.7000324E0, 11.96052486E0, 248.454 ,
            2E0,2E0,  -1.90E0, 14.9869335E0, 12.11031632E0, 253.755 ,
            2E0,2E0,  -2.18E0, 11.4831564E0, 12.12363121E0, 254.556 ,
            2E0,2E0, -23.58E0,  4.3105437E0, 12.13990896E0, 255.545 ,
            2E0,2E0, 631.92E0,  7.6060827E0, 12.14083318E0, 255.555 ,
            2E0,2E0,   1.92E0,  3.7290090E0, 12.15803515E0, 256.554 ,
            2E0,2E0,  -4.66E0, 10.6350594E0, 12.33834347E0, 263.655 ,
            2E0,2E0, -17.86E0,  3.2542086E0, 12.36886033E0, 265.455 ,
            2E0,2E0,   4.47E0, 12.7336164E0, 12.37274905E0, 265.655 ,
            2E0,2E0,   1.97E0, 16.0291555E0, 12.37367327E0, 265.665 ,
            2E0,2E0,  17.20E0, 10.1602590E0, 12.54916865E0, 272.556 ,
            2E0,2E0, 294.00E0,  6.2831853E0, 12.56637061E0, 273.555 ,
            2E0,2E0,  -2.46E0,  2.4061116E0, 12.58357258E0, 274.554 ,
            2E0,2E0,  -1.02E0,  5.0862033E0, 12.59985198E0, 275.545 ,
            2E0,2E0,  79.96E0,  8.3817423E0, 12.60077620E0, 275.555 ,
            2E0,2E0,  23.83E0, 11.6772814E0, 12.60170041E0, 275.565 ,
            2E0,2E0,   2.59E0, 14.9728205E0, 12.60262463E0, 275.575 ,
            2E0,2E0,   4.47E0,  4.0298682E0, 12.82880334E0, 285.455 ,
            2E0,2E0,   1.95E0,  7.3254073E0, 12.82972756E0, 285.465 ,
            2E0,2E0,   1.17E0,  9.1574019E0, 13.06071921E0, 295.555
        };
        
        for(int i = 0 ; i< 71; i++ )
        {
            NJ[i] = data[i][0];
            MJ[i] = data[i][1];
            HS[i] = data[i][2];
            PHASE[i] = data[i][3];
            FREQ[i] = data[i][4];
        }
        
        /* Compute the time dependent potential matrix*/
        for(K = -1; K<=1; K++ )
        {
            DT60 = (DMJD - K*DT) - D1960;
            
            ANM[K+1][1][0] = 0.0;ANM[K+1][2][0] = 0.0;
            BNM[K+1][1][0] = 0.0;BNM[K+1][2][0] = 0.0;
            for( J=1; J<=NLINES; J++ )
            {
                N = NJ[J-1];
                M = MJ[J-1];
                PINM =  fmod(N+M, 2.0)*TWOPI/4.0;
                ALPHA = fmod(PHASE[J-1]-PINM, TWOPI) + fmod(FREQ[J-1]*DT60, TWOPI);
                
                ANM[K+1][M][N-2] += HS[J-1]*cos(ALPHA);
                BNM[K+1][M][N-2] -= HS[J-1]*sin(ALPHA);
                
            }
            
            
        }
        
        
        /* orthogonalize the response terms*/
        
        for( M=1;M<=2; M++ )
        {
            AP = ANM[2][M][0] + ANM[0][M][0];
            AM = ANM[2][M][0] - ANM[0][M][0];
            BP = BNM[2][M][0] + BNM[0][M][0];
            BM = BNM[2][M][0] - BNM[0][M][0];
            
            P[M-1][0] = SP[M-1][0]*ANM[1][M][0];
            P[M-1][1] = SP[M-1][1]*ANM[1][M][0] - SP[M-1][2]*AP;
            P[M-1][2] = SP[M-1][3]*ANM[1][M][0] - SP[M-1][4]*AP + SP[M-1][5]*BM;
            Q[M-1][0] = SP[M-1][0]*BNM[1][M][0];
            Q[M-1][1] = SP[M-1][1]*BNM[1][M][0] - SP[M-1][2]*BP;
            Q[M-1][2] = SP[M-1][3]*BNM[1][M][0] - SP[M-1][4]*BP - SP[M-1][5]*AM;
            
            ANM[0][M][0] = P[M-1][0];
            ANM[1][M][0] = P[M-1][1];
            ANM[2][M][0] = P[M-1][2];
            BNM[0][M][0] = Q[M-1][0];
            BNM[1][M][0] = Q[M-1][1];
            BNM[2][M][0] = Q[M-1][2];
        }
        
        
        /* fill partials vector*/
        
        J  = 0 ;
        for(N =2; N<= NMAX; N++ )
        {
            for(M=1; M<=N; M++)
            {
                for(K = -1; K<=1; K++)
                {
                    H[J] = ANM[K+1][M][N-2];
                    H[J+1] = BNM[K+1][M][N-2];
                    J = J + 2;
                }
            }
        }
    
    } // end of function CNMTX
    
    
    
    /*
     *+
     *  - - - - - - - - - -
     *   O R T H O _ E O P
     *  - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  The purpose of the subroutine is to compute the diurnal and semi-
     *  diurnal variations in Earth Orientation Parameters (x,y, UT1) from
     *  ocean tides.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Class 1 model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     TIME           d     Modified Julian Date
     *
     *  Returned:
     *     EOP            d     delta_x, in microarcseconds
     *                          delta_y, in microarcseconds
     *                          delta_UT1, in microseconds
     *
     *  Notes:
     *
     *  1) The diurnal and semidiurnal orthoweights fit to the 8 constituents
     *     are listed in Reference 1.
     *
     *  Called:
     *     CNMTX                Compute time dependent part of second degree
     *                          diurnal and semidiurnal tidal potential from
     *                          the dominant spectral lines in the Cartwright-
     *                          Tayler-Edden harmonic decomposition                  *
     *  Test case:
     *     given input: MJD = 47100D0
     *
     *     expected output: delta_x = -162.8386373279636530D0 microarcseconds
     *                      delta_y = 117.7907525842668974D0 microarcseconds
     *                      delta_UT1 = -23.39092370609808214D0 microseconds
     *
     *  References:
     *
     *     Ray, R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
     *     "Diurnal and Semidiurnal Variations in the Earth's Rotation
     *     Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  1997 February    R. Eanes         Original code
     *  2008 November 07 B.E. Stetzler    Added header and copyright
     *  2008 November 07 B.E. Stetzler    Provided test case
     *  2009 May      12 B.E. Stetzler    Replaced ENDDO statements with
     *                                    CONTINUE statements
     *  2009 June     09 B.E. Stetzler    Updated validated test case values
     *                                    based on changes to CNMTX.F
     *  2010 March    19 B.E. Stetzler    Capitalized variables for FORTRAN
     *                                    77 backwards compatibility
     *-----------------------------------------------------------------------
     */
    void GIERS::ORTHO_EOP(double MJD_TT,double* EOP_correction)
    {
        int K = 0, J = 0;
        double ORTHOW[3][12] = {  -6.77832E0,-14.86323E0,0.47884E0,-1.45303E0,0.16406E0,  0.42030E0,0.09398E0,25.73054E0,-4.77974E0, 0.28080E0,1.94539E0, -0.73089E0,
            14.86283E0,-6.77846E0, 1.45234E0, 0.47888E0,-0.42056E0, 0.16469E0,15.30276E0,-4.30615E0, 0.07564E0, 2.28321E0,-0.45717E0,-1.62010E0,
            -1.76335E0,1.03364E0, -0.27553E0, 0.34569E0,-0.12343E0,-0.10146E0,-0.47119E0,1.28997E0, -0.19336E0, 0.02724E0, 0.08955E0, 0.04726E0
        };
        
        double H[12]={0.0};
        CNMTX(MJD_TT, H);
        
        for( K = 1; K<=3; K++ )
        {
            EOP_correction[K-1] = 0.0;
            for(J =1; J<=12; J++)
            {
                EOP_correction[K-1] += H[J-1]*ORTHOW[K-1][J-1];
            }
        }
        
    } // end of function ORTHO_EOP
    
    
    /*
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine evaluates the model of subdiurnal libration
     *  in the axial component of rotation, expressed by UT1 and LOD.
     *  This effect is due to the influence of tidal gravitation on the
     *  departures of the Earth's mass distribution from the rotational
     *  symmetry, expressed by the non-zonal components of geopotential.
     *  The amplitudes have been computed for an elastic Earth with liquid
     *  core. The adopted truncation level is 0.033 microseconds in UT1
     *  corresponding to the angular displacement of 0.5 microarcseconds
     *  or to 0.015 mm at the planet surface. With this truncation level
     *  the model contains 11 semidiurnal terms. The coefficients of
     *  the model are given in Table 5.1b of the IERS Conventions (2010).
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status:  Class 3 model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as
     *     a Class 1, 2, or 3 model.
     *
     *  Given:
     *     rmjd        d      Time expressed as modified Julian date
     *
     *  Returned:
     *     dUT1        d      Incremental UT1 in microseconds
     *     dLOD        d      Incremental LOD in microseconds per day
     *
     *  Notes:
     *  1) The procedure FUNDARG.F is the same as used by the program PMSDNUT2.F
     *     which implements the corresponding model of the lunisolar libration in
     *     polar motion.
     *
     *  Called:
     *     FUNDARG             Compute the angular fundamental arguments
     *
     *  Test cases:
     *     given input:  rmjd_a = 44239.1 ( January 1, 1980 2:24.00 )
     *                   rmjd_b = 55227.4 ( January 31, 2010 9:35.59 )
     *
     *     expected output: dUT1_a =   2.441143834386761746D0 mus;
     *                      dLOD_a = -14.78971247349449492D0 mus / day
     *                      dUT1_b = - 2.655705844335680244D0 mus;
     *                      dLOD_b =  27.39445826599846967D0 mus / day
     *
     *  References:
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  2010 May       A.Brzezinski   Original code
     *  2010 June  1   B.E.Stetzler   Initial changes to code
     *  2010 June  2   B.E.Stetzler   Provided test case
     *  2010 June  2   B.E.Stetzler   Capitalized all variables for FORTRAN
     *                                77 compatibility
     *  2010 June  2   B.E.Stetzler   Replaced call to PMARGS to FUNDARG
     *                                for universal fundamental argument
     *                                subroutine
     *  2010 June  2   B.E.Stetzler   Validated test case using internally
     *                                computed GMST and call to FUNDARG
     *                                matched previous external call to
     *                                PMARGS for all six parameters
     *  2010 June  23  B.E.Stetzler   Modified coefficients of semi-diurnal
     *                                variations in UT1 and LOD due to
     *                                libration for a non-rigid Earth to
     *                                coincide with Table 5.1b
     *  http://iers-conventions.obspm.fr/2010/2010_official/chapter5/software/UTLIBR.F
     */
    void GIERS::UTLIBR(double RMJD, double &dut1, double &dlod)
    {
        int i =0,j = 0;
        
        double T, GMST, L, LP, F,D,OM, ANGLE;
        double ARG[6];
        double DAS2R = 4.848136811095359935899141E-6;
        double TURNAS = 1296000.0;
        double RMJD0 = 51544.5;
        double PI = 3.141592653589793238462643;
        double TWOPI = 6.283185307179586476925287;
        double RAD2SEC = 86400.0/TWOPI;
        int IARG[11][6]={
                            {2, -2,  0, -2,  0, -2},
                            {2,  0,  0, -2, -2, -2},
                            {2, -1,  0, -2,  0, -2},
                            {2,  1,  0, -2, -2, -2},
                            {2,  0,  0, -2,  0, -1},
                            {2,  0,  0, -2,  0, -2},
                            {2,  1,  0, -2,  0, -2},
                            {2,  0, -1, -2,  2, -2},
                            {2,  0,  0, -2,  2, -2},
                            {2,  0,  0,  0,  0,  0},
                            {2,  0,  0,  0,  0, -1}
                        };
        double PER[11] = { 0.5377239, 0.5363232,0.5274312,0.5260835,0.5175645,0.5175251,0.5079842,0.5006854,0.5000000,0.4986348,0.4985982  };
        double DUT1S[11]={ 0.05,0.06,0.35,0.07,-0.07,1.75,-0.05, 0.04,0.76,0.21,0.06 };
        double DUT1C[11]={ -0.03,-0.03,-0.20,-0.04,0.04,-1.01,0.03, -0.03,-0.44,-0.12,-0.04};
        double DLODS[11]={ -0.3,-0.4,-2.4,-0.5,0.5,-12.2,0.3,-0.3,-5.5,-1.5,-0.4};
        double DLODC[11]={ -0.6,-0.7,-4.1,-0.8, 0.8,-21.3, 0.6,-0.6,-9.6,-2.6, -0.8};
        
        dut1 = 0.0;
        dlod = 0.0;
        T = (RMJD-RMJD0)/36525;
        GMST = fmod (   67310.54841 +
                      T*( (8640184.812866 + 3155760000.0) +
                      T*( 0.093104 +
                         T*( -0.0000062 ))), 86400.0 );
        GIERS::FUNDARG(T, L, LP, F, D, OM);
        ARG[0] = GMST/RAD2SEC + PI; ARG[0] = fmod(ARG[0], TWOPI);
        ARG[1] = L; ARG[2] = LP; ARG[3] = F; ARG[4] = D; ARG[5]= OM;
        for( j= 0 ; j < 11; j++ )
        {
            ANGLE = 0.0;
            for(i= 0; i< 6; i++ )
            {
                ANGLE+= IARG[j][i]*ARG[i];
            }
            ANGLE = fmod(ANGLE, TWOPI);
            dut1 += DUT1S[j]*sin(ANGLE) + DUT1C[j]*cos(ANGLE);
            dlod += DLODS[j]*sin(ANGLE) + DLODC[j]*cos(ANGLE);
        }
        
    } // end of function UTLIBR
    
    
    
    
    
    
    /*
     
     !*+
     !*  - - - - - - - - - - -
     !*   P M S D N U T 2
     !*  - - - - - - - - - - -
     !*
     !*  This routine is part of the International Earth Rotation and
     !*  Reference Systems Service (IERS) Conventions software collection.
     !*
     !*  This subroutine evaluates the model of polar motion for
     !*  a nonrigid Earth due to tidal gravitation. This polar motion
     !*  is equivalent to the so-called "subdiurnal nutation." The model
     !*  is a sum of a first order polynomial and 25 trigonometric terms
     !*  (15 long periodic and 10 quasi diurnal) with coefficients given
     !*  in Table 5.1a of the IERS Conventions (2010).
     !*  Subroutine to compute the diurnal lunisolar effect on polar motion
     !*
     !*     :------------------------------------------:
     !*     :                                          :
     !*     :                 IMPORTANT                :
     !*     :                                          :
     !*     : In the present version this subroutine   :
     !*     : neglects the linear trend and the long   :
     !*     : periodic terms of the expansion, for the :
     !*     : reasons explained in Section 5.x.x of    :
     !*     : the IERS Conventions (2010), last para-  :
     !*     : graph before Table 5.1. If the full      :
     !*     : expansion is needed, set the parameter   :
     !*     : iband to 0 instead of 1, that is replace :
     !*     : the statement                            :
     !*     :     PARAMETER ( iband = 1 )              :
     !*     : to  PARAMETER ( iband = 0 )              :
     !*     :                                          :
     !*     :__________________________________________:
     !*
     !*  In general, Class 1, 2, and 3 models represent physical effects that
     !*  act on geodetic parameters while canonical models provide lower-level
     !*  representations or basic computations that are used by Class 1, 2, or
     !*  3 models.
     !*
     !*  Status:  Class 3 model
     !*
     !*     Class 1 models are those recommended to be used a priori in the
     !*     reduction of raw space geodetic data in order to determine
     !*     geodetic parameter estimates.
     !*     Class 2 models are those that eliminate an observational
     !*     singularity and are purely conventional in nature.
     !*     Class 3 models are those that are not required as either Class
     !*     1 or 2.
     !*     Canonical models are accepted as is and cannot be classified as
     !*     a Class 1, 2, or 3 model.
     !*
     !*  Given:
     !*     rmjd        d      Time expressed as modified Julian date
     !*
     !*  Returned:
     !*     pm          d(2)      Vector of length 2 (Note 1)
     !*
     !*  Notes:
     !*
     !*  1) The polar motion coordinates (dx, dy) are expressed in
     !*     microarcseconds.
     !*
     !*  Called:
     !*     FUNDARG             Compute the angular fundamental arguments
     !*
     !*  Test case:
     !*     given input: rmjd = 54335D0 ( August 23, 2007 )
     !*
     !*     expected output: (dx) pm(1)  = 24.65518398386097942D0 microarcsecond
     !*                      (dy) pm(2) = -14.11070254891893327D0 microarcsecond
     !*
     !*  References:
     !*
     !*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     !*     IERS Technical Note No. 36, BKG (2010)
     !*
     !*  Revisions:
     !*  2005 March       A.Brzezinski   Original code
     !*  2008 November 26 B.E.Stetzler   Initial changes to code
     !*  2008 December 01 B.E.Stetzler   Provided test case
     !*  2009 August   18 B.E.Stetzler   Capitalized all variables for FORTRAN
     !*                                  77 compatibility
     !*  2010 May      14 B.E.Stetzler   Replaced call to PMARGS to FUNDARG
     !*                                  for universal fundamental argument
     !*                                  subroutine
     !*  2010 May      17 B.E.Stetzler   Validated test case using internally
     !*                                  computed GMST and call to FUNDARG
     !*                                  matched previous external call to
     !*                                  PMARGS for all six parameters
     !*  2010 June     23 B.E.Stetzler   Modified coefficients of the long
     !*                                  and short period terms in polar
     !*                                  motion and secular polar motion
     !*                                  rate to coincide with Table 5.1a
     !*-----------------------------------------------------------------------
     */
    // http://iers-conventions.obspm.fr/2010/2010_official/chapter5/software/PMSDNUT2.F
    void GIERS::PMSDNUT2(double RMJD, double* pm )
    {
        int IBAND = 1, I =0,J =0, JSTART =0;
        double ARG[6] = {0.0}, ANGLE = 0.0;
        static double DAS2R = 4.848136811095359935899141E-6;
        static double TURNAS = 1296000E0;
        static double  RMJD0   = 51544.5E0;
        static double PI = 3.141592653589793238462643E0 ;
        static double TWOPI = 6.283185307179586476925287E0;
        static double RAD2SEC = 86400E0/TWOPI;
        
        static int IARG[25][6] = {
            0,  0, 0,  0,  0, -1,
            0, -1, 0,  1,  0,  2,
            0, -1, 0,  1,  0,  1,
            0, -1, 0,  1,  0,  0,
            0,  1, 1, -1,  0,  0,
            0,  1, 1, -1,  0, -1,
            0,  0, 0,  1, -1,  1,
            0,  1, 0,  1, -2,  1,
            0,  0, 0,  1,  0,  2,
            0,  0, 0,  1,  0,  1,
            0,  0, 0,  1,  0,  0,
            0, -1, 0,  1,  2,  1,
            0,  1, 0,  1,  0,  1,
            0,  0, 0,  3,  0,  3,
            0,  0, 0,  3,  0,  2,
            1, -1, 0, -2,  0, -1,
            1, -1, 0, -2,  0, -2,
            1,  1, 0, -2, -2, -2,
            1,  0, 0, -2,  0, -1,
            1,  0, 0, -2,  0, -2,
            1, -1, 0,  0,  0,  0,
            1,  0, 0, -2,  2, -2,
            1,  0, 0,  0,  0,  0,
            1,  0, 0,  0,  0, -1,
            1,  1, 0,  0,  0,  0
        };
        
        static double PER[25] = {6798.3837,6159.1355,3231.4956,2190.3501,438.35990,411.80661,365.24219,193.55971,27.431826,27.321582,27.212221,14.698136,13.718786,9.1071941,9.0950103,1.1196992,1.1195149,1.1134606,1.0759762,1.0758059,1.0347187,1.0027454,0.9972696,0.9971233,0.9624365};
        static double XS[25] = {  0.0,1.5,-28.5,-4.7,-0.7,1.0,1.2,1.3,-0.1,0.9,0.1,0.0,-0.1,-0.1,-0.1,-0.4,-2.3,-0.4,-2.1,-11.4,0.8,-4.8,14.3,1.9,0.8};
        static double XC[25] = {  0.6,0.0,-0.2,-0.1,0.2,0.3,0.2,0.4,-0.2,4.0,0.6,0.1,0.3,0.1,0.1,0.3,1.3,0.3,1.2,6.5,-0.5,2.7,-8.2,-1.1,-0.4};
        static double YS[25] = {  -0.1,-0.2,3.4,0.6,-0.2,-0.3,-0.2,-0.2,0.0,-0.1,0.0,0.0,0.0,0.0,0.0,-0.3,-1.3,-0.3,-1.2,-6.5,0.5,-2.7,8.2,1.1,0.4};
        static double YC[25] = { -0.1,0.1,-3.9,-0.9,-0.7,1.0,1.4,2.9,-1.7,32.4,5.1,0.6,2.7,0.9,0.6,-0.4,-2.3,-0.4,-2.1,-11.4,0.8,-4.8,14.3,1.9,0.8};
        
        double XRATE = -3.8, YRATE = -4.3;
        
        pm[0] = 0.0; pm[1] = 0.0;
        
        double T  = (RMJD-RMJD0)/36525E0;
        //!*  Compute GMST + pi
        double GMST = fmod(   67310.54841E0 + T*( (8640184.812866E0 + 3155760000E0) + T*( 0.093104E0 + T*( -0.0000062 ))), 86400.0 );
        
        double L=0.0,LP=0.0,F=0.0,D=0.0,OM=0.0;
        FUNDARG(T, L, LP, F, D, OM);
        
        ARG[0] = GMST/ RAD2SEC + PI;
        ARG[0] = fmod(ARG[0], TWOPI);
        ARG[1] = L;
        ARG[2] = LP;
        ARG[3] = F;
        ARG[4] = D;
        ARG[5] = OM;
        
        
        if(IBAND == 1 ) {JSTART = 16;}
        else {JSTART = 1;}
        
        for( J = JSTART; J<= 25; J++)
        {
            ANGLE = 0.0;
            for(I = 1; I<= 6; I++)
            {
                ANGLE += IARG[J-1][I-1]*ARG[I-1];
            }
            
            ANGLE = fmod(ANGLE, TWOPI);
            double s = sin(ANGLE);
            double c = cos(ANGLE);
            
            pm[0] +=  (XS[J-1]*s + XC[J-1]*c);
            pm[1] +=  (YS[J-1]*s + YC[J-1]*c);
        }
        
        if(IBAND == 1 ) return;
        
        pm[0] += XRATE*(RMJD - RMJD0)/365.25;
        pm[1] += YRATE*(RMJD - RMJD0)/365.25;
        
        return;
        
    } // end of function
    

    /*  IAU2000 model
     *   get the Greenwich Mean Sidereal Time
     *  JD_UT1 : Julian Date in UT1 since J2000
     *  JD_TT  : Julian Date in TT  since J2000
     */
        double GIERS::getGMST( double JD_UT1_I,double JD_UT1_F, double JD_TT)
        {
            double gmst = 0.0;
            double PI = GCONST("PI");
            double D2PI = PI+PI;
            double DAS2R = GCONST("AS2R");
            double T = JD_TT/36525.0; // Julian Century in TT since J2000
            double era = getEarthRotationAngle(JD_UT1_I, JD_UT1_F);
            double tmp =  era + ((0.014506 +(4612.15739966+(1.39667721+(-0.00009344+(0.00001882)*T)*T)*T)*T))*DAS2R;
            
            double w = fmod(tmp,D2PI);
            if(w<0) w+=D2PI;
            gmst = w;
            return gmst;
        }
    
    /*
     * ref: ERA2000.f era00.c in sofa
     * the earth rotation angle
     * input : time in Julia date, JD_UT1_I is the integer part while JD_UT1_F is the fractional part
     * output: ear in radians(0-2pi)
     */
    double GIERS::getEarthRotationAngle(double JD_UT1_I, double JD_UT1_F)
    {
        /* Reference epoch (J2000.0), Julian Date */
        // ref: sofa sofam.h
        double DJ00 = 2451545.0;
        
        double era2000 = 0.0;
        double PI = GCONST("PI");
        double D2PI = PI*2.0;
        
        double  T =0.0, F =0.0;
        T = JD_UT1_I + JD_UT1_F - DJ00;
        //F = fmod(T, 1.0);  // the fractional part
        F = fmod(JD_UT1_I,1.0) + fmod(JD_UT1_F,1.0);
        
        era2000 = D2PI * (F + 0.7790572732640 + 0.00273781191135448 * T);
        //transfer the era2000 to the range 0 to 2PI
        double w = fmod(era2000,D2PI);
        if(w<0) w+=D2PI;
        era2000 = w;
        return era2000;
    }
    
    /*Earth rotation angle first order rate.
     *  output: d(GAST)/d(t) in radians
     */
    double GIERS::getEarthRotationAngleRate(double JC_TT)
    {
        double T = JC_TT;  // Julian Centuries in TT since J2000.
        double PI = GCONST("PI");
        double D2PI = PI + PI;
        double dera = (1.002737909350795+5.9006e-11*T-5.9E-15*T*T)* D2PI / 86400.0;
        return dera;
    }

        /*
         * Greenwich Apparent Sidereal Time (model consistent with IAU 2000 resolutions).
         * inputs: JD_UT1 : Julian Date in UT1 since J2000, JD_TT : Julian Date in TT since J2000,
         *  psi: nutation in longitude (radians)
         *  return: gast(0-2pi)
         */
        double GIERS::getGAST( double JD_UT1_I, double JD_UT1_F ,double JD_TT, double dpsi )
        {
            double D2PI = GCONST("2PI");
            double gast = 0.0;
            
            double gmst = getGMST(JD_UT1_I,JD_UT1_F, JD_TT);
            
            gast  = gmst + getEE2000(JD_TT, dpsi);
            
            double w = fmod(gast, D2PI);
            if(w<0) w+= D2PI;
            gast = w;
            return gast;
        }

    // Normalize angle into the range -pi <= a < +pi.
    double GIERS::normalizeAngle(double a)
    {
        double D2PI = 2.0*GCONST("PI");
        double w = fmod(a, D2PI);
        if (fabs(w) >= (D2PI*0.5))
        {
            w-= ((a<0.0)?-D2PI:D2PI);
        }
        
        return w;
    }
    
    /*get the mean obliquity of the ecliptic in radians
     * **  Mean obliquity of the ecliptic, IAU 2006 precession model.
     * The result is the angle between the ecliptic and mean equator of
     **     date JD_TT.
     * ref: sofa iauObl06.c
     */
    double GIERS::getMeanObliquity(double JC_TT)
    {
        double DAS2R = GCONST("AS2R");
        double T = JC_TT; // Julian Centuries in TT
        
        /* Mean obliquity. */
        double eps0 = (84381.406     +
                       (-46.836769    +
                        ( -0.0001831   +
                         (  0.00200340  +
                          ( -0.000000576 +
                           ( -0.0000000434) * T) * T) * T) * T) * T) * DAS2R;
        
        return eps0;
    }

    /*
     *  Equation of the equinoxes complementary terms, consistent with
     *  IAU 2000 resolutions.
     *
     *  Annexe to IERS Conventions 2000, Chapter 5
     *  ref: EECT2000.f
     *  input: Julian Date in TT since J2000
     *
     */
    double GIERS::getEECT(double JC_TT)
    {
        double EECT2000 = 0.0;
        double PI = GCONST("PI");
        double D2PI = PI+PI;
        double DAS2R = GCONST("AS2R");
        double A =0.0, S0 =0.0, S1 =0.0, FA[14]={0.0};
        
        double T = JC_TT;  //Julian Centuries
        
        /*  -----------------------------------------
         *  The series for the EE complementary terms
         *  -----------------------------------------
         */
        int NE0 = 33, NE1 = 1;
        // the order of members in this structure is very important
        //        struct COEF
        //        {
        //            int nl,nlp,nf,nd,nom; // coefficients of l,l',F,D,Om
        //            int nlme, nlve,nle,nlma, nlju,nlsa,nlu,nln,npa;
        //        };
        
        /*  Argument coefficients for t^0*/
        double KE0[][14] =
        {
            {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1},
            {  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {  1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}
        };
        
        /*Argument coefficients for t^1*/
        double KE1[14] = {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
        
        /*Sine and cosine coefficients for t^0*/
        double SE0[][2] =
        {
            {            +2640.96E-6,          -0.39E-6},
            {              +63.52E-6,          -0.02E-6},
            {              +11.75E-6,          +0.01E-6},
            {              +11.21E-6,          +0.01E-6},
            {               -4.55E-6,          +0.00E-6},
            {               +2.02E-6,          +0.00E-6},
            {               +1.98E-6,          +0.00E-6},
            {               -1.72E-6,          +0.00E-6},
            {               -1.41E-6,          -0.01E-6},
            {               -1.26E-6,          -0.01E-6},
            {               -0.63E-6,          +0.00E-6},
            {               -0.63E-6,          +0.00E-6},
            {               +0.46E-6,          +0.00E-6},
            {               +0.45E-6,          +0.00E-6},
            {               +0.36E-6,          +0.00E-6},
            {               -0.24E-6,          -0.12E-6},
            {               +0.32E-6,          +0.00E-6},
            {               +0.28E-6,          +0.00E-6},
            {               +0.27E-6,          +0.00E-6},
            {               +0.26E-6,          +0.00E-6},
            {               -0.21E-6,          +0.00E-6},
            {               +0.19E-6,          +0.00E-6},
            {               +0.18E-6,          +0.00E-6},
            {               -0.10E-6,          +0.05E-6},
            {               +0.15E-6,          +0.00E-6},
            {               -0.14E-6,          +0.00E-6},
            {               +0.14E-6,          +0.00E-6},
            {               -0.14E-6,          +0.00E-6},
            {               +0.14E-6,          +0.00E-6},
            {               +0.13E-6,          +0.00E-6},
            {               -0.11E-6,          +0.00E-6},
            {               +0.11E-6,          +0.00E-6},
            {               +0.11E-6,          +0.00E-6}
        };
        
        /*Sine and cosine coefficients for t^1*/
        double SE1[2] ={-0.87E-6,          +0.00E-6};
        
        /*Mean Anomaly of the Moon*/
        FA[0] = normalizeAngle((485868.249036E0 +(715923.2178E0 +(31.8792E0 +( 0.051635E0 +( -0.00024470E0 )*T)*T)*T)*T)*DAS2R+fmod(1325.0*T, 1.0)*D2PI);
        
        /*Mean Anomaly of the Sun*/
        FA[1] = normalizeAngle((1287104.793048E0+( 1292581.0481E0+(-0.5532E0+(+0.000136E0+(-0.00001149E0)*T)*T)*T)*T)*DAS2R+fmod(99E0*T, 1.0E0)*D2PI);
        
        /*Mean Longitude of the Moon minus Mean Longitude of the Ascending*/
        FA[2] = normalizeAngle((335779.526232E0+(295262.8478E0 +(-12.7512E0 +(-0.001037E0 +(0.00000417E0)*T)*T)*T)*T)*DAS2R+ fmod( 1342E0*T, 1E0 )*D2PI);
        
        /*Mean Elongation of the Moon from the Sun*/
        FA[3]=normalizeAngle((1072260.703692E0 +( 1105601.2090E0+(-6.3706E0 +( 0.006593E0 +(-0.00003169E0 )*T)*T)*T)*T)*DAS2R+fmod(1236E0*T,1E0)*D2PI);
        
        /* Mean Longitude of the Ascending Node of the Moon*/
        FA[4]=normalizeAngle((450160.398036E0 +( -482890.5431E0 +(7.4722E0+(0.007702E0 +(-0.00005939E0)*T)*T)*T)*T)*DAS2R+fmod( -5E0*T,1E0)*D2PI);
        
        FA[5] = normalizeAngle ( 4.402608842E0 + 2608.7903141574E0 * T );
        FA[6] = normalizeAngle ( 3.176146697E0 + 1021.3285546211E0 * T );
        FA[7] = normalizeAngle ( 1.753470314E0 +  628.3075849991E0 * T );
        FA[8] = normalizeAngle ( 6.203480913E0 +  334.0612426700E0 * T );
        FA[9] = normalizeAngle ( 0.599546497E0 +   52.9690962641E0 * T );
        FA[10] = normalizeAngle ( 0.874016757E0 +   21.3299104960E0 * T );
        FA[11] = normalizeAngle ( 5.481293872E0 +    7.4781598567E0 * T );
        FA[12] = normalizeAngle ( 5.311886287E0 +    3.8133035638E0 * T );
        FA[13] =          ( 0.024381750E0 +    0.00000538691E0 * T ) * T;
        
        S0 = 0.0;
        S1 = 0.0;
        for( int i = NE0-1 ; i>0; i--)
        {
            A = 0.0;
            for( int j =0; j< 14; j++ )
            {
                A = A + KE0[i][j]*FA[j];
            }
            S0 = S0 + ( SE0[i][0]*sin(A) + SE0[i][1]*cos(A) );
        }
        
        A =0.0;
        for( int j =0; j< 14; j++ )
        {
            A = A + KE1[j]*FA[j];
        }
        S1 = S1 + ( SE1[0]*sin(A) + SE1[1]*cos(A) );
        
        EECT2000 = ( S0 + S1 * T ) * DAS2R;
        
        return EECT2000;
    }

    /*
     * inputs: dpsi: nutation in longitude (radians) , JD_TT: Julian Date in TT since J2000
     * The equation of the equinoxes, compatible with IAU 2000 resolutions,
     * ref: EE2000.f
     */
    double GIERS::getEE2000( double JC_TT, double dpsi)
    {
        double ee2000 = 0.0;
        //double PI = GCONST("PI");
        //double DAS2R = GCONST("AS2R");
        //double T = JD_TT/36525.0; // Julian Centuries in TT
        
        double EPSA = getMeanObliquity(JC_TT);
        //double EPS0 = 84381.448E0 * DAS2R; // J2000 obliquity(Lieske et al. 1977)
        //Mean obliquity from Chapter 5, expression (32)
        //double EPSA = EPS0 + (-46.8402E0 +(-0.00059E0 +(0.001813E0)*T)*T)*T*DAS2R;
        
        ee2000 = dpsi*cos(EPSA) + getEECT(JC_TT);
        return ee2000;
    }
    
    /*
     ref: IERS ARG2.f
     *+
     *  - - - - - - - - -
     *   A R G 2
     *  - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  The purpose of the subroutine is to compute the angular astronomical
     *  argument, which depends on time, for 11 tidal argument calculations.
     *  The order of the 11 angular quantities in vector angle are given below:
     *
     *  01-M2, 02-S2, 03-N2, 04-K2, 05-K1, 06-O1, 07-P1, 08-Q1, 09-Mf,
     *  10-Mm, 11-Ssa (See Reference 1)
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Canonical model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     IYEAR          i      Four digit year (Note 1)
     *     DAY            d      Day of Year Greenwich Time (Note 2)
     *
     *  Returned:
     *     ANGLE(K)       d      Angular argument for Schwiderski
     *                           computation, in radians (Notes 3 and 4)
     *
     *  Notes:
     *
     *  1) This subroutine is valid only after 1973 CE.  A validation
     *     test has been added to stop the subroutine if an invalid
     *     year is used.
     *
     *  2) Example: 32.5 for February 1 12 Noon
     *     Example: 1.25 for January 1 6 AM
     *
     *  3) Ocean loading phases computed from Schwiderski's models
     *     refer to the phase of the associated solid Earth tide
     *     generating potential at the zero meridian according to
     *
     *      OL_DR = OL_AMP ' COS (SE_PHASE" - OL_PHASE)
     *
     *     where OL = OCEAN LOADING TIDE,
     *           SE = SOLID EARTH TIDE GENERATING POTENTIAL.
     *
     *     If the harmonic tide development of Cartwright, et al.
     *     (CTE) (1971, 1973) is used, make sure that SE_PHASE"
     *     take into account:
     *
     *     (1) the sign of SE_AMP in the tables of Cartwright et al.
     *     (2) that CTE'S SE_PHASE refers to a sine rather than a
     *     cosine function if (N+M) = (DEGREE + ORDER) of the tide
     *     spherical harmonic is odd.
     *
     *     i.e. SE_PHASE" = TAU(T) ' N1 + S(T) ' N2 + H(T) ' N3
     *                 + P(T) ' N4 + N'(T) ' N5 + PS(T) ' N6
     *                 + PI   If CTE'S amplitude coefficient < 0
     *                 + PI/2 If (DEGREE + N1) is odd
     *
     *     where TAU ... PS = astronomical arguments,
     *           N1 ... N6 = CTE'S argument numbers.
     *
     *     Most tide generating software compute SE_PHASE" (for use
     *     with cosines).
     *
     *  4) The double precision change from the original routine ARG.f
     *     to ARG2.F yields output differences on the order of 10^-9 radians.
     *
     *  Called:
     *     None
     *
     *  Test case:
     *     given input: IYEAR = 2008
     *                  DAY = 311.5 (November 6 Noon)
     *     expected output: ANGLE(1)  = 2.849663065753787805D0  rad
     *                      ANGLE(2)  = 6.28318080000000023D0   rad
     *                      ANGLE(3)  = 4.926040134021299366D0  rad
     *                      ANGLE(4)  = 1.608450491115348768D0  rad
     *                      ANGLE(5)  = 2.375021572352622456D0  rad
     *                      ANGLE(6)  = 0.4746414933980958040D0 rad
     *                      ANGLE(7)  = 3.908159227647345801D0  rad
     *                      ANGLE(8)  = 2.551018561669245344D0  rad
     *                      ANGLE(9)  = 5.041990012540757959D0  rad
     *                      ANGLE(10) = 4.206816878908014701D0  rad
     *                      ANGLE(11) = 1.608463638294885811D0  rad
     *
     *  References:
     *
     *     Schwiderski, E., 1983, "Atlas of Ocean Tidal Charts and Maps, Part I:
     *     The Semidiurnal Principal Lunar Tide M2," Marine Geodesy, 6, pp. 219-256.
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  2008 November  05 B.E. Stetzler    Added header and copyright
     *  2008 November  07 B.E. Stetzler    Provided test case
     *  2008 November  13 B.E. Stetzler    Re-defined variables
     *  2008 November  14 B.E. Stetzler    Corrected test case results
     *  2010 April     20 B.E. Stetzler    New version of routine without the
     *                                     use of EQUIVALENCE
     *  2011 September 20 B.E. Stetzler    Changed IYEAR variable from two
     *                                     digit value to four digit value,
     *                                     added year validation and IYMIN
     *                                     variable, and updated test case
     *  2011 October   07 B.E. Stetzler    Removed IYEAR-1900 from code and
     *                                     modified ICAPD calculation 
     *-----------------------------------------------------------------------
     */
    void GIERS::ARG2( int IYEAR, double DAY, double *ANGLE)
    {
        int I, ID, ICAPD, K = 11, IYMIN = 1974;
        double FDAY, CAPT,H0,S0,P0,DTR;
        double SIGM2, SIGS2, SIGN2, SIGK2, SIGK1, SIGO1, SIGP1, SIGQ1, SIGMF, SIGMM, SIGSSA;
        
        double SPEED[11] = { 1.40519E-4, 1.45444E-4,1.37880E-4,1.45842E-4,0.72921E-4,
                            0.67598E-4,0.72523E-4,0.64959E-4 ,0.053234E-4,0.026392E-4,0.003982E-4};
        
        SIGM2 = 1.40519E-4; SIGS2 = 1.45444E-4; SIGN2 = 1.37880E-4; SIGK2 = 1.45842E-4; SIGK1 = .72921E-4;
        SIGO1 = 0.67598E-4 ; SIGP1 = 0.72523E-4; SIGQ1 = 0.64959E-4; SIGMF = 0.053234E-4; SIGMM = 0.026392E-4;
        SIGSSA = 0.003982E-4;
        
        double ANGFAC[11][4] = {
            {2.E0,-2.E0,0.E0,0.E0},
            {0.E0,0.E0,0.E0,0.E0},
            {2.E0,-3.E0,1.E0,0.E0},
            {2.E0,0.E0,0.E0,0.E0},
            {1.E0,0.E0,0.E0,.25E0},
            {1.E0,-2.E0,0.E0,-.25E0},
            {-1.E0,0.E0,0.E0,-.25E0},
            {1.E0,-3.E0,1.E0,-.25E0},
            {0.E0,2.E0,0.E0,0.E0},
            {0.E0,1.E0,-1.E0,0.E0},
            {2.E0,0.E0,0.E0,0.E0}
            };
        
        double TWOPI = 6.283185307179586476925287;
        DTR = 0.174532925199E-1;
        
        /*  Validate year */
        if(IYEAR < IYMIN )
        {
            printf("ERROR: IYEAR must be larger than 1974\n");
        }
        
        /*  Initialize day of year(integer part) */
        ID = (int)DAY;
        
        /* ------------------------------------------
        *  Compute fractional part of day in seconds
        * ------------------------------------------
        */
        FDAY = (DAY-ID)*86400.0;
        
        /* Revision 07 October 2011: ICAPD modified*/
        ICAPD = ID+365*(IYEAR-1975)+((IYEAR-1973)/4);
        CAPT = (27392.500528+1.000000035*ICAPD)/36525.0;
        
        /*  Compute mean longitude of Sun at beginning of day*/
        H0=(279.69668+(36000.768930485+3.03E-4*CAPT)*CAPT)*DTR;
        /*  Compute mean longitude of Moon at beginning of day */
        S0=(((1.9E-6*CAPT-0.001133)*CAPT+481267.88314137)*CAPT +270.434358)*DTR;
        
        /*  Compute mean longitude of lunar perigee at beginning of day */
        P0=(((-1.2E-5*CAPT-.010325)*CAPT+4069.0340329577)*CAPT+334.329653)*DTR;
        
        /* Compute the tidal angle arguments*/
        
        for(I = 1; I<=K; I++ )
        {
            ANGLE[I-1] = SPEED[I-1]*FDAY + ANGFAC[I-1][0]*H0 + ANGFAC[I-1][1]*S0
            + ANGFAC[I-1][2]*P0 + ANGFAC[I-1][3]*TWOPI;
            
            ANGLE[I-1] = fmod(ANGLE[I-1],TWOPI);
            if(ANGLE[I-1] < 0.0)
            {
                ANGLE[I-1] += TWOPI;
            }
        }
        
        
        
    } // end of function ARG2
    
    
    /*
     
     *+
     *  - - - - - - - - - - -
     *   R G _ Z O N T 2
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine evaluates the effects of zonal Earth tides on the
     *  rotation of the Earth.  The model used is a combination of Yoder
     *  et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
     *  body tide, and Kantha et al. (1998) ocean tide models
     *  as recommended by the IERS Conventions (2010).  Refer to
     *  Chapter 8 pp. xx - xx.  The latest version of the model is located
     *  at http://tai.bipm.org/iers/convupdt/convupdt_c8.html.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status:  Class 3 model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as
     *     a Class 1, 2, or 3 model.
     *
     *  Given:
     *     T           d      TT, Julian centuries since J2000 (Note 1)
     *
     *  Returned:
     *     DUT         d      Effect on UT1 (Note 2)
     *     DLOD        d      Effect on excess length of day (LOD) (Note 3)
     *     DOMEGA      d      Effect on rotational speed (Note 4)
     *
     *  Notes:
     *
     *  1) Though T is strictly TDB, it is usually more convenient to use
     *     TT, which makes no significant difference.  Julian centuries since
     *     J2000 is (JD - 2451545.0)/36525.
     *
     *  2) The expression used is as adopted in IERS Conventions (2010).
     *     DUT is expressed in seconds and is double precision.
     *
     *  3) The expression used is as adopted in IERS Conventions (2010).
     *     DLOD is the excess in LOD and is expressed in seconds per day
     *     and is double precision.  The phrase 'per day' is generally
     *     understood, so it has been omitted commonly in speech and
     *     literature.
     *     See: Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
     *     "Long-Term Changes in the Rotation of the Earth: 700 B. C. to
     *     A. D. 1980 [and Discussion]", Phil. Trans. Roy. Soc. of London.
     *     Series A, 313, pp. 47 - 70.
     *
     *  4) The expression used is as adopted in IERS Conventions (2010).
     *     Rotational speed is expressed in radians per second and is
     *     double precision.
     *
     *  Called:
     *     FUNDARG      Computation of the fundamental lunisolar arguments
     *
     *  Test case:
     *     given input: T = .07995893223819302 Julian centuries since J2000
     *                  (MJD = 54465)
     *     expected output: DUT    =  7.983287678576557467E-002 seconds
     *                      DLOD   =  5.035303035410713729E-005 seconds / day
     *                      DOMEGA = -4.249711616463017E-014 radians / second
     *
     *  References:
     *
     *     Yoder, C. F., Williams, J. G., and Parke, M. E., (1981),
     *     "Tidal Variations of Earth Rotation," J. Geophys. Res., 86,
     *     pp. 881 - 891.
     *
     *     Wahr, J. and Bergen, Z., (1986), "The effects of mantle
     *     anelasticity on nutations, Earth tides, and tidal variations
     *     in rotation rate," Geophys. J. Roy. astr. Soc., 87, pp. 633 - 668.
     *
     *     Kantha, L. H., Stewart, J. S., and Desai, S. D., (1998), "Long-
     *     period lunar fortnightly and monthly ocean tides," J. Geophys.
     *     Res., 103, pp. 12639 - 12647.
     *
     *     Gross, R. S., (2009), "Ocean tidal effects on Earth rotation,"
     *     J. Geodyn., 48(3-5), pp. 219 - 225.
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. xx, BKG (to be issued 2010)
     *
     *  Revisions:
     *  2008 January 18 B.E. Stetzler  Initial changes to header
     *               and used 2PI instead of PI as parameter
     *  2008 January 25 B.E. Stetzler Additional changes to header
     *  2008 February 21 B.E. Stetzler Definition of (excess) LOD clarified
     *  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
     *  2008 March   14 B.E. Stetzler Further changes applied to code.
     *  2008 April   03 B.E. Stetzler Provided example test case
     *  2009 February 11 B.E. Stetzler Updated test case due to changes made
     *                                 to FUNDARG.F subroutine
     *  2009 April   10 B.E. Stetzler DLOD corrected to say it is expressed
     *                                in seconds per day
     *  2009 May     04 B.E. Stetzler Code formatting changes based on
     *                                client recommendations
     *  2009 May     07 B.E. Stetzler Updated test case due to above changes
     *  2010 February 19 B.E. Stetzler Replaced Conventions 2003 recommended
     *                                 model with Conventions 2010 model
     *  2010 February 22 B.E. Stetzler Provided example test case
     *  2010 February 23 B.E. Stetzler Updated values to two decimal places
     *  2010 February 23 B.E. Stetzler Split fundamental arguments and
     *                                 coefficients for four decimal place
     *                                 precision
     *  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
     *  2010 March    01 B.E. Stetzler Updated table values to four decimal
     *                                 places and double precision
     *  2010 March    12 B.E. Stetzler Applied changes to wording of notes.
     *  2010 March    22 B.E. Stetzler Corrected DOMEGA output for test case
     
    */
    void GIERS::RG_ZONT2(double JC_TT, double &DUT, double &DLOD, double &DOMEGA)
    {
        
        double T = JC_TT;
        
        int I,J;
        double L,LP,F,D,OM, ARG;
        double D2PI = 6.283185307179586476925287;
        double TURNAS = 1296000;
        double DAS2R = 4.848136811095359935899141E-6 ;
        /*  Number of terms in the zonal Earth tide model */
        int NZONT = 62;
        double NFUND[62][5] =
        {
            { 1,  0,  2,  2,  2},
            { 2,  0,  2,  0,  1},
            { 2,  0,  2,  0,  2},
            { 0,  0,  2,  2,  1},
            { 0,  0,  2,  2,  2},
            { 1,  0,  2,  0,  0},
            { 1,  0,  2,  0,  1},
            { 1,  0,  2,  0,  2},
            { 3,  0,  0,  0,  0},
            {-1,  0,  2,  2,  1},
            {-1,  0,  2,  2,  2},
            { 1,  0,  0,  2,  0},
            { 2,  0,  2, -2,  2},
            { 0,  1,  2,  0,  2},
            { 0,  0,  2,  0,  0},
            { 0,  0,  2,  0,  1},
            { 0,  0,  2,  0,  2},
            { 2,  0,  0,  0, -1},
            { 2,  0,  0,  0,  0},
            { 2,  0,  0,  0,  1},
            { 0, -1,  2,  0,  2},
            { 0,  0,  0,  2, -1},
            { 0,  0,  0,  2,  0},
            { 0,  0,  0,  2,  1},
            { 0, -1,  0,  2,  0},
            { 1,  0,  2, -2,  1},
            { 1,  0,  2, -2,  2},
            { 1,  1,  0,  0,  0},
            {-1,  0,  2,  0,  0},
            {-1,  0,  2,  0,  1},
            {-1,  0,  2,  0,  2},
            { 1,  0,  0,  0, -1},
            { 1,  0,  0,  0,  0},
            { 1,  0,  0,  0,  1},
            { 0,  0,  0,  1,  0},
            { 1, -1,  0,  0,  0},
            {-1,  0,  0,  2, -1},
            {-1,  0,  0,  2,  0},
            {-1,  0,  0,  2,  1},
            { 1,  0, -2,  2, -1},
            {-1, -1,  0,  2,  0},
            { 0,  2,  2, -2,  2},
            { 0,  1,  2, -2,  1},
            { 0,  1,  2, -2,  2},
            { 0,  0,  2, -2,  0},
            { 0,  0,  2, -2,  1},
            { 0,  0,  2, -2,  2},
            { 0,  2,  0,  0,  0},
            { 2,  0,  0, -2, -1},
            { 2,  0,  0, -2,  0},
            { 2,  0,  0, -2,  1}, 
            { 0, -1,  2, -2,  1}, 
            { 0,  1,  0,  0, -1}, 
            { 0, -1,  2, -2,  2}, 
            { 0,  1,  0,  0,  0}, 
            { 0,  1,  0,  0,  1}, 
            { 1,  0,  0, -1,  0}, 
            { 2,  0, -2,  0,  0}, 
            {-2,  0,  2,  0,  1}, 
            {-1,  1,  0,  1,  0},
            { 0,  0,  0,  0,  2}, 
            { 0,  0,  0,  0,  1}
        };
        
        double TIDE[62][6] =
        {
            {    -0.0235E0,0.0000E0, 0.2617E0, 0.0000E0, -0.2209E0, 0.0000E0},
            {    -0.0404E0,0.0000E0, 0.3706E0, 0.0000E0, -0.3128E0, 0.0000E0},
            {    -0.0987E0,0.0000E0, 0.9041E0, 0.0000E0, -0.7630E0, 0.0000E0},
            {    -0.0508E0,0.0000E0, 0.4499E0, 0.0000E0, -0.3797E0, 0.0000E0},
            {    -0.1231E0,0.0000E0, 1.0904E0, 0.0000E0, -0.9203E0, 0.0000E0},
            {    -0.0385E0,0.0000E0, 0.2659E0, 0.0000E0, -0.2244E0, 0.0000E0},
            {    -0.4108E0,0.0000E0, 2.8298E0, 0.0000E0, -2.3884E0, 0.0000E0},
            {    -0.9926E0,0.0000E0, 6.8291E0, 0.0000E0, -5.7637E0, 0.0000E0},
            {    -0.0179E0,0.0000E0, 0.1222E0, 0.0000E0, -0.1031E0, 0.0000E0},
            {    -0.0818E0,0.0000E0, 0.5384E0, 0.0000E0, -0.4544E0, 0.0000E0},
            {    -0.1974E0,0.0000E0, 1.2978E0, 0.0000E0, -1.0953E0, 0.0000E0},
            {    -0.0761E0,0.0000E0, 0.4976E0, 0.0000E0, -0.4200E0, 0.0000E0},
            {     0.0216E0,0.0000E0,-0.1060E0, 0.0000E0,  0.0895E0, 0.0000E0},
            {     0.0254E0,0.0000E0,-0.1211E0, 0.0000E0,  0.1022E0, 0.0000E0},
            {    -0.2989E0,0.0000E0, 1.3804E0, 0.0000E0, -1.1650E0, 0.0000E0},
            {    -3.1873E0,0.2010E0,14.6890E0, 0.9266E0,-12.3974E0,-0.7820E0},
            {    -7.8468E0,0.5320E0,36.0910E0, 2.4469E0,-30.4606E0,-2.0652E0},
            {     0.0216E0,0.0000E0,-0.0988E0, 0.0000E0,  0.0834E0, 0.0000E0},
            {    -0.3384E0,0.0000E0, 1.5433E0, 0.0000E0, -1.3025E0, 0.0000E0},
            {     0.0179E0,0.0000E0,-0.0813E0, 0.0000E0,  0.0686E0, 0.0000E0},
            {    -0.0244E0,0.0000E0, 0.1082E0, 0.0000E0, -0.0913E0, 0.0000E0},
            {     0.0470E0,0.0000E0,-0.2004E0, 0.0000E0,  0.1692E0, 0.0000E0},
            {    -0.7341E0,0.0000E0, 3.1240E0, 0.0000E0, -2.6367E0, 0.0000E0},
            {    -0.0526E0,0.0000E0, 0.2235E0, 0.0000E0, -0.1886E0, 0.0000E0},
            {    -0.0508E0,0.0000E0, 0.2073E0, 0.0000E0, -0.1749E0, 0.0000E0},
            {     0.0498E0,0.0000E0,-0.1312E0, 0.0000E0,  0.1107E0, 0.0000E0},
            {     0.1006E0,0.0000E0,-0.2640E0, 0.0000E0,  0.2228E0, 0.0000E0},
            {     0.0395E0,0.0000E0,-0.0968E0, 0.0000E0,  0.0817E0, 0.0000E0},
            {     0.0470E0,0.0000E0,-0.1099E0, 0.0000E0,  0.0927E0, 0.0000E0},
            {     0.1767E0,0.0000E0,-0.4115E0, 0.0000E0,  0.3473E0, 0.0000E0},
            {     0.4352E0,0.0000E0,-1.0093E0, 0.0000E0,  0.8519E0, 0.0000E0},
            {     0.5339E0,0.0000E0,-1.2224E0, 0.0000E0,  1.0317E0, 0.0000E0},
            {    -8.4046E0,0.2500E0,19.1647E0, 0.5701E0,-16.1749E0,-0.4811E0},
            {     0.5443E0,0.0000E0,-1.2360E0, 0.0000E0,  1.0432E0, 0.0000E0},
            {     0.0470E0,0.0000E0,-0.1000E0, 0.0000E0,  0.0844E0, 0.0000E0},
            {    -0.0555E0,0.0000E0, 0.1169E0, 0.0000E0, -0.0987E0, 0.0000E0},
            {     0.1175E0,0.0000E0,-0.2332E0, 0.0000E0,  0.1968E0, 0.0000E0},
            {    -1.8236E0,0.0000E0, 3.6018E0, 0.0000E0, -3.0399E0, 0.0000E0},
            {     0.1316E0,0.0000E0,-0.2587E0, 0.0000E0,  0.2183E0, 0.0000E0},
            {     0.0179E0,0.0000E0,-0.0344E0, 0.0000E0,  0.0290E0, 0.0000E0},
            {    -0.0855E0,0.0000E0, 0.1542E0, 0.0000E0, -0.1302E0, 0.0000E0},
            {    -0.0573E0,0.0000E0, 0.0395E0, 0.0000E0, -0.0333E0, 0.0000E0},
            {     0.0329E0,0.0000E0,-0.0173E0, 0.0000E0,  0.0146E0, 0.0000E0},
            {    -1.8847E0,0.0000E0, 0.9726E0, 0.0000E0, -0.8209E0, 0.0000E0},
            {     0.2510E0,0.0000E0,-0.0910E0, 0.0000E0,  0.0768E0, 0.0000E0},
            {     1.1703E0,0.0000E0,-0.4135E0, 0.0000E0,  0.3490E0, 0.0000E0},
            {   -49.7174E0,0.4330E0,17.1056E0, 0.1490E0,-14.4370E0,-0.1257E0},
            {    -0.1936E0,0.0000E0, 0.0666E0, 0.0000E0, -0.0562E0, 0.0000E0},
            {     0.0489E0,0.0000E0,-0.0154E0, 0.0000E0,  0.0130E0, 0.0000E0},
            {    -0.5471E0,0.0000E0, 0.1670E0, 0.0000E0, -0.1409E0, 0.0000E0},
            {     0.0367E0,0.0000E0,-0.0108E0, 0.0000E0,  0.0092E0, 0.0000E0},
            {    -0.0451E0,0.0000E0, 0.0082E0, 0.0000E0, -0.0069E0, 0.0000E0},
            {     0.0921E0,0.0000E0,-0.0167E0, 0.0000E0,  0.0141E0, 0.0000E0},
            {     0.8281E0,0.0000E0,-0.1425E0, 0.0000E0,  0.1202E0, 0.0000E0},
            {   -15.8887E0,0.1530E0, 2.7332E0, 0.0267E0, -2.3068E0,-0.0222E0},
            {    -0.1382E0,0.0000E0, 0.0225E0, 0.0000E0, -0.0190E0, 0.0000E0},
            {     0.0348E0,0.0000E0,-0.0053E0, 0.0000E0,  0.0045E0, 0.0000E0},
            {    -0.1372E0,0.0000E0,-0.0079E0, 0.0000E0,  0.0066E0, 0.0000E0},
            {     0.4211E0,0.0000E0,-0.0203E0, 0.0000E0,  0.0171E0, 0.0000E0},
            {    -0.0404E0,0.0000E0, 0.0008E0, 0.0000E0, -0.0007E0, 0.0000E0},
            {     7.8998E0,0.0000E0, 0.1460E0, 0.0000E0, -0.1232E0, 0.0000E0},
            { -1617.2681E0,0.0000E0,-14.9471E0,0.0000E0, 12.6153E0, 0.0000E0}
        };
        
        /*   Computation of fundamental arguments */
        FUNDARG(T, L, LP, F, D, OM);
        /*  Set initial values to zero. */
        DUT    = 0.0;
        DLOD   = 0.0;
        DOMEGA = 0.0;
        
        /*  Sum zonal tide terms. */
        for(I =1 ; I <= NZONT; I++ )
        {
            
            double temp = NFUND[I-1][0] * L+  NFUND[I-1][1] * LP+  NFUND[I-1][2]  * F+  NFUND[I-1][3]  * D+  NFUND[I-1][4]  * OM;
            /*     Formation of multiples of arguments. */
            ARG = fmod( temp, D2PI );
            
            if (ARG < 0.0) ARG = ARG + D2PI;
            
            /*     Evaluate zonal tidal terms. */
            DUT    = DUT    + TIDE[I-1][0] *sin(ARG) + TIDE[I-1][1] *cos(ARG);
            DLOD   = DLOD   + TIDE[I-1][2] *cos(ARG) + TIDE[I-1][3] *sin(ARG);
            DOMEGA = DOMEGA + TIDE[I-1][4] *cos(ARG) + TIDE[I-1][5] *sin(ARG);
            
        }
        
        /*  Rescale corrections so that they are in units involving seconds. */
        
        DUT    = DUT    * 1.0E-4;
        DLOD   = DLOD   * 1.0E-5;
        DOMEGA = DOMEGA * 1.0E-14;
        
    } // end of function RG_ZONT2
    
    
    /*
     
     *  - - - - - - - - - - -
     *   M D A Y
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This function finds the day number of days before start of month m,
     *  of year iy, in Gregorian intercalation.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Canonical model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     iy           i      a year
     *     m            i      a month
     *
     *  Returned:
     *     mday         i      day number of day before start of a month
     *
     *  Notes:
     *
     *  1)  This function needs to test for a leap year.
     *
     *  Called:
     *     None
     *
     *  Test case:  This is a support function of the main program HARDISP.F.
     *     given input: iy = 2009
     *                   m = 5
     *     expected output: mday = 120
     *
     *  References:
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  2009 April 22 B.E. Stetzler Initial standardization of function
     *                              and provided a test case 
     *  2009 July  29 B.E. Stetzler Capitalized all variables for FORTRAN 77
     *                              compatibility 
     *-----------------------------------------------------------------------
     
     */
    int GIERS::MDAY( int iy, int im)
    {
        int dayN = -1 , LEAP = -1;
        LEAP = 1 - (iy%4+3)/4;
        if( iy%100 == 0 && iy%400 != 0) LEAP=0;
        dayN =  (367*(im-2-12*((im-14)/12)))/12+29%365 + LEAP*((9+im)/12);
        return dayN;
    } // end of function MDAY
    
    
    /*  - - - - - - - - - - -
     *   L E A P
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This function determines whether a given integer year is a leap
     *  year.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Canonical model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     iy           i      a year (Note 1)
     *
     *  Returned:
     *     0            i      if year is not a leap year
     *     1            i      if year is a leap year
     *
     *  Notes:
     *
     *  1) The year is a Gregorian year.
     *
     *  Called:
     *     None
     *
     *  Test case:  This is a support function of the main program HARDISP.F.
     *     given input: IY = 2009
     *
     *     expected output: 0
     *
     *  References:
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  2009 April  20 B.E.Stetzler  Initial standardization of function
     *                               and provided a test case 
     *  2009 August 19 B.E.Stetzler  Capitalized all variables for FORTRAN
     *                               77 compatibility
     *-----------------------------------------------------------------------
     */
    bool GIERS::LEAP(int iy)
    {
        int leap = 1 - ( iy%4 +3)/4;
        if( iy%100 == 0&& iy%400 != 0 ) leap = 0;
        
        return leap;
    } // end of function LEAP
    
    
    
    /*  - - - - - - - - - - -
    *   T O Y M D
    *  - - - - - - - - - - -
    *
    *  This routine is part of the International Earth Rotation and
    *  Reference Systems Service (IERS) Conventions software collection.
    *
    *  This subroutine converts times given in it1 expressed in year and
    *  day of year to year-month-day in it2.
    *
    *  In general, Class 1, 2, and 3 models represent physical effects that
    *  act on geodetic parameters while canonical models provide lower-level
        *  representations or basic computations that are used by Class 1, 2, or
        *  3 models.
        *
        *  Status: Canonical model
        *
        *     Class 1 models are those recommended to be used a priori in the
        *     reduction of raw space geodetic data in order to determine
        *     geodetic parameter estimates.
        *     Class 2 models are those that eliminate an observational
        *     singularity and are purely conventional in nature.
        *     Class 3 models are those that are not required as either Class
        *     1 or 2.
        *     Canonical models are accepted as is and cannot be classified as a
        *     Class 1, 2, or 3 model.
        *
        *  Given:
        *     it1           i(2)      time given in year and day of year (Note 1)
        *
        *  Returned:
        *     it2           i(3)      time given in year-month-day format
        *
        *  Notes:
        *
        *  1) The time is split into a year, given as it1(1) and the day of the
        *     year, given as it1(2).
        *
        *  Called:
        *    LEAP
        *
        *  Test case:
        *    Given input:  it1(1) = 2008
        *                  it1(2) = 120
        *
        *    Expected output: it2(1) = 2008
        *                     it2(2) = 4
        *                     it2(3) = 29
        *
        *  References:
        *
        *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
        *     IERS Technical Note No. 36, BKG (2010)
        *
        *  Revisions:
        *  2009 April   20 B.E.Stetzler Initial standardization of subroutine 
        *  2009 April   23 B.E.Stetzler Provided test case
        *  2009 August  19 B.E.Stetzler Capitalized all variables for FORTRAN
        *                               77 compatibility
        *  ref: ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/hardisp/TOYMD.F
        */
    void GIERS::TOYMD( int IT1[2], int IT2[3])
    {
        
        /*
         IMPLICIT NONE
         INTEGER IDN,IT1,IT2,JJ,M,MON,LEAP
         DIMENSION IT1(*),IT2(*)
         
         IDN(M) = MOD((367*(M-2-12*((M-14)/12)))/12+29,365)
         MON(JJ,M) = (12*(JJ-29-M))/367 + 2 + (JJ-200)/169
         IT2(1) = IT1(1)
         IT2(2) = MON(IT1(2),LEAP(IT1(1)))
         IT2(3) = IT1(2) - IDN(IT2(2)) - LEAP(IT2(1))*((9+IT2(2))/12)
         
         RETURN
         */
        int IDN, JJ,M,MON,LEAP;
        
        IDN = (367*(M-2-12*((M-14)/12)))/12+29%365;
        
        
    } // end of function TOYMD
    
    /*
    *+
    *  - - - - - - - - - - - -
    *   F C U L Z D _ H P A
    *  - - - - - - - - - - - -
    *
    *  This routine is part of the International Earth Rotation and
    *  Reference Systems Service (IERS) Conventions software collection.
    *
    *  This subroutine determines the total zenith delay following (Mendes and Pavlis, 2004).
    *
    *  In general, Class 1, 2, and 3 models represent physical effects that
    *  act on geodetic parameters while canonical models provide lower-level
        *  representations or basic computations that are used by Class 1, 2, or
        *  3 models.
        *
        *  Status: Class 1 model
        *
        *     Class 1 models are those recommended to be used a priori in the
        *     reduction of raw space geodetic data in order to determine
        *     geodetic parameter estimates.
        *     Class 2 models are those that eliminate an observational
        *     singularity and are purely conventional in nature.
        *     Class 3 models are those that are not required as either Class
        *     1 or 2.
        *     Canonical models are accepted as is and cannot be classified as a
        *     Class 1, 2, or 3 model.
        *
        *  Given:
        *     LATITUDE       d      Geodetic Latitude given in degrees (North Latitude)
        *     ELLIP_HT       d      Height above ellipsoid given in meters
        *     PRESSURE       d      Surface pressure given in hPa (mbars) (Note 1)
        *     WVP            d      Water vapor pressure in hPa (mbars) (Note 1)
        *     LAMBDA_UM      d      Laser wavelength (micrometers)
        *
        *  Returned:
        *     FCUL_ZTD       d      Zenith total delay in meters
        *     FCUL_ZHD       d      Zenith hydrostatic delay in meters
        *     FCUL_ZWD       d      Zenith non-hydrostatic delay in meters
        *
        *  Notes:
        *
        *  1) The surface pressure provided was converted from inches Hg.
        *     The water vapor pressure was calculated from the surface
        *     temperature (Celsius) and Relative Humidity (% R.H.) at the station.
        *
        *  Test case:
        *     given input: LATITUDE  = 30.67166667D0 degrees (McDonald Observatory)
        *                  ELLIP_HT  = 2010.344D0 meters
        *                  PRESSURE  = 798.4188D0 hPa (August 14, 2009)
        *                  WVP       = 14.322D0 hPa (August 14, 2009)
        *                  LAMBDA_UM = 0.532D0 micrometers (See Mendes et al.)
        *     expected output: FCUL_ZTD = 1.935225924846803114D0 m
        *                      FCUL_ZHD = 1.932992176591644462D0 m
        *                      FCUL_ZWD = 0.2233748255158703871D-02 m
        *
        *  References:
        *     Mendes, V.B. and E.C. Pavlis, 2004,
        *     "High-accuracy zenith delay prediction at optical wavelengths,"
        *     Geophysical Res. Lett., 31, L14602, doi:10.1029/2004GL020308, 2004
        *
        *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
        *     IERS Technical Note No. 36, BKG (2010)
        *
        *  Revisions:
        *  2000 March  31 V.B. Mendes   Original code
        *  2009 August 14 B.E. Stetzler Added header and copyright
        *  2009 August 14 B.E. Stetzler Use of DOUBLE PRECISION
        *  2009 August 14 B.E. Stetzler Provided test case
        *  2009 August 14 B.E. Stetzler Capitalized all variables for FORTRAN 77
        *                               compatibility and provided more comments
        *-----------------------------------------------------------------------
        */
    void GIERS::FCULZD_HPA(double latitude, double ellip_ht, double pressure, double wvp, double lambda_um, double &fcul_ztd, double &fcul_zhd, double &fcul_zwd)
    {
        static double C = 2.99792458E8;
        static double PI = 3.1415926535897932384626433;
        /*  CO2 content in ppm */
        double XC = 375.0;
        /*         constant values to be used in Equation (20)
        *         k1 and k3 are k1* and k3*
         */
        static double K0 = 238.0185;
        static double K1 = 19990.975;
        static double K2 = 57.362;
        static double K3 = 579.55174;
        
        /*         constant values to be used in Equation (32) */
        static double W0 = 295.235;
        static double W1 = 2.6422;
        static double W2 = -0.032380;
        static double W3 = 0.004028;
        
        /*  Wave number */
        double SIGMA = 1/lambda_um;
        
        /*     correction factor - Equation (24) */
        double  F = 1 - 0.00266*cos(2*PI/180*latitude) - 0.00028E-3*ellip_ht;
        
        /*     correction for CO2 content */
        double CORR = 1.0 + 0.534E-6*(XC-450);
            
        /*     dispersion equation for the hydrostatic component - Equation (20) */
        double  FH = 0.01*CORR*((K1*(K0+SIGMA*SIGMA))/( ( K0-SIGMA*SIGMA )*( K0-SIGMA*SIGMA ) )
                                + K3*(K2+SIGMA*SIGMA)/( (K2-SIGMA*SIGMA)*(K2-SIGMA*SIGMA) ) );
            
        /*     computation of the hydrostatic component - Equation (26)
        *     caution: pressure in hectoPascal units
        */
        fcul_zhd = 2.416579E-3*FH*pressure/F;
        
        /*     dispersion equation for the non-hydrostatic component - Equation (32) */
        double FNH = 0.003101*(W0+3.0*W1*SIGMA*SIGMA + 5.0*W2*pow(SIGMA,4.0)+7.0*W3*pow(SIGMA,6.0) );
        
        /*     computation of the non-hydrostatic component - Equation (38)
        *     caution: pressure in hectoPascal units
         */
        fcul_zwd = 1.E-4*(5.316E0*FNH-3.759*FH)*wvp/F;
        
        /*      compute the zenith total delay */
        fcul_ztd = fcul_zhd + fcul_zwd;
        
    } // end of function FCUL_ZD_HPA
    
    
    /*
     *  - - - - - - - - -
     *   F C U L _ A
     *  - - - - - - - - -
     *
     *  This function is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This function computes the global total FCULa mapping function (Mendes et al. 2002).
     *  It is dependent on latitude, height, and surface temperature.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Class 1 model
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     LATITUDE       d      Latitude given in degrees (North Latitude)
     *     HEIGHT_M       d      Height given in meters (mean sea level)
     *     T_K            d      Surface temperature given in Kelvin
     *     ELEV_DEG       d      Elevation angle given in degrees (See references)
     
     *  Returned:
     *     FCUL_A         d      Mapping function to scale total delay (Note 1)
     *
     *  Notes:
     *
     *  1) These coefficients are based on a LS adjustment of 87766 (cleaned)
     *     set of traces, based on Ciddor routines to compute refractivity,
     *     according to IUGG recommendations (1999).
     *
     *  Test case:
     *     given input: LATITUDE = 30.67166667D0 degrees (McDonald Observatory)
     *                  HEIGHT_M = 2075D0 meters (mean sea level)
     *                  T_K      = 300.15D0 Kelvin (August 12, 2009)
     *                  ELEV_DEG = 15D0 degrees (See Mendes et al.)
     *     expected output: FCUL_A = 3.800243667312344087D0
     *
     *  References:
     *     Mendes, V.B., G. Prates, E.C. Pavlis, D.E. Pavlis,
     *     and R.B. Langley (2002). "Improved mapping functions for
     *     atmospheric refraction correction in SLR", Geophysical
     *     Res. Lett., 29(10), 1414, doi:10.1029/2001GL014394, 2002
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  2000 March  31 V.B. Mendes   Original code
     *  2009 August 12 B.E. Stetzler Added header and copyright
     *  2009 August 13 B.E. Stetzler Use of DOUBLE PRECISION
     *  2009 August 13 B.E. Stetzler Provided test case
     *  2009 August 13 B.E. Stetzler Capitalized all variables for FORTRAN 77
     *                               compatibility and provided more comments
     *-----------------------------------------------------------------------
     */
    double GIERS::FCUL_A(double LATITUDE, double HEIGHT_M, double T_K, double ELEV_DEG)
    {
        static double PI = 3.1415926535897932384626433;
        /* Convert elevation angle to radians */
        double EPSILON = ELEV_DEG * (PI/180);
        double SINE    = sin(EPSILON);
        /* Convert temperature to Celsius */
        double T_C     = T_K - 273.15 ;
        double COSPHI  = cos (LATITUDE*(PI/180));
        
        /* Define coefficients used in the model */
        
        double A10 =  0.121008E-02;
        double A11 =  0.17295E-05;
        double A12 =  0.3191E-04;
        double A13 = -0.18478E-07;
        
        double A20 =  0.304965E-02;
        double A21 =  0.2346E-05;
        double A22 = -0.1035E-03;
        double A23 = -0.1856E-07;
        
        
        double A30 =  0.68777E-01;
        double A31 =  0.1972E-04;
        double A32 = -0.3458E-02;
        double A33 =  0.1060E-06;
        
        /*     a, b, and c in Marini continued fraction (Eq. 5)*/
        double A1 = A10+A11*T_C+A12*COSPHI+A13*HEIGHT_M;
        double A2 = A20+A21*T_C+A22*COSPHI+A23*HEIGHT_M;
        double A3 = A30+A31*T_C+A32*COSPHI+A33*HEIGHT_M;
        
        /*     numerator in continued fraction */
        double MAP_ZEN   = (1.0 + A1/(1.0 + A2/(1.0+A3)));
        
        double mapping = MAP_ZEN/(SINE+A1/(SINE+A2/(SINE+A3)));
        
        return mapping;
    } // end of function FCUL_A
    
    
    /*
     
     *+
     *  - - - - - - - - - - -
     *   S T 1 I D I U
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine gives the out-of-phase corrections induced by
     *  mantle anelasticity in the diurnal band.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Class 1
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
     *     XSUN          d(3)   Geocentric position of the Sun (Note 2)
     *     XMON          d(3)   Geocentric position of the Moon (Note 2)
     *     FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)
     *     FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3)
     *
     *  Returned:
     *     XCORSTA       d(3)   Out of phase station corrections for diurnal band
     *
     *  Notes:
     *
     *  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
     *     expressed in meters.
     *
     *  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
     *     coordinates are expressed in meters.
     *
     *  3) The expressions are computed in the main program.  TGP is the tide
     *     generated potential.  The units are inverse meters.
     *
     *  Test case:
     *     given input: XSTA(1) = 4075578.385D0 meters
     *                  XSTA(2) =  931852.890D0 meters
     *                  XSTA(3) = 4801570.154D0 meters
     *                  XSUN(1) = 137859926952.015D0 meters
     *                  XSUN(2) = 54228127881.4350D0 meters
     *                  XSUN(3) = 23509422341.6960D0 meters
     *                  XMON(1) = -179996231.920342D0 meters
     *                  XMON(2) = -312468450.131567D0 meters
     *                  XMON(3) = -169288918.592160D0 meters
     *                  FAC2SUN =  0.163271964478954D0 1/meters
     *                  FAC2MON =  0.321989090026845D0 1/meters
     *
     *     expected output:  XCORSTA(1) = -0.2836337012840008001D-03 meters
     *                       XCORSTA(2) =  0.1125342324347507444D-03 meters
     *                       XCORSTA(3) = -0.2471186224343683169D-03 meters
     *
     *  References:
     *
     *     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
     *     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  1996 March    23 V. Dehant      Original code
     *  2009 July     30 B.E. Stetzler  Initial standardization of code 
     *  2009 July     31 B.E. Stetzler  Provided a test case
     *-----------------------------------------------------------------------
     
     */
    
    void GIERS::ST1IDIU( double* XSTA, double* XSUN, double* XMON, double FAC2SUN,double FAC2MON, double* XCORSTA)
    {
        double NORM8, RSTA, SINPHI, COSPHI, COS2PHI, SINLA,COSLA,RMON;
        double RSUN, DRSUN, DRMON, DNSUN, DNMON, DESUN, DEMON, DR, DN, DE;
        double DHI = -0.0025, DLI = -0.0007;
        
        /* Compute the normalized position vector of the IGS station */
        RSTA = sqrt( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] + XSTA[2]*XSTA[2] );
        SINPHI = XSTA[2] / RSTA;
        COSPHI = sqrt( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] ) / RSTA;
        COS2PHI = COSPHI*COSPHI - SINPHI*SINPHI;
        SINLA = XSTA[1]/COSPHI/RSTA;
        COSLA = XSTA[0]/COSPHI/RSTA;
        
        /* Compute the normalized position vector of the Moon. */
        RMON = sqrt( XMON[0]*XMON[0] + XMON[1]*XMON[1] + XMON[2]*XMON[2] ) ;
        /* Compute the normalized position vector of the Sun */
        RSUN = sqrt( XSUN[0]*XSUN[0] + XSUN[1]*XSUN[1] + XSUN[2]*XSUN[2] ) ;
        
        DRSUN=-3E0*DHI*SINPHI*COSPHI*FAC2SUN*XSUN[2]*(XSUN[0]*SINLA-XSUN[1]*COSLA)/pow(RSUN, 2.0);
        
        DRMON=-3E0*DHI*SINPHI*COSPHI*FAC2MON*XMON[2]*(XMON[0]*SINLA-XMON[1]*COSLA)/pow(RMON,2.0);
        
        DNSUN=-3E0*DLI*COS2PHI*FAC2SUN*XSUN[2]*(XSUN[0]*SINLA-XSUN[1]*COSLA)/pow(RSUN,2.0);
        
        DNMON=-3E0*DLI*COS2PHI*FAC2MON*XMON[2]*(XMON[0]*SINLA-XMON[1]*COSLA)/pow(RMON,2.0);
        
        DESUN=-3E0*DLI*SINPHI*FAC2SUN*XSUN[2]*(XSUN[0]*COSLA+XSUN[1]*SINLA)/pow(RSUN,2.0);
        
        DEMON=-3E0*DLI*SINPHI*FAC2MON*XMON[2]*(XMON[0]*COSLA+XMON[1]*SINLA)/pow(RMON,2.0);
        
        DR = DRSUN+DRMON ;
        DN = DNSUN+DNMON ;
        DE = DESUN+DEMON ;
        
        /*  Compute the corrections for the station. */
        XCORSTA[0]=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
        XCORSTA[1]=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
        XCORSTA[2]=DR*SINPHI+DN*COSPHI;
        
        
    } // end of function ST1IDIU
    
    
    /*
    *+
    *  - - - - - - - - - - -
    *   S T 1 L 1
    *  - - - - - - - - - - -
    *
    *  This routine is part of the International Earth Rotation and
    *  Reference Systems Service (IERS) Conventions software collection.
    *
    *  This subroutine gives the corrections induced by the latitude
    *  dependence given by L^1 in Mathews et al. 1991 (See References).
    *
    *  In general, Class 1, 2, and 3 models represent physical effects that
    *  act on geodetic parameters while canonical models provide lower-level
        *  representations or basic computations that are used by Class 1, 2, or
        *  3 models.
        *
        *  Status: Class 1
        *
        *     Class 1 models are those recommended to be used a priori in the
        *     reduction of raw space geodetic data in order to determine
        *     geodetic parameter estimates.
        *     Class 2 models are those that eliminate an observational
        *     singularity and are purely conventional in nature.
        *     Class 3 models are those that are not required as either Class
        *     1 or 2.
        *     Canonical models are accepted as is and cannot be classified as a
        *     Class 1, 2, or 3 model.
        *
        *  Given:
        *     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
        *     XSUN          d(3)   Geocentric position of the Sun (Note 2)
        *     XMON          d(3)   Geocentric position of the Moon (Note 2)
        *     FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)
            *     FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3)
                *
                *  Returned:
                *     XCORSTA       d(3)   Out of phase station corrections for
                    *                          semi-diurnal band
                    *
                    *  Notes:
                    *
                    *  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
                    *     expressed in meters.
                    *
                    *  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
                    *     coordinates are expressed in meters.
                    *
                    *  3) The expressions are computed in the main program. TGP is the tide
                    *     generated potential.  The units are inverse meters.
                    *
                    *  Test case:
                    *     given input: XSTA(1) = 4075578.385D0 meters
                    *                  XSTA(2) =  931852.890D0 meters
                    *                  XSTA(3) = 4801570.154D0 meters
                    *                  XSUN(1) = 137859926952.015D0 meters
                    *                  XSUN(2) = 54228127881.4350D0 meters
                    *                  XSUN(3) = 23509422341.6960D0 meters
                    *                  XMON(1) = -179996231.920342D0 meters
                    *                  XMON(2) = -312468450.131567D0 meters
                    *                  XMON(3) = -169288918.592160D0 meters
                    *                  FAC2SUN =  0.163271964478954D0 1/meters
                    *                  FAC2MON =  0.321989090026845D0 1/meters
                    *
                    *     expected output:  XCORSTA(1) = 0.2367189532359759044D-03 meters
                    *                       XCORSTA(2) = 0.5181609907284959182D-03 meters
                    *                       XCORSTA(3) = -0.3014881422940427977D-03 meters
                    *
                    *  References:
                    *
                    *     Mathews, P. M., Buffett, B. A., Herring, T. A., Shapiro, I. I.,
                    *     1991b, Forced nutations of the Earth: Influence of inner core
                    *     Dynamics 2. Numerical results and comparisons, J. Geophys. Res.,
                    *     96, 8243-8257
                    *
                    *     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
                    *     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
                    *
                    *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
                    *     IERS Technical Note No. 36, BKG (2010)
                    *
                    *  Revisions:
                    *  1996 March    23 V. Dehant      Original code
                    *  2009 July     31 B.E. Stetzler  Initial standardization of code 
                    *  2009 July     31 B.E. Stetzler  Provided a test case and Mathews
                    *                                  reference
                    *-----------------------------------------------------------------------
                    */
    void GIERS::ST1L1(double* XSTA, double* XSUN, double* XMON, double FAC2SUN,double FAC2MON, double* XCORSTA)
    {
        double RSTA, SINPHI, COSPHI, COSTWOLA, SINLA, COSLA, RMON, RSUN, DRSUN, DRMON, DNSUN, DNMON, DESUN, DEMON, DR, DN, DE;
        double  DHI, DLI, SINTWOLA, L1;
        double L1D = 0.0012, L1SD = 0.0024;
        
        /* Compute the normalized position vector of the IGS station */
        RSTA = sqrt( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] + XSTA[2]*XSTA[2] );
        SINPHI = XSTA[2] / RSTA;
        COSPHI = sqrt( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] ) / RSTA;
        SINLA = XSTA[1]/COSPHI/RSTA;
        COSLA = XSTA[0]/COSPHI/RSTA;
        
        /* Compute the normalized position vector of the Moon. */
        RMON = sqrt( XMON[0]*XMON[0] + XMON[1]*XMON[1] + XMON[2]*XMON[2] ) ;
        /* Compute the normalized position vector of the Sun */
        RSUN = sqrt( XSUN[0]*XSUN[0] + XSUN[1]*XSUN[1] + XSUN[2]*XSUN[2] ) ;
        
        /* Compute the station corrections for the diurnal band */
            
        L1=L1D;
        DNSUN=-L1*SINPHI*SINPHI*FAC2SUN*XSUN[2]*(XSUN[0]*COSLA+XSUN[1]*SINLA)/pow(RSUN,2.0);
        
        DNMON=-L1*SINPHI*SINPHI*FAC2MON*XMON[2]*(XMON[0]*COSLA+XMON[1]*SINLA)/pow(RMON,2.0);
        
        DESUN=L1*SINPHI*(COSPHI*COSPHI-SINPHI*SINPHI)*FAC2SUN*XSUN[2]*(XSUN[0]*SINLA-XSUN[1]*COSLA)/pow(RSUN,2.0);
        
        DEMON=L1*SINPHI*(COSPHI*COSPHI-SINPHI*SINPHI)*FAC2MON*XMON[2]*(XMON[0]*SINLA-XMON[1]*COSLA)/pow(RMON,2.0);
            
        DE = 3.0*(DESUN+DEMON);
        DN = 3.0*(DNSUN+DNMON);
            
        XCORSTA[0] = -DE*SINLA-DN*SINPHI*COSLA;
        XCORSTA[1] = DE*COSLA-DN*SINPHI*SINLA;
        XCORSTA[2] = DN*COSPHI;
        
       /* Compute the station corrections for the semi-diurnal band. */
            
        L1=L1SD;
        COSTWOLA=COSLA*COSLA-SINLA*SINLA;
        SINTWOLA=2.0*COSLA*SINLA;
            
        DNSUN=-L1/2.0*SINPHI*COSPHI*FAC2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*COSTWOLA+2.0*XSUN[0]*XSUN[1]*SINTWOLA)/pow(RSUN,2.0);
            
        DNMON=-L1/2.0*SINPHI*COSPHI*FAC2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*COSTWOLA+2.0*XMON[0]*XMON[1]*SINTWOLA)/pow(RMON,2.0);
            
        DESUN=-L1/2.0*SINPHI*SINPHI*COSPHI*FAC2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*SINTWOLA-2.0*XSUN[0]*XSUN[1]*COSTWOLA)/pow(RSUN,2.0);
            
        DEMON=-L1/2.0*SINPHI*SINPHI*COSPHI*FAC2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*SINTWOLA-2.0*XMON[0]*XMON[1]*COSTWOLA)/pow(RMON,2.0);
            
        DE = 3.0*(DESUN+DEMON);
        DN = 3.0*(DNSUN+DNMON);
            
        XCORSTA[0]=XCORSTA[0]-DE*SINLA-DN*SINPHI*COSLA;
        XCORSTA[1]=XCORSTA[1]+DE*COSLA-DN*SINPHI*SINLA;
        XCORSTA[2]=XCORSTA[2]+DN*COSPHI;
        
        
    }  // end of function ST1L1
    
    
    /*
     *  - - - - - - - - - - -
     *   S T 1 I S E M
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine gives the out-of-phase corrections induced by
     *  mantle anelasticity in the semi-diurnal band.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Class 1
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
     *     XSUN          d(3)   Geocentric position of the Sun (Note 2)
     *     XMON          d(3)   Geocentric position of the Moon (Note 2)
     *     FAC2SUN       d      Degree 2 TGP factor for the Sun (Note 3)
     *     FAC2MON       d      Degree 2 TGP factor for the Moon (Note 3)
     *
     *  Returned:
     *     XCORSTA       d(3)   Out of phase station corrections for
     *                          semi-diurnal band
     *
     *  Notes:
     *
     *  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
     *     expressed in meters.
     *
     *  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
     *     coordinates are expressed in meters.
     *
     *  3) The expressions are computed in the main program.  TGP is the tide
     *     generated potential.  The units are inverse meters.
     *
     *  Test case:
     *     given input: XSTA(1) = 4075578.385D0 meters
     *                  XSTA(2) =  931852.890D0 meters
     *                  XSTA(3) = 4801570.154D0 meters
     *                  XSUN(1) = 137859926952.015D0 meters
     *                  XSUN(2) = 54228127881.4350D0 meters
     *                  XSUN(3) = 23509422341.6960D0 meters
     *                  XMON(1) = -179996231.920342D0 meters
     *                  XMON(2) = -312468450.131567D0 meters
     *                  XMON(3) = -169288918.592160D0 meters
     *                  FAC2SUN =  0.163271964478954D0 1/meters
     *                  FAC2MON =  0.321989090026845D0 1/meters
     *
     *     expected output:  XCORSTA(1) = -0.2801334805106874015D-03 meters
     *                       XCORSTA(2) =  0.2939522229284325029D-04 meters
     *                       XCORSTA(3) = -0.6051677912316721561D-04 meters
     *
     *  References:
     *
     *     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
     *     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  1996 March    23 V. Dehant      Original code
     *  2009 July     31 B.E. Stetzler  Initial standardization of code 
     *  2009 July     31 B.E. Stetzler  Provided a test case
     *-----------------------------------------------------------------------
     */
    void GIERS::ST1ISEM ( double* XSTA, double* XSUN, double* XMON, double FAC2SUN,double FAC2MON, double* XCORSTA)
    {
        double RSTA, SINPHI, COSPHI, COSTWOLA, SINLA, COSLA, RMON, RSUN, DRSUN, DRMON, DNSUN, DNMON, DESUN, DEMON, DR, DN, DE,SINTWOLA;
        double DHI = -0.0022, DLI = -0.0007;
        
        /* Compute the normalized position vector of the IGS station */
        RSTA = sqrt( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] + XSTA[2]*XSTA[2] );
        SINPHI = XSTA[2] / RSTA;
        COSPHI = sqrt( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] ) / RSTA;
        SINLA = XSTA[1]/COSPHI/RSTA;
        COSLA = XSTA[0]/COSPHI/RSTA;
        COSTWOLA = COSLA*COSLA - SINLA*SINLA;
        SINTWOLA = 2.0*SINLA*COSLA;
        
        /* Compute the normalized position vector of the Moon. */
        RMON = sqrt( XMON[0]*XMON[0] + XMON[1]*XMON[1] + XMON[2]*XMON[2] ) ;
        /* Compute the normalized position vector of the Sun */
        RSUN = sqrt( XSUN[0]*XSUN[0] + XSUN[1]*XSUN[1] + XSUN[2]*XSUN[2] ) ;

        DRSUN=-3E0/4E0*DHI*COSPHI*COSPHI*FAC2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*SINTWOLA-2E0*XSUN[0]*XSUN[1]*COSTWOLA)/pow(RSUN, 2.0);
        
        DRMON=-3E0/4E0*DHI*COSPHI*COSPHI*FAC2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*SINTWOLA-2E0*XMON[0]*XMON[1]*COSTWOLA)/pow(RMON,2.0);
        
        DNSUN=3E0/2E0*DLI*SINPHI*COSPHI*FAC2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*SINTWOLA-2E0*XSUN[0]*XSUN[1]*COSTWOLA)/pow(RSUN,2.0);
        
        DNMON=3E0/2E0*DLI*SINPHI*COSPHI*FAC2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*SINTWOLA-2E0*XMON[0]*XMON[1]*COSTWOLA)/pow(RMON,2.0);
        
        DESUN=-3E0/2E0*DLI*COSPHI*FAC2SUN*((XSUN[0]*XSUN[0]-XSUN[1]*XSUN[1])*COSTWOLA+2E0*XSUN[0]*XSUN[1]*SINTWOLA)/pow(RSUN,2.0);
        
        DEMON=-3E0/2E0*DLI*COSPHI*FAC2MON*((XMON[0]*XMON[0]-XMON[1]*XMON[1])*COSTWOLA+2E0*XMON[0]*XMON[1]*SINTWOLA)/pow(RMON,2.0);
        
        DR=DRSUN+DRMON ;
        DN=DNSUN+DNMON ;
        DE=DESUN+DEMON ;
        
        XCORSTA[0]=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA ;
        XCORSTA[1]=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA ;
        XCORSTA[2]=DR*SINPHI+DN*COSPHI ;
        
    } // end of function ST1ISEM
    
    /*
     
     *  - - - - - - - - - - -
     *   S T E P 2 D I U
     *  - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine gives the in-phase and out-of-phase corrections
     *  induced by mantle anelasticity in the diurnal band.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Class 1
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
     *     FHR           d      Fractional hours in the day (Note 2)
     *     T             d      Centuries since J2000
     *
     *  Returned:
     *     XCORSTA       d(3)   In phase and out of phase station corrections
     *                          for diurnal band (Note 4)
     *
     *  Notes:
     *
     *  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
     *     expressed in meters.
     *
     *  2) The fractional hours in the day is computed as the hour + minutes/60.0
     *     + sec/3600.0.  The unit is expressed in Universal Time (UT).
     *
     *  4) All coordinates are expressed in meters.
     *
     *  Test case:
     *     given input: XSTA(1) = 4075578.385D0 meters
     *                  XSTA(2) =  931852.890D0 meters
     *                  XSTA(3) = 4801570.154D0 meters
     *                  FHR     = 0.00D0 hours
     *                  T       = 0.1059411362080767D0 Julian centuries
     *
     *     expected output:  XCORSTA(1) = 0.4193085327321284701D-02 meters
     *                       XCORSTA(2) = 0.1456681241014607395D-02 meters
     *                       XCORSTA(3) = 0.5123366597450316508D-02 meters
     *
     *  References:
     *
     *     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
     *     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *  Revisions:
     *  1996 March    23 V. Dehant      Original code
     *  2009 July     31 B.E. Stetzler  Initial standardization of code
     *  2009 August   06 B.E. Stetzler  Provided a test case
     *  2009 August   06 B.E. Stetzler  Capitalized all variables for 
     *                                  Fortran 77 compatibility
     *  2010 October  20 B.E. Stetzler  Input T corrected to be number of
     *                                  centuries since J2000
     *-----------------------------------------------------------------------
     */
    void GIERS::STEP2DIU( double* XSTA, double FHR, double T, double* XCORSTA)
    {
        
        int i, j;
        double DEG2RAD, S, TAU, PR, H,P, ZNS, PS, RSTA, SINPHI, COSPHI, COSLA, SINLA,ZLA, THETAF, DR, DN, DE;
        double D2PI = 6.283185307179586476925287;
        
        double DATDI[31][9] = {
            {-3E0, 0E0, 2E0, 0E0, 0E0,-0.01E0, 0E0, 0E0, 0E0},
            {-3E0, 2E0, 0E0, 0E0, 0E0,-0.01E0, 0E0, 0E0, 0E0},
            {-2E0, 0E0, 1E0,-1E0, 0E0,-0.02E0, 0E0, 0E0, 0E0},
            {-2E0, 0E0, 1E0, 0E0, 0E0,-0.08E0, 0E0,-0.01E0, 0.01E0},
            {-2E0, 2E0,-1E0, 0E0, 0E0,-0.02E0, 0E0, 0E0, 0E0},
            {-1E0, 0E0, 0E0,-1E0, 0E0,-0.10E0, 0E0, 0E0, 0E0},
            {-1E0, 0E0, 0E0, 0E0, 0E0,-0.51E0, 0E0,-0.02E0, 0.03E0},
            {-1E0, 2E0, 0E0, 0E0, 0E0, 0.01E0, 0E0, 0E0, 0E0},
            {0E0,-2E0, 1E0, 0E0, 0E0, 0.01E0, 0E0, 0E0, 0E0},
            {0E0, 0E0,-1E0, 0E0, 0E0, 0.02E0, 0E0, 0E0, 0E0},
            {0E0, 0E0, 1E0, 0E0, 0E0, 0.06E0, 0E0, 0E0, 0E0},
            {0E0, 0E0, 1E0, 1E0, 0E0, 0.01E0, 0E0, 0E0, 0E0},
            {0E0, 2E0,-1E0, 0E0, 0E0, 0.01E0, 0E0, 0E0, 0E0},
            {1E0,-3E0, 0E0, 0E0, 1E0,-0.06E0, 0E0, 0E0, 0E0},
            {1E0,-2E0, 0E0,-1E0, 0E0, 0.01E0, 0E0, 0E0, 0E0},
            {1E0,-2E0, 0E0, 0E0, 0E0,-1.23E0,-0.07E0, 0.06E0, 0.01E0},
            {1E0,-1E0, 0E0, 0E0,-1E0, 0.02E0, 0E0, 0E0, 0E0},
            {1E0,-1E0, 0E0, 0E0, 1E0, 0.04E0, 0E0, 0E0, 0E0},
            {1E0, 0E0, 0E0,-1E0, 0E0,-0.22E0, 0.01E0, 0.01E0, 0E0},
            {1E0, 0E0, 0E0, 0E0, 0E0,12.00E0,-0.80E0,-0.67E0,-0.03E0},
            {1E0, 0E0, 0E0, 1E0, 0E0, 1.73E0,-0.12E0,-0.10E0, 0E0},
            {1E0, 0E0, 0E0, 2E0, 0E0,-0.04E0, 0E0, 0E0, 0E0},
            {1E0, 1E0, 0E0, 0E0,-1E0,-0.50E0,-0.01E0, 0.03E0, 0E0},
            {1E0, 1E0, 0E0, 0E0, 1E0, 0.01E0, 0E0, 0E0, 0E0},
            {0E0, 1E0, 0E0, 1E0,-1E0,-0.01E0, 0E0, 0E0, 0E0},
            {1E0, 2E0,-2E0, 0E0, 0E0,-0.01E0, 0E0, 0E0, 0E0},
            {1E0, 2E0, 0E0, 0E0, 0E0,-0.11E0, 0.01E0, 0.01E0, 0E0},
            {2E0,-2E0, 1E0, 0E0, 0E0,-0.01E0, 0E0, 0E0, 0E0},
            {2E0, 0E0,-1E0, 0E0, 0E0,-0.02E0, 0E0, 0E0, 0E0},
            {3E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E0, 0E0},
            {3E0, 0E0, 0E0, 1E0, 0E0, 0E0, 0E0, 0E0, 0E0}
                                };
        
        DEG2RAD = D2PI / 360.0;
        
        /*  Compute the phase angles in degrees.*/
        S = 218.31664563 + (481267.88194 + (-0.0014663889 + (0.00000185139)*T)*T)*T;
        TAU = FHR*15 + 280.4606184 + (36000.7700536 + (0.00038793 + (-0.0000000258)*T)*T)*T + (-S);
        PR = (1.396971278 + (0.000308889 + (0.000000021 + (0.000000007)*T)*T)*T)*T;
        
        S = S + PR;
        H = 280.46645 + (36000.7697489 + (0.00030322222 + (0.000000020 + (-0.00000000654)*T)*T)*T)*T;
        P = 83.35324312 + (4069.01363525+ (-0.01032172222+ (-0.0000124991 + (0.00000005263)*T)*T)*T)*T;
        ZNS = 234.95544499 + (1934.13626197 + (-0.00207561111 + (-0.00000213944+ (0.00000001650)*T)*T)*T)*T;
        PS = 282.93734098+ (1.71945766667+ (0.00045688889+ (-0.00000001778+ (-0.00000000334)*T)*T)*T)*T;
        
        /* Reduce angles to between the range 0 and 360.*/
        S =  fmod(S,360.0);
        TAU = fmod(TAU,360.0);
        H =  fmod(H,360.0);
        P =  fmod(P,360.0);
        ZNS = fmod(ZNS,360.0);
        PS = fmod(PS,360.0);
        
        RSTA = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1]+XSTA[2]*XSTA[2]);
        SINPHI = XSTA[2]/RSTA;
        COSPHI = sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
        
        COSLA = XSTA[0]/COSPHI/RSTA;
        SINLA = XSTA[1]/COSPHI/RSTA;
        ZLA = atan2(XSTA[1],XSTA[0]);
        
        memset(XCORSTA,0,sizeof(double)*3);
        
        for(int J = 1 ; J<= 31; J++)
        {
            THETAF=(TAU+DATDI[J-1][0]*S+DATDI[J-1][1]*H+DATDI[J-1][2]*P+ DATDI[J-1][3]*ZNS+DATDI[J-1][4]*PS)*DEG2RAD;
            
            DR=DATDI[J-1][5]*2.0*SINPHI*COSPHI*sin(THETAF+ZLA)+DATDI[J-1][6]*2.0*SINPHI*COSPHI*cos(THETAF+ZLA);
            
            DN=DATDI[J-1][7]*(COSPHI*COSPHI-SINPHI*SINPHI)*sin(THETAF+ZLA)+DATDI[J-1][8]*(COSPHI*COSPHI-SINPHI*SINPHI)*cos(THETAF+ZLA);
            
            DE=DATDI[J-1][7]*SINPHI*cos(THETAF+ZLA)-DATDI[J-1][8]*SINPHI*sin(THETAF+ZLA);
            
            XCORSTA[0]=XCORSTA[0]+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
            XCORSTA[1]=XCORSTA[1]+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
            XCORSTA[2]=XCORSTA[2]+DR*SINPHI+DN*COSPHI;
            
        }
        
        for(int J = 0 ; J< 3; J++)
        {
            XCORSTA[J] = XCORSTA[J]/1000.0;
        }
        
    } // end of function STEP2DIU
    
    void GIERS::STEP2LON(double *XSTA, double T, double *XCORSTA)
    {
        int I =0, J =0;
        double DEG2RAD, S, PR, H, P, ZNS, PS, RSTA, SINPHI,COSPHI, COSLA, SINLA, THETAF, DR, DN, DE,DR_TOT, DN_TOT;
        double D2PI = 6.283185307179586476925287;
        double DATDI[5][9] = {
            {0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07},
            {0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05},
            {1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04},
            {2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07},
            {2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03}
                                };
        
        DEG2RAD = D2PI/360.0;
        
        /*  Compute the phase angles in degrees.*/
        S = 218.31664563+ (481267.88194+ (-0.0014663889+ (0.00000185139)*T)*T)*T ;
        PR = (1.396971278 + (0.000308889 + (0.000000021 + (0.000000007)*T)*T)*T)*T ;
        S = S + PR;
        H = 280.46645+ (36000.7697489+ (0.00030322222+ (0.000000020+ (-0.00000000654)*T)*T)*T)*T;
        P = 83.35324312+ (4069.01363525+ (-0.01032172222+ (-0.0000124991+ (0.00000005263)*T)*T)*T)*T;
        ZNS = 234.95544499+ (1934.13626197+ (-0.00207561111+ (-0.00000213944+ (0.00000001650)*T)*T)*T)*T;
        PS = 282.93734098+ (1.71945766667+ (0.00045688889+ (-0.00000001778+ (-0.00000000334)*T)*T)*T)*T;
        
        RSTA=sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1]+XSTA[2]*XSTA[2]);
        SINPHI=XSTA[2]/RSTA;
        COSPHI=sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
        
        COSLA=XSTA[0]/COSPHI/RSTA;
        SINLA=XSTA[1]/COSPHI/RSTA;
        
        /* Reduce angles to between the range 0 and 360.*/
        S =  fmod(S,360.0);
        H =  fmod(H,360.0);
        P =  fmod(P,360.0);
        ZNS = fmod(ZNS,360.0);
        PS = fmod(PS,360.0);
        
        DR_TOT = 0.0;
        DN_TOT = 0.0;
        
        memset(XCORSTA,0,sizeof(double)*3);
        
        for(int J =1 ; J <=5; J++)
        {
            THETAF=(DATDI[J-1][0]*S+DATDI[J-1][1]*H+DATDI[J-1][2]*P+ DATDI[J-1][3]*ZNS+DATDI[J-1][4]*PS)*DEG2RAD;
            
            DR=DATDI[J-1][5]*(3.0*SINPHI*SINPHI-1.0)/2.0*cos(THETAF)+DATDI[J-1][7]*(3.0*SINPHI*SINPHI-1.0)/2.0*sin(THETAF);
            
            DN=DATDI[J-1][6]*(COSPHI*SINPHI*2.0)*cos(THETAF)+DATDI[J-1][8]*(COSPHI*SINPHI*2.0)*sin(THETAF);
            
            DE = 0.0;
            DR_TOT = DR_TOT+DR;
            DN_TOT = DN_TOT+DN;
            
            XCORSTA[0]=XCORSTA[0]+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
            XCORSTA[1]=XCORSTA[1]+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
            XCORSTA[2]=XCORSTA[2]+DR*SINPHI+DN*COSPHI;
            
        }
        
        for(int J = 0; J < 3; J++)
        {
            XCORSTA[J] = XCORSTA[J]/1000.0;
        }
        
        
    } // end of function STEP2LON
    
    
    /*
     
     *  - - - - - - - - - - -
     *   C A L 2 J D
     *  - - - - - - - - - - -
     *
     *  Gregorian Calendar to Julian Date.
     *
     *  This routine is part of the International Astronomical Union's
     *  SOFA (Standards of Fundamental Astronomy) software collection.
     *
     *  Status:  support routine.
     *
     *  Given:
     *     IY,IM,ID    i     year, month, day in Gregorian calendar (Note 1)
     *
     *  Returned:
     *     DJM0        d     MJD zero-point: always 2400000.5
     *     DJM         d     Modified Julian Date for 0 hrs
     *     J           i     status:
     *                           0 = OK
     *                          -1 = bad year   (Note 3: JD not computed)
     *                          -2 = bad month  (JD not computed)
     *                          -3 = bad day    (JD computed)
     *
     *  Notes:
     *
     *  1) The algorithm used is valid from -4800 March 1, but this
     *     implementation rejects dates before -4799 January 1.
     *
     *  2) The Julian Date is returned in two pieces, in the usual SOFA
     *     manner, which is designed to preserve time resolution.  The
     *     Julian Date is available as a single number by adding DJM0 and
     *     DJM.
     *
     *  3) In early eras the conversion is from the "Proleptic Gregorian
     *     Calendar";  no account is taken of the date(s) of adoption of
     *     the Gregorian Calendar, nor is the AD/BC numbering convention
     *     observed.
     *
     *  Reference:
     *
     *     Explanatory Supplement to the Astronomical Almanac,
     *     P. Kenneth Seidelmann (ed), University Science Books (1992),
     *     Section 12.92 (p604).
     *
     *  This revision:  2001 September 16
     *
     *  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
     *
     *-----------------------------------------------------------------------
     
     */
    int GIERS::CAL2JD(int IY, int IM, int ID, double &DJM0, double &DJM)
    {
        int J, MY, IYPMY, IYMIN;
        /*  Earliest year allowed (4800BC) */
        IYMIN = -4799 ;
        /*  Month lengths in days*/
        int MTAB[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        /*  Preset status. */
        J = 0;
        /*  Validate year. */
        if( IY < IYMIN)
        {
            J = -1;
        }
        else
        {
           /*     Validate month. */
            if( IM >= 1 && IM <= 12)
            {
                /* Allow for leap year. */
                if ( IY%4 == 0 )
                {
                    MTAB[1] = 29;
                }
                else
                {
                    MTAB[1] = 28;
                }
                
                if ( IY%100 == 0 && IY%400 != 0 ) MTAB[1] = 28;
                /*        Validate day. */
                if ( ID < 1 || ID > MTAB[IM-1] ) J = -3;
                
                /*        Result. */
                MY = ( IM - 14 ) / 12;
                IYPMY = IY + MY;
                DJM0 = 2400000.5;
                DJM = double( ( 1461 * ( IYPMY + 4800 ) ) / 4 + (  367 * ( IM-2 - 12*MY ) ) / 12 - (    3 * ( ( IYPMY + 4900 ) / 100 ) ) / 4+ ID - 2432076 );
               
            }
            else
            {
                J = -2;
            }
            
        }
        
        return J;
        
    } // end of CAL2JD
    
    
    
    
    /*
     *  - - - - - - - - - - - - - - -
     *   D E H A N T T I D E I N E L
     *  - - - - - - - - - - - - - - -
     *
     *  This routine is part of the International Earth Rotation and
     *  Reference Systems Service (IERS) Conventions software collection.
     *
     *  This subroutine computes the tidal corrections of station displacements
     *  caused by lunar and solar gravitational attraction (see References).
     *  The computations are calculated by the following steps:
     *
     *  Step 1): General degree 2 and degree 3 corrections + CALL ST1IDIU
     *  + CALL ST1ISEM + CALL ST1L1.
     *
     *  Step 2): CALL STEP2DIU + CALL STEP2LON
     *
     *  It has been decided that the Step 3 non-correction for permanent tide
     *  would not be applied in order to avoid a jump in the reference frame.
     *  This Step 3 must be added in order to get the non-tidal station position
     *  and to conform with the IAG Resolution.
     *
     *  In general, Class 1, 2, and 3 models represent physical effects that
     *  act on geodetic parameters while canonical models provide lower-level
     *  representations or basic computations that are used by Class 1, 2, or
     *  3 models.
     *
     *  Status: Class 1
     *
     *     Class 1 models are those recommended to be used a priori in the
     *     reduction of raw space geodetic data in order to determine
     *     geodetic parameter estimates.
     *     Class 2 models are those that eliminate an observational
     *     singularity and are purely conventional in nature.
     *     Class 3 models are those that are not required as either Class
     *     1 or 2.
     *     Canonical models are accepted as is and cannot be classified as a
     *     Class 1, 2, or 3 model.
     *
     *  Given:
     *     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
     *     XSUN          d(3)   Geocentric position of the Sun (Note 2)
     *     XMON          d(3)   Geocentric position of the Moon (Note 2)
     *     YR            i      Year (Note 3)
     *     MONTH         i      Month (Note 3)
     *     DAY           i      Day of Month (Note 3)
     *     FHR           d      Hour in the day (Note 4)
     *
     *  Returned:
     *     DXTIDE        d(3)   Displacement vector (Note 5)
     *
     *  Notes:
     *
     *  1) The IGS station is in ITRF co-rotating frame.  All coordinates,
     *     X, Y, and Z, are expressed in meters.
     *
     *  2) The position is in Earth Centered Earth Fixed (ECEF) frame.  All
     *     coordinates are expressed in meters.
     *
     *  3) The values are expressed in Coordinated Universal Time (UTC).
     *
     *  4) The fractional hours in the day is computed as the hour + minutes/60.0
     *     + sec/3600.0.  The unit is expressed in Universal Time (UT).
     *
     *  5) The displacement vector is in the geocentric ITRF.  All components are
     *     expressed in meters.
     *
     *  Called:
     *     SPROD             Finds the scalar product and unit vector of two vectors
     *     ZERO_VEC8         Returns the zero vector
     *     ST1IDIU           Corrects for the out-of-phase part of Love numbers
     *                       for the diurnal band
     *     ST1ISEM           Same as above for the semi-diurnal band
     *     ST1L1             Corrects for the latitude dependence of Love numbers
     *     CAL2JD            Computes Julian Date from Gregorian calendar date
     *     DAT               Computes the difference TAI-UTC
     *     STEP2DIU          Computes in-phase and out-of-phase corrections in
     *                       the diurnal band
     *     STEP2LON          Same as above for the long period band
     *
     *  Test case:
     *     given input: XSTA(1) = 4075578.385D0 meters
     *                  XSTA(2) =  931852.890D0 meters
     *                  XSTA(3) = 4801570.154D0 meters
     *                  XSUN(1) = 137859926952.015D0 meters
     *                  XSUN(2) = 54228127881.4350D0 meters
     *                  XSUN(3) = 23509422341.6960D0 meters
     *                  XMON(1) = -179996231.920342D0 meters
     *                  XMON(2) = -312468450.131567D0 meters
     *                  XMON(3) = -169288918.592160D0 meters
     *                  YR      = 2009
     *                  MONTH   = 4
     *                  DAY     = 13
     *                  FHR     = 0.00D0 seconds
     *
     *     expected output:  DXTIDE(1) = 0.7700420357108125891D-01 meters
     *                       DXTIDE(2) = 0.6304056321824967613D-01 meters
     *                       DXTIDE(3) = 0.5516568152597246810D-01 meters
     *
     *  References:
     *
     *     Groten, E., 2000, Geodesists Handbook 2000, Part 4,
     *     http://www.gfy.ku.dk/~iag/HB2000/part4/groten.htm. See also
     *     ''Parameters of Common Relevance of Astronomy, Geodesy, and
     *     Geodynamics," J. Geod., 74, pp. 134-140
     *
     *     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
     *     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
     *
     *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
     *     IERS Technical Note No. 36, BKG (2010)
     *
     *     Pitjeva, E. and Standish, E. M., 2009, ''Proposals for the masses
     *     of the three largest asteroids, the Moon-Earth mass ratio and the
     *     Astronomical Unit," Celest. Mech. Dyn. Astr., 103, pp. 365-372
     *
     *     Ries, J. C., Eanes, R. J., Shum, C. K. and Watkins, M. M., 1992,
     *     ''Progress in the Determination of the Gravitational Coefficient
     *     of the Earth," Geophys. Res. Lett., 19(6), pp. 529-531
     *
     *  Revisions:
     *  1996 March    23 V. Dehant      Original code
     *                   P. M. Mathews
     *                   J. Gipson
     *  2000 May      17 V. Dehant      Last modifications
     *                   P. M. Mathews
     *  2006 February 06 J. Ray         Header comments modified to clarify
     *                                  input/output units and systems
     *  2006 February 06 J. Ray         Subroutine DUTC modified for leap
     *                                  second on 2006.0 and to correct
     *                                  do 5 i=1,87 from 84 to 87
     *  2006 August   31 G. Petit       Correct DUTC for dates after 2007
     *  2007 June     20 H. Manche      Modified DUTC to correct past mistake
     *                                  and corrected DE line in subroutine
     *                                  STEP2DIU
     *  2007 October  23 H. Manche      Replace subroutines DUTC and FJLDY with
     *                   G. Petit       SOFA subroutines iau_CAL2JD and iau_DAT
     *                                  and correct time arguments of subroutine
     *                                  STEP2DIU
     *  2009 February 19 G. Petit       Update routine iau_DAT for 2009.0 leap
     *                                  second
     *  2009 August   06 B.E. Stetzler  Initial standardization of code 
     *  2009 August   07 B.E. Stetzler  Updated MASS_RATIO_SUN, 
     *                                  MASS_RATIO_MOON and RE to CBEs and
     *                                  provided a test case
     *  2009 August  07  B.E. Stetzler  Capitalized all variables for Fortran
     *                                  77 compatibility
     *  2009 September 01 B.E. Stetzler Removed 'iau_' from redistributed SOFA
     *                                  subroutines
     *-----------------------------------------------------------------------
     */
    //void GIERS::DEHANTTIDEINEL(double *XSTA, double JD_UTC,double leadsec, double *XSUN, double *XMON, double *DXTIDE)
    void GIERS::DEHANTTIDEINEL(double *XSTA, int YR,int MONTH, int DAY, double FHR, double leadsec, double *XSUN, double *XMON, double *DXTIDE)
    {
        
        int I =0, J =0;
        double XCORSTA[3],H20,L20,H3,L3,H2,L2,SCS,RSTA,SCM,RSUN,RMON,SCSUN,SCMON,COSPHI,P2SUN,P2MON,P3SUN,P3MON;
        double X2SUN,X2MON,X3SUN,X3MON,MASS_RATIO_SUN,MASS_RATIO_MOON,RE,FAC2SUN,FAC2MON,FAC3SUN,FAC3MON,JJM0,JJM1,DTT,T,PI,SINPHI,COSLA,SINLA,DR,DN;
        PI = 3.1415926535897932384626433;
        
        /*----------------------------------------------------------------------
        * NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS
        *----------------------------------------------------------------------
        */
        H20 = 0.6078; L20 = 0.0847; H3 = 0.292; L3 = 0.015;
        /*----------------------------------------------------------------------
        * SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR
        *----------------------------------------------------------------------
        */
        RSTA = sqrt ( XSTA[0]*XSTA[0] + XSTA[1]*XSTA[1] + XSTA[2]*XSTA[2] );
        RSUN = sqrt ( XSUN[0]*XSUN[0] + XSUN[1]*XSUN[1] + XSUN[2]*XSUN[2] );
        RMON = sqrt ( XMON[0]*XMON[0] + XMON[1]*XMON[1] + XMON[2]*XMON[2] );
        SCS  =  XSTA[0]*XSUN[0] + XSTA[1]*XSUN[1] + XSTA[2]*XSUN[2] ;
        SCM =   XSTA[0]*XMON[0] + XSTA[1]*XMON[1] + XSTA[2]*XMON[2] ;
        
        SCSUN=SCS/RSTA/RSUN;
        SCMON=SCM/RSTA/RMON;
        /*----------------------------------------------------------------------
        * COMPUTATION OF NEW H2 AND L2
        *----------------------------------------------------------------------
        */
        COSPHI=sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
        H2=H20-0.0006*(1.0-3.0/2.0*COSPHI*COSPHI);
        L2=L20+0.0002*(1.0-3.0/2.0*COSPHI*COSPHI);
        
        /* P2 term */
        P2SUN=3.0*(H2/2.0-L2)*SCSUN*SCSUN-H2/2.0;
        P2MON=3.0*(H2/2.0-L2)*SCMON*SCMON-H2/2.0;
        
        /* P3 term */
        P3SUN=5.0/2.0*(H3-3.0*L3)*pow(SCSUN,3.0)+3.0/2.0*(L3-H3)*SCSUN;
        P3MON=5.0/2.0*(H3-3.0*L3)*pow(SCMON,3.0)+3.0/2.0*(L3-H3)*SCMON;
        
        /*----------------------------------------------------------------------
        * TERM IN DIRECTION OF SUN/MOON VECTOR
        *----------------------------------------------------------------------
        */
        X2SUN=3.0*L2*SCSUN;
        X2MON=3.0*L2*SCMON;
        X3SUN=3.0*L3/2.0*(5.0*SCSUN*SCSUN-1.0);
        X3MON=3.0*L3/2.0*(5.0*SCMON*SCMON-1.0);
        
        /*----------------------------------------------------------------------
        * FACTORS FOR SUN/MOON USING IAU CURRENT BEST ESTIMATES (SEE REFERENCES)
        *----------------------------------------------------------------------
        */
        MASS_RATIO_SUN=332946.0482;
        MASS_RATIO_MOON=0.0123000371;
        RE=6378136.6;
        FAC2SUN=MASS_RATIO_SUN*RE*pow((RE/RSUN),3.0);
        FAC2MON=MASS_RATIO_MOON*RE*pow((RE/RMON),3.0);
        FAC3SUN=FAC2SUN*(RE/RSUN);
        FAC3MON=FAC2MON*(RE/RMON);
        
        /* TOTAL DISPLACEMENT  */
        for(I = 0; I< 3; I++)
        {
            DXTIDE[I]= FAC2SUN*( X2SUN*XSUN[I]/RSUN + P2SUN*XSTA[I]/RSTA ) +
                       FAC2MON*( X2MON*XMON[I]/RMON + P2MON*XSTA[I]/RSTA ) +
                       FAC3SUN*( X3SUN*XSUN[I]/RSUN + P3SUN*XSTA[I]/RSTA ) +
                       FAC3MON*( X3MON*XMON[I]/RMON + P3MON*XSTA[I]/RSTA );
        }
        
        memset(XCORSTA,0,sizeof(double)*3);
        
        /*+---------------------------------------------------------------------
        * CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I
                                                                 * AND L_2^(0)I )
        *----------------------------------------------------------------------
        */
        
        /* FIRST, FOR THE DIURNAL BAND */
        GIERS::ST1IDIU(XSTA, XSUN, XMON, FAC2SUN, FAC2MON, XCORSTA);
        for( I = 0; I< 3; I++)
        {
            DXTIDE[I] += XCORSTA[I];
        }
        
        /* SECOND, FOR THE SEMI-DIURNAL BAND */
        GIERS::ST1ISEM(XSTA, XSUN, XMON, FAC2SUN, FAC2MON, XCORSTA);
        for( I = 0; I< 3; I++)
        {
            DXTIDE[I] += XCORSTA[I];
        }
        
        /*+---------------------------------------------------------------------
        * CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )
        *----------------------------------------------------------------------
        */
        GIERS::ST1L1(XSTA, XSUN, XMON, FAC2SUN, FAC2MON, XCORSTA);
        for( I = 0; I< 3; I++)
        {
            DXTIDE[I] += XCORSTA[I];
        }
        
        /* CONSIDER CORRECTIONS FOR STEP 2 */
        
        /*+---------------------------------------------------------------------
        * CORRECTIONS FOR THE DIURNAL BAND:
        *
        *  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES
        *
        *   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE
        *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        */
        int STATUS = GIERS::CAL2JD(YR, MONTH, DAY, JJM0, JJM1);
        T=((JJM0-2451545.0)+JJM1+FHR/24.0)/36525.0;
        T = T + ( leadsec + 32.184 )/(3600.0*24.0*36525.0);
        
//        double t1 = JD_UTC + 0.5;
//        FHR = (t1 - int(t1))*24 ;
//        T = (JD_UTC - 2451545.0)/36525.0;
//        T = T + ( leadsec + 32.184 )/(3600.0*24.0*36525);
        
        /*  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND
        *  CORRECTIONS, (in-phase and out-of-phase frequency dependence):
        */
        
        GIERS::STEP2DIU(XSTA, FHR, T, XCORSTA);
        for( I = 0; I< 3; I++)
        {
            DXTIDE[I] += XCORSTA[I];
        }
        
        
        /*  CORRECTIONS FOR THE LONG-PERIOD BAND,
        *  (in-phase and out-of-phase frequency dependence):
        */
        
        GIERS::STEP2LON(XSTA, T, XCORSTA);
        for( I = 0; I< 3; I++)
        {
            DXTIDE[I] += XCORSTA[I];
        }
        
        
        /* CONSIDER CORRECTIONS FOR STEP 3 */
         
         /*----------------------------------------------------------------------
             //UNCORRECT FOR THE PERMANENT TIDE
         
        SINPHI=XSTA[2]/RSTA;
        COSPHI=sqrt(XSTA[0]*XSTA[0]+XSTA[1]*XSTA[1])/RSTA;
        COSLA=XSTA[0]/COSPHI/RSTA;
        SINLA=XSTA[1]/COSPHI/RSTA;
        DR=-sqrt(5.0/4.0/PI)*H2*0.31460*(3.0/2.0*SINPHI*SINPHI-0.5);
        DN=-sqrt(5.0/4.0/PI)*L2*0.31460*3.0*COSPHI*SINPHI;
        DXTIDE[0]=DXTIDE[0]-DR*COSLA*COSPHI+DN*COSLA*SINPHI;
        DXTIDE[1]=DXTIDE[1]-DR*SINLA*COSPHI+DN*SINLA*SINPHI;
        DXTIDE[2]=DXTIDE[2]-DR*SINPHI      -DN*COSPHI;
         // *-----------------------------------------------------------------------
        */
        
        
        
    } // end of function DEHANTTIDEINEL
    
    
    
    
} // end of namespace gfc




