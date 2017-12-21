//
//  GIAU.cpp
//  GFC
//
//  Created by lizhen on 09/09/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#include "GIAU.hpp"
namespace gfc
{
    
    double GIAU::DJ00 = 2451545.0;
    double GIAU::DJC  = 36525.0;
    double GIAU::DAS2R =4.848136811095359935899141e-6;
    double GIAU::TURNAS = 1296000.0;
    double GIAU::D2PI = 6.283185307179586476925287;
    double GIAU::DPI  = 3.141592653589793238462643;
    
    
    /*
     **  - - - - - - - - -
     **   i a u O b l 0 6
     **  - - - - - - - - -
     **
     **  Mean obliquity of the ecliptic, IAU 2006 precession model.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     date1,date2  double   TT as a 2-part Julian Date (Note 1)
     **
     **  Returned (function value):
     **                  double   obliquity of the ecliptic (radians, Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The result is the angle between the ecliptic and mean equator of
     **     date date1+date2.
     **
     **  Reference:
     **
     **     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauObl06(double date1, double date2)
    {
        double t, eps0;
        
        /* Interval between fundamental date J2000.0 and given date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* Mean obliquity. */
        eps0 = (84381.406     +
                (-46.836769    +
                 ( -0.0001831   +
                  (  0.00200340  +
                   ( -0.000000576 +
                    ( -0.0000000434) * t) * t) * t) * t) * t) * DAS2R;
        
        return eps0;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - -
     **   i a u P f w 0 6
     **  - - - - - - - - -
     **
     **  Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     date1,date2  double   TT as a 2-part Julian Date (Note 1)
     **
     **  Returned:
     **     gamb         double   F-W angle gamma_bar (radians)
     **     phib         double   F-W angle phi_bar (radians)
     **     psib         double   F-W angle psi_bar (radians)
     **     epsa         double   F-W angle epsilon_A (radians)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) Naming the following points:
     **
     **           e = J2000.0 ecliptic pole,
     **           p = GCRS pole,
     **           E = mean ecliptic pole of date,
     **     and   P = mean pole of date,
     **
     **     the four Fukushima-Williams angles are as follows:
     **
     **        gamb = gamma_bar = epE
     **        phib = phi_bar = pE
     **        psib = psi_bar = pEP
     **        epsa = epsilon_A = EP
     **
     **  3) The matrix representing the combined effects of frame bias and
     **     precession is:
     **
     **        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)
     **
     **  4) The matrix representing the combined effects of frame bias,
     **     precession and nutation is simply:
     **
     **        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
     **
     **     where dP and dE are the nutation components with respect to the
     **     ecliptic of date.
     **
     **  Reference:
     **
     **     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
     **
     **  Called:
     **     iauObl06     mean obliquity, IAU 2006
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauPfw06(double date1, double date2, double *gamb, double *phib, double *psib, double *epsa)
    {
        double t;
        
        /* Interval between fundamental date J2000.0 and given date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* P03 bias+precession angles. */
        *gamb = (    -0.052928     +
                 (    10.556378     +
                  (     0.4932044    +
                   (    -0.00031238   +
                    (    -0.000002788  +
                     (     0.0000000260 )
                     * t) * t) * t) * t) * t) * DAS2R;
        *phib = ( 84381.412819     +
                 (   -46.811016     +
                  (     0.0511268    +
                   (     0.00053289   +
                    (    -0.000000440  +
                     (    -0.0000000176 )
                     * t) * t) * t) * t) * t) * DAS2R;
        *psib = (    -0.041775     +
                 (  5038.481484     +
                  (     1.5584175    +
                   (    -0.00018522   +
                    (    -0.000026452  +
                     (    -0.0000000148 )
                     * t) * t) * t) * t) * t) * DAS2R;
        *epsa =  iauObl06(date1, date2);
        
        return;
    } // end of function
    
    
    
    /*
     **  - - - - - - - - - -
     **   i a u N u t 0 0 a
     **  - - - - - - - - - -
     **
     **  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
     **  with free core nutation omitted).
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     date1,date2   double   TT as a 2-part Julian Date (Note 1)
     **
     **  Returned:
     **     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The nutation components in longitude and obliquity are in radians
     **     and with respect to the equinox and ecliptic of date.  The
     **     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
     **     value of 84381.448 arcsec.
     **
     **     Both the luni-solar and planetary nutations are included.  The
     **     latter are due to direct planetary nutations and the
     **     perturbations of the lunar and terrestrial orbits.
     **
     **  3) The function computes the MHB2000 nutation series with the
     **     associated corrections for planetary nutations.  It is an
     **     implementation of the nutation part of the IAU 2000A precession-
     **     nutation model, formally adopted by the IAU General Assembly in
     **     2000, namely MHB2000 (Mathews et al. 2002), but with the free
     **     core nutation (FCN - see Note 4) omitted.
     **
     **  4) The full MHB2000 model also contains contributions to the
     **     nutations in longitude and obliquity due to the free-excitation
     **     of the free-core-nutation during the period 1979-2000.  These FCN
     **     terms, which are time-dependent and unpredictable, are NOT
     **     included in the present function and, if required, must be
     **     independently computed.  With the FCN corrections included, the
     **     present function delivers a pole which is at current epochs
     **     accurate to a few hundred microarcseconds.  The omission of FCN
     **     introduces further errors of about that size.
     **
     **  5) The present function provides classical nutation.  The MHB2000
     **     algorithm, from which it is adapted, deals also with (i) the
     **     offsets between the GCRS and mean poles and (ii) the adjustments
     **     in longitude and obliquity due to the changed precession rates.
     **     These additional functions, namely frame bias and precession
     **     adjustments, are supported by the SOFA functions iauBi00  and
     **     iauPr00.
     **
     **  6) The MHB2000 algorithm also provides "total" nutations, comprising
     **     the arithmetic sum of the frame bias, precession adjustments,
     **     luni-solar nutation and planetary nutation.  These total
     **     nutations can be used in combination with an existing IAU 1976
     **     precession implementation, such as iauPmat76,  to deliver GCRS-
     **     to-true predictions of sub-mas accuracy at current dates.
     **     However, there are three shortcomings in the MHB2000 model that
     **     must be taken into account if more accurate or definitive results
     **     are required (see Wallace 2002):
     **
     **       (i) The MHB2000 total nutations are simply arithmetic sums,
     **           yet in reality the various components are successive Euler
     **           rotations.  This slight lack of rigor leads to cross terms
     **           that exceed 1 mas after a century.  The rigorous procedure
     **           is to form the GCRS-to-true rotation matrix by applying the
     **           bias, precession and nutation in that order.
     **
     **      (ii) Although the precession adjustments are stated to be with
     **           respect to Lieske et al. (1977), the MHB2000 model does
     **           not specify which set of Euler angles are to be used and
     **           how the adjustments are to be applied.  The most literal
     **           and straightforward procedure is to adopt the 4-rotation
     **           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
     **           to psi_A and DEPSPR to both omega_A and eps_A.
     **
     **     (iii) The MHB2000 model predates the determination by Chapront
     **           et al. (2002) of a 14.6 mas displacement between the
     **           J2000.0 mean equinox and the origin of the ICRS frame.  It
     **           should, however, be noted that neglecting this displacement
     **           when calculating star coordinates does not lead to a
     **           14.6 mas change in right ascension, only a small second-
     **           order distortion in the pattern of the precession-nutation
     **           effect.
     **
     **     For these reasons, the SOFA functions do not generate the "total
     **     nutations" directly, though they can of course easily be
     **     generated by calling iauBi00, iauPr00 and the present function
     **     and adding the results.
     **
     **  7) The MHB2000 model contains 41 instances where the same frequency
     **     appears multiple times, of which 38 are duplicates and three are
     **     triplicates.  To keep the present code close to the original MHB
     **     algorithm, this small inefficiency has not been corrected.
     **
     **  Called:
     **     iauFal03     mean anomaly of the Moon
     **     iauFaf03     mean argument of the latitude of the Moon
     **     iauFaom03    mean longitude of the Moon's ascending node
     **     iauFame03    mean longitude of Mercury
     **     iauFave03    mean longitude of Venus
     **     iauFae03     mean longitude of Earth
     **     iauFama03    mean longitude of Mars
     **     iauFaju03    mean longitude of Jupiter
     **     iauFasa03    mean longitude of Saturn
     **     iauFaur03    mean longitude of Uranus
     **     iauFapa03    general accumulated precession in longitude
     **
     **  References:
     **
     **     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
     **     Astron.Astrophys. 387, 700
     **
     **     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
     **     Astron.Astrophys. 58, 1-16
     **
     **     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
     **     107, B4.  The MHB_2000 code itself was obtained on 9th September
     **     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **     Wallace, P.T., "Software for Implementing the IAU 2000
     **     Resolutions", in IERS Workshop 5.1 (2002)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauNut00a(double date1, double date2, double *dpsi, double *deps)
    {
        int i;
        double t, el, elp, f, d, om, arg, dp, de, sarg, carg,
        al, af, ad, aom, alme, alve, alea, alma,
        alju, alsa, alur, alne, apa, dpsils, depsls,
        dpsipl, depspl;
        
        /* Units of 0.1 microarcsecond to radians */
        const double U2R = DAS2R / 1e7;
        
        /* ------------------------- */
        /* Luni-Solar nutation model */
        /* ------------------------- */
        
        /* The units for the sine and cosine coefficients are */
        /* 0.1 microarcsecond and the same per Julian century */
        
        static const struct {
            int nl,nlp,nf,nd,nom; /* coefficients of l,l',F,D,Om */
            double sp,spt,cp;     /* longitude sin, t*sin, cos coefficients */
            double ce,cet,se;     /* obliquity cos, t*cos, sin coefficients */
        } xls[] = {
            
            /* 1- 10 */
            { 0, 0, 0, 0, 1,
                -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
            { 0, 0, 2,-2, 2,
                -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0},
            { 0, 0, 2, 0, 2,-2276413.0,-234.0,2796.0,978459.0,-485.0, 1374.0},
            { 0, 0, 0, 0, 2,2074554.0, 207.0, -698.0,-897492.0,470.0, -291.0},
            { 0, 1, 0, 0, 0,1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0},
            { 0, 1, 2,-2, 2,-516821.0,1226.0, -524.0,224386.0,-677.0, -174.0},
            { 1, 0, 0, 0, 0, 711159.0,  73.0, -872.0,  -6750.0,  0.0,  358.0},
            { 0, 0, 2, 0, 1,-387298.0,-367.0,  380.0, 200728.0, 18.0,  318.0},
            { 1, 0, 2, 0, 2,-301461.0, -36.0,  816.0, 129025.0,-63.0,  367.0},
            { 0,-1, 2,-2, 2, 215829.0,-494.0,  111.0, -95929.0,299.0,  132.0},
            
            /* 11-20 */
            { 0, 0, 2,-2, 1, 128227.0, 137.0,  181.0, -68982.0, -9.0,   39.0},
            {-1, 0, 2, 0, 2, 123457.0,  11.0,   19.0, -53311.0, 32.0,   -4.0},
            {-1, 0, 0, 2, 0, 156994.0,  10.0, -168.0,  -1235.0,  0.0,   82.0},
            { 1, 0, 0, 0, 1,  63110.0,  63.0,   27.0, -33228.0,  0.0,   -9.0},
            {-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0,  31429.0,  0.0,  -75.0},
            {-1, 0, 2, 2, 2, -59641.0, -11.0,  149.0,  25543.0,-11.0,   66.0},
            { 1, 0, 2, 0, 1, -51613.0, -42.0,  129.0,  26366.0,  0.0,   78.0},
            {-2, 0, 2, 0, 1,  45893.0,  50.0,   31.0, -24236.0,-10.0,   20.0},
            { 0, 0, 0, 2, 0,  63384.0,  11.0, -150.0,  -1220.0,  0.0,   29.0},
            { 0, 0, 2, 2, 2, -38571.0,  -1.0,  158.0,  16452.0,-11.0,   68.0},
            
            /* 21-30 */
            { 0,-2, 2,-2, 2,  32481.0,   0.0,    0.0, -13870.0,  0.0,    0.0},
            {-2, 0, 0, 2, 0, -47722.0,   0.0,  -18.0,    477.0,  0.0,  -25.0},
            { 2, 0, 2, 0, 2, -31046.0,  -1.0,  131.0,  13238.0,-11.0,   59.0},
            { 1, 0, 2,-2, 2,  28593.0,   0.0,   -1.0, -12338.0, 10.0,   -3.0},
            {-1, 0, 2, 0, 1,  20441.0,  21.0,   10.0, -10758.0,  0.0,   -3.0},
            { 2, 0, 0, 0, 0,  29243.0,   0.0,  -74.0,   -609.0,  0.0,   13.0},
            { 0, 0, 2, 0, 0,  25887.0,   0.0,  -66.0,   -550.0,  0.0,   11.0},
            { 0, 1, 0, 0, 1, -14053.0, -25.0,   79.0,   8551.0, -2.0,  -45.0},
            {-1, 0, 0, 2, 1,  15164.0,  10.0,   11.0,  -8001.0,  0.0,   -1.0},
            { 0, 2, 2,-2, 2, -15794.0,  72.0,  -16.0,   6850.0,-42.0,   -5.0},
            
            /* 31-40 */
            { 0, 0,-2, 2, 0,  21783.0,   0.0,   13.0,   -167.0,  0.0,   13.0},
            { 1, 0, 0,-2, 1, -12873.0, -10.0,  -37.0,   6953.0,  0.0,  -14.0},
            { 0,-1, 0, 0, 1, -12654.0,  11.0,   63.0,   6415.0,  0.0,   26.0},
            {-1, 0, 2, 2, 1, -10204.0,   0.0,   25.0,   5222.0,  0.0,   15.0},
            { 0, 2, 0, 0, 0,  16707.0, -85.0,  -10.0,    168.0, -1.0,   10.0},
            { 1, 0, 2, 2, 2,  -7691.0,   0.0,   44.0,   3268.0,  0.0,   19.0},
            {-2, 0, 2, 0, 0, -11024.0,   0.0,  -14.0,    104.0,  0.0,    2.0},
            { 0, 1, 2, 0, 2,   7566.0, -21.0,  -11.0,  -3250.0,  0.0,   -5.0},
            { 0, 0, 2, 2, 1,  -6637.0, -11.0,   25.0,   3353.0,  0.0,   14.0},
            { 0,-1, 2, 0, 2,  -7141.0,  21.0,    8.0,   3070.0,  0.0,    4.0},
            
            /* 41-50 */
            { 0, 0, 0, 2, 1,  -6302.0, -11.0,    2.0,   3272.0,  0.0,    4.0},
            { 1, 0, 2,-2, 1,   5800.0,  10.0,    2.0,  -3045.0,  0.0,   -1.0},
            { 2, 0, 2,-2, 2,   6443.0,   0.0,   -7.0,  -2768.0,  0.0,   -4.0},
            {-2, 0, 0, 2, 1,  -5774.0, -11.0,  -15.0,   3041.0,  0.0,   -5.0},
            { 2, 0, 2, 0, 1,  -5350.0,   0.0,   21.0,   2695.0,  0.0,   12.0},
            { 0,-1, 2,-2, 1,  -4752.0, -11.0,   -3.0,   2719.0,  0.0,   -3.0},
            { 0, 0, 0,-2, 1,  -4940.0, -11.0,  -21.0,   2720.0,  0.0,   -9.0},
            {-1,-1, 0, 2, 0,   7350.0,   0.0,   -8.0,    -51.0,  0.0,    4.0},
            { 2, 0, 0,-2, 1,   4065.0,   0.0,    6.0,  -2206.0,  0.0,    1.0},
            { 1, 0, 0, 2, 0,   6579.0,   0.0,  -24.0,   -199.0,  0.0,    2.0},
            
            /* 51-60 */
            { 0, 1, 2,-2, 1,   3579.0,   0.0,    5.0,  -1900.0,  0.0,    1.0},
            { 1,-1, 0, 0, 0,   4725.0,   0.0,   -6.0,    -41.0,  0.0,    3.0},
            {-2, 0, 2, 0, 2,  -3075.0,   0.0,   -2.0,   1313.0,  0.0,   -1.0},
            { 3, 0, 2, 0, 2,  -2904.0,   0.0,   15.0,   1233.0,  0.0,    7.0},
            { 0,-1, 0, 2, 0,   4348.0,   0.0,  -10.0,    -81.0,  0.0,    2.0},
            { 1,-1, 2, 0, 2,  -2878.0,   0.0,    8.0,   1232.0,  0.0,    4.0},
            { 0, 0, 0, 1, 0,  -4230.0,   0.0,    5.0,    -20.0,  0.0,   -2.0},
            {-1,-1, 2, 2, 2,  -2819.0,   0.0,    7.0,   1207.0,  0.0,    3.0},
            {-1, 0, 2, 0, 0,  -4056.0,   0.0,    5.0,     40.0,  0.0,   -2.0},
            { 0,-1, 2, 2, 2,  -2647.0,   0.0,   11.0,   1129.0,  0.0,    5.0},
            
            /* 61-70 */
            {-2, 0, 0, 0, 1,  -2294.0,   0.0,  -10.0,   1266.0,  0.0,   -4.0},
            { 1, 1, 2, 0, 2,   2481.0,   0.0,   -7.0,  -1062.0,  0.0,   -3.0},
            { 2, 0, 0, 0, 1,   2179.0,   0.0,   -2.0,  -1129.0,  0.0,   -2.0},
            {-1, 1, 0, 1, 0,   3276.0,   0.0,    1.0,     -9.0,  0.0,    0.0},
            { 1, 1, 0, 0, 0,  -3389.0,   0.0,    5.0,     35.0,  0.0,   -2.0},
            { 1, 0, 2, 0, 0,   3339.0,   0.0,  -13.0,   -107.0,  0.0,    1.0},
            {-1, 0, 2,-2, 1,  -1987.0,   0.0,   -6.0,   1073.0,  0.0,   -2.0},
            { 1, 0, 0, 0, 2,  -1981.0,   0.0,    0.0,    854.0,  0.0,    0.0},
            {-1, 0, 0, 1, 0,   4026.0,   0.0, -353.0,   -553.0,  0.0, -139.0},
            { 0, 0, 2, 1, 2,   1660.0,   0.0,   -5.0,   -710.0,  0.0,   -2.0},
            
            /* 71-80 */
            {-1, 0, 2, 4, 2,  -1521.0,   0.0,    9.0,    647.0,  0.0,    4.0},
            {-1, 1, 0, 1, 1,   1314.0,   0.0,    0.0,   -700.0,  0.0,    0.0},
            { 0,-2, 2,-2, 1,  -1283.0,   0.0,    0.0,    672.0,  0.0,    0.0},
            { 1, 0, 2, 2, 1,  -1331.0,   0.0,    8.0,    663.0,  0.0,    4.0},
            {-2, 0, 2, 2, 2,   1383.0,   0.0,   -2.0,   -594.0,  0.0,   -2.0},
            {-1, 0, 0, 0, 2,   1405.0,   0.0,    4.0,   -610.0,  0.0,    2.0},
            { 1, 1, 2,-2, 2,   1290.0,   0.0,    0.0,   -556.0,  0.0,    0.0},
            {-2, 0, 2, 4, 2,  -1214.0,   0.0,    5.0,    518.0,  0.0,    2.0},
            {-1, 0, 4, 0, 2,   1146.0,   0.0,   -3.0,   -490.0,  0.0,   -1.0},
            { 2, 0, 2,-2, 1,   1019.0,   0.0,   -1.0,   -527.0,  0.0,   -1.0},
            
            /* 81-90 */
            { 2, 0, 2, 2, 2,  -1100.0,   0.0,    9.0,    465.0,  0.0,    4.0},
            { 1, 0, 0, 2, 1,   -970.0,   0.0,    2.0,    496.0,  0.0,    1.0},
            { 3, 0, 0, 0, 0,   1575.0,   0.0,   -6.0,    -50.0,  0.0,    0.0},
            { 3, 0, 2,-2, 2,    934.0,   0.0,   -3.0,   -399.0,  0.0,   -1.0},
            { 0, 0, 4,-2, 2,    922.0,   0.0,   -1.0,   -395.0,  0.0,   -1.0},
            { 0, 1, 2, 0, 1,    815.0,   0.0,   -1.0,   -422.0,  0.0,   -1.0},
            { 0, 0,-2, 2, 1,    834.0,   0.0,    2.0,   -440.0,  0.0,    1.0},
            { 0, 0, 2,-2, 3,   1248.0,   0.0,    0.0,   -170.0,  0.0,    1.0},
            {-1, 0, 0, 4, 0,   1338.0,   0.0,   -5.0,    -39.0,  0.0,    0.0},
            { 2, 0,-2, 0, 1,    716.0,   0.0,   -2.0,   -389.0,  0.0,   -1.0},
            
            /* 91-100 */
            {-2, 0, 0, 4, 0,   1282.0,   0.0,   -3.0,    -23.0,  0.0,    1.0},
            {-1,-1, 0, 2, 1,    742.0,   0.0,    1.0,   -391.0,  0.0,    0.0},
            {-1, 0, 0, 1, 1,   1020.0,   0.0,  -25.0,   -495.0,  0.0,  -10.0},
            { 0, 1, 0, 0, 2,    715.0,   0.0,   -4.0,   -326.0,  0.0,    2.0},
            { 0, 0,-2, 0, 1,   -666.0,   0.0,   -3.0,    369.0,  0.0,   -1.0},
            { 0,-1, 2, 0, 1,   -667.0,   0.0,    1.0,    346.0,  0.0,    1.0},
            { 0, 0, 2,-1, 2,   -704.0,   0.0,    0.0,    304.0,  0.0,    0.0},
            { 0, 0, 2, 4, 2,   -694.0,   0.0,    5.0,    294.0,  0.0,    2.0},
            {-2,-1, 0, 2, 0,  -1014.0,   0.0,   -1.0,      4.0,  0.0,   -1.0},
            { 1, 1, 0,-2, 1,   -585.0,   0.0,   -2.0,    316.0,  0.0,   -1.0},
            
            /* 101-110 */
            {-1, 1, 0, 2, 0,   -949.0,   0.0,    1.0,      8.0,  0.0,   -1.0},
            {-1, 1, 0, 1, 2,   -595.0,   0.0,    0.0,    258.0,  0.0,    0.0},
            { 1,-1, 0, 0, 1,    528.0,   0.0,    0.0,   -279.0,  0.0,    0.0},
            { 1,-1, 2, 2, 2,   -590.0,   0.0,    4.0,    252.0,  0.0,    2.0},
            {-1, 1, 2, 2, 2,    570.0,   0.0,   -2.0,   -244.0,  0.0,   -1.0},
            { 3, 0, 2, 0, 1,   -502.0,   0.0,    3.0,    250.0,  0.0,    2.0},
            { 0, 1,-2, 2, 0,   -875.0,   0.0,    1.0,     29.0,  0.0,    0.0},
            {-1, 0, 0,-2, 1,   -492.0,   0.0,   -3.0,    275.0,  0.0,   -1.0},
            { 0, 1, 2, 2, 2,    535.0,   0.0,   -2.0,   -228.0,  0.0,   -1.0},
            {-1,-1, 2, 2, 1,   -467.0,   0.0,    1.0,    240.0,  0.0,    1.0},
            
            /* 111-120 */
            { 0,-1, 0, 0, 2,    591.0,   0.0,    0.0,   -253.0,  0.0,    0.0},
            { 1, 0, 2,-4, 1,   -453.0,   0.0,   -1.0,    244.0,  0.0,   -1.0},
            {-1, 0,-2, 2, 0,    766.0,   0.0,    1.0,      9.0,  0.0,    0.0},
            { 0,-1, 2, 2, 1,   -446.0,   0.0,    2.0,    225.0,  0.0,    1.0},
            { 2,-1, 2, 0, 2,   -488.0,   0.0,    2.0,    207.0,  0.0,    1.0},
            { 0, 0, 0, 2, 2,   -468.0,   0.0,    0.0,    201.0,  0.0,    0.0},
            { 1,-1, 2, 0, 1,   -421.0,   0.0,    1.0,    216.0,  0.0,    1.0},
            {-1, 1, 2, 0, 2,    463.0,   0.0,    0.0,   -200.0,  0.0,    0.0},
            { 0, 1, 0, 2, 0,   -673.0,   0.0,    2.0,     14.0,  0.0,    0.0},
            { 0,-1,-2, 2, 0,    658.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 121-130 */
            { 0, 3, 2,-2, 2,   -438.0,   0.0,    0.0,    188.0,  0.0,    0.0},
            { 0, 0, 0, 1, 1,   -390.0,   0.0,    0.0,    205.0,  0.0,    0.0},
            {-1, 0, 2, 2, 0,    639.0, -11.0,   -2.0,    -19.0,  0.0,    0.0},
            { 2, 1, 2, 0, 2,    412.0,   0.0,   -2.0,   -176.0,  0.0,   -1.0},
            { 1, 1, 0, 0, 1,   -361.0,   0.0,    0.0,    189.0,  0.0,    0.0},
            { 1, 1, 2, 0, 1,    360.0,   0.0,   -1.0,   -185.0,  0.0,   -1.0},
            { 2, 0, 0, 2, 0,    588.0,   0.0,   -3.0,    -24.0,  0.0,    0.0},
            { 1, 0,-2, 2, 0,   -578.0,   0.0,    1.0,      5.0,  0.0,    0.0},
            {-1, 0, 0, 2, 2,   -396.0,   0.0,    0.0,    171.0,  0.0,    0.0},
            { 0, 1, 0, 1, 0,    565.0,   0.0,   -1.0,     -6.0,  0.0,    0.0},
            
            /* 131-140 */
            { 0, 1, 0,-2, 1,   -335.0,   0.0,   -1.0,    184.0,  0.0,   -1.0},
            {-1, 0, 2,-2, 2,    357.0,   0.0,    1.0,   -154.0,  0.0,    0.0},
            { 0, 0, 0,-1, 1,    321.0,   0.0,    1.0,   -174.0,  0.0,    0.0},
            {-1, 1, 0, 0, 1,   -301.0,   0.0,   -1.0,    162.0,  0.0,    0.0},
            { 1, 0, 2,-1, 2,   -334.0,   0.0,    0.0,    144.0,  0.0,    0.0},
            { 1,-1, 0, 2, 0,    493.0,   0.0,   -2.0,    -15.0,  0.0,    0.0},
            { 0, 0, 0, 4, 0,    494.0,   0.0,   -2.0,    -19.0,  0.0,    0.0},
            { 1, 0, 2, 1, 2,    337.0,   0.0,   -1.0,   -143.0,  0.0,   -1.0},
            { 0, 0, 2, 1, 1,    280.0,   0.0,   -1.0,   -144.0,  0.0,    0.0},
            { 1, 0, 0,-2, 2,    309.0,   0.0,    1.0,   -134.0,  0.0,    0.0},
            
            /* 141-150 */
            {-1, 0, 2, 4, 1,   -263.0,   0.0,    2.0,    131.0,  0.0,    1.0},
            { 1, 0,-2, 0, 1,    253.0,   0.0,    1.0,   -138.0,  0.0,    0.0},
            { 1, 1, 2,-2, 1,    245.0,   0.0,    0.0,   -128.0,  0.0,    0.0},
            { 0, 0, 2, 2, 0,    416.0,   0.0,   -2.0,    -17.0,  0.0,    0.0},
            {-1, 0, 2,-1, 1,   -229.0,   0.0,    0.0,    128.0,  0.0,    0.0},
            {-2, 0, 2, 2, 1,    231.0,   0.0,    0.0,   -120.0,  0.0,    0.0},
            { 4, 0, 2, 0, 2,   -259.0,   0.0,    2.0,    109.0,  0.0,    1.0},
            { 2,-1, 0, 0, 0,    375.0,   0.0,   -1.0,     -8.0,  0.0,    0.0},
            { 2, 1, 2,-2, 2,    252.0,   0.0,    0.0,   -108.0,  0.0,    0.0},
            { 0, 1, 2, 1, 2,   -245.0,   0.0,    1.0,    104.0,  0.0,    0.0},
            
            /* 151-160 */
            { 1, 0, 4,-2, 2,    243.0,   0.0,   -1.0,   -104.0,  0.0,    0.0},
            {-1,-1, 0, 0, 1,    208.0,   0.0,    1.0,   -112.0,  0.0,    0.0},
            { 0, 1, 0, 2, 1,    199.0,   0.0,    0.0,   -102.0,  0.0,    0.0},
            {-2, 0, 2, 4, 1,   -208.0,   0.0,    1.0,    105.0,  0.0,    0.0},
            { 2, 0, 2, 0, 0,    335.0,   0.0,   -2.0,    -14.0,  0.0,    0.0},
            { 1, 0, 0, 1, 0,   -325.0,   0.0,    1.0,      7.0,  0.0,    0.0},
            {-1, 0, 0, 4, 1,   -187.0,   0.0,    0.0,     96.0,  0.0,    0.0},
            {-1, 0, 4, 0, 1,    197.0,   0.0,   -1.0,   -100.0,  0.0,    0.0},
            { 2, 0, 2, 2, 1,   -192.0,   0.0,    2.0,     94.0,  0.0,    1.0},
            { 0, 0, 2,-3, 2,   -188.0,   0.0,    0.0,     83.0,  0.0,    0.0},
            
            /* 161-170 */
            {-1,-2, 0, 2, 0,    276.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2, 1, 0, 0, 0,   -286.0,   0.0,    1.0,      6.0,  0.0,    0.0},
            { 0, 0, 4, 0, 2,    186.0,   0.0,   -1.0,    -79.0,  0.0,    0.0},
            { 0, 0, 0, 0, 3,   -219.0,   0.0,    0.0,     43.0,  0.0,    0.0},
            { 0, 3, 0, 0, 0,    276.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0, 2,-4, 1,   -153.0,   0.0,   -1.0,     84.0,  0.0,    0.0},
            { 0,-1, 0, 2, 1,   -156.0,   0.0,    0.0,     81.0,  0.0,    0.0},
            { 0, 0, 0, 4, 1,   -154.0,   0.0,    1.0,     78.0,  0.0,    0.0},
            {-1,-1, 2, 4, 2,   -174.0,   0.0,    1.0,     75.0,  0.0,    0.0},
            { 1, 0, 2, 4, 2,   -163.0,   0.0,    2.0,     69.0,  0.0,    1.0},
            
            /* 171-180 */
            {-2, 2, 0, 2, 0,   -228.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-2,-1, 2, 0, 1,     91.0,   0.0,   -4.0,    -54.0,  0.0,   -2.0},
            {-2, 0, 0, 2, 2,    175.0,   0.0,    0.0,    -75.0,  0.0,    0.0},
            {-1,-1, 2, 0, 2,   -159.0,   0.0,    0.0,     69.0,  0.0,    0.0},
            { 0, 0, 4,-2, 1,    141.0,   0.0,    0.0,    -72.0,  0.0,    0.0},
            { 3, 0, 2,-2, 1,    147.0,   0.0,    0.0,    -75.0,  0.0,    0.0},
            {-2,-1, 0, 2, 1,   -132.0,   0.0,    0.0,     69.0,  0.0,    0.0},
            { 1, 0, 0,-1, 1,    159.0,   0.0,  -28.0,    -54.0,  0.0,   11.0},
            { 0,-2, 0, 2, 0,    213.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            {-2, 0, 0, 4, 1,    123.0,   0.0,    0.0,    -64.0,  0.0,    0.0},
            
            /* 181-190 */
            {-3, 0, 0, 0, 1,   -118.0,   0.0,   -1.0,     66.0,  0.0,    0.0},
            { 1, 1, 2, 2, 2,    144.0,   0.0,   -1.0,    -61.0,  0.0,    0.0},
            { 0, 0, 2, 4, 1,   -121.0,   0.0,    1.0,     60.0,  0.0,    0.0},
            { 3, 0, 2, 2, 2,   -134.0,   0.0,    1.0,     56.0,  0.0,    1.0},
            {-1, 1, 2,-2, 1,   -105.0,   0.0,    0.0,     57.0,  0.0,    0.0},
            { 2, 0, 0,-4, 1,   -102.0,   0.0,    0.0,     56.0,  0.0,    0.0},
            { 0, 0, 0,-2, 2,    120.0,   0.0,    0.0,    -52.0,  0.0,    0.0},
            { 2, 0, 2,-4, 1,    101.0,   0.0,    0.0,    -54.0,  0.0,    0.0},
            {-1, 1, 0, 2, 1,   -113.0,   0.0,    0.0,     59.0,  0.0,    0.0},
            { 0, 0, 2,-1, 1,   -106.0,   0.0,    0.0,     61.0,  0.0,    0.0},
            
            /* 191-200 */
            { 0,-2, 2, 2, 2,   -129.0,   0.0,    1.0,     55.0,  0.0,    0.0},
            { 2, 0, 0, 2, 1,   -114.0,   0.0,    0.0,     57.0,  0.0,    0.0},
            { 4, 0, 2,-2, 2,    113.0,   0.0,   -1.0,    -49.0,  0.0,    0.0},
            { 2, 0, 0,-2, 2,   -102.0,   0.0,    0.0,     44.0,  0.0,    0.0},
            { 0, 2, 0, 0, 1,    -94.0,   0.0,    0.0,     51.0,  0.0,    0.0},
            { 1, 0, 0,-4, 1,   -100.0,   0.0,   -1.0,     56.0,  0.0,    0.0},
            { 0, 2, 2,-2, 1,     87.0,   0.0,    0.0,    -47.0,  0.0,    0.0},
            {-3, 0, 0, 4, 0,    161.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-1, 1, 2, 0, 1,     96.0,   0.0,    0.0,    -50.0,  0.0,    0.0},
            {-1,-1, 0, 4, 0,    151.0,   0.0,   -1.0,     -5.0,  0.0,    0.0},
            
            /* 201-210 */
            {-1,-2, 2, 2, 2,   -104.0,   0.0,    0.0,     44.0,  0.0,    0.0},
            {-2,-1, 2, 4, 2,   -110.0,   0.0,    0.0,     48.0,  0.0,    0.0},
            { 1,-1, 2, 2, 1,   -100.0,   0.0,    1.0,     50.0,  0.0,    0.0},
            {-2, 1, 0, 2, 0,     92.0,   0.0,   -5.0,     12.0,  0.0,   -2.0},
            {-2, 1, 2, 0, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0},
            { 2, 1, 0,-2, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0},
            {-3, 0, 2, 0, 1,    -78.0,   0.0,    0.0,     41.0,  0.0,    0.0},
            {-2, 0, 2,-2, 1,    -77.0,   0.0,    0.0,     43.0,  0.0,    0.0},
            {-1, 1, 0, 2, 2,      2.0,   0.0,    0.0,     54.0,  0.0,    0.0},
            { 0,-1, 2,-1, 2,     94.0,   0.0,    0.0,    -40.0,  0.0,    0.0},
            
            /* 211-220 */
            {-1, 0, 4,-2, 2,    -93.0,   0.0,    0.0,     40.0,  0.0,    0.0},
            { 0,-2, 2, 0, 2,    -83.0,   0.0,   10.0,     40.0,  0.0,   -2.0},
            {-1, 0, 2, 1, 2,     83.0,   0.0,    0.0,    -36.0,  0.0,    0.0},
            { 2, 0, 0, 0, 2,    -91.0,   0.0,    0.0,     39.0,  0.0,    0.0},
            { 0, 0, 2, 0, 3,    128.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-2, 0, 4, 0, 2,    -79.0,   0.0,    0.0,     34.0,  0.0,    0.0},
            {-1, 0,-2, 0, 1,    -83.0,   0.0,    0.0,     47.0,  0.0,    0.0},
            {-1, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -44.0,  0.0,    0.0},
            { 3, 0, 0, 0, 1,     83.0,   0.0,    0.0,    -43.0,  0.0,    0.0},
            {-1, 0, 2, 3, 2,     91.0,   0.0,    0.0,    -39.0,  0.0,    0.0},
            
            /* 221-230 */
            { 2,-1, 2, 0, 1,    -77.0,   0.0,    0.0,     39.0,  0.0,    0.0},
            { 0, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -43.0,  0.0,    0.0},
            { 0,-1, 2, 4, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0},
            { 2,-1, 2, 2, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0},
            { 0, 2,-2, 2, 0,    -94.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 2,-1, 1,     68.0,   0.0,    0.0,    -36.0,  0.0,    0.0},
            { 0,-2, 0, 0, 1,    -61.0,   0.0,    0.0,     32.0,  0.0,    0.0},
            { 1, 0, 2,-4, 2,     71.0,   0.0,    0.0,    -31.0,  0.0,    0.0},
            { 1,-1, 0,-2, 1,     62.0,   0.0,    0.0,    -34.0,  0.0,    0.0},
            {-1,-1, 2, 0, 1,    -63.0,   0.0,    0.0,     33.0,  0.0,    0.0},
            
            /* 231-240 */
            { 1,-1, 2,-2, 2,    -73.0,   0.0,    0.0,     32.0,  0.0,    0.0},
            {-2,-1, 0, 4, 0,    115.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 0, 0, 3, 0,   -103.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2,-1, 2, 2, 2,     63.0,   0.0,    0.0,    -28.0,  0.0,    0.0},
            { 0, 2, 2, 0, 2,     74.0,   0.0,    0.0,    -32.0,  0.0,    0.0},
            { 1, 1, 0, 2, 0,   -103.0,   0.0,   -3.0,      3.0,  0.0,   -1.0},
            { 2, 0, 2,-1, 2,    -69.0,   0.0,    0.0,     30.0,  0.0,    0.0},
            { 1, 0, 2, 1, 1,     57.0,   0.0,    0.0,    -29.0,  0.0,    0.0},
            { 4, 0, 0, 0, 0,     94.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 2, 1, 2, 0, 1,     64.0,   0.0,    0.0,    -33.0,  0.0,    0.0},
            
            /* 241-250 */
            { 3,-1, 2, 0, 2,    -63.0,   0.0,    0.0,     26.0,  0.0,    0.0},
            {-2, 2, 0, 2, 1,    -38.0,   0.0,    0.0,     20.0,  0.0,    0.0},
            { 1, 0, 2,-3, 1,    -43.0,   0.0,    0.0,     24.0,  0.0,    0.0},
            { 1, 1, 2,-4, 1,    -45.0,   0.0,    0.0,     23.0,  0.0,    0.0},
            {-1,-1, 2,-2, 1,     47.0,   0.0,    0.0,    -24.0,  0.0,    0.0},
            { 0,-1, 0,-1, 1,    -48.0,   0.0,    0.0,     25.0,  0.0,    0.0},
            { 0,-1, 0,-2, 1,     45.0,   0.0,    0.0,    -26.0,  0.0,    0.0},
            {-2, 0, 0, 0, 2,     56.0,   0.0,    0.0,    -25.0,  0.0,    0.0},
            {-2, 0,-2, 2, 0,     88.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1, 0,-2, 4, 0,    -75.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 251-260 */
            { 1,-2, 0, 0, 0,     85.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 0, 1, 1,     49.0,   0.0,    0.0,    -26.0,  0.0,    0.0},
            {-1, 2, 0, 2, 0,    -74.0,   0.0,   -3.0,     -1.0,  0.0,   -1.0},
            { 1,-1, 2,-2, 1,    -39.0,   0.0,    0.0,     21.0,  0.0,    0.0},
            { 1, 2, 2,-2, 2,     45.0,   0.0,    0.0,    -20.0,  0.0,    0.0},
            { 2,-1, 2,-2, 2,     51.0,   0.0,    0.0,    -22.0,  0.0,    0.0},
            { 1, 0, 2,-1, 1,    -40.0,   0.0,    0.0,     21.0,  0.0,    0.0},
            { 2, 1, 2,-2, 1,     41.0,   0.0,    0.0,    -21.0,  0.0,    0.0},
            {-2, 0, 0,-2, 1,    -42.0,   0.0,    0.0,     24.0,  0.0,    0.0},
            { 1,-2, 2, 0, 2,    -51.0,   0.0,    0.0,     22.0,  0.0,    0.0},
            
            /* 261-270 */
            { 0, 1, 2, 1, 1,    -42.0,   0.0,    0.0,     22.0,  0.0,    0.0},
            { 1, 0, 4,-2, 1,     39.0,   0.0,    0.0,    -21.0,  0.0,    0.0},
            {-2, 0, 4, 2, 2,     46.0,   0.0,    0.0,    -18.0,  0.0,    0.0},
            { 1, 1, 2, 1, 2,    -53.0,   0.0,    0.0,     22.0,  0.0,    0.0},
            { 1, 0, 0, 4, 0,     82.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 1, 0, 2, 2, 0,     81.0,   0.0,   -1.0,     -4.0,  0.0,    0.0},
            { 2, 0, 2, 1, 2,     47.0,   0.0,    0.0,    -19.0,  0.0,    0.0},
            { 3, 1, 2, 0, 2,     53.0,   0.0,    0.0,    -23.0,  0.0,    0.0},
            { 4, 0, 2, 0, 1,    -45.0,   0.0,    0.0,     22.0,  0.0,    0.0},
            {-2,-1, 2, 0, 0,    -44.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 271-280 */
            { 0, 1,-2, 2, 1,    -33.0,   0.0,    0.0,     16.0,  0.0,    0.0},
            { 1, 0,-2, 1, 0,    -61.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 0,-1,-2, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
            { 2,-1, 0,-2, 1,    -38.0,   0.0,    0.0,     19.0,  0.0,    0.0},
            {-1, 0, 2,-1, 2,    -33.0,   0.0,    0.0,     21.0,  0.0,    0.0},
            { 1, 0, 2,-3, 2,    -60.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 2,-2, 3,     48.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
            { 0, 0, 2,-3, 1,     27.0,   0.0,    0.0,    -14.0,  0.0,    0.0},
            {-1, 0,-2, 2, 1,     38.0,   0.0,    0.0,    -20.0,  0.0,    0.0},
            { 0, 0, 2,-4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0},
            
            /* 281-290 */
            {-2, 1, 0, 0, 1,    -29.0,   0.0,    0.0,     15.0,  0.0,    0.0},
            {-1, 0, 0,-1, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
            { 2, 0, 2,-4, 2,    -32.0,   0.0,    0.0,     15.0,  0.0,    0.0},
            { 0, 0, 4,-4, 4,     45.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
            { 0, 0, 4,-4, 2,    -44.0,   0.0,    0.0,     19.0,  0.0,    0.0},
            {-1,-2, 0, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
            {-2, 0, 0, 3, 0,    -51.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0,-2, 2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0},
            {-3, 0, 2, 2, 2,     44.0,   0.0,    0.0,    -19.0,  0.0,    0.0},
            {-3, 0, 2, 2, 1,     26.0,   0.0,    0.0,    -14.0,  0.0,    0.0},
            
            /* 291-300 */
            {-2, 0, 2, 2, 0,    -60.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2,-1, 0, 0, 1,     35.0,   0.0,    0.0,    -18.0,  0.0,    0.0},
            {-2, 1, 2, 2, 2,    -27.0,   0.0,    0.0,     11.0,  0.0,    0.0},
            { 1, 1, 0, 1, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 0, 1, 4,-2, 2,     36.0,   0.0,    0.0,    -15.0,  0.0,    0.0},
            {-1, 1, 0,-2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0},
            { 0, 0, 0,-4, 1,    -35.0,   0.0,    0.0,     19.0,  0.0,    0.0},
            { 1,-1, 0, 2, 1,    -37.0,   0.0,    0.0,     19.0,  0.0,    0.0},
            { 1, 1, 0, 2, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
            {-1, 2, 2, 2, 2,     35.0,   0.0,    0.0,    -14.0,  0.0,    0.0},
            
            /* 301-310 */
            { 3, 1, 2,-2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0},
            { 0,-1, 0, 4, 0,     65.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2,-1, 0, 2, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 0, 0, 4, 0, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
            { 2, 0, 4,-2, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
            {-1,-1, 2, 4, 1,    -30.0,   0.0,    0.0,     15.0,  0.0,    0.0},
            { 1, 0, 0, 4, 1,    -32.0,   0.0,    0.0,     16.0,  0.0,    0.0},
            { 1,-2, 2, 2, 2,    -31.0,   0.0,    0.0,     13.0,  0.0,    0.0},
            { 0, 0, 2, 3, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0},
            {-1, 1, 2, 4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0},
            
            /* 311-320 */
            { 3, 0, 0, 2, 0,     49.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 0, 4, 2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0},
            { 1, 1, 2, 2, 1,     23.0,   0.0,    0.0,    -12.0,  0.0,    0.0},
            {-2, 0, 2, 6, 2,    -43.0,   0.0,    0.0,     18.0,  0.0,    0.0},
            { 2, 1, 2, 2, 2,     26.0,   0.0,    0.0,    -11.0,  0.0,    0.0},
            {-1, 0, 2, 6, 2,    -32.0,   0.0,    0.0,     14.0,  0.0,    0.0},
            { 1, 0, 2, 4, 1,    -29.0,   0.0,    0.0,     14.0,  0.0,    0.0},
            { 2, 0, 2, 4, 2,    -27.0,   0.0,    0.0,     12.0,  0.0,    0.0},
            { 1, 1,-2, 1, 0,     30.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3, 1, 2, 1, 2,    -11.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            
            /* 321-330 */
            { 2, 0,-2, 0, 2,    -21.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            {-1, 0, 0, 1, 2,    -34.0,   0.0,    0.0,     15.0,  0.0,    0.0},
            {-4, 0, 2, 2, 1,    -10.0,   0.0,    0.0,      6.0,  0.0,    0.0},
            {-1,-1, 0, 1, 0,    -36.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0,-2, 2, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 1, 0, 0,-1, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            { 0,-1, 2,-2, 3,    -21.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            {-2, 1, 2, 0, 0,    -29.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 0, 0, 2,-2, 4,    -15.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-2,-2, 0, 2, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 331-340 */
            {-2, 0,-2, 4, 0,     28.0,   0.0,    0.0,      0.0,  0.0,   -2.0},
            { 0,-2,-2, 2, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 2, 0,-2, 1,    -22.0,   0.0,    0.0,     12.0,  0.0,    0.0},
            { 3, 0, 0,-4, 1,    -14.0,   0.0,    0.0,      7.0,  0.0,    0.0},
            {-1, 1, 2,-2, 2,     24.0,   0.0,    0.0,    -11.0,  0.0,    0.0},
            { 1,-1, 2,-4, 1,     11.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
            { 1, 1, 0,-2, 2,     14.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
            {-3, 0, 2, 0, 0,     24.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3, 0, 2, 0, 2,     18.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
            {-2, 0, 0, 1, 0,    -38.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 341-350 */
            { 0, 0,-2, 1, 0,    -31.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3, 0, 0, 2, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0},
            {-1,-1,-2, 2, 0,     29.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 2,-4, 1,    -18.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            { 2, 1, 0,-4, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            { 0, 2, 0,-2, 1,    -17.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            { 1, 0, 0,-3, 1,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            {-2, 0, 2,-2, 2,     16.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
            {-2,-1, 0, 0, 1,     22.0,   0.0,    0.0,    -12.0,  0.0,    0.0},
            {-4, 0, 0, 2, 0,     20.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 351-360 */
            { 1, 1, 0,-4, 1,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0},
            {-1, 0, 2,-4, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0},
            { 0, 0, 4,-4, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0},
            { 0, 3, 2,-2, 2,      0.0,   0.0,    0.0,     -7.0,  0.0,    0.0},
            {-3,-1, 0, 4, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3, 0, 0, 4, 1,     19.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
            { 1,-1,-2, 2, 0,    -34.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 0, 2, 2,    -20.0,   0.0,    0.0,      8.0,  0.0,    0.0},
            { 1,-2, 0, 0, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            { 1,-1, 0, 0, 2,    -18.0,   0.0,    0.0,      7.0,  0.0,    0.0},
            
            /* 361-370 */
            { 0, 0, 0, 1, 2,     13.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
            {-1,-1, 2, 0, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1,-2, 2,-2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            { 0,-1, 2,-1, 1,     15.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
            {-1, 0, 2, 0, 3,    -11.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 1, 1, 0, 0, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            {-1, 1, 2, 0, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 2, 0, 0, 0,    -35.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 2, 2, 0, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            {-1, 0, 4,-2, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            
            /* 371-380 */
            { 3, 0, 2,-4, 2,    -26.0,   0.0,    0.0,     11.0,  0.0,    0.0},
            { 1, 2, 2,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 1, 0, 4,-4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            {-2,-1, 0, 4, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
            { 0,-1, 0, 2, 2,    -21.0,   0.0,    0.0,      9.0,  0.0,    0.0},
            {-2, 1, 0, 4, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2,-1, 2, 2, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            { 2, 0,-2, 2, 0,    -29.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0, 0, 1, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            { 0, 1, 0, 2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            
            /* 381-390 */
            { 1,-1, 2,-1, 2,     22.0,   0.0,    0.0,     -9.0,  0.0,    0.0},
            {-2, 0, 4, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            { 2, 1, 0, 0, 1,    -20.0,   0.0,    0.0,     11.0,  0.0,    0.0},
            { 0, 1, 2, 0, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0,-1, 4,-2, 2,    -17.0,   0.0,    0.0,      7.0,  0.0,    0.0},
            { 0, 0, 4,-2, 4,     15.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 0, 2, 2, 0, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            {-3, 0, 0, 6, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 0, 4, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0},
            { 1,-2, 0, 2, 0,     25.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 391-400 */
            {-1, 0, 0, 4, 2,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0},
            {-1,-2, 2, 2, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0},
            {-1, 0, 0,-2, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            { 1, 0,-2,-2, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0},
            { 0, 0,-2,-2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0},
            {-2, 0,-2, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            { 0, 0, 0, 3, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0},
            { 0, 0, 0, 3, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 1, 0, 4, 0,    -22.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 2, 2, 0,     28.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            
            /* 401-410 */
            {-2, 0, 2, 3, 2,     15.0,   0.0,    0.0,     -7.0,  0.0,    0.0},
            { 1, 0, 0, 2, 2,     23.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
            { 0,-1, 2, 1, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            { 3,-1, 0, 0, 0,     29.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 2, 0, 0, 1, 0,    -25.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 1,-1, 2, 0, 0,     22.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0, 2, 1, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0, 2, 0, 3,     15.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 3, 1, 0, 0, 0,    -23.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 3,-1, 2,-2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0},
            
            /* 411-420 */
            { 2, 0, 2,-1, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 1, 1, 2, 0, 0,    -19.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0, 4,-1, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 1, 2, 2, 0, 2,     21.0,   0.0,    0.0,     -9.0,  0.0,    0.0},
            {-2, 0, 0, 6, 0,     23.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 0,-1, 0, 4, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0},
            {-2,-1, 2, 4, 1,    -19.0,   0.0,    0.0,      9.0,  0.0,    0.0},
            { 0,-2, 2, 2, 1,    -22.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            { 0,-1, 2, 2, 0,     27.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-1, 0, 2, 3, 1,     16.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
            
            /* 421-430 */
            {-2, 1, 2, 4, 2,     19.0,   0.0,    0.0,     -8.0,  0.0,    0.0},
            { 2, 0, 0, 2, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 2,-2, 2, 0, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            {-1, 1, 2, 3, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 3, 0, 2,-1, 2,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 4, 0, 2,-2, 1,     18.0,   0.0,    0.0,     -9.0,  0.0,    0.0},
            {-1, 0, 0, 6, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-1,-2, 2, 4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            {-3, 0, 2, 6, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0},
            {-1, 0, 2, 4, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            
            /* 431-440 */
            { 3, 0, 0, 2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0},
            { 3,-1, 2, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 3, 0, 2, 0, 0,     30.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 1, 0, 4, 0, 2,     24.0,   0.0,    0.0,    -10.0,  0.0,    0.0},
            { 5, 0, 2,-2, 2,     10.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 0,-1, 2, 4, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0},
            { 2,-1, 2, 2, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0},
            { 0, 1, 2, 4, 2,     17.0,   0.0,    0.0,     -7.0,  0.0,    0.0},
            { 1,-1, 2, 4, 2,    -24.0,   0.0,    0.0,     10.0,  0.0,    0.0},
            { 3,-1, 2, 2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            
            /* 441-450 */
            { 3, 0, 2, 2, 1,    -24.0,   0.0,    0.0,     11.0,  0.0,    0.0},
            { 5, 0, 2, 0, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0},
            { 0, 0, 2, 6, 2,    -13.0,   0.0,    0.0,      5.0,  0.0,    0.0},
            { 4, 0, 2, 2, 2,    -15.0,   0.0,    0.0,      7.0,  0.0,    0.0},
            { 0,-1, 1,-1, 1,      0.0,   0.0,-1988.0,      0.0,  0.0,-1679.0},
            {-1, 0, 1, 0, 3,      0.0,   0.0,  -63.0,      0.0,  0.0,  -27.0},
            { 0,-2, 2,-2, 3,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0,-1, 0, 1,      0.0,   0.0,    5.0,      0.0,  0.0,    4.0},
            { 2,-2, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            {-1, 0, 1, 0, 2,      0.0,   0.0,  364.0,      0.0,  0.0,  176.0},
            
            /* 451-460 */
            {-1, 0, 1, 0, 1,      0.0,   0.0,-1044.0,      0.0,  0.0, -891.0},
            {-1,-1, 2,-1, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-2, 2, 0, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 0, 1, 0, 0,      0.0,   0.0,  330.0,      0.0,  0.0,    0.0},
            {-4, 1, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-3, 0, 2, 1, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-2,-1, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 1, 0,-2, 1, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2,-1,-2, 0, 1,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-4, 0, 2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 461-470 */
            {-3, 1, 0, 3, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 0,-1, 2, 0,      0.0,   0.0,    5.0,      0.0,  0.0,    0.0},
            { 0,-2, 0, 0, 2,      0.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 0,-2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-3, 0, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2,-1, 0, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 0,-2, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-4, 0, 0, 4, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2, 1,-2, 0, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 2,-1, 0,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            
            /* 471-480 */
            { 0, 0, 1,-1, 0,     -5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 2, 0, 1, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2, 1, 2, 0, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 1, 1, 0,-1, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 1, 0, 1,-2, 1,      0.0,   0.0,  -12.0,      0.0,  0.0,  -10.0},
            { 0, 2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 1,-1, 2,-3, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 1, 2,-1, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2, 0, 4,-2, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-2, 0, 4,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            
            /* 481-490 */
            {-2,-2, 0, 2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-2, 0,-2, 4, 0,      0.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 2, 2,-4, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 1, 1, 2,-4, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            {-1, 2, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2, 0, 0,-3, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 2, 0, 0, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 0, 0, 0,-2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 2,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1, 1, 0, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 491-500 */
            { 0, 0, 0,-1, 2,     -8.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-2, 1, 0, 1, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1,-2, 0,-2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 1, 0,-2, 0, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-3, 1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 1,-2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 0, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-3, 0, 0, 2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3,-1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            
            /* 501-510 */
            { 0, 1, 2,-4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-2, 1, 2,-2, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0,-1, 2,-4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 0, 1, 0,-2, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            {-1, 0, 0,-2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2, 0,-2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-4, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1,-1, 0,-1, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0,-2, 0, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            
            /* 511-520 */
            {-3, 0, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 0,-2, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2, 0,-2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 0, 0,-4, 2, 0,      8.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2,-1,-2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1, 0, 2,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 1, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 2, 1, 2,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 2, 1, 2,-4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            
            /* 521-530 */
            { 0, 1, 4,-4, 4,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 4,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-1,-1,-2, 4, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 0,-2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2,-1, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0,-2, 3, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2, 0, 0, 3, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 0,-1, 0, 1, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3, 0, 2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 531-540 */
            { 1, 1,-2, 2, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 1, 0, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 1,-2, 2,-2, 1,     10.0,   0.0,   13.0,      6.0,  0.0,   -5.0},
            { 0, 0, 1, 0, 2,      0.0,   0.0,   30.0,      0.0,  0.0,   14.0},
            { 0, 0, 1, 0, 1,      0.0,   0.0, -162.0,      0.0,  0.0, -138.0},
            { 0, 0, 1, 0, 0,      0.0,   0.0,   75.0,      0.0,  0.0,    0.0},
            {-1, 2, 0, 2, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            { 0, 0, 2, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2, 0, 2, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2, 0, 0,-1, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 541-550 */
            { 3, 0, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 1, 0, 2,-2, 3,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 2, 0, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2, 0, 2,-3, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1, 1, 4,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2,-2, 0, 4, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0,-2, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 0, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2, 0, 0, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            
            /* 551-560 */
            {-1, 0, 0, 3, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2,-2, 0, 0, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1,-1, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0,-2, 2, 0, 1,     -6.0,   0.0,   -3.0,      3.0,  0.0,    1.0},
            {-1, 0, 1, 2, 1,      0.0,   0.0,   -3.0,      0.0,  0.0,   -2.0},
            {-1, 1, 0, 3, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-1, 2, 1, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 0,-1, 2, 0, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2, 1, 2, 2, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            
            /* 561-570 */
            { 2,-2, 2,-2, 2,     -1.0,   0.0,    3.0,      3.0,  0.0,   -1.0},
            { 1, 1, 0, 1, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 1, 0, 1, 0, 1,      0.0,   0.0,  -13.0,      0.0,  0.0,  -11.0},
            { 1, 0, 1, 0, 0,      3.0,   0.0,    6.0,      0.0,  0.0,    0.0},
            { 0, 2, 0, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2,-1, 2,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 0,-1, 4,-2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 0, 0, 4,-2, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 4, 0, 2,-4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            
            /* 571-580 */
            { 2, 2, 2,-2, 2,      8.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 2, 0, 4,-4, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1,-2, 0, 4, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1,-3, 2, 2, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-3, 0, 2, 4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-3, 0, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1,-1, 0,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            {-3, 0, 0, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-3, 0,-2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 0,-4, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            
            /* 581-590 */
            {-2, 1, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-4, 0, 0, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            {-1, 0, 0,-4, 1,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-3, 0, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0, 0, 3, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-1, 1, 0, 4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 1,-2, 2, 0, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 0, 1, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-1, 0, 2, 2, 3,      6.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 0, 0, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 591-600 */
            {-2, 0, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 3, 0, 0, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2, 1, 0, 1, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2,-1, 2,-1, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 0, 0, 2, 0, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0, 3, 0, 3,      0.0,   0.0,  -26.0,      0.0,  0.0,  -11.0},
            { 0, 0, 3, 0, 2,      0.0,   0.0,  -10.0,      0.0,  0.0,   -5.0},
            {-1, 2, 2, 2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            {-1, 0, 4, 0, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 601-610 */
            { 1, 2, 2, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 3, 1, 2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 1, 1, 4,-2, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            {-2,-1, 0, 6, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0,-2, 0, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-2, 0, 0, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0,-3, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0, 0, 4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-1,-1, 2, 3, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 611-620 */
            {-2, 0, 2, 4, 0,     13.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2,-1, 0, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 1, 0, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 0, 4, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 0, 1, 0, 4, 0,    -11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1,-1, 2, 1, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 0, 0, 2, 2, 3,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0, 2, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-1, 0, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-2, 0, 4, 2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            
            /* 621-630 */
            { 2, 1, 0, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2, 1, 0, 2, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2,-1, 2, 0, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0, 2, 1, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2, 0, 2, 0, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 3, 0, 2, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 1, 0, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 1, 0, 3, 0, 3,      0.0,   0.0,   -5.0,      0.0,  0.0,   -2.0},
            { 1, 1, 2, 1, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0},
            
            /* 631-640 */
            { 0, 2, 2, 2, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 2, 1, 2, 0, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2, 0, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 4, 1, 2,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            {-1,-1, 0, 6, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            {-3,-1, 2, 6, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            {-1, 0, 0, 6, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-3, 0, 2, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 1,-1, 0, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 1,-1, 0, 4, 0,     12.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 641-650 */
            {-2, 0, 2, 5, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 1,-2, 2, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 3,-1, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1,-1, 2, 2, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0, 2, 3, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            {-1, 1, 2, 4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 0, 1, 2, 3, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-1, 0, 4, 2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2, 0, 2, 1, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 5, 0, 0, 0, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 651-660 */
            { 2, 1, 2, 1, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 1, 0, 4, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 3, 1, 2, 0, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0},
            { 3, 0, 4,-2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            {-2,-1, 2, 6, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0, 0, 6, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            {-2, 0, 2, 6, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0},
            { 2, 0, 0, 4, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 2, 0, 0, 4, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            
            /* 661-670 */
            { 2,-2, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 0, 0, 2, 4, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 1, 0, 2, 3, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0},
            { 4, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 2, 0, 2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0},
            { 0, 0, 4, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 4,-1, 2, 0, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 3, 0, 2, 1, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 2, 1, 2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 4, 1, 2, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            
            /* 671-678 */
            {-1,-1, 2, 6, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            {-1, 0, 2, 6, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 1,-1, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0},
            { 1, 1, 2, 4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0},
            { 3, 1, 2, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0},
            { 5, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 2,-1, 2, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0},
            { 2, 0, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0}
        };
        
        /* Number of terms in the luni-solar nutation model */
        const int NLS = (int) (sizeof xls / sizeof xls[0]);
        
        /* ------------------------ */
        /* Planetary nutation model */
        /* ------------------------ */
        
        /* The units for the sine and cosine coefficients are */
        /* 0.1 microarcsecond                                 */
        
        static const struct {
            int nl,               /* coefficients of l, F, D and Omega */
            nf,
            nd,
            nom,
            nme,              /* coefficients of planetary longitudes */
            nve,
            nea,
            nma,
            nju,
            nsa,
            nur,
            nne,
            npa;              /* coefficient of general precession */
            int sp,cp;            /* longitude sin, cos coefficients */
            int se,ce;            /* obliquity sin, cos coefficients */
        } xpl[] = {
            
            /* 1-10 */
            { 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 0, 1440,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0, -8, 16,-4,-5, 0, 0, 2,   56,-117,  -42, -40},
            { 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 2,  125, -43,    0, -54},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0,-1, 2, 2,    0,   5,    0,   0},
            { 0, 0, 0, 0, 0,  0, -4,  8,-1,-5, 0, 0, 2,    3,  -7,   -3,   0},
            { 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 1,    3,   0,    0,  -2},
            { 0, 1,-1, 1, 0,  0,  3, -8, 3, 0, 0, 0, 0, -114,   0,    0,  61},
            {-1, 0, 0, 0, 0, 10, -3,  0, 0, 0, 0, 0, 0, -219,  89,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0,-2, 6,-3, 0, 2,   -3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0, -462,1604,    0,   0},
            
            /* 11-20 */
            { 0, 1,-1, 1, 0,  0, -5,  8,-3, 0, 0, 0, 0,   99,   0,    0, -53},
            { 0, 0, 0, 0, 0,  0, -4,  8,-3, 0, 0, 0, 1,   -3,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  4, -8, 1, 5, 0, 0, 2,    0,   6,    2,   0},
            { 0, 0, 0, 0, 0, -5,  6,  4, 0, 0, 0, 0, 2,    3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 2,  -12,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 1,   14,-218,  117,   8},
            { 0, 1,-1, 1, 0,  0, -1,  0, 2,-5, 0, 0, 0,   31,-481, -257, -17},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 0, -491, 128,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0,-2, 5, 0, 0, 0,-3084,5123, 2735,1647},
            { 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 1,-1444,2409,-1286,-771},
            
            /* 21-30 */
            { 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 2,   11, -24,  -11,  -9},
            { 2,-1,-1, 0, 0,  0,  3, -7, 0, 0, 0, 0, 0,   26,  -9,    0,   0},
            { 1, 0,-2, 0, 0, 19,-21,  3, 0, 0, 0, 0, 0,  103, -60,    0,   0},
            { 0, 1,-1, 1, 0,  2, -4,  0,-3, 0, 0, 0, 0,    0, -13,   -7,   0},
            { 1, 0,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -26, -29,  -16,  14},
            { 0, 1,-1, 1, 0,  0, -1,  0,-4,10, 0, 0, 0,    9, -27,  -14,  -5},
            {-2, 0, 2, 1, 0,  0,  2,  0, 0,-5, 0, 0, 0,   12,   0,    0,  -6},
            { 0, 0, 0, 0, 0,  3, -7,  4, 0, 0, 0, 0, 0,   -7,   0,    0,   0},
            { 0,-1, 1, 0, 0,  0,  1,  0, 1,-1, 0, 0, 0,    0,  24,    0,   0},
            {-2, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,  284,   0,    0,-151},
            
            /* 31-40 */
            {-1, 0, 0, 0, 0, 18,-16,  0, 0, 0, 0, 0, 0,  226, 101,    0,   0},
            {-2, 1, 1, 2, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -8,   -2,   0},
            {-1, 1,-1, 1, 0, 18,-17,  0, 0, 0, 0, 0, 0,    0,  -6,   -3,   0},
            {-1, 0, 1, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -3},
            { 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 2,  -41, 175,   76,  17},
            { 0, 2,-2, 2, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,  15,    6,   0},
            { 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 1,  425, 212, -133, 269},
            { 0, 1,-1, 1, 0, -8, 12,  0, 0, 0, 0, 0, 0, 1200, 598,  319,-641},
            { 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 0,  235, 334,    0,   0},
            { 0, 1,-1, 1, 0,  8,-14,  0, 0, 0, 0, 0, 0,   11, -12,   -7,  -6},
            
            /* 41-50 */
            { 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 1,    5,  -6,    3,   3},
            {-2, 0, 2, 1, 0,  0,  2,  0,-4, 5, 0, 0, 0,   -5,   0,    0,   3},
            {-2, 0, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3},
            {-2, 0, 2, 0, 0,  0,  2,  0,-3, 1, 0, 0, 0,   15,   0,    0,   0},
            { 0, 0, 0, 1, 0,  3, -5,  0, 2, 0, 0, 0, 0,   13,   0,    0,  -7},
            {-2, 0, 2, 0, 0,  0,  2,  0,-4, 3, 0, 0, 0,   -6,  -9,    0,   0},
            { 0,-1, 1, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,  266, -78,    0,   0},
            { 0, 0, 0, 1, 0,  0, -1,  2, 0, 0, 0, 0, 0, -460,-435, -232, 246},
            { 0, 1,-1, 2, 0,  0, -2,  2, 0, 0, 0, 0, 0,    0,  15,    7,   0},
            {-1, 1, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
            
            /* 51-60 */
            {-1, 0, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0, 131,    0,   0},
            {-2, 0, 2, 0, 0,  0,  2,  0,-2,-2, 0, 0, 0,    4,   0,    0,   0},
            {-2, 2, 0, 2, 0,  0, -5,  9, 0, 0, 0, 0, 0,    0,   3,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0, 0,-1, 0, 0,    0,   4,    2,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 1, 0, 0,    0,   3,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 2, 0,  -17, -19,  -10,   9},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 1,   -9, -11,    6,  -5},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 2,   -6,   0,    0,   3},
            {-1, 0, 1, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -16,   8,    0,   0},
            { 0,-1, 1, 0, 0,  0,  1,  0, 0, 2, 0, 0, 0,    0,   3,    0,   0},
            
            /* 61-70 */
            { 0, 1,-1, 2, 0,  0, -1,  0, 0, 2, 0, 0, 0,   11,  24,   11,  -5},
            { 0, 0, 0, 1, 0,  0, -9, 17, 0, 0, 0, 0, 0,   -3,  -4,   -2,   1},
            { 0, 0, 0, 2, 0, -3,  5,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 0, 1,-1, 1, 0,  0, -1,  0,-1, 2, 0, 0, 0,    0,  -8,   -4,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 1,-2, 0, 0, 0,    0,   3,    0,   0},
            { 1, 0,-2, 0, 0, 17,-16,  0,-2, 0, 0, 0, 0,    0,   5,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0, 1,-3, 0, 0, 0,    0,   3,    2,   0},
            {-2, 0, 2, 1, 0,  0,  5, -6, 0, 0, 0, 0, 0,   -6,   4,    2,   3},
            { 0,-2, 2, 0, 0,  0,  9,-13, 0, 0, 0, 0, 0,   -3,  -5,    0,   0},
            { 0, 1,-1, 2, 0,  0, -1,  0, 0, 1, 0, 0, 0,   -5,   0,    0,   2},
            
            /* 71-80 */
            { 0, 0, 0, 1, 0,  0,  0,  0, 0, 1, 0, 0, 0,    4,  24,   13,  -2},
            { 0,-1, 1, 0, 0,  0,  1,  0, 0, 1, 0, 0, 0,  -42,  20,    0,   0},
            { 0,-2, 2, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,  -10, 233,    0,   0},
            { 0,-1, 1, 1, 0,  5, -7,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
            {-2, 0, 2, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   78, -18,    0,   0},
            { 2, 1,-3, 1, 0, -6,  7,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
            { 0, 0, 0, 2, 0,  0,  0,  0, 1, 0, 0, 0, 0,    0,  -3,   -1,   0},
            { 0,-1, 1, 1, 0,  0,  1,  0, 1, 0, 0, 0, 0,    0,  -4,   -2,   1},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 2, 0, 0,    0,  -8,   -4,  -1},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 1,    0,  -5,    3,   0},
            
            /* 81-90 */
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 2,   -7,   0,    0,   3},
            { 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 2,  -14,   8,    3,   6},
            { 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 1,    0,   8,   -4,   0},
            { 0, 1,-1, 1, 0,  0, -9, 15, 0, 0, 0, 0, 0,    0,  19,   10,   0},
            { 0, 0, 0, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   45, -22,    0,   0},
            { 1,-1,-1, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 2, 0,-2, 0, 0,  2, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
            {-2, 0, 2, 0, 0,  0,  2,  0,-5, 5, 0, 0, 0,    0,   3,    0,   0},
            { 2, 0,-2, 1, 0,  0, -6,  8, 0, 0, 0, 0, 0,    3,   5,    3,  -2},
            { 2, 0,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   89, -16,   -9, -48},
            
            /* 91-100 */
            {-2, 1, 1, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,    0,   3,    0,   0},
            {-2, 1, 1, 1, 0,  0,  1,  0,-3, 0, 0, 0, 0,   -3,   7,    4,   2},
            {-2, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0, -349, -62,    0,   0},
            {-2, 0, 2, 0, 0,  0,  6, -8, 0, 0, 0, 0, 0,  -15,  22,    0,   0},
            {-2, 0, 2, 0, 0,  0,  2,  0,-1,-5, 0, 0, 0,   -3,   0,    0,   0},
            {-1, 0, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,  -53,   0,    0,   0},
            {-1, 1, 1, 1, 0,-20, 20,  0, 0, 0, 0, 0, 0,    5,   0,    0,  -3},
            { 1, 0,-2, 0, 0, 20,-21,  0, 0, 0, 0, 0, 0,    0,  -8,    0,   0},
            { 0, 0, 0, 1, 0,  0,  8,-15, 0, 0, 0, 0, 0,   15,  -7,   -4,  -8},
            { 0, 2,-2, 1, 0,  0,-10, 15, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
            
            /* 101-110 */
            { 0,-1, 1, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,  -21, -78,    0,   0},
            { 0, 0, 0, 1, 0,  0,  0,  0, 1, 0, 0, 0, 0,   20, -70,  -37, -11},
            { 0, 1,-1, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,   6,    3,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0,-2, 4, 0, 0, 0,    5,   3,    2,  -2},
            { 2, 0,-2, 1, 0, -6,  8,  0, 0, 0, 0, 0, 0,  -17,  -4,   -2,   9},
            { 0,-2, 2, 1, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   6,    3,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0,-1, 0, 0, 1,   32,  15,   -8,  17},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0,-1, 0, 0, 0,  174,  84,   45, -93},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 0,   11,  56,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0, 1, 0, 0, 0,  -66, -12,   -6,  35},
            
            /* 111-120 */
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 1,   47,   8,    4, -25},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 2,    0,   8,    4,   0},
            { 0, 2,-2, 1, 0,  0, -9, 13, 0, 0, 0, 0, 0,   10, -22,  -12,  -5},
            { 0, 0, 0, 1, 0,  0,  7,-13, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
            {-2, 0, 2, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,  -24,  12,    0,   0},
            { 0, 0, 0, 0, 0,  0,  9,-17, 0, 0, 0, 0, 0,    5,  -6,    0,   0},
            { 0, 0, 0, 0, 0,  0, -9, 17, 0, 0, 0, 0, 2,    3,   0,    0,  -2},
            { 1, 0,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    4,   3,    1,  -2},
            { 1, 0,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  29,   15,   0},
            { 0, 0, 0, 2, 0,  0, -1,  2, 0, 0, 0, 0, 0,   -5,  -4,   -2,   2},
            
            /* 121-130 */
            { 0,-1, 1, 1, 0,  0,  0,  2, 0, 0, 0, 0, 0,    8,  -3,   -1,  -5},
            { 0,-2, 2, 0, 1,  0, -2,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
            { 0, 0, 0, 0, 0,  3, -5,  0, 2, 0, 0, 0, 0,   10,   0,    0,   0},
            {-2, 0, 2, 1, 0,  0,  2,  0,-3, 1, 0, 0, 0,    3,   0,    0,  -2},
            {-2, 0, 2, 1, 0,  3, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   3},
            { 0, 0, 0, 1, 0,  8,-13,  0, 0, 0, 0, 0, 0,   46,  66,   35, -25},
            { 0,-1, 1, 0, 0,  8,-12,  0, 0, 0, 0, 0, 0,  -14,   7,    0,   0},
            { 0, 2,-2, 1, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0},
            {-1, 0, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
            {-1, 0, 0, 1, 0, 18,-16,  0, 0, 0, 0, 0, 0,  -68, -34,  -18,  36},
            
            /* 131-140 */
            { 0, 1,-1, 1, 0,  0, -1,  0,-1, 1, 0, 0, 0,    0,  14,    7,   0},
            { 0, 0, 0, 1, 0,  3, -7,  4, 0, 0, 0, 0, 0,   10,  -6,   -3,  -5},
            {-2, 1, 1, 1, 0,  0, -3,  7, 0, 0, 0, 0, 0,   -5,  -4,   -2,   3},
            { 0, 1,-1, 2, 0,  0, -1,  0,-2, 5, 0, 0, 0,   -3,   5,    2,   1},
            { 0, 0, 0, 1, 0,  0,  0,  0,-2, 5, 0, 0, 0,   76,  17,    9, -41},
            { 0, 0, 0, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,   84, 298,  159, -45},
            { 1, 0, 0, 1, 0,-10,  3,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 0, 2,-2, 1, 0,  0, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
            {-1, 0, 0, 1, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
            { 0, 0, 0, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,  -82, 292,  156,  44},
            
            /* 141-150 */
            { 0, 0, 0, 1, 0,  0,  0,  0, 2,-5, 0, 0, 0,  -73,  17,    9,  39},
            { 0,-1, 1, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,   -9, -16,    0,   0},
            { 2,-1,-1, 1, 0,  0,  3, -7, 0, 0, 0, 0, 0,    3,   0,   -1,  -2},
            {-2, 0, 2, 0, 0,  0,  2,  0, 0,-5, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 1, 0, -3,  7, -4, 0, 0, 0, 0, 0,   -9,  -5,   -3,   5},
            {-2, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0, -439,   0,    0,   0},
            { 1, 0, 0, 1, 0,-18, 16,  0, 0, 0, 0, 0, 0,   57, -28,  -15, -30},
            {-2, 1, 1, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -6,   -3,   0},
            { 0, 1,-1, 2, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2},
            { 0, 0, 0, 1, 0, -8, 13,  0, 0, 0, 0, 0, 0,  -40,  57,   30,  21},
            
            /* 151-160 */
            { 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 1,   23,   7,    3, -13},
            { 0, 1,-1, 1, 0,  0,  0, -2, 0, 0, 0, 0, 0,  273,  80,   43,-146},
            { 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 0, -449, 430,    0,   0},
            { 0, 1,-1, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,   -8, -47,  -25,   4},
            { 0, 0, 0, 0, 0,  0, -1,  2, 0, 0, 0, 0, 1,    6,  47,   25,  -3},
            {-1, 0, 1, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  23,   13,   0},
            {-1, 0, 1, 1, 0,  0,  3, -4, 0, 0, 0, 0, 0,   -3,   0,    0,   2},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0,-2, 0, 0, 0,    3,  -4,   -2,  -2},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0, 2, 0, 0, 0,  -48,-110,  -59,  26},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 1,   51, 114,   61, -27},
            
            /* 161-170 */
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 2, -133,   0,    0,  57},
            { 0, 1,-1, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0},
            { 0, 0, 0, 1, 0, -3,  5,  0, 0, 0, 0, 0, 0,  -21,  -6,   -3,  11},
            { 0, 1,-1, 2, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  -3,   -1,   0},
            { 0, 0, 0, 1, 0,  0, -2,  4, 0, 0, 0, 0, 0,  -11, -21,  -11,   6},
            { 0, 2,-2, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,  -18,-436, -233,   9},
            { 0,-1, 1, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   35,  -7,    0,   0},
            { 0, 0, 0, 1, 0,  5, -8,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
            {-2, 0, 2, 1, 0,  6, -8,  0, 0, 0, 0, 0, 0,   11,  -3,   -1,  -6},
            { 0, 0, 0, 1, 0,  0, -8, 15, 0, 0, 0, 0, 0,   -5,  -3,   -1,   3},
            
            /* 171-180 */
            {-2, 0, 2, 1, 0,  0,  2,  0,-3, 0, 0, 0, 0,  -53,  -9,   -5,  28},
            {-2, 0, 2, 1, 0,  0,  6, -8, 0, 0, 0, 0, 0,    0,   3,    2,   1},
            { 1, 0,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    4,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  0,  0, 3,-5, 0, 0, 0,    0,  -4,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0,-1, 0, 0, 0, 0,  -50, 194,  103,  27},
            { 0, 0, 0, 0, 0,  0,  0,  0,-1, 0, 0, 0, 1,  -13,  52,   28,   7},
            { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 0,  -91, 248,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    6,  49,   26,  -3},
            { 0, 1,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   -6, -47,  -25,   3},
            { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    0,   5,    3,   0},
            
            /* 181-190 */
            { 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 2,   52,  23,   10, -23},
            { 0, 1,-1, 2, 0,  0, -1,  0, 0,-1, 0, 0, 0,   -3,   0,    0,   1},
            { 0, 0, 0, 1, 0,  0,  0,  0, 0,-1, 0, 0, 0,    0,   5,    3,   0},
            { 0,-1, 1, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,   -4,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0, -7, 13, 0, 0, 0, 0, 2,   -4,   8,    3,   2},
            { 0, 0, 0, 0, 0,  0,  7,-13, 0, 0, 0, 0, 0,   10,   0,    0,   0},
            { 2, 0,-2, 1, 0,  0, -5,  6, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
            { 0, 2,-2, 1, 0,  0, -8, 11, 0, 0, 0, 0, 0,    0,   8,    4,   0},
            { 0, 2,-2, 1,-1,  0,  2,  0, 0, 0, 0, 0, 0,    0,   8,    4,   1},
            {-2, 0, 2, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
            
            /* 191-200 */
            { 0, 0, 0, 0, 0,  0,  0,  0, 2,-2, 0, 0, 0,   -4,   0,    0,   0},
            { 0, 1,-1, 1, 0,  0, -1,  0, 0, 3, 0, 0, 0,   -8,   4,    2,   4},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 1,    8,  -4,   -2,  -4},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 2,    0,  15,    7,   0},
            {-2, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0, -138,   0,    0,   0},
            { 0, 0, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,  -7,   -3,   0},
            { 0, 0, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -7,   -3,   0},
            { 2, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   54,   0,    0, -29},
            { 0, 1,-1, 2, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,  10,    4,   0},
            { 0, 1,-1, 2, 0,  0,  0, -2, 0, 0, 0, 0, 0,   -7,   0,    0,   3},
            
            /* 201-210 */
            { 0, 0, 0, 1, 0,  0,  1, -2, 0, 0, 0, 0, 0,  -37,  35,   19,  20},
            { 0,-1, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,    0,   4,    0,   0},
            { 0,-1, 1, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -4,   9,    0,   0},
            { 0, 2,-2, 1, 0,  0, -2,  0, 0, 2, 0, 0, 0,    8,   0,    0,  -4},
            { 0, 1,-1, 1, 0,  3, -6,  0, 0, 0, 0, 0, 0,   -9, -14,   -8,   5},
            { 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 1,   -3,  -9,   -5,   3},
            { 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 0, -145,  47,    0,   0},
            { 0, 1,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,  -10,  40,   21,   5},
            { 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 1,   11, -49,  -26,  -7},
            { 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,-2150,   0,    0, 932},
            
            /* 211-220 */
            { 0, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,  -12,   0,    0,   5},
            { 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,   85,   0,    0, -37},
            { 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 1,    4,   0,    0,  -2},
            { 0, 1,-1, 1, 0,  0,  1, -4, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 0,  -86, 153,    0,   0},
            { 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -6,   9,    5,   3},
            { 0, 1,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    9, -13,   -7,  -5},
            { 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -8,  12,    6,   4},
            { 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 2,  -51,   0,    0,  22},
            { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,  -11,-268, -116,   5},
            
            /* 221-230 */
            { 0, 2,-2, 2, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,  12,    5,   0},
            { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,    0,   7,    3,   0},
            { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   31,   6,    3, -17},
            { 0, 1,-1, 1, 0, -5,  7,  0, 0, 0, 0, 0, 0,  140,  27,   14, -75},
            { 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   57,  11,    6, -30},
            { 0, 0, 0, 0, 0,  5, -8,  0, 0, 0, 0, 0, 0,  -14, -39,    0,   0},
            { 0, 1,-1, 2, 0,  0, -1,  0,-1, 0, 0, 0, 0,    0,  -6,   -2,   0},
            { 0, 0, 0, 1, 0,  0,  0,  0,-1, 0, 0, 0, 0,    4,  15,    8,  -2},
            { 0,-1, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    0,   4,    0,   0},
            { 0, 2,-2, 1, 0,  0, -2,  0, 1, 0, 0, 0, 0,   -3,   0,    0,   1},
            
            /* 231-240 */
            { 0, 0, 0, 0, 0,  0, -6, 11, 0, 0, 0, 0, 2,    0,  11,    5,   0},
            { 0, 0, 0, 0, 0,  0,  6,-11, 0, 0, 0, 0, 0,    9,   6,    0,   0},
            { 0, 0, 0, 0,-1,  0,  4,  0, 0, 0, 0, 0, 2,   -4,  10,    4,   2},
            { 0, 0, 0, 0, 1,  0, -4,  0, 0, 0, 0, 0, 0,    5,   3,    0,   0},
            { 2, 0,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,   16,   0,    0,  -9},
            {-2, 0, 2, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 2,-2, 1, 0,  0, -7,  9, 0, 0, 0, 0, 0,    0,   3,    2,  -1},
            { 0, 0, 0, 0, 0,  0,  0,  0, 4,-5, 0, 0, 2,    7,   0,    0,  -3},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 0,  -25,  22,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,   42, 223,  119, -22},
            
            /* 241-250 */
            { 0, 1,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -27,-143,  -77,  14},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,    9,  49,   26,  -5},
            { 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 2,-1166,   0,    0, 505},
            { 0, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -5,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 5, 0, 0, 2,   -6,   0,    0,   3},
            { 0, 0, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -8,   0,    1,   4},
            { 0,-1, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
            { 0, 2,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,  117,   0,    0, -63},
            { 0, 0, 0, 1, 0,  0,  2, -4, 0, 0, 0, 0, 0,   -4,   8,    4,   2},
            { 0, 2,-2, 1, 0,  0, -4,  4, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
            
            /* 251-260 */
            { 0, 1,-1, 2, 0, -5,  7,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  3, -6, 0, 0, 0, 0, 0,    0,  31,    0,   0},
            { 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -5,   0,    1,   3},
            { 0, 1,-1, 1, 0,  0, -4,  6, 0, 0, 0, 0, 0,    4,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 2,  -24, -13,   -6,  10},
            { 0,-1, 1, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0},
            { 0, 0, 0, 1, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0, -32,  -17,   0},
            { 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 2,    8,  12,    5,  -3},
            { 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 1,    3,   0,    0,  -1},
            
            /* 261-270 */
            { 0, 0, 0, 0, 0,  0,  5, -9, 0, 0, 0, 0, 0,    7,  13,    0,   0},
            { 0,-1, 1, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0,   -3,  16,    0,   0},
            { 0, 2,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   50,   0,    0, -27},
            {-2, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -5,   -3,   0},
            { 0,-2, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0},
            { 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 1,    0,   5,    3,   1},
            { 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 2,   24,   5,    2, -11},
            { 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 2,    5, -11,   -5,  -2},
            { 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 1,   30,  -3,   -2, -16},
            { 0, 1,-1, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   18,   0,    0,  -9},
            
            /* 271-280 */
            { 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 0,    8, 614,    0,   0},
            { 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 1,    3,  -3,   -1,  -2},
            { 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    6,  17,    9,  -3},
            { 0, 1,-1, 1, 0,  0, -1,  0, 3, 0, 0, 0, 0,   -3,  -9,   -5,   2},
            { 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    0,   6,    3,  -1},
            { 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 2, -127,  21,    9,  55},
            { 0, 0, 0, 0, 0,  0,  4, -8, 0, 0, 0, 0, 0,    3,   5,    0,   0},
            { 0, 0, 0, 0, 0,  0, -4,  8, 0, 0, 0, 0, 2,   -6, -10,   -4,   3},
            { 0,-2, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    5,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 2,   16,   9,    4,  -7},
            
            /* 281-290 */
            { 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 1,    3,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  4, -7, 0, 0, 0, 0, 0,    0,  22,    0,   0},
            { 0, 0, 0, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,  19,   10,   0},
            { 0, 2,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,    7,   0,    0,  -4},
            { 0, 0, 0, 0, 0,  0, -5, 10, 0, 0, 0, 0, 2,    0,  -5,   -2,   0},
            { 0, 0, 0, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0,  0,  0,  0, 4, 0, 0, 0, 2,   -9,   3,    1,   4},
            { 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 2,   17,   0,    0,  -7},
            { 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 1,    0,  -3,   -2,  -1},
            { 0, 0, 0, 0, 0,  0,  3, -5, 0, 0, 0, 0, 0,  -20,  34,    0,   0},
            
            /* 291-300 */
            { 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 1,  -10,   0,    1,   5},
            { 0, 1,-1, 1, 0,  1, -3,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 0,   22, -87,    0,   0},
            { 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 2,   -3,  -6,   -2,   1},
            { 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 2,  -16,  -3,   -1,   7},
            { 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0},
            { 0,-2, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  2, -3, 0, 0, 0, 0, 0,  -68,  39,    0,   0},
            { 0, 2,-2, 1, 0, -4,  4,  0, 0, 0, 0, 0, 0,   27,   0,    0, -14},
            
            /* 301-310 */
            { 0,-1, 1, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1, -1, 0, 0, 0, 0, 0,  -25,   0,    0,   0},
            { 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 1,  -12,  -3,   -2,   6},
            { 0, 1,-1, 1, 0, -4,  6,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 2,    3,  66,   29,  -1},
            { 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 2,  490,   0,    0,-213},
            { 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,  -22,  93,   49,  12},
            { 0, 1,-1, 1, 0, -4,  5,  0, 0, 0, 0, 0, 0,   -7,  28,   15,   4},
            { 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,   -3,  13,    7,   2},
            { 0, 0, 0, 0, 0,  4, -6,  0, 0, 0, 0, 0, 0,  -46,  14,    0,   0},
            
            /* 311-320 */
            {-2, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  1, 0, 0, 0, 0, 0,    2,   1,    0,   0},
            { 0,-1, 1, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
            { 0, 0, 0, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,  -28,   0,    0,  15},
            { 0, 0, 0, 0, 0,  0, -1,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  1, -3, 0, 0, 0, 0, 0,    0,   3,    0,   0},
            { 0, 0, 0, 0, 0,  0, -1,  3, 0, 0, 0, 0, 2,  -11,   0,    0,   5},
            { 0, 0, 0, 0, 0,  0, -7, 12, 0, 0, 0, 0, 2,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 1,   25, 106,   57, -13},
            
            /* 321-330 */
            { 0, 1,-1, 1, 0, -1,  0,  0, 0, 0, 0, 0, 0,    5,  21,   11,  -3},
            { 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0, 1485,   0,    0,   0},
            { 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 1,   -7, -32,  -17,   4},
            { 0, 1,-1, 1, 0,  1, -2,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
            { 0, 0, 0, 0, 0,  0, -2,  5, 0, 0, 0, 0, 2,   -6,  -3,   -2,   3},
            { 0, 0, 0, 0, 0,  0, -1,  0, 4, 0, 0, 0, 2,   30,  -6,   -2, -13},
            { 0, 0, 0, 0, 0,  0,  1,  0,-4, 0, 0, 0, 0,   -4,   4,    0,   0},
            { 0, 0, 0, 1, 0, -1,  1,  0, 0, 0, 0, 0, 0,  -19,   0,    0,  10},
            { 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 2,    0,   4,    2,  -1},
            { 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 0,    0,   3,    0,   0},
            
            /* 331-340 */
            { 0, 2,-2, 1, 0,  0, -3,  0, 3, 0, 0, 0, 0,    4,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0, -3,  7, 0, 0, 0, 0, 2,    0,  -3,   -1,   0},
            {-2, 0, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0, -5,  8, 0, 0, 0, 0, 2,    5,   3,    1,  -2},
            { 0, 0, 0, 0, 0,  0,  5, -8, 0, 0, 0, 0, 0,    0,  11,    0,   0},
            { 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 2,  118,   0,    0, -52},
            { 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 1,    0,  -5,   -3,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,  -28,  36,    0,   0},
            { 0, 0, 0, 0, 0,  2, -4,  0, 0, 0, 0, 0, 0,    5,  -5,    0,   0},
            { 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 1,   14, -59,  -31,  -8},
            
            /* 341-350 */
            { 0, 1,-1, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   9,    5,   1},
            { 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 2, -458,   0,    0, 198},
            { 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 2,    0, -45,  -20,   0},
            { 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 1,    9,   0,    0,  -5},
            { 0, 0, 0, 0, 0,  6, -9,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
            { 0, 0, 0, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -4,   -2,  -1},
            { 0, 2,-2, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   11,   0,    0,  -6},
            { 0, 0, 0, 0, 0,  0, -4,  6, 0, 0, 0, 0, 2,    6,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  4, -6, 0, 0, 0, 0, 0,  -16,  23,    0,   0},
            { 0, 0, 0, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,   -2,   0},
            
            /* 351-360 */
            { 0, 0, 0, 0, 0,  0, -1,  0, 2, 0, 0, 0, 2,   -5,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0, -166, 269,    0,   0},
            { 0, 0, 0, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   15,   0,    0,  -8},
            { 0, 0, 0, 0, 0, -5,  9,  0, 0, 0, 0, 0, 2,   10,   0,    0,  -4},
            { 0, 0, 0, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -78,  45,    0,   0},
            { 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 2,    0,  -5,   -2,   0},
            { 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 1,    7,   0,    0,  -4},
            { 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,   -5, 328,    0,   0},
            { 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2},
            { 0, 0, 0, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -2},
            
            /* 361-370 */
            { 0, 0, 0, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 0,-3, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 1,-5, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 1,    0,  -4,   -2,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,-1223, -26,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 1,    0,   7,    3,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0,-3, 5, 0, 0, 0,    3,   0,    0,   0},
            { 0, 0, 0, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -6,  20,    0,   0},
            { 0, 0, 0, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0, -368,   0,    0,   0},
            
            /* 371-380 */
            { 0, 0, 0, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,  -75,   0,    0,   0},
            { 0, 0, 0, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   11,   0,    0,  -6},
            { 0, 0, 0, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,    3,   0,    0,  -2},
            { 0, 0, 0, 0, 0, -8, 14,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 0, 0, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,  -13, -30,    0,   0},
            { 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 0,   21,   3,    0,   0},
            { 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 0, 0, 0, 0,  0, -1,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    8, -27,    0,   0},
            { 0, 0, 0, 0, 0,  0,  3, -8, 3, 0, 0, 0, 0,  -19, -11,    0,   0},
            
            /* 381-390 */
            { 0, 0, 0, 0, 0,  0, -3,  8,-3, 0, 0, 0, 2,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  1,  0,-2, 5, 0, 0, 2,    0,   5,    2,   0},
            { 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   2},
            { 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 1,-2, 0, 0, 0,   -1,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 0, 1, 0, 0, 2,  -14,   0,    0,   6},
            { 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,    6,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 2,  -74,   0,    0,  32},
            { 0, 0, 0, 0, 0,  0,  1,  0, 0, 2, 0, 0, 2,    0,  -3,   -1,   0},
            { 0, 2,-2, 1, 0, -5,  5,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2},
            
            /* 391-400 */
            { 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,    8,  11,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 1,    0,   3,    2,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 2, -262,   0,    0, 114},
            { 0, 0, 0, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
            { 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 1,   -7,   0,    0,   4},
            { 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 2,    0, -27,  -12,   0},
            { 0, 0, 0, 0, 0,  0, -1,  4, 0, 0, 0, 0, 2,  -19,  -8,   -4,   8},
            { 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 2,  202,   0,    0, -87},
            { 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 1,   -8,  35,   19,   5},
            { 0, 1,-1, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,   4,    2,   0},
            
            /* 401-410 */
            { 0, 0, 0, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   16,  -5,    0,   0},
            { 0, 2,-2, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    5,   0,    0,  -3},
            { 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,  -3,    0,   0},
            { 0, 0, 0, 0,-1,  0,  3,  0, 0, 0, 0, 0, 2,    1,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  0, 2, 0, 0, 0, 2,  -35, -48,  -21,  15},
            { 0, 0, 0, 0, 0,  0, -2,  6, 0, 0, 0, 0, 2,   -3,  -5,   -2,   1},
            { 0, 0, 0, 1, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3},
            { 0, 0, 0, 0, 0,  0, -6,  9, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  0,  6, -9, 0, 0, 0, 0, 0,    0,  -5,    0,   0},
            { 0, 0, 0, 0, 0, -2,  2,  0, 0, 0, 0, 0, 1,   12,  55,   29,  -6},
            
            /* 411-420 */
            { 0, 1,-1, 1, 0, -2,  1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
            { 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0, -598,   0,    0,   0},
            { 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 1,   -3, -13,   -7,   1},
            { 0, 0, 0, 0, 0,  0,  1,  0, 3, 0, 0, 0, 2,   -5,  -7,   -3,   2},
            { 0, 0, 0, 0, 0,  0, -5,  7, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  0,  5, -7, 0, 0, 0, 0, 0,    5,  -7,    0,   0},
            { 0, 0, 0, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  4, -5, 0, 0, 0, 0, 0,   16,  -6,    0,   0},
            { 0, 0, 0, 0, 0,  1, -3,  0, 0, 0, 0, 0, 0,    8,  -3,    0,   0},
            { 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 1,    8, -31,  -16,  -4},
            
            /* 421-430 */
            { 0, 1,-1, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 2,  113,   0,    0, -49},
            { 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 2,    0, -24,  -10,   0},
            { 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 1,    4,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  3, -3, 0, 0, 0, 0, 0,   27,   0,    0,   0},
            { 0, 0, 0, 0, 0, -4,  8,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 2,    0,  -4,   -2,   0},
            { 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 1,    5,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0},
            { 0, 0, 0, 0, 0,  0,  1,  1, 0, 0, 0, 0, 2,  -13,   0,    0,   6},
            
            /* 431-440 */
            { 0, 0, 0, 0, 0,  0, -2,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  0,  0,  3, 0, 0, 0, 0, 2,  -18, -10,   -4,   8},
            { 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,   -4, -28,    0,   0},
            { 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 2,   -5,   6,    3,   2},
            { 0, 0, 0, 0, 0, -9, 13,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 0, 0, 0, 0,  0, -1,  5, 0, 0, 0, 0, 2,   -5,  -9,   -4,   2},
            { 0, 0, 0, 0, 0,  0, -2,  0, 4, 0, 0, 0, 2,   17,   0,    0,  -7},
            { 0, 0, 0, 0, 0,  0,  2,  0,-4, 0, 0, 0, 0,   11,   4,    0,   0},
            { 0, 0, 0, 0, 0,  0, -2,  7, 0, 0, 0, 0, 2,    0,  -6,   -2,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   83,  15,    0,   0},
            
            /* 441-450 */
            { 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 2,    0,-114,  -49,   0},
            { 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 2,  117,   0,    0, -51},
            { 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 1,   -5,  19,   10,   2},
            { 0, 0, 0, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0, -3,  9, 0, 0, 0, 0, 2,    0,  -3,   -1,   0},
            { 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,    3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 2,    0,  -6,   -2,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  393,   3,    0,   0},
            
            /* 451-460 */
            { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 1,   -4,  21,   11,   2},
            { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 2,   -6,   0,   -1,   3},
            { 0, 0, 0, 0, 0, -5, 10,  0, 0, 0, 0, 0, 2,   -3,   8,    4,   1},
            { 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,    8,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 2,   18, -29,  -13,  -8},
            { 0, 0, 0, 0, 0, -3,  3,  0, 0, 0, 0, 0, 1,    8,  34,   18,  -4},
            { 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   89,   0,    0,   0},
            { 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 1,    3,  12,    6,  -1},
            { 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 2,   54, -15,   -7, -24},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0,-3, 0, 0, 0,    0,   3,    0,   0},
            
            /* 461-470 */
            { 0, 0, 0, 0, 0,  0, -5, 13, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 0,    0,  35,    0,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 2, -154, -30,  -13,  67},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   15,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 1,    0,   4,    2,   0},
            { 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 0,    0,   9,    0,   0},
            { 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 2,   80, -71,  -31, -35},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0,-1, 0, 0, 2,    0, -20,   -9,   0},
            { 0, 0, 0, 0, 0,  0, -6, 15, 0, 0, 0, 0, 2,   11,   5,    2,  -5},
            { 0, 0, 0, 0, 0, -8, 15,  0, 0, 0, 0, 0, 2,   61, -96,  -42, -27},
            
            /* 471-480 */
            { 0, 0, 0, 0, 0, -3,  9, -4, 0, 0, 0, 0, 2,   14,   9,    4,  -6},
            { 0, 0, 0, 0, 0,  0,  2,  0, 2,-5, 0, 0, 2,  -11,  -6,   -3,   5},
            { 0, 0, 0, 0, 0,  0, -2,  8,-1,-5, 0, 0, 2,    0,  -3,   -1,   0},
            { 0, 0, 0, 0, 0,  0,  6, -8, 3, 0, 0, 0, 2,  123,-415, -180, -53},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,    0,   0,    0, -35},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    7, -32,  -17,  -4},
            { 0, 1,-1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -9,   -5,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    0,  -4,    2,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 2,  -89,   0,    0,  38},
            
            /* 481-490 */
            { 0, 0, 0, 0, 0,  0, -6, 16,-4,-5, 0, 0, 2,    0, -86,  -19,  -6},
            { 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2,    0,   0,  -19,   6},
            { 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2, -123,-416, -180,  53},
            { 0, 0, 0, 0, 0,  0,  6, -8, 1, 5, 0, 0, 2,    0,  -3,   -1,   0},
            { 0, 0, 0, 0, 0,  0,  2,  0,-2, 5, 0, 0, 2,   12,  -6,   -3,  -5},
            { 0, 0, 0, 0, 0,  3, -5,  4, 0, 0, 0, 0, 2,  -13,   9,    4,   6},
            { 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,    0, -15,   -7,   0},
            { 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,  -62, -97,  -42,  27},
            { 0, 0, 0, 0, 0,  0, 11,  0, 0, 0, 0, 0, 2,  -11,   5,    2,   5},
            
            /* 491-500 */
            { 0, 0, 0, 0, 0,  0,  2,  0, 0, 1, 0, 0, 2,    0, -19,   -8,   0},
            { 0, 0, 0, 0, 0,  3, -3,  0, 2, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 2,-2, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   4,    2,   0},
            { 0, 1,-1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0},
            { 0, 2,-2, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   4,    2,   0},
            { 0, 0, 0, 0, 0,  0,  1,  2, 0, 0, 0, 0, 2,  -85, -70,  -31,  37},
            { 0, 0, 0, 0, 0,  0,  2,  0, 1, 0, 0, 0, 2,  163, -12,   -5, -72},
            { 0, 0, 0, 0, 0, -3,  7,  0, 0, 0, 0, 0, 2,  -63, -16,   -7,  28},
            { 0, 0, 0, 0, 0,  0,  0,  4, 0, 0, 0, 0, 2,  -21, -32,  -14,   9},
            { 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 2,    0,  -3,   -1,   0},
            
            /* 501-510 */
            { 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2},
            { 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   8,    0,   0},
            { 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 2,    3,  10,    4,  -1},
            { 0, 0, 0, 0, 0,  0,  2,  0, 2, 0, 0, 0, 2,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  0, -1,  6, 0, 0, 0, 0, 2,    0,  -7,   -3,   0},
            { 0, 0, 0, 0, 0,  0,  7, -9, 0, 0, 0, 0, 2,    0,  -4,   -2,   0},
            { 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 0,    6,  19,    0,   0},
            { 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 2,    5,-173,  -75,  -2},
            { 0, 0, 0, 0, 0,  0,  6, -7, 0, 0, 0, 0, 2,    0,  -7,   -3,   0},
            { 0, 0, 0, 0, 0,  0,  5, -5, 0, 0, 0, 0, 2,    7, -12,   -5,  -3},
            
            /* 511-520 */
            { 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 1,   -3,   0,    0,   2},
            { 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 2,    3,  -4,   -2,  -1},
            { 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 2,   74,   0,    0, -32},
            { 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 1,   -3,  12,    6,   2},
            { 0, 0, 0, 0, 0,  0,  4, -3, 0, 0, 0, 0, 2,   26, -14,   -6, -11},
            { 0, 0, 0, 0, 0,  0,  3, -1, 0, 0, 0, 0, 2,   19,   0,    0,  -8},
            { 0, 0, 0, 0, 0, -4,  4,  0, 0, 0, 0, 0, 1,    6,  24,   13,  -3},
            { 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   83,   0,    0,   0},
            { 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 1,    0, -10,   -5,   0},
            { 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 2,   11,  -3,   -1,  -5},
            
            /* 521-530 */
            { 0, 0, 0, 0, 0,  0,  2,  1, 0, 0, 0, 0, 2,    3,   0,    1,  -1},
            { 0, 0, 0, 0, 0,  0, -3,  0, 5, 0, 0, 0, 2,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
            { 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 1,    5, -23,  -12,  -3},
            { 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 2, -339,   0,    0, 147},
            { 0, 0, 0, 0, 0, -9, 12,  0, 0, 0, 0, 0, 2,    0, -10,   -5,   0},
            { 0, 0, 0, 0, 0,  0,  3,  0,-4, 0, 0, 0, 0,    5,   0,    0,   0},
            { 0, 2,-2, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  0,  7, -8, 0, 0, 0, 0, 2,    0,  -4,   -2,   0},
            { 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 0,   18,  -3,    0,   0},
            
            /* 531-540 */
            { 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 2,    9, -11,   -5,  -4},
            { 0, 0, 0, 0, 0, -2,  6,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4},
            { 0, 0, 0, 0, 0, -6,  7,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1},
            { 0, 0, 0, 0, 0,  6, -7,  0, 0, 0, 0, 0, 0,    0,   9,    0,   0},
            { 0, 0, 0, 0, 0,  0,  6, -6, 0, 0, 0, 0, 2,    6,  -9,   -4,  -2},
            { 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 0,   -4, -12,    0,   0},
            { 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 2,   67, -91,  -39, -29},
            { 0, 0, 0, 0, 0,  0,  5, -4, 0, 0, 0, 0, 2,   30, -18,   -8, -13},
            { 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 0,    0,   0,    0,   0},
            { 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 2,    0,-114,  -50,   0},
            
            /* 541-550 */
            { 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,    0,   0,    0,  23},
            { 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,  517,  16,    7,-224},
            { 0, 0, 0, 0, 0,  0,  3,  0, 0,-2, 0, 0, 2,    0,  -7,   -3,   0},
            { 0, 0, 0, 0, 0,  0,  4, -2, 0, 0, 0, 0, 2,  143,  -3,   -1, -62},
            { 0, 0, 0, 0, 0,  0,  3,  0, 0,-1, 0, 0, 2,   29,   0,    0, -13},
            { 0, 2,-2, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -4,   0,    0,   2},
            { 0, 0, 0, 0, 0, -8, 16,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   3},
            { 0, 0, 0, 0, 0,  0,  3,  0, 2,-5, 0, 0, 2,    5,  12,    5,  -2},
            { 0, 0, 0, 0, 0,  0,  7, -8, 3, 0, 0, 0, 2,  -25,   0,    0,  11},
            { 0, 0, 0, 0, 0,  0, -5, 16,-4,-5, 0, 0, 2,   -3,   0,    0,   1},
            
            /* 551-560 */
            { 0, 0, 0, 0, 0,  0,  3,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0},
            { 0, 0, 0, 0, 0,  0, -1,  8,-3, 0, 0, 0, 2,  -22,  12,    5,  10},
            { 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,   50,   0,    0, -22},
            { 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0},
            { 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0,  0,  2,  2, 0, 0, 0, 0, 2,   -4,   4,    2,   2},
            { 0, 0, 0, 0, 0,  0,  3,  0, 1, 0, 0, 0, 2,   -5, -11,   -5,   2},
            { 0, 0, 0, 0, 0, -3,  8,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0},
            { 0, 0, 0, 0, 0, -5,  5,  0, 0, 0, 0, 0, 1,    4,  17,    9,  -2},
            { 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 0,   59,   0,    0,   0},
            
            /* 561-570 */
            { 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 1,    0,  -4,   -2,   0},
            { 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4},
            { 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 1,    4, -15,   -8,  -2},
            { 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 2,  370,  -8,    0,-160},
            { 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   0,   -3,   0},
            { 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0,  0,  6, -5, 0, 0, 0, 0, 2,   -6,   3,    1,   3},
            { 0, 0, 0, 0, 0,  7, -8,  0, 0, 0, 0, 0, 0,    0,   6,    0,   0},
            { 0, 0, 0, 0, 0,  0,  5, -3, 0, 0, 0, 0, 2,  -10,   0,    0,   4},
            
            /* 571-580 */
            { 0, 0, 0, 0, 0,  4, -3,  0, 0, 0, 0, 0, 2,    0,   9,    4,   0},
            { 0, 0, 0, 0, 0,  1,  2,  0, 0, 0, 0, 0, 2,    4,  17,    7,  -2},
            { 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 2,   34,   0,    0, -15},
            { 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0},
            { 0, 0, 0, 0, 0,  0,  4,  0,-4, 0, 0, 0, 2,   -5,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  4,  0,-3, 0, 0, 0, 2,  -37,  -7,   -3,  16},
            { 0, 0, 0, 0, 0, -6,  6,  0, 0, 0, 0, 0, 1,    3,  13,    7,  -2},
            { 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 0,   40,   0,    0,   0},
            { 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0},
            { 0, 0, 0, 0, 0,  0,  4,  0,-2, 0, 0, 0, 2, -184,  -3,   -1,  80},
            
            /* 581-590 */
            { 0, 0, 0, 0, 0,  0,  6, -4, 0, 0, 0, 0, 2,   -3,   0,    0,   1},
            { 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 1,    0, -10,   -6,  -1},
            { 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 2,   31,  -6,    0, -13},
            { 0, 0, 0, 0, 0,  0,  4,  0,-1, 0, 0, 0, 2,   -3, -32,  -14,   1},
            { 0, 0, 0, 0, 0,  0,  4,  0, 0,-2, 0, 0, 2,   -7,   0,    0,   3},
            { 0, 0, 0, 0, 0,  0,  5, -2, 0, 0, 0, 0, 2,    0,  -8,   -4,   0},
            { 0, 0, 0, 0, 0,  0,  4,  0, 0, 0, 0, 0, 0,    3,  -4,    0,   0},
            { 0, 0, 0, 0, 0,  8, -9,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0},
            { 0, 0, 0, 0, 0,  5, -4,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0},
            
            /* 591-600 */
            { 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 2,   19, -23,  -10,   2},
            { 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   0,    0, -10},
            { 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   3,    2,   0},
            { 0, 0, 0, 0, 0, -7,  7,  0, 0, 0, 0, 0, 1,    0,   9,    5,  -1},
            { 0, 0, 0, 0, 0,  7, -7,  0, 0, 0, 0, 0, 0,   28,   0,    0,   0},
            { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 1,    0,  -7,   -4,   0},
            { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 2,    8,  -4,    0,  -4},
            { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   0,   -2,   0},
            { 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0},
            { 0, 0, 0, 0, 0,  0,  5,  0,-4, 0, 0, 0, 2,   -3,   0,    0,   1},
            
            /* 601-610 */
            { 0, 0, 0, 0, 0,  0,  5,  0,-3, 0, 0, 0, 2,   -9,   0,    1,   4},
            { 0, 0, 0, 0, 0,  0,  5,  0,-2, 0, 0, 0, 2,    3,  12,    5,  -1},
            { 0, 0, 0, 0, 0,  3,  0,  0, 0, 0, 0, 0, 2,   17,  -3,   -1,   0},
            { 0, 0, 0, 0, 0, -8,  8,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0},
            { 0, 0, 0, 0, 0,  8, -8,  0, 0, 0, 0, 0, 0,   19,   0,    0,   0},
            { 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 1,    0,  -5,   -3,   0},
            { 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 2,   14,  -3,    0,  -1},
            { 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,   -1,   0},
            { 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,    0,  -5},
            { 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0},
            
            /* 611-620 */
            { 0, 0, 0, 0, 0,  9, -9,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0},
            { 0, 0, 0, 0, 0,  6, -4,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    2,   9,    4,   3},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    0,   0,    0,  -4},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    8,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   4,    2,   0},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    6,   0,    0,  -3},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   3,    1,   0},
            { 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    5,   0,    0,  -2},
            
            /* 621-630 */
            { 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 2,    3,   0,    0,  -1},
            { 1, 0,-2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 1, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0},
            { 1, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    7,   0,    0,   0},
            { 1, 0,-2, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
            {-1, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0},
            {-1, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    6,   0,    0,   0},
            {-1, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0},
            { 1, 0,-2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0},
            {-2, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    5,   0,    0,   0},
            
            /* 631-640 */
            {-1, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -3,   0,    0,   0},
            {-1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0},
            {-1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0},
            {-1, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0},
            { 1,-1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0},
            {-1, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   13,   0,    0,   0},
            {-2, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   21,  11,    0,   0},
            { 1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0},
            {-1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,  -5,   -2,   0},
            { 1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0},
            
            /* 641-650 */
            {-1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0},
            {-1, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2},
            { 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   20,  10,    0,   0},
            {-1, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  -34,   0,    0,   0},
            {-1, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,  -19,   0,    0,   0},
            { 1, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,    3,   0,    0,  -2},
            { 1, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
            { 1, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -6,   0,    0,   3},
            { 1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0},
            { 1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    3,   0,    0,   0},
            
            /* 651-660 */
            { 0, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0},
            { 0, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0},
            { 0, 2, 0, 2, 0, -2,  2,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 0, 2, 0, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    6,   0,    0,  -3},
            { 0, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   3},
            { 0, 2, 0, 2, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
            { 0, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0},
            { 0, 1, 1, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -3,   -2,   0},
            { 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  126, -63,  -27, -55},
            {-1, 2, 0, 2, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    1,   2},
            
            /* 661-670 */
            { 0, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,   -3,  28,   15,   2},
            { 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    5,   0,    1,  -2},
            { 0, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   9,    4,   1},
            { 0, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   9,    4,  -1},
            {-1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0, -126, -63,  -27,  55},
            { 2, 2,-2, 2, 0,  0, -2,  0, 3, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 1, 2, 0, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   21, -11,   -6, -11},
            { 0, 1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0},
            {-1, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -21, -11,   -6,  11},
            {-2, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   1},
            
            /* 671-680 */
            { 0, 2, 0, 2, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0},
            { 0, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    8,   0,    0,  -4},
            { 0, 2, 0, 2, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -6,   0,    0,   3},
            { 0, 2, 0, 2, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
            {-1, 2, 2, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 1, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1},
            {-1, 2, 2, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -5,   0,    0,   2},
            { 2, 2, 0, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   24, -12,   -5, -11},
            { 1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   3,    1,   0},
            { 1, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   3,    1,   0},
            
            /* 681-687 */
            { 1, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0},
            { 0, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -24, -12,   -5,  10},
            { 2, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    4,   0,   -1,  -2},
            {-1, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   13,   0,    0,  -6},
            {-1, 2, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    7,   0,    0,  -3},
            { 1, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1},
            { 0, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,    3,   0,    0,  -1}
        };
        
        /* Number of terms in the planetary nutation model */
        const int NPL = (int) (sizeof xpl / sizeof xpl[0]);
        
        /*--------------------------------------------------------------------*/
        
        /* Interval between fundamental date J2000.0 and given date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* ------------------- */
        /* LUNI-SOLAR NUTATION */
        /* ------------------- */
        
        /* Fundamental (Delaunay) arguments */
        
        /* Mean anomaly of the Moon (IERS 2003). */
        el = iauFal03(t);
        
        /* Mean anomaly of the Sun (MHB2000). */
        elp = fmod(1287104.79305  +
                   t * (129596581.0481  +
                        t * (-0.5532  +
                             t * (0.000136  +
                                  t * (-0.00001149)))), TURNAS) * DAS2R;
        
        /* Mean longitude of the Moon minus that of the ascending node */
        /* (IERS 2003. */
        f = iauFaf03(t);
        
        /* Mean elongation of the Moon from the Sun (MHB2000). */
        d = fmod(1072260.70369  +
                 t * (1602961601.2090  +
                      t * (-6.3706  +
                           t * (0.006593  +
                                t * (-0.00003169)))), TURNAS) * DAS2R;
        
        /* Mean longitude of the ascending node of the Moon (IERS 2003). */
        om = iauFaom03(t);
        
        /* Initialize the nutation values. */
        dp = 0.0;
        de = 0.0;
        
        /* Summation of luni-solar nutation series (in reverse order). */
        for (i = NLS-1; i >= 0; i--) {
            
            /* Argument and functions. */
            arg = fmod((double)xls[i].nl  * el +
                       (double)xls[i].nlp * elp +
                       (double)xls[i].nf  * f +
                       (double)xls[i].nd  * d +
                       (double)xls[i].nom * om, D2PI);
            sarg = sin(arg);
            carg = cos(arg);
            
            /* Term. */
            dp += (xls[i].sp + xls[i].spt * t) * sarg + xls[i].cp * carg;
            de += (xls[i].ce + xls[i].cet * t) * carg + xls[i].se * sarg;
        }
        
        /* Convert from 0.1 microarcsec units to radians. */
        dpsils = dp * U2R;
        depsls = de * U2R;
        
        /* ------------------ */
        /* PLANETARY NUTATION */
        /* ------------------ */
        
        /* n.b.  The MHB2000 code computes the luni-solar and planetary nutation */
        /* in different functions, using slightly different Delaunay */
        /* arguments in the two cases.  This behaviour is faithfully */
        /* reproduced here.  Use of the IERS 2003 expressions for both */
        /* cases leads to negligible changes, well below */
        /* 0.1 microarcsecond. */
        
        /* Mean anomaly of the Moon (MHB2000). */
        al = fmod(2.35555598 + 8328.6914269554 * t, D2PI);
        
        /* Mean longitude of the Moon minus that of the ascending node */
        /*(MHB2000). */
        af = fmod(1.627905234 + 8433.466158131 * t, D2PI);
        
        /* Mean elongation of the Moon from the Sun (MHB2000). */
        ad = fmod(5.198466741 + 7771.3771468121 * t, D2PI);
        
        /* Mean longitude of the ascending node of the Moon (MHB2000). */
        aom = fmod(2.18243920 - 33.757045 * t, D2PI);
        
        /* General accumulated precession in longitude (IERS 2003). */
        apa = iauFapa03(t);
        
        /* Planetary longitudes, Mercury through Uranus (IERS 2003). */
        alme = iauFame03(t);
        alve = iauFave03(t);
        alea = iauFae03(t);
        alma = iauFama03(t);
        alju = iauFaju03(t);
        alsa = iauFasa03(t);
        alur = iauFaur03(t);
        
        /* Neptune longitude (MHB2000). */
        alne = fmod(5.321159000 + 3.8127774000 * t, D2PI);
        
        /* Initialize the nutation values. */
        dp = 0.0;
        de = 0.0;
        
        /* Summation of planetary nutation series (in reverse order). */
        for (i = NPL-1; i >= 0; i--) {
            
            /* Argument and functions. */
            arg = fmod((double)xpl[i].nl  * al   +
                       (double)xpl[i].nf  * af   +
                       (double)xpl[i].nd  * ad   +
                       (double)xpl[i].nom * aom  +
                       (double)xpl[i].nme * alme +
                       (double)xpl[i].nve * alve +
                       (double)xpl[i].nea * alea +
                       (double)xpl[i].nma * alma +
                       (double)xpl[i].nju * alju +
                       (double)xpl[i].nsa * alsa +
                       (double)xpl[i].nur * alur +
                       (double)xpl[i].nne * alne +
                       (double)xpl[i].npa * apa, D2PI);
            sarg = sin(arg);
            carg = cos(arg);
            
            /* Term. */
            dp += (double)xpl[i].sp * sarg + (double)xpl[i].cp * carg;
            de += (double)xpl[i].se * sarg + (double)xpl[i].ce * carg;
            
        }
        
        /* Convert from 0.1 microarcsec units to radians. */
        dpsipl = dp * U2R;
        depspl = de * U2R;
        
        /* ------- */
        /* RESULTS */
        /* ------- */
        
        /* Add luni-solar and planetary components. */
        *dpsi = dpsils + dpsipl;
        *deps = depsls + depspl;
        
        return;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - -
     **   i a u F a l 0 3
     **  - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean anomaly of the Moon.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    l, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     is from Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFal03(double t)
    {
        double a;
        
        /* Mean anomaly of the Moon (IERS Conventions 2003). */
        a = fmod(           485868.249036  +
                 t * ( 1717915923.2178 +
                      t * (         31.8792 +
                           t * (          0.051635 +
                                t * (        - 0.00024470 ) ) ) ), TURNAS ) * DAS2R;
        
        return a;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - -
     **   i a u F a f 0 3
     **  - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of the Moon minus mean longitude of the ascending
     **  node.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    F, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     is from Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFaf03(double t)
    {
        double a;
        
        /* Mean longitude of the Moon minus that of the ascending node */
        /* (IERS Conventions 2003).                                    */
        a = fmod(           335779.526232 +
                 t * ( 1739527262.8478 +
                      t * (       - 12.7512 +
                           t * (        - 0.001037 +
                                t * (          0.00000417 ) ) ) ), TURNAS ) * DAS2R;
        
        return a;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u F a o m 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of the Moon's ascending node.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    Omega, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     is from Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFaom03(double t)
    {
        double a;
        
        /* Mean longitude of the Moon's ascending node */
        /* (IERS Conventions 2003).                    */
        a = fmod(          450160.398036 +
                 t * ( - 6962890.5431 +
                      t * (         7.4722 +
                           t * (         0.007702 +
                                t * (       - 0.00005939 ) ) ) ), TURNAS ) * DAS2R;
        
        return a;
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u F a p a 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  general accumulated precession in longitude.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    general precession in longitude, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003).  It
     **     is taken from Kinoshita & Souchay (1990) and comes originally
     **     from Lieske et al. (1977).
     **
     **  References:
     **
     **     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
     **     48, 187
     **
     **     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
     **     Astron.Astrophys. 58, 1-16
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFapa03(double t)
    {
        double a;
        
        /* General accumulated precession in longitude. */
        a = (0.024381750 + 0.00000538691 * t) * t;
        
        return a;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u F a v e 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Venus.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    mean longitude of Venus, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     comes from Souchay et al. (1999) after Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFave03(double t)
    {
        double a;
        
        /* Mean longitude of Venus (IERS Conventions 2003). */
        a = fmod(3.176146697 + 1021.3285546211 * t, D2PI);
        
        return a;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - -
     **   i a u F a e 0 3
     **  - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Earth.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    mean longitude of Earth, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     comes from Souchay et al. (1999) after Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFae03(double t)
    {
        double a;
        
        /* Mean longitude of Earth (IERS Conventions 2003). */
        a = fmod(1.753470314 + 628.3075849991 * t, D2PI);
        
        return a;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u F a m a 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Mars.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    mean longitude of Mars, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     comes from Souchay et al. (1999) after Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFama03(double t)
    {
        double a;
        
        /* Mean longitude of Mars (IERS Conventions 2003). */
        a = fmod(6.203480913 + 334.0612426700 * t, D2PI);
        
        return a;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u F a j u 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Jupiter.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    mean longitude of Jupiter, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     comes from Souchay et al. (1999) after Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFaju03(double t)
    {
        double a;
        
        /* Mean longitude of Jupiter (IERS Conventions 2003). */
        a = fmod(0.599546497 + 52.9690962641 * t, D2PI);
        
        return a;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u F a s a 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Saturn.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    mean longitude of Saturn, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     comes from Souchay et al. (1999) after Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFasa03(double t)
    {
        
        double a;
        
        /* Mean longitude of Saturn (IERS Conventions 2003). */
        a = fmod(0.874016757 + 21.3299104960 * t, D2PI);
        
        return a;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u F a u r 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Uranus.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned  (function value):
     **           double    mean longitude of Uranus, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     is adapted from Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFaur03(double t)
    {
        double a;
        
        /* Mean longitude of Uranus (IERS Conventions 2003). */
        a = fmod(5.481293872 + 7.4781598567 * t, D2PI);
        
        return a;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u F a m e 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean longitude of Mercury.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    mean longitude of Mercury, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     comes from Souchay et al. (1999) after Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFame03(double t)
    {
        double a;
        
        /* Mean longitude of Mercury (IERS Conventions 2003). */
        a = fmod(4.402608842 + 2608.7903141574 * t, D2PI);
        
        return a;
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u N u t 0 6 a
     **  - - - - - - - - - -
     **
     **  IAU 2000A nutation with adjustments to match the IAU 2006
     **  precession.
     **
     **  Given:
     **     date1,date2   double   TT as a 2-part Julian Date (Note 1)
     **
     **  Returned:
     **     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
     **
     **  Status:  canonical model.
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The nutation components in longitude and obliquity are in radians
     **     and with respect to the mean equinox and ecliptic of date,
     **     IAU 2006 precession model (Hilton et al. 2006, Capitaine et al.
     **     2005).
     **
     **  3) The function first computes the IAU 2000A nutation, then applies
     **     adjustments for (i) the consequences of the change in obliquity
     **     from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
     **     secular variation in the Earth's dynamical form factor J2.
     **
     **  4) The present function provides classical nutation, complementing
     **     the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
     **     pole which is at current epochs accurate to a few tens of
     **     microarcseconds, apart from the free core nutation.
     **
     **  Called:
     **     iauNut00a    nutation, IAU 2000A
     **
     **  References:
     **
     **     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
     **     Astron.Astrophys. 387, 700
     **
     **     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
     **     Astron.Astrophys. 58, 1-16
     **
     **     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
     **     107, B4.  The MHB_2000 code itself was obtained on 9th September
     **     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     **     Astron.Astrophys.Supp.Ser. 135, 111
     **
     **     Wallace, P.T., "Software for Implementing the IAU 2000
     **     Resolutions", in IERS Workshop 5.1 (2002)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauNut06a(double date1, double date2, double *dpsi, double *deps)
    {
        double t, fj2, dp, de;
        
        /* Interval between fundamental date J2000.0 and given date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* Factor correcting for secular variation of J2. */
        fj2 = -2.7774e-6 * t;
        
        /* Obtain IAU 2000A nutation. */
        iauNut00a(date1, date2, &dp, &de);
        
        /* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). */
        *dpsi = dp + dp * (0.4697e-6 + fj2);
        *deps = de + de * fj2;
        
        return;
        
    } // end of function
    
    
    /*
     **  - - - - - - - -
     **   i a u F w 2 m
     **  - - - - - - - -
     **
     **  Form rotation matrix given the Fukushima-Williams angles.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     gamb     double         F-W angle gamma_bar (radians)
     **     phib     double         F-W angle phi_bar (radians)
     **     psi      double         F-W angle psi (radians)
     **     eps      double         F-W angle epsilon (radians)
     **
     **  Returned:
     **     r        double[3][3]   rotation matrix
     **
     **  Notes:
     **
     **  1) Naming the following points:
     **
     **           e = J2000.0 ecliptic pole,
     **           p = GCRS pole,
     **           E = ecliptic pole of date,
     **     and   P = CIP,
     **
     **     the four Fukushima-Williams angles are as follows:
     **
     **        gamb = gamma = epE
     **        phib = phi = pE
     **        psi = psi = pEP
     **        eps = epsilon = EP
     **
     **  2) The matrix representing the combined effects of frame bias,
     **     precession and nutation is:
     **
     **        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)
     **
     **  3) Three different matrices can be constructed, depending on the
     **     supplied angles:
     **
     **     o  To obtain the nutation x precession x frame bias matrix,
     **        generate the four precession angles, generate the nutation
     **        components and add them to the psi_bar and epsilon_A angles,
     **        and call the present function.
     **
     **     o  To obtain the precession x frame bias matrix, generate the
     **        four precession angles and call the present function.
     **
     **     o  To obtain the frame bias matrix, generate the four precession
     **        angles for date J2000.0 and call the present function.
     **
     **     The nutation-only and precession-only matrices can if necessary
     **     be obtained by combining these three appropriately.
     **
     **  Called:
     **     iauIr        initialize r-matrix to identity
     **     iauRz        rotate around Z-axis
     **     iauRx        rotate around X-axis
     **
     **  Reference:
     **
     **     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauFw2m(double gamb, double phib, double psi, double eps, double (*r)[3])
    {
        /* Construct the matrix. */
        iauIr(r);
        iauRz(gamb, r);
        iauRx(phib, r);
        iauRz(-psi, r);
        iauRx(-eps, r);
        
        return;
    } // end of function
    
    void GIAU::iauIr(double (*r)[3])
    {
        r[0][0] = 1.0;
        r[0][1] = 0.0;
        r[0][2] = 0.0;
        r[1][0] = 0.0;
        r[1][1] = 1.0;
        r[1][2] = 0.0;
        r[2][0] = 0.0;
        r[2][1] = 0.0;
        r[2][2] = 1.0;
        
        return;
    }
    
    void GIAU::iauRx(double phi, double (*r)[3])
    {
        double s, c, a10, a11, a12, a20, a21, a22;
        
        s = sin(phi);
        c = cos(phi);
        
        a10 =   c*r[1][0] + s*r[2][0];
        a11 =   c*r[1][1] + s*r[2][1];
        a12 =   c*r[1][2] + s*r[2][2];
        a20 = - s*r[1][0] + c*r[2][0];
        a21 = - s*r[1][1] + c*r[2][1];
        a22 = - s*r[1][2] + c*r[2][2];
        
        r[1][0] = a10;
        r[1][1] = a11;
        r[1][2] = a12;
        r[2][0] = a20;
        r[2][1] = a21;
        r[2][2] = a22;
        
        return;
    }
    
    void GIAU::iauRy(double theta, double (*r)[3])
    {
        double s, c, a00, a01, a02, a20, a21, a22;
        
        s = sin(theta);
        c = cos(theta);
        
        a00 = c*r[0][0] - s*r[2][0];
        a01 = c*r[0][1] - s*r[2][1];
        a02 = c*r[0][2] - s*r[2][2];
        a20 = s*r[0][0] + c*r[2][0];
        a21 = s*r[0][1] + c*r[2][1];
        a22 = s*r[0][2] + c*r[2][2];
        
        r[0][0] = a00;
        r[0][1] = a01;
        r[0][2] = a02;
        r[2][0] = a20;
        r[2][1] = a21;
        r[2][2] = a22;
    }
    
    void GIAU::iauRz(double psi, double (*r)[3])
    {
        double s, c, a00, a01, a02, a10, a11, a12;
        
        s = sin(psi);
        c = cos(psi);
        
        a00 =   c*r[0][0] + s*r[1][0];
        a01 =   c*r[0][1] + s*r[1][1];
        a02 =   c*r[0][2] + s*r[1][2];
        a10 = - s*r[0][0] + c*r[1][0];
        a11 = - s*r[0][1] + c*r[1][1];
        a12 = - s*r[0][2] + c*r[1][2];
        
        r[0][0] = a00;
        r[0][1] = a01;
        r[0][2] = a02;
        r[1][0] = a10;
        r[1][1] = a11;
        r[1][2] = a12;
        
        return;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u P n m 0 6 a
     **  - - - - - - - - - -
     **
     **  Form the matrix of precession-nutation for a given date (including
     **  frame bias), IAU 2006 precession and IAU 2000A nutation models.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     date1,date2 double       TT as a 2-part Julian Date (Note 1)
     **
     **  Returned:
     **     rnpb        double[3][3] bias-precession-nutation matrix (Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The matrix operates in the sense V(date) = rnpb * V(GCRS), where
     **     the p-vector V(date) is with respect to the true equatorial triad
     **     of date date1+date2 and the p-vector V(GCRS) is with respect to
     **     the Geocentric Celestial Reference System (IAU, 2000).
     **
     **  Called:
     **     iauPfw06     bias-precession F-W angles, IAU 2006
     **     iauNut06a    nutation, IAU 2006/2000A
     **     iauFw2m      F-W angles to r-matrix
     **
     **  Reference:
     **
     **     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauPnm06a(double date1, double date2, double (*rnpb)[3])
    {
        double gamb, phib, psib, epsa, dp, de;
        
        /* Fukushima-Williams angles for frame bias and precession. */
        iauPfw06(date1, date2, &gamb, &phib, &psib, &epsa);
        
        /* Nutation components. */
        iauNut06a(date1, date2, &dp, &de);
        
        /* Equinox based nutation x precession x bias matrix. */
        iauFw2m(gamb, phib, psib + dp, epsa + de, rnpb);
        
        return;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u B p n 2 x y
     **  - - - - - - - - - -
     **
     **  Extract from the bias-precession-nutation matrix the X,Y coordinates
     **  of the Celestial Intermediate Pole.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     rbpn      double[3][3]  celestial-to-true matrix (Note 1)
     **
     **  Returned:
     **     x,y       double        Celestial Intermediate Pole (Note 2)
     **
     **  Notes:
     **
     **  1) The matrix rbpn transforms vectors from GCRS to true equator (and
     **     CIO or equinox) of date, and therefore the Celestial Intermediate
     **     Pole unit vector is the bottom row of the matrix.
     **
     **  2) The arguments x,y are components of the Celestial Intermediate
     **     Pole unit vector in the Geocentric Celestial Reference System.
     **
     **  Reference:
     **
     **     "Expressions for the Celestial Intermediate Pole and Celestial
     **     Ephemeris Origin consistent with the IAU 2000A precession-
     **     nutation model", Astron.Astrophys. 400, 1145-1154
     **     (2003)
     **
     **     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     **          intermediate origin" (CIO) by IAU 2006 Resolution 2.
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauBpn2xy(double (*rbpn)[3], double *x, double *y)
    {
        /* Extract the X,Y coordinates. */
        *x = rbpn[2][0];
        *y = rbpn[2][1];
        
        return;
    }
    
    
    /*
     **  - - - - - - -
     **   i a u S 0 0
     **  - - - - - - -
     **
     **  The CIO locator s, positioning the Celestial Intermediate Origin on
     **  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
     **  coordinates.  Compatible with IAU 2000A precession-nutation.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     date1,date2   double    TT as a 2-part Julian Date (Note 1)
     **     x,y           double    CIP coordinates (Note 3)
     **
     **  Returned (function value):
     **                   double    the CIO locator s in radians (Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The CIO locator s is the difference between the right ascensions
     **     of the same point in two systems:  the two systems are the GCRS
     **     and the CIP,CIO, and the point is the ascending node of the
     **     CIP equator.  The quantity s remains below 0.1 arcsecond
     **     throughout 1900-2100.
     **
     **  3) The series used to compute s is in fact for s+XY/2, where X and Y
     **     are the x and y components of the CIP unit vector;  this series
     **     is more compact than a direct series for s would be.  This
     **     function requires X,Y to be supplied by the caller, who is
     **     responsible for providing values that are consistent with the
     **     supplied date.
     **
     **  4) The model is consistent with the IAU 2000A precession-nutation.
     **
     **  Called:
     **     iauFal03     mean anomaly of the Moon
     **     iauFalp03    mean anomaly of the Sun
     **     iauFaf03     mean argument of the latitude of the Moon
     **     iauFad03     mean elongation of the Moon from the Sun
     **     iauFaom03    mean longitude of the Moon's ascending node
     **     iauFave03    mean longitude of Venus
     **     iauFae03     mean longitude of Earth
     **     iauFapa03    general accumulated precession in longitude
     **
     **  References:
     **
     **     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
     **     "Expressions for the Celestial Intermediate Pole and Celestial
     **     Ephemeris Origin consistent with the IAU 2000A precession-
     **     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
     **
     **     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     **          intermediate origin" (CIO) by IAU 2006 Resolution 2.
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauS00(double date1, double date2, double x, double y)
    {
        /* Time since J2000.0, in Julian centuries */
        double t;
        
        /* Miscellaneous */
        int i, j;
        double a, w0, w1, w2, w3, w4, w5;
        
        /* Fundamental arguments */
        double fa[8];
        
        /* Returned value */
        double s;
        
        /* --------------------- */
        /* The series for s+XY/2 */
        /* --------------------- */
        
        typedef struct {
            int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
            double s, c;     /* sine and cosine coefficients */
        } TERM;
        
        /* Polynomial coefficients */
        static const double sp[] = {
            
            /* 1-6 */
            94.00e-6,
            3808.35e-6,
            -119.94e-6,
            -72574.09e-6,
            27.70e-6,
            15.61e-6
        };
        
        /* Terms of order t^0 */
        static const TERM s0[] = {
            
            /* 1-10 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 },
            {{ 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 },
            {{ 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 },
            {{ 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 },
            {{ 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 },
            {{ 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 },
            
            /* 11-20 */
            {{ 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 },
            {{ 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 },
            {{ 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 },
            {{ 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 },
            {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 },
            
            /* 21-30 */
            {{ 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 },
            {{ 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 },
            {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 },
            {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 },
            {{ 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
            {{ 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
            {{ 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 },
            
            /* 31-33 */
            {{ 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 }
        };
        
        /* Terms of order t^1 */
        static const TERM s1[] ={
            
            /* 1-3 */
            {{ 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 },
            {{ 0,  0,  0,  0,  1,  0,  0,  0},     1.71e-6,  -0.03e-6 },
            {{ 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 }
        };
        
        /* Terms of order t^2 */
        static const TERM s2[] ={
            
            /* 1-10 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},   743.53e-6,  -0.17e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 },
            {{ 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 },
            {{ 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 },
            
            /* 11-20 */
            {{ 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 },
            {{ 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 },
            {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 },
            {{ 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 },
            {{ 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 },
            {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 },
            
            /* 21-25 */
            {{ 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 },
            {{ 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 },
            {{ 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 },
            {{ 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 }
        };
        
        /* Terms of order t^3 */
        static const TERM s3[] ={
            
            /* 1-4 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.51e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.39e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.24e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.22e-6 }
        };
        
        /* Terms of order t^4 */
        static const TERM s4[] ={
            
            /* 1-1 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 }
        };
        
        /* Number of terms in the series */
        const int NS0 = (int) (sizeof s0 / sizeof (TERM));
        const int NS1 = (int) (sizeof s1 / sizeof (TERM));
        const int NS2 = (int) (sizeof s2 / sizeof (TERM));
        const int NS3 = (int) (sizeof s3 / sizeof (TERM));
        const int NS4 = (int) (sizeof s4 / sizeof (TERM));
        
        /*--------------------------------------------------------------------*/
        
        /* Interval between fundamental epoch J2000.0 and current date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* Fundamental Arguments (from IERS Conventions 2003) */
        
        /* Mean anomaly of the Moon. */
        fa[0] = iauFal03(t);
        
        /* Mean anomaly of the Sun. */
        fa[1] = iauFalp03(t);
        
        /* Mean longitude of the Moon minus that of the ascending node. */
        fa[2] = iauFaf03(t);
        
        /* Mean elongation of the Moon from the Sun. */
        fa[3] = iauFad03(t);
        
        /* Mean longitude of the ascending node of the Moon. */
        fa[4] = iauFaom03(t);
        
        /* Mean longitude of Venus. */
        fa[5] = iauFave03(t);
        
        /* Mean longitude of Earth. */
        fa[6] = iauFae03(t);
        
        /* General precession in longitude. */
        fa[7] = iauFapa03(t);
        
        /* Evaluate s. */
        w0 = sp[0];
        w1 = sp[1];
        w2 = sp[2];
        w3 = sp[3];
        w4 = sp[4];
        w5 = sp[5];
        
        for (i = NS0-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s0[i].nfa[j] * fa[j];
            }
            w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
        }
        
        for (i = NS1-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s1[i].nfa[j] * fa[j];
            }
            w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
        }
        
        for (i = NS2-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s2[i].nfa[j] * fa[j];
            }
            w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
        }
        
        for (i = NS3-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s3[i].nfa[j] * fa[j];
            }
            w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
        }
        
        for (i = NS4-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s4[i].nfa[j] * fa[j];
            }
            w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
        }
        
        s = (w0 +
             (w1 +
              (w2 +
               (w3 +
                (w4 +
                 w5 * t) * t) * t) * t) * t) * DAS2R - x*y/2.0;
        
        return s;
        
    }  // end of function
    
    
    
    /*
     **  - - - - - - -
     **   i a u S 0 6
     **  - - - - - - -
     **
     **  The CIO locator s, positioning the Celestial Intermediate Origin on
     **  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
     **  coordinates.  Compatible with IAU 2006/2000A precession-nutation.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     date1,date2   double    TT as a 2-part Julian Date (Note 1)
     **     x,y           double    CIP coordinates (Note 3)
     **
     **  Returned (function value):
     **                   double    the CIO locator s in radians (Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The CIO locator s is the difference between the right ascensions
     **     of the same point in two systems:  the two systems are the GCRS
     **     and the CIP,CIO, and the point is the ascending node of the
     **     CIP equator.  The quantity s remains below 0.1 arcsecond
     **     throughout 1900-2100.
     **
     **  3) The series used to compute s is in fact for s+XY/2, where X and Y
     **     are the x and y components of the CIP unit vector;  this series
     **     is more compact than a direct series for s would be.  This
     **     function requires X,Y to be supplied by the caller, who is
     **     responsible for providing values that are consistent with the
     **     supplied date.
     **
     **  4) The model is consistent with the "P03" precession (Capitaine et
     **     al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
     **     IAU 2000A nutation (with P03 adjustments).
     **
     **  Called:
     **     iauFal03     mean anomaly of the Moon
     **     iauFalp03    mean anomaly of the Sun
     **     iauFaf03     mean argument of the latitude of the Moon
     **     iauFad03     mean elongation of the Moon from the Sun
     **     iauFaom03    mean longitude of the Moon's ascending node
     **     iauFave03    mean longitude of Venus
     **     iauFae03     mean longitude of Earth
     **     iauFapa03    general accumulated precession in longitude
     **
     **  References:
     **
     **     Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
     **     Astrophys. 432, 355
     **
     **     McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauS06(double date1, double date2, double x, double y)
    {
        /* Time since J2000.0, in Julian centuries */
        double t;
        
        /* Miscellaneous */
        int i, j;
        double a, w0, w1, w2, w3, w4, w5;
        
        /* Fundamental arguments */
        double fa[8];
        
        /* Returned value */
        double s;
        
        /* --------------------- */
        /* The series for s+XY/2 */
        /* --------------------- */
        
        typedef struct {
            int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
            double s, c;     /* sine and cosine coefficients */
        } TERM;
        
        /* Polynomial coefficients */
        static const double sp[] = {
            
            /* 1-6 */
            94.00e-6,
            3808.65e-6,
            -122.68e-6,
            -72574.11e-6,
            27.98e-6,
            15.62e-6
        };
        
        /* Terms of order t^0 */
        static const TERM s0[] = {
            
            /* 1-10 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 },
            {{ 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 },
            {{ 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 },
            {{ 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 },
            {{ 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 },
            {{ 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 },
            
            /* 11-20 */
            {{ 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 },
            {{ 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 },
            {{ 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 },
            {{ 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 },
            {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 },
            
            /* 21-30 */
            {{ 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 },
            {{ 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 },
            {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 },
            {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 },
            {{ 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
            {{ 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
            {{ 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 },
            
            /* 31-33 */
            {{ 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 }
        };
        
        /* Terms of order t^1 */
        static const TERM s1[] = {
            
            /* 1 - 3 */
            {{ 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 },
            {{ 0,  0,  0,  0,  1,  0,  0,  0},     1.73e-6,  -0.03e-6 },
            {{ 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 }
        };
        
        /* Terms of order t^2 */
        static const TERM s2[] = {
            
            /* 1-10 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},   743.52e-6,  -0.17e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 },
            {{ 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 },
            {{ 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 },
            
            /* 11-20 */
            {{ 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 },
            {{ 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 },
            {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 },
            {{ 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 },
            {{ 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 },
            {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 },
            
            /* 21-25 */
            {{ 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 },
            {{ 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 },
            {{ 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 },
            {{ 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 }
        };
        
        /* Terms of order t^3 */
        static const TERM s3[] = {
            
            /* 1-4 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.42e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.46e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.25e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.23e-6 }
        };
        
        /* Terms of order t^4 */
        static const TERM s4[] = {
            
            /* 1-1 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 }
        };
        
        /* Number of terms in the series */
        static const int NS0 = (int) (sizeof s0 / sizeof (TERM));
        static const int NS1 = (int) (sizeof s1 / sizeof (TERM));
        static const int NS2 = (int) (sizeof s2 / sizeof (TERM));
        static const int NS3 = (int) (sizeof s3 / sizeof (TERM));
        static const int NS4 = (int) (sizeof s4 / sizeof (TERM));
        
        /*--------------------------------------------------------------------*/
        
        /* Interval between fundamental epoch J2000.0 and current date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* Fundamental Arguments (from IERS Conventions 2003) */
        
        /* Mean anomaly of the Moon. */
        fa[0] = iauFal03(t);
        
        /* Mean anomaly of the Sun. */
        fa[1] = iauFalp03(t);
        
        /* Mean longitude of the Moon minus that of the ascending node. */
        fa[2] = iauFaf03(t);
        
        /* Mean elongation of the Moon from the Sun. */
        fa[3] = iauFad03(t);
        
        /* Mean longitude of the ascending node of the Moon. */
        fa[4] = iauFaom03(t);
        
        /* Mean longitude of Venus. */
        fa[5] = iauFave03(t);
        
        /* Mean longitude of Earth. */
        fa[6] = iauFae03(t);
        
        /* General precession in longitude. */
        fa[7] = iauFapa03(t);
        
        /* Evaluate s. */
        w0 = sp[0];
        w1 = sp[1];
        w2 = sp[2];
        w3 = sp[3];
        w4 = sp[4];
        w5 = sp[5];
        
        for (i = NS0-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s0[i].nfa[j] * fa[j];
            }
            w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
        }
        
        for (i = NS1-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s1[i].nfa[j] * fa[j];
            }
            w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
        }
        
        for (i = NS2-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s2[i].nfa[j] * fa[j];
            }
            w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
        }
        
        for (i = NS3-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s3[i].nfa[j] * fa[j];
            }
            w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
        }
        
        for (i = NS4-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s4[i].nfa[j] * fa[j];
            }
            w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
        }
        
        s = (w0 +
             (w1 +
              (w2 +
               (w3 +
                (w4 +
                 w5 * t) * t) * t) * t) * t) * DAS2R - x*y/2.0;
        
        return s;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u F a l p 0 3
     **  - - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean anomaly of the Sun.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    l', radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     is from Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFalp03(double t)
    {
        double a;
        
        /* Mean anomaly of the Sun (IERS Conventions 2003). */
        a = fmod(         1287104.793048 +
                 t * ( 129596581.0481 +
                      t * (       - 0.5532 +
                           t * (         0.000136 +
                                t * (       - 0.00001149 ) ) ) ), TURNAS ) * DAS2R;
        
        return a;
    } // end of function
  
    /*
     **  - - - - - - - - -
     **   i a u F a d 0 3
     **  - - - - - - - - -
     **
     **  Fundamental argument, IERS Conventions (2003):
     **  mean elongation of the Moon from the Sun.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     t     double    TDB, Julian centuries since J2000.0 (Note 1)
     **
     **  Returned (function value):
     **           double    D, radians (Note 2)
     **
     **  Notes:
     **
     **  1) Though t is strictly TDB, it is usually more convenient to use
     **     TT, which makes no significant difference.
     **
     **  2) The expression used is as adopted in IERS Conventions (2003) and
     **     is from Simon et al. (1994).
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     **     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauFad03(double t)
    {
        double a;
        
        /* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
        a = fmod(          1072260.703692 +
                 t * ( 1602961601.2090 +
                      t * (        - 6.3706 +
                           t * (          0.006593 +
                                t * (        - 0.00003169 ) ) ) ), TURNAS ) * DAS2R;
        
        return a;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u C 2 i 0 6 a
     **  - - - - - - - - - -
     **
     **  Form the celestial-to-intermediate matrix for a given date using the
     **  IAU 2006 precession and IAU 2000A nutation models.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     date1,date2 double       TT as a 2-part Julian Date (Note 1)
     **
     **  Returned:
     **     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The matrix rc2i is the first stage in the transformation from
     **     celestial to terrestrial coordinates:
     **
     **        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
     **
     **               =  RC2T * [CRS]
     **
     **     where [CRS] is a vector in the Geocentric Celestial Reference
     **     System and [TRS] is a vector in the International Terrestrial
     **     Reference System (see IERS Conventions 2003), ERA is the Earth
     **     Rotation Angle and RPOM is the polar motion matrix.
     **
     **  Called:
     **     iauPnm06a    classical NPB matrix, IAU 2006/2000A
     **     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
     **     iauS06       the CIO locator s, given X,Y, IAU 2006
     **     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
     **
     **  References:
     **
     **     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauC2i06a(double date1, double date2, double (*rc2i)[3])
    {
        double rbpn[3][3], x, y, s;
        
        /* Obtain the celestial-to-true matrix (IAU 2006/2000A). */
        iauPnm06a(date1, date2, rbpn);
        
        /* Extract the X,Y coordinates. */
        iauBpn2xy(rbpn, &x, &y);
        
        /* Obtain the CIO locator. */
        s = iauS06(date1, date2, x, y);
        
        /* Form the celestial-to-intermediate matrix. */
        iauC2ixys(x, y, s, rc2i);
        
        return;
    } // end of function
    
    
    /*
     **  - - - - - - - - -
     **   i a u C 2 i x y
     **  - - - - - - - - -
     **
     **  Form the celestial to intermediate-frame-of-date matrix for a given
     **  date when the CIP X,Y coordinates are known.  IAU 2000.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     date1,date2 double       TT as a 2-part Julian Date (Note 1)
     **     x,y         double       Celestial Intermediate Pole (Note 2)
     **
     **  Returned:
     **     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The Celestial Intermediate Pole coordinates are the x,y components
     **     of the unit vector in the Geocentric Celestial Reference System.
     **
     **  3) The matrix rc2i is the first stage in the transformation from
     **     celestial to terrestrial coordinates:
     **
     **        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
     **
     **              = RC2T * [CRS]
     **
     **     where [CRS] is a vector in the Geocentric Celestial Reference
     **     System and [TRS] is a vector in the International Terrestrial
     **     Reference System (see IERS Conventions 2003), ERA is the Earth
     **     Rotation Angle and RPOM is the polar motion matrix.
     **
     **  4) Although its name does not include "00", This function is in fact
     **     specific to the IAU 2000 models.
     **
     **  Called:
     **     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
     **     iauS00       the CIO locator s, given X,Y, IAU 2000A
     **
     **  Reference:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauC2ixy(double date1, double date2, double x, double y, double (*rc2i)[3])
    {
        /* Compute s and then the matrix. */
        iauC2ixys(x, y, iauS00(date1, date2, x, y), rc2i);
        
        return;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u C 2 i x y s
     **  - - - - - - - - - -
     **
     **  Form the celestial to intermediate-frame-of-date matrix given the CIP
     **  X,Y and the CIO locator s.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     x,y      double         Celestial Intermediate Pole (Note 1)
     **     s        double         the CIO locator s (Note 2)
     **
     **  Returned:
     **     rc2i     double[3][3]   celestial-to-intermediate matrix (Note 3)
     **
     **  Notes:
     **
     **  1) The Celestial Intermediate Pole coordinates are the x,y
     **     components of the unit vector in the Geocentric Celestial
     **     Reference System.
     **
     **  2) The CIO locator s (in radians) positions the Celestial
     **     Intermediate Origin on the equator of the CIP.
     **
     **  3) The matrix rc2i is the first stage in the transformation from
     **     celestial to terrestrial coordinates:
     **
     **        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
     **
     **              = RC2T * [CRS]
     **
     **     where [CRS] is a vector in the Geocentric Celestial Reference
     **     System and [TRS] is a vector in the International Terrestrial
     **     Reference System (see IERS Conventions 2003), ERA is the Earth
     **     Rotation Angle and RPOM is the polar motion matrix.
     **
     **  Called:
     **     iauIr        initialize r-matrix to identity
     **     iauRz        rotate around Z-axis
     **     iauRy        rotate around Y-axis
     **
     **  Reference:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2014 November 7
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauC2ixys(double x, double y, double s, double (*rc2i)[3])
    {
        double r2, e, d;
        
        
        /* Obtain the spherical angles E and d. */
        r2 = x*x + y*y;
        e = (r2 > 0.0) ? atan2(y, x) : 0.0;
        d = atan(sqrt(r2 / (1.0 - r2)));
        
        /* Form the matrix. */
        iauIr(rc2i);
        iauRz(e, rc2i);
        iauRy(d, rc2i);
        iauRz(-(e+s), rc2i);
        
        return;
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u C 2 t 0 6 a
     **  - - - - - - - - - -
     **
     **  Form the celestial to terrestrial matrix given the date, the UT1 and
     **  the polar motion, using the IAU 2006 precession and IAU 2000A
     **  nutation models.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
     **     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
     **     xp,yp    double         coordinates of the pole (radians, Note 2)
     **
     **  Returned:
     **     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
     **
     **  Notes:
     **
     **  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
     **     apportioned in any convenient way between the arguments uta and
     **     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
     **     these ways, among others:
     **
     **             uta            utb
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution is
     **     acceptable.  The J2000 and MJD methods are good compromises
     **     between resolution and convenience.  In the case of uta,utb, the
     **     date & time method is best matched to the Earth rotation angle
     **     algorithm used:  maximum precision is delivered when the uta
     **     argument is for 0hrs UT1 on the day in question and the utb
     **     argument lies in the range 0 to 1, or vice versa.
     **
     **  2) The arguments xp and yp are the coordinates (in radians) of the
     **     Celestial Intermediate Pole with respect to the International
     **     Terrestrial Reference System (see IERS Conventions 2003),
     **     measured along the meridians to 0 and 90 deg west respectively.
     **
     **  3) The matrix rc2t transforms from celestial to terrestrial
     **     coordinates:
     **
     **        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
     **
     **              = rc2t * [CRS]
     **
     **     where [CRS] is a vector in the Geocentric Celestial Reference
     **     System and [TRS] is a vector in the International Terrestrial
     **     Reference System (see IERS Conventions 2003), RC2I is the
     **     celestial-to-intermediate matrix, ERA is the Earth rotation
     **     angle and RPOM is the polar motion matrix.
     **
     **  Called:
     **     iauC2i06a    celestial-to-intermediate matrix, IAU 2006/2000A
     **     iauEra00     Earth rotation angle, IAU 2000
     **     iauSp00      the TIO locator s', IERS 2000
     **     iauPom00     polar motion matrix
     **     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
     **
     **  Reference:
     **
     **     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauC2t06a(double tta, double ttb, double uta, double utb, double xp, double yp, double (*rc2t)[3])
    {
        double rc2i[3][3], era, sp, rpom[3][3];
        
        /* Form the celestial-to-intermediate matrix for this TT. */
        iauC2i06a(tta, ttb, rc2i);
        
        /* Predict the Earth rotation angle for this UT1. */
        era = iauEra00(uta, utb);
        
        /* Estimate s'. */
        sp = iauSp00(tta, ttb);
        
        /* Form the polar motion matrix. */
        iauPom00(xp, yp, sp, rpom);
        
        /* Combine to form the celestial-to-terrestrial matrix. */
        iauC2tcio(rc2i, era, rpom, rc2t);
        
        return;
        
    } // end of function
    
    
    /*
     **  - - - - - - - - - -
     **   i a u C 2 t c i o
     **  - - - - - - - - - -
     **
     **  Assemble the celestial to terrestrial matrix from CIO-based
     **  components (the celestial-to-intermediate matrix, the Earth Rotation
     **  Angle and the polar motion matrix).
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     rc2i     double[3][3]    celestial-to-intermediate matrix
     **     era      double          Earth rotation angle (radians)
     **     rpom     double[3][3]    polar-motion matrix
     **
     **  Returned:
     **     rc2t     double[3][3]    celestial-to-terrestrial matrix
     **
     **  Notes:
     **
     **  1) This function constructs the rotation matrix that transforms
     **     vectors in the celestial system into vectors in the terrestrial
     **     system.  It does so starting from precomputed components, namely
     **     the matrix which rotates from celestial coordinates to the
     **     intermediate frame, the Earth rotation angle and the polar motion
     **     matrix.  One use of the present function is when generating a
     **     series of celestial-to-terrestrial matrices where only the Earth
     **     Rotation Angle changes, avoiding the considerable overhead of
     **     recomputing the precession-nutation more often than necessary to
     **     achieve given accuracy objectives.
     **
     **  2) The relationship between the arguments is as follows:
     **
     **        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
     **
     **              = rc2t * [CRS]
     **
     **     where [CRS] is a vector in the Geocentric Celestial Reference
     **     System and [TRS] is a vector in the International Terrestrial
     **     Reference System (see IERS Conventions 2003).
     **
     **  Called:
     **     iauCr        copy r-matrix
     **     iauRz        rotate around Z-axis
     **     iauRxr       product of two r-matrices
     **
     **  Reference:
     **
     **     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG
     **
     **  This revision:  2013 August 24
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauC2tcio(double (*rc2i)[3], double era, double (*rpom)[3], double (*rc2t)[3])
    {
        double r[3][3];
        
        /* Construct the matrix. */
        iauCr(rc2i, r);
        iauRz(era, r);
        iauRxr(rpom, r, rc2t);
        
        return;
        
    } // end of function
    
    void GIAU::iauCp(double *p, double *c)
    {
        c[0] = p[0];
        c[1] = p[1];
        c[2] = p[2];
        
        return;
        
    } // end of function
    void GIAU::iauCr(double (*r)[3], double (*c)[3])
    {
        iauCp(r[0], c[0]);
        iauCp(r[1], c[1]);
        iauCp(r[2], c[2]);
        
        return;
    } // end of function
    
    
    /*
     **  - - - - - - -
     **   i a u R x r
     **  - - - - - - -
     **
     **  Multiply two r-matrices.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  vector/matrix support function.
     **
     **  Given:
     **     a        double[3][3]    first r-matrix
     **     b        double[3][3]    second r-matrix
     **
     **  Returned:
     **     atb      double[3][3]    a * b
     **
     **  Note:
     **     It is permissible to re-use the same array for any of the
     **     arguments.
     **
     **  Called:
     **     iauCr        copy r-matrix
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauRxr(double (*a)[3], double (*b)[3], double (*atb)[3])
    {
        int i, j, k;
        double w, wm[3][3];
        
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                w = 0.0;
                for (k = 0; k < 3; k++) {
                    w +=  a[i][k] * b[k][j];
                }
                wm[i][j] = w;
            }
        }
        iauCr(wm, atb);
        
        return;
        
    } // end of function
    
    /*
     **  - - - - - - - - - -
     **   i a u P o m 0 0
     **  - - - - - - - - - -
     **
     **  Form the matrix of polar motion for a given date, IAU 2000.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  support function.
     **
     **  Given:
     **     xp,yp    double    coordinates of the pole (radians, Note 1)
     **     sp       double    the TIO locator s' (radians, Note 2)
     **
     **  Returned:
     **     rpom     double[3][3]   polar-motion matrix (Note 3)
     **
     **  Notes:
     **
     **  1) The arguments xp and yp are the coordinates (in radians) of the
     **     Celestial Intermediate Pole with respect to the International
     **     Terrestrial Reference System (see IERS Conventions 2003),
     **     measured along the meridians to 0 and 90 deg west respectively.
     **
     **  2) The argument sp is the TIO locator s', in radians, which
     **     positions the Terrestrial Intermediate Origin on the equator.  It
     **     is obtained from polar motion observations by numerical
     **     integration, and so is in essence unpredictable.  However, it is
     **     dominated by a secular drift of about 47 microarcseconds per
     **     century, and so can be taken into account by using s' = -47*t,
     **     where t is centuries since J2000.0.  The function iauSp00
     **     implements this approximation.
     **
     **  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
     **     that it is the final rotation when computing the pointing
     **     direction to a celestial source.
     **
     **  Called:
     **     iauIr        initialize r-matrix to identity
     **     iauRz        rotate around Z-axis
     **     iauRy        rotate around Y-axis
     **     iauRx        rotate around X-axis
     **
     **  Reference:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    void GIAU::iauPom00(double xp, double yp, double sp, double (*rpom)[3])
    {
        /* Construct the matrix. */
        iauIr(rpom);
        iauRz(sp, rpom);
        iauRy(-xp, rpom);
        iauRx(-yp, rpom);
        
        return;
        
    } // end of function
    
    /*
     **  - - - - - - - -
     **   i a u S p 0 0
     **  - - - - - - - -
     **
     **  The TIO locator s', positioning the Terrestrial Intermediate Origin
     **  on the equator of the Celestial Intermediate Pole.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     date1,date2  double    TT as a 2-part Julian Date (Note 1)
     **
     **  Returned (function value):
     **                  double    the TIO locator s' in radians (Note 2)
     **
     **  Notes:
     **
     **  1) The TT date date1+date2 is a Julian Date, apportioned in any
     **     convenient way between the two arguments.  For example,
     **     JD(TT)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **            date1          date2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 method is best matched to the way
     **     the argument is handled internally and will deliver the
     **     optimum resolution.  The MJD method and the date & time methods
     **     are both good compromises between resolution and convenience.
     **
     **  2) The TIO locator s' is obtained from polar motion observations by
     **     numerical integration, and so is in essence unpredictable.
     **     However, it is dominated by a secular drift of about
     **     47 microarcseconds per century, which is the approximation
     **     evaluated by the present function.
     **
     **  Reference:
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauSp00(double date1, double date2)
    {
        double t, sp;
        
        /* Interval between fundamental epoch J2000.0 and current date (JC). */
        t = ((date1 - DJ00) + date2) / DJC;
        
        /* Approximate s'. */
        sp = -47e-6 * t * DAS2R;
        
        return sp;
        
    } // end of function
    
    /*
     **  - - - - - - - - -
     **   i a u E r a 0 0
     **  - - - - - - - - -
     **
     **  Earth rotation angle (IAU 2000 model).
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  canonical model.
     **
     **  Given:
     **     dj1,dj2   double    UT1 as a 2-part Julian Date (see note)
     **
     **  Returned (function value):
     **               double    Earth rotation angle (radians), range 0-2pi
     **
     **  Notes:
     **
     **  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
     **     convenient way between the arguments dj1 and dj2.  For example,
     **     JD(UT1)=2450123.7 could be expressed in any of these ways,
     **     among others:
     **
     **             dj1            dj2
     **
     **         2450123.7           0.0       (JD method)
     **         2451545.0       -1421.3       (J2000 method)
     **         2400000.5       50123.2       (MJD method)
     **         2450123.5           0.2       (date & time method)
     **
     **     The JD method is the most natural and convenient to use in
     **     cases where the loss of several decimal digits of resolution
     **     is acceptable.  The J2000 and MJD methods are good compromises
     **     between resolution and convenience.  The date & time method is
     **     best matched to the algorithm used:  maximum precision is
     **     delivered when the dj1 argument is for 0hrs UT1 on the day in
     **     question and the dj2 argument lies in the range 0 to 1, or vice
     **     versa.
     **
     **  2) The algorithm is adapted from Expression 22 of Capitaine et al.
     **     2000.  The time argument has been expressed in days directly,
     **     and, to retain precision, integer contributions have been
     **     eliminated.  The same formulation is given in IERS Conventions
     **     (2003), Chap. 5, Eq. 14.
     **
     **  Called:
     **     iauAnp       normalize angle into range 0 to 2pi
     **
     **  References:
     **
     **     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
     **     Astrophys., 355, 398-405.
     **
     **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     **     IERS Technical Note No. 32, BKG (2004)
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauEra00(double dj1, double dj2)
    {
        double d1, d2, t, f, theta;
        
        /* Days since fundamental epoch. */
        if (dj1 < dj2) {
            d1 = dj1;
            d2 = dj2;
        } else {
            d1 = dj2;
            d2 = dj1;
        }
        t = d1 + (d2- DJ00);
        
        /* Fractional part of T (days). */
        f = fmod(d1, 1.0) + fmod(d2, 1.0);
        
        /* Earth rotation angle at this UT1. */
        theta = iauAnp(D2PI * (f + 0.7790572732640
                               + 0.00273781191135448 * t));
        
        return theta;
        
    } // end of function
    
    
    /*
     **  - - - - - - -
     **   i a u A n p
     **  - - - - - - -
     **
     **  Normalize angle into the range 0 <= a < 2pi.
     **
     **  This function is part of the International Astronomical Union's
     **  SOFA (Standards Of Fundamental Astronomy) software collection.
     **
     **  Status:  vector/matrix support function.
     **
     **  Given:
     **     a        double     angle (radians)
     **
     **  Returned (function value):
     **              double     angle in range 0-2pi
     **
     **  This revision:  2013 June 18
     **
     **  SOFA release 2015-02-09
     **
     **  Copyright (C) 2015 IAU SOFA Board.  See notes at end.
     */
    double GIAU::iauAnp(double a)
    {
        double w;
        
        w = fmod(a, D2PI);
        if (w < 0) w += D2PI;
        
        return w;
    }
    
    
    
} // end of namespace

