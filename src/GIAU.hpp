//
//  GIAU.hpp
//  GFC
//
//  Created by lizhen on 09/09/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#ifndef GIAU_hpp
#define GIAU_hpp

#include <stdio.h>
#include <math.h>

namespace gfc
{
    // this is a wrape class of sofa software
    class GIAU
    {
        
    public:
        static double DJ00 ;
        static double DJC  ;
        static double DAS2R ;
        static double TURNAS;
        static double D2PI;
        static double DPI;
        
        double iauObl06(double date1, double date2);
        void iauPfw06(double date1, double date2,
                      double *gamb, double *phib, double *psib, double *epsa);
        void iauNut00a(double date1, double date2, double *dpsi, double *deps);
        void iauNut06a(double date1, double date2, double *dpsi, double *deps);
        void iauPnm06a(double date1, double date2, double rnpb[3][3]);
        
        void iauFw2m(double gamb, double phib, double psi, double eps,
                     double r[3][3]);
        
        void iauBpn2xy(double rbpn[3][3], double *x, double *y);
        double iauS06(double date1, double date2, double x, double y);
        double iauS00(double date1, double date2, double x, double y);
        
        void iauC2i06a(double date1, double date2, double rc2i[3][3]);
        
        void iauC2t06a(double tta, double ttb, double uta, double utb,
                       double xp, double yp, double rc2t[3][3]);
        
        double iauEra00(double dj1, double dj2);
        
        void iauC2ixy(double date1, double date2, double x, double y,
                      double rc2i[3][3]);
        
        void iauC2ixys(double x, double y, double s, double rc2i[3][3]);
        
        double iauSp00(double date1, double date2);
        void iauPom00(double xp, double yp, double sp, double rpom[3][3]);
        void iauC2tcio(double rc2i[3][3], double era, double rpom[3][3],
                       double rc2t[3][3]);
        
        void iauIr(double r[3][3]);
        void iauRx(double phi, double r[3][3]);
        void iauRy(double theta, double r[3][3]);
        void iauRz(double psi, double r[3][3]);
        double iauAnp(double a);
        void iauCr(double r[3][3], double c[3][3]);
        void iauCp(double p[3], double c[3]);
        void iauRxr(double a[3][3], double b[3][3], double atb[3][3]);
        
        double iauFame03(double t);
        double iauFal03(double t);
        double iauFalp03(double t);
        double iauFad03(double t);
        double iauFaf03(double t);
        double iauFaom03(double t);
        double iauFapa03(double t);
        double iauFave03(double t);
        double iauFae03(double t);
        double iauFama03(double t);
        double iauFaju03(double t);
        double iauFasa03(double t);
        double iauFaur03(double t);
        
        
    };

} // end of namespace



#endif /* GIAU_hpp */
