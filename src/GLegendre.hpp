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
//  GLegendre.hpp
//  GFC
//
//  Created by lizhen on 13/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GLegendre_hpp
#define GLegendre_hpp

#include <stdio.h>
#include <math.h>
#include "GFCCONST.h"
namespace gfc {
/*

 勒让德多项式计算类
 */

class GLegendre {
 private:
  double NPNM_sin(int n, int m, double sx);  //以sinx展开的归一化勒让德多项式的值
 public:
  double factorial(float n);     //计算阶乘
  double P(double n, double x);  // 计算非归一化勒让德多项式的值Pn(x)
  double PNM(double n, double m, double x);  //计算非归一化的连带勒让德多项式的值Pnm(x)
  double NPNM(int n, int m, double x);  //以角度进行展开的归一化勒让德多项式的值
  double YNM(int n, int m, double xx);  //计算归一化球谐函数值
  void  Sph_Harmonics(double e, double z,double *P);  //利用归一化的连带勒让德多项式计算球谐展开系数
  void  SH_Fit(int countD, double *data, double *x,double *dta);  ///利用球谐系数，采用最小二乘计算对球谐系数进行拟合
  
  double SPH_Function( double e, double z,double *X);  //利用最小二乘得到的系数，计算球谐函数展开的函数值
  
  GLegendre(int n, int m);
  
  ~GLegendre(void);
    
  static double PI;
    
 private:
  int N;
  int M;
};

}  // end of namespace gfc

#endif /* GLegendre_hpp */
