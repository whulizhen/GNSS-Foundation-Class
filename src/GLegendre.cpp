
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
//  GLegendre.cpp
//  GFC
//
//  Created by lizhen on 13/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GLegendre.hpp"

namespace gfc
{
   
    double GLegendre::PI = GCONST("PI");
    
    GLegendre::GLegendre(int n , int m)
    {
        N = n ;
        M = m ;
    }
    
    GLegendre::~GLegendre(void)
    {
        
    }
    
    
    //计算n的阶乘
    double  GLegendre::factorial(float n)
    {
        double fac = 0 ;
        if(n==0)
            fac =1 ;
        else
        {
            fac = n*factorial(n-1);
        }
        return fac ;
    }
    
    /*
     //参数y是一个以弧度为单位的角度值,没有归一化；
     已验证，与matlab结果一致，无需修改
     */
    double GLegendre::P(double n , double x)
    {
        double a = 0 ;
        if(n==0)
            a =    1;
        else if(n==1)
            a =   x;
        else if(n>1)
        {
            double m= (n-1)*1.0;
            a= (2*m+1)/(m+1)*x*P(m,x)- m/(m+1)*P((m-1),x);
        }
        return a ;
    }
    
    /*
     超高阶缔合勒让德多项式计算（约2160阶）
     参考文献[1]：常用超高阶次缔合勒让德函数计算方法对比分析；王贱强 赵国强 朱广彬
     武汉大学测绘学院；大地测量与地球动力学；2009年4月第29卷第2期
     参考文献[2] ：数学物理方法与仿真 杨华军编著 电子工业出版社 2011年第2版 P310 16.4.6 勒让德函数的递推公式
     参考[3]：http://mitgcm.org/~mlosch/geoidcookbook/node11.html
     本函数中的x为cos(phi)
     没有归一化的计算结果和matlab中函数legendre一致；但是文献3归一化后结果不同
     */
    
    double GLegendre::PNM( double n , double m , double x)
    {
        double y =-1.0;
        //double k = 2.0 ;
        if(n<m)  y = 0 ;
        if(m==0)
        {
            if(n==0)  y = 1 ;
            else if(n==1)  y =x;
            else if(n==2)  y =0.5*(3*x*x-1);
            else y=((2*n-1)*x*PNM(n-1 , m , x)-(n+m-1)*PNM(n-2, m, x))/(n-m);
        }
        else if(m>0)
        {
            y= ( -(n+m-1)*(n+m)*PNM(n-1,m-1,x) + (n-m+1)*(n-m+2)*PNM(n+1,m-1,x) )/((2*n+1)*sqrt(1.0-x*x));
        }
        else if(m<0)
        {
            y  = ( -PNM(n+1,m+1, x) + PNM(n-1,m+1,x))/ ((2*n+1)*sqrt(1.0-x*x));
        }
        return y;       //无归一化 已验证正确
    }
    
    /*
     直接递推计算归一化的连带勒让德多项式
     参考文献：A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre funcitons
     S.A.Holmes, W.E.Featherstone  ; Journal of Geodesy
     参考：http://mitgcm.org/~mlosch/geoidcookbook/node11.html
     本函数所采用的方法是标准向前列递推法
     xx=sin(phi)   ;展开为sin(phi)的函数
     已验证正确
     */
    double GLegendre::NPNM_sin(int n , int m  , double sx)
    {
        double y = -1;
        double anm =0;
        double bnm = 0;
        if (n==0&&m==0)
            y = 1 ;
        else if(n==1&&m==1)
            y = sqrt(3.0)*sqrt(1-sx*sx);
        else
        {
            if(n>m)
            {
                anm = sqrt(  (2.0*n-1)*(2*n+1)/(n-m)/(n+m) );
                bnm = sqrt(  (2.0*n+1)*(n+m-1)*(n-m-1)/(n-m)/(n+m)/(2*n-3) );
                y = anm*sx*NPNM_sin(n-1, m, sx) - bnm*NPNM_sin(n-2 , m , sx );
            }
            else if(n==m)
            {
                y = sqrt(1-sx*sx) *sqrt( (2.0*m+1)/(2*m) )*NPNM_sin(m-1 , m -1 , sx);
            }
            else if(n<m)  //m不能大于n
            {
                y = 0 ;
            }
        }
        return y ;
    }
    
    /*
     以弧度展开的勒让德多项式的值
     */
    double GLegendre::NPNM( int n , int  m , double x)
    {
        double y = -1 ;
        double sx = sin(x);
        y = NPNM_sin( n ,  m  ,  sx) ;
        return y ;
    }
    
    /*
     //计算x的Pnm球谐函数（已经归一化）
     e是一个余弦值
     */
    double GLegendre::YNM( int n, int m, double xx)
    {
        double x = cos(xx);
        int i;
        double y, yp, ypp, f, fp, s = 1 - x*x;
        if (m < 0) m = -m;
        if (m > n|| s < 0)
            return 0;
        for (yp = 1, i = 1; i <= m; i++)
            yp *= (1 + .5/i)*s;
        yp = sqrt(yp/(4*PI)) * ((m % 2 && m > 0) ? -1: 1);
        for (fp = 1, ypp = 0, i = m + 1; i <= n; i++, fp = f, ypp = yp, yp = y)
        {
            f = sqrt((4.*i*i-1)/((i-m)*(i+m)));
            y = f*(x*yp - ypp/fp);
        }
        return yp;
    }
    
    
    /*
     //球谐函数的参数拟合spherical harmonics
     将某一未知函数f(e , z)展开到N阶M级的球谐函数
     并用最小二乘法求其系数
     对于PCV标定, e为卫星高度角,z为天线的方位角,单位是弧度
     其中M<=N
     返回值为P ;其中的元素排列顺序为：A00 A10 A11 B11 A20 A21 B21 A22 B22......
     */
    void GLegendre::Sph_Harmonics( double e, double z , double *P)
    {
        //int total = N*N+2*M+1 ;  //总共的参数个数
        //int cntA = N+ N*(N-1)/2 + M +1 ;
        for(int n=0;n<N+1;n++ )
        {
            if( n<N )   //前N项中M都是以N为准
            {
                for( int t=0; t<2*n+1 ;t++)
                {
                    if(t==0)
                    {
                        // P[n*n] = NPNM(n , 0 , PI/2.0-e);
                        P[n*n] = NPNM(n , 0 , e);
                    }
                    else if((t%2==0)&&(t!=0))    //项内偶数
                    {
                        int m = t/2;
                        P[n*n+t] =NPNM(n, m, e)*sin(m*z);
                        //P[n*n+t] =NPNM(n, m, PI/2.0-e)*sin(m*z);
                    }
                    if((t+1)%2==0)            //项内奇数
                    {
                        int m = int(t/2) +1 ;
                        //P[n*n +t] = NPNM(n, m , PI/2.0-e)*cos(m*z);
                        P[n*n +t] = NPNM(n, m , e)*cos(m*z);
                    }
                }
            }
            if(n==N)  //最后一个N才是以M为准
            {
                for(int t=0; t<2*M+1 ;t++)
                {
                    if(t==0)
                    {
                        //P[n*n] = NPNM(n, 0 , PI/2.0-e);
                        P[n*n] = NPNM(n, 0 , e);
                    }
                    else if((t%2==0)&&(t!=0))
                    {
                        int m = t/2 ;
                        //P[n*n+t] = NPNM(n, m, PI/2.0-e)*sin(m*z);
                        P[n*n+t] = NPNM(n, m, e)*sin(m*z);
                    }
                    if((t+1)%2==0)
                    {
                        int m =int(t/2) +1 ;
                        //P[n*n +t] = NPNM(n, m , PI/2.0-e)*cos(m*z);
                        P[n*n +t] = NPNM(n, m , e)*cos(m*z);
                    }
                }
            }
        }
        //cout<<endl<<"球谐系数展开完成!"<<endl;
    }
    
    /*
     球谐函数拟合
     计算出系数
     data中包含进行拟合的系数
     进行电离层TEC拟合时，每一行包含经度，纬度，TEC值,即n行3列
     countX是需要解算的未知数个数
     countD是观测值的个数
     M和N勒让德级数展开的阶数
     x是输出的拟合系数
     */
    
    void GLegendre::SH_Fit( int countD  , double *data , double *x , double* dta)
    {
        int countX = N*N+2*M+1 ; 
        double *b = new double[countX*countD];   //设计矩阵B
        memset(b , 0 , sizeof(double)*countX*countD);
        double *l = new double[countD];    //L=B*X
        memset(l , 0 , sizeof(double)*countD) ;
        double *tempX = new double[countX];     //某一行的勒让德多项式系数
        memset(tempX , 0 , sizeof(double)*countX);
        for(int i=0; i<countD ; i++ )    //行
        {
            l[i] = data[i*3+2];   
            Sph_Harmonics(data[i*3+1], data[i*3+0] , tempX);   //对于电离层数据,纬度是高度角，经度是方位角
            for(int j=0 ; j<countX ; j++)  //列
            {
                b[i*countX+j] = tempX[j];
            }
        }
        delete[] tempX ;   //清楚掉临时变量，释放内存
        
//        CMatrix B(b , countD , countX);
//        CMatrix L(l , countD ,1 );
//        CMatrix X , P , V , Dta;
//        P.SetIdentity(countD);
//        X = (B.GetTransposedMatrix()*B).GetInverseMatrix() *B.GetTransposedMatrix()*L;   //解算得到各个系数
//        X.GetData(x);
//        
//        V = B*X -L;
//        Dta = V.GetTransposedMatrix()*P*V;
//        Dta.GetData(dta);
//        (*dta)=sqrt((*dta)/(countD-countX));
//        //cout<<"拟合中误差："<<sqrt(dta/(countD -countX))<<endl;
//        delete[] b;
//        delete[] l;
//        //cout<<endl<<"球谐系数拟合完成！"<<endl;
    }
    
    /*
     球谐函数值的计算
     */
    double GLegendre::SPH_Function( double e , double z , double *X )
    {
        int total = N*N+2*M+1 ;
        double *P = new double[total];
        Sph_Harmonics( e, z ,P);
        double sum = 0.0 ;
        for(int i = 0 ; i< total ; i++)
        {
            sum = sum + P[i]*X[i];
        }
        delete []P ;
        return sum ;
    }
    
} // end of namespace gfc