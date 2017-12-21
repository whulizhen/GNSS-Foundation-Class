#ifndef  GFCMATRIX_H
#define  GFCMATRIX_H


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

#include "Platform.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <valarray>

//#include "../lapacklib/lapacke.h"
//#ifdef _WIN32_
//动态链接
//#pragma comment(lib,"../lapacklib/libblas.lib")
//#pragma comment(lib,"../lapacklib/liblapacke.lib")
//#endif

using namespace std;

namespace gfc
{
	class GMatrix
	{
	  	
        
	public:
        
        friend double dotproduct(GMatrix a, GMatrix b);
        friend GMatrix crossproduct(GMatrix a, GMatrix b);
        
		GMatrix();  //构造函数
		GMatrix( double* data,int nrow,int ncol);
        GMatrix(int nrow,int ncol);
		GMatrix( const GMatrix& matrix);   //拷贝构造函数
		GMatrix& operator= (const GMatrix& right);  //赋值重载
        
		GMatrix  operator/ (const GMatrix& right);  //矩阵右除运算,对于右除式A/B，相当于A*inv(B)
		//gfcMatrix  operator\ (const gfcMatrix& right);  //矩阵左除运算
		
		
        void  operator *=( GMatrix right);
        void  operator *=(double right);
        void operator/=(double right);
        void operator+=(GMatrix right);
        void operator-=(gfc::GMatrix right);
        GMatrix operator-();  // negative
        
        double& operator[](int index);
        double& operator()(int index_row, int index_col);
        
        void normalise();
        void resize(int row, int col);
		double det();             //矩阵的行列式  
        double norm();   // get the norm of matrix
        int getRowNO() { return m_rowno;}
        int getColNO() { return m_colno;}
        
        void clear()
        {
            m_data.resize(m_rowno*m_colno);
        }
        
        void getData(double* data)
        {
            memcpy(data,&m_data[0],sizeof(double)*m_rowno*m_colno);
        }
        void setData(double* data,int row,int col);
        
		//流输入运算符重载
		friend inline ostream& operator>> (ostream& os, GMatrix& m)
		{
			return os;
		}
        
		//流输出运算符重载
		friend inline ostream &operator<< (ostream& os,GMatrix m)
		{
			//os<<"matrix object"<<std::endl;
            os.setf(ios::fixed);//固定小数位数
            os<<setw(20);
			os<<setprecision(16);
			for( int i = 0 ; i< m.m_rowno ; i++ )
			{
				for( int j = 0 ; j < m.m_colno ; j++ )
				{
					os<<m.m_data[i*m.m_colno+j]<<"  ";
				}
				os<<std::endl;
			}
			return os;
		}
		
		void dump();  //输出
		//用于测试调用lapack的函数
		void testLapack();
        
		//void print_matrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );
		//void print_int_vector( char* desc, lapack_int n, lapack_int* a );
		
		~GMatrix(); //析构函数
	  
      //LU分解相关的函数
      static  int  ludcom( double *A, int n, int *indx, double& d);
      static  void lubksb( const double *A, int n, const int *indx, double *b);
      static  void multiply( double* A, int rA, int cA, double* B, int cB, double* C);
      
      
    private:
		
        std::valarray<double> m_data;
        
		int      m_rowno; //行数
		int      m_colno; //列数
        
	};
    
    
    GMatrix operator+(GMatrix a, GMatrix b); //!< returns a+b
    GMatrix operator-(GMatrix a, GMatrix b); //!< returns a-b
    
    GMatrix operator*(GMatrix a, double b);
    GMatrix operator*(double b, GMatrix a);
    GMatrix operator*(GMatrix a, GMatrix b);
    
    GMatrix operator/(GMatrix a, double b);
    
    GMatrix  operator~ (GMatrix a); //转置
    GMatrix  operator! (GMatrix a); //求逆
    
    GMatrix normalise(GMatrix a);  // normalise
    
    double dotproduct(GMatrix a, GMatrix b);
    
    //GMatrix crossproduct(GMatrix a, GMatrix b);
    
}

#endif
