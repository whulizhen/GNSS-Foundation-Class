//
//  GQuaternion.cpp
//  GFC
//
//  Created by 李桢 on 2017/12/28.
//

#include "GQuaternion.hpp"
namespace gfc
{
    GQuaternion::GQuaternion()
    {
        r = 0.0;
    }
    
    GQuaternion::GQuaternion(double in_r, double in_i, double in_j, double in_k)
    {
        r = in_r;
        v.set(in_i, in_j, in_k);
    }
    
    GQuaternion::GQuaternion(double real, GVector vector)
    {
        r = real;
        v = vector;
    }
    
    
    void GQuaternion::set(double half_angle, GVector u)
    {
        u.normalise();
        r = cos(half_angle);
        v = u*sin(half_angle);
    }
    
    // tranfrom a vector using this quaternion
    // this gets the same results as matlab quatrotate
    GVector GQuaternion::rotate(gfc::GVector a)
    {
        //normalise itself firstly
        this->normalise();
        GVector b;
        
        GQuaternion q= *this;
        GQuaternion q_star = q;
        q_star.conjugate();
        GQuaternion c =  q_star *  a * q;
        b = c.v;
        return b;
    }
    
    
    void GQuaternion::normalise()
    {
        double n =  this->norm();
        r = r/n;
        v = v/n;
    }
    
    void GQuaternion::conjugate()
    {
        v = -v;
    }
    
    void GQuaternion::inverse()
    {
        GQuaternion q = *this;
        GQuaternion q_star=q;
        q_star.conjugate();
        double n2 = q.norm2();
        
        r = q_star.r/n2;
        v = q_star.v/n2;
    }
    
    double GQuaternion::norm()
    {
        return sqrt(this->norm2());
    }
    
    double GQuaternion::norm2()
    {
        return v.norm2() + r*r;
    }
    
    GQuaternion GQuaternion::operator-() const
    {
        return GQuaternion(-r, -v);
    }
    
    void GQuaternion::operator+=(GQuaternion a)
    {
        r += a.r;
        v += a.v;
    }
    
    void GQuaternion::operator-=(GQuaternion a)
    {
        r -= a.r;
        v -= a.v;
    }
    
    void GQuaternion::operator*=(GQuaternion a)
    {
        double r1 = r*a.r - dotproduct(v,a.v);
        v = r*a.v + a.r*v + crossproduct(v, a.v);
        r = r1;
        
    }
    
    void GQuaternion::operator*=(GVector a)
    {
        GQuaternion q(0.0,a);
        *this *= q;
    }
    
    
    void GQuaternion::operator/=(GQuaternion a)
    {
        double n = a.norm2();
        double t0 = a.r*r + a.v.x*v.x + a.v.y*v.y + a.v.z*v.z;
        double t1 = a.r*v.x - a.v.x*r - a.v.y*v.z + a.v.z*v.y;
        double t2 = a.r*v.y + a.v.x*v.z - a.v.y*r - a.v.z*v.x;
        double t3 = a.r*v.z - a.v.x*v.y + a.v.y*v.x - a.v.z*r;
        
        r = t0/n;
        v.set(t1/n, t2/n, t3/n);
    }
    
    void GQuaternion::operator*=(double k)
    {
        r = r*k;
        v = v*k;
    }
    
    void GQuaternion::operator/=(double k)
    {
        r = r/k;
        v = v/k;
    }
    
    
    
    GQuaternion operator+(GQuaternion a, GQuaternion b)
    {
        a += b;
        return a;
    }
    
    
    GQuaternion operator-(GQuaternion a, GQuaternion b)
    {
        a -= b;
        return a;
    }
    GQuaternion operator*(GQuaternion a, GQuaternion b)
    {
        a *= b;
        return a;
    }
    
    GQuaternion operator*(GQuaternion a, GVector b)  //!< Hamilton product
    {
         a*=b;
        return a;
    }
    
    GQuaternion operator*(GVector a, GQuaternion b)  //!< Hamilton product
    {
        b*=a;
        return b;
    }
    
    GQuaternion operator/(GQuaternion a, GQuaternion b)  //!< Hamilton product
    {
        a/=b;
        return a;
    }
    
}
