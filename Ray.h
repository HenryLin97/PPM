/*#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>*/
#include<stdio.h>
#include<vector>
#include<iostream>
#include<cmath>

#ifndef RAY_H
#define RAY_H

struct vector3d{
    double x, y, z;
    vector3d(double xx = 0, double yy = 0, double zz = 0)
    {
        x = xx; y = yy; z = zz;
    }
    
    inline vector3d operator + (const vector3d& r) const
    {
        vector3d temp(x + r.x, y + r.y, z + r.z);
        return temp;
    }
    inline vector3d operator - (const vector3d& r) const
    {
        vector3d temp(x - r.x, y - r.y, z - r.z);
        return temp;
    }
    inline vector3d operator + (double t) const
    {
        vector3d temp(x + t, y + t, z + t);
        return temp;
    }
    inline vector3d operator - (double t) const
    {
        vector3d temp(x - t, y - t, z - t);
        return temp;
    }
    inline vector3d operator * (double t) const
    {
        vector3d temp(x * t, y * t, z * t);
        return temp;
    }
    vector3d operator / (double t) const
    {
        vector3d temp(x / t, y / t, z / t);
        return temp;
    }
    void normalize()
    {
        double l2norm = std::sqrt(x*x + y*y + z*z);
        x /= l2norm;
        y /= l2norm;
        z /= l2norm;
    }
    void multiply(vector3d b)
    {
        x *= b.x;
        y *= b.y;
        z *= b.z;
    }
    inline vector3d mul(const vector3d& b)
    {
        return vector3d(x*b.x, y*b.y, z*b.z);
    }
    void threshold(double thres, int sign)
    {
        if(sign == 0) //lower bound
        {
            x = x-thres > 0?x : thres;
            y = y-thres > 0?y : thres;
            z = z-thres > 0?z : thres; 
        }
        else if(sign == 1) //upper bound
        {
            x = x-thres < 0?x : thres;
            y = y-thres < 0?y : thres;
            z = z-thres < 0?z : thres;
        }
    }
    double l2norm()
    {
        return std::sqrt(x*x + y*y + z*z);
    }
    inline vector3d norm()
    {
        double temp = std::sqrt(x*x + y*y + z*z);
        return vector3d(x/temp, y/temp, z/temp);
    }
    vector3d operator% (vector3d& b)
    {
        return vector3d(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
    }
    inline double dot(const vector3d &b) const {return x * b.x + y * b.y + z * b.z;}
};

std::ostream& operator<<(std::ostream& out, const vector3d& src)
    {
        out<<"("<<src.x<<" "<<src.y<<" "<<src.z<<")\n";
        //out<<'\n';
    }

double dot3d(vector3d a, vector3d b)
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

vector3d dot3d(double a, vector3d b)
{
    return vector3d(a*b.x, a*b.y, a*b.z);
}

vector3d neg(vector3d a)
{
    return vector3d(-a.x, -a.y, -a.z);
}
vector3d reflect(vector3d I, vector3d normal)
{
    I.normalize();
    normal.normalize();
    return I - dot3d(2 * dot3d(I, normal), normal);
}
vector3d refract(vector3d I_in, vector3d N, double ratio_iT)
{
    vector3d I = neg(I_in);
    I.normalize();
    N.normalize();
    double cosi = dot3d(I, N);
    if(1 - ratio_iT * ratio_iT * (1 - cosi * cosi) < 0)
    {
        return vector3d();
    }
    double cosT = std::sqrt(1 - ratio_iT * ratio_iT * (1 - cosi * cosi));
    vector3d T = dot3d(ratio_iT * cosi - cosT, N) - dot3d(ratio_iT, I);
    T.normalize();
    return T;
}
double det(vector3d a, vector3d b, vector3d c)
{
    return (a.x * b.y * c.z + b.x * c.y * a.z + c.x * a.y * b.z - a.x * c.y * b.z - b.x * a.y * c.z - c.x * b.y * a.z);
}
vector3d crossProd(vector3d a, vector3d b)
{
    return vector3d(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
vector3d element_max(vector3d a, vector3d b)
{
    vector3d result;
    result.x = a.x > b.x? a.x:b.x;
    result.y = a.y > b.y? a.y:b.y;
    result.z = a.z > b.z? a.z:b.z;
    return result;
}
vector3d element_min(vector3d a, vector3d b)
{
    vector3d result;
    result.x = a.x < b.x? a.x:b.x;
    result.y = a.y < b.y? a.y:b.y;
    result.z = a.z < b.z? a.z:b.z;
    return result;
}
vector3d axis_rotate(int axis, vector3d a, double theta)
{
    vector3d result;
    double Cos = cos(theta);
    double Sin = sin(theta);
    if(axis == 0)//x axis
    {
        result.x = a.x;
        result.y = Cos * a.y - Sin * a.z;
        result.z = Sin * a.y + Cos * a.z;
    }
    else if(axis == 1)//y axis
    {
        result.x = Cos * a.x + Sin * a.z;
        result.y = a.y;
        result.z = -Sin * a.x + Cos * a.z;
    }
    else if(axis == 2)//z axis
    {
        result.x = Cos * a.x - Sin * a.y;
        result.y = Sin * a.x + Cos * a.y;
        result.z = a.z;
    }
    else
    {
        result = a;
    }
    return result;
}

vector3d axis_rotate(int axis, vector3d a, double theta, double Cos, double Sin)
{
    vector3d result;
    if(axis == 0)//x axis
    {
        result.x = a.x;
        result.y = Cos * a.y - Sin * a.z;
        result.z = Sin * a.y + Cos * a.z;
    }
    else if(axis == 1)//y axis
    {
        result.x = Cos * a.x + Sin * a.z;
        result.y = a.y;
        result.z = -Sin * a.x + Cos * a.z;
    }
    else if(axis == 2)//z axis
    {
        result.x = Cos * a.x - Sin * a.y;
        result.y = Sin * a.x + Cos * a.y;
        result.z = a.z;
    }
    else
    {
        result = a;
    }
    return result;
}

class Ray{
public:
    vector3d R0;
    vector3d Rd;
    Ray()
    {
        R0 = vector3d(0, 0, 0);
        Rd = vector3d(0, 0, 0);
    }
    Ray(vector3d r0, vector3d rd)
    {
        R0 = r0;
        Rd = rd;
    }
    void NormalizeDirect()
    {
        double l2norm = std::sqrt(Rd.x * Rd.x + Rd.y * Rd.y + Rd.z * Rd.z);
        Rd.x /= l2norm;
        Rd.y /= l2norm;
        Rd.z /= l2norm;
    }
};

#endif