#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include "Object.h"
#include<iostream>

class Sphere: public Object{
public:
    vector3d Center;
    double Radius;
    Sphere(vector3d center, double radius)
    {
        Center = center;
        Radius = radius;
        Kd = vector3d(0, 0, 0);
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
    }
    Sphere(double radius, vector3d center, vector3d KD, OBJMT mt)
    {
        Center = center;
        Radius = radius;
        MT = mt;
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
    }
    bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
    {
        vector3d R0 = ray.R0;
        vector3d Rd = ray.Rd;
        vector3d l = Center - R0;
        double tp = dot3d(l, Rd);
        //std::cout<<"here1!\n";
        //std::cout<<"l:"<<l.x<<" "<<l.y<<" "<<l.z<<"\n";
        //std::cout<<"Rd:"<<Rd.x<<" "<<Rd.y<<" "<<Rd.z<<"\n";
        //std::cout<<"tp:"<<tp<<"\n";
        if(tp < 0 && dot3d(l, l) - Radius * Radius >= 0)
        {
            return false;
        }
        double d2 = dot3d(l, l) - tp * tp;
        //std::cout<<"here2!\n";
        if(d2 - Radius * Radius > 0)
        {
            return false;
        }
        double tb = std::sqrt(Radius * Radius - d2); 
        //std::cout<<"here3!\n";
        if(dot3d(l, l) - Radius * Radius <= min_t) //inside the ball
        {
            //std::cout<<"here4!\n";
            if(tp + tb < min_t)
            {
                return false;
            }
            t = tp + tb;
            vector3d intersect_point = R0 + dot3d(t, Rd);
            //n = Center - intersect_point;
            n = intersect_point - Center;
            n.normalize();
            return true;
        }
        else if(dot3d(l, l) - Radius * Radius > min_t) //outside the ball
        {
            if(tp - tb < min_t) //then we should consider the light source as inside the ball
            {
                t = tp + tb;
                vector3d intersect_point = R0 + dot3d(t, Rd);
                n = Center - intersect_point;
                n.normalize();
                return true;
            }
            t = tp - tb;
            vector3d intersect_point = R0 + dot3d(t, Rd);
            n = intersect_point - Center;
            n.normalize();
            return true;
        }
        else
        {
            t = tp + tb;
            vector3d intersect_point = R0 + dot3d(t, Rd);
            n = Center - intersect_point;
            n.normalize();
            return true;
        }
    }
    bool GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const
    {
        ka = Ka;
        kd = Kd;
        ks = Ks;
        return true;
    }
    bool IsInside(vector3d point) const
    {
        vector3d l = Center - point;
        if(l.l2norm() > Radius)
        {
            return false;
        }
        return true;
    }
};

#endif