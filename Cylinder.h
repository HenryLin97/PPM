#ifndef CYLINDER_H
#define CYLINDER_H

#include"Ray.h"
#include"Object.h"
#include"Plane.h"
#include"math.h"

class Cylinder: public Object{
public:
    vector3d Center;
    double Height;
    double Radius;
    Plane* bottom;
    Plane* top;
    Cylinder(vector3d center, double height, double r, vector3d KD, OBJMT mt)
    {
        Center = center;
        Height = height;
        Radius = r;
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
        MT = mt;
        bottom = new Plane(Center, vector3d(0, -1, 0), KD);
        top = new Plane(Center + vector3d(0, height, 0), vector3d(0, 1, 0), KD);
    }
    bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
    {
        double result_t = 2*1e5;
        vector3d result_n = vector3d(0, 0, 0);
        double temp_t;
        vector3d temp_n;
        
        //try bottom
        if(bottom->intersect(ray, min_t, temp_t, temp_n))
        {
            if(temp_t < result_t)
            {
                vector3d inter_point = ray.R0 + ray.Rd * temp_t;
                if((inter_point - Center).l2norm() < Radius)
                {
                    result_t = temp_t;
                    result_n = vector3d(0, -1, 0);
                }
            }
        }
        //try top
        if(top->intersect(ray, min_t, temp_t, temp_n))
        {
            if(temp_t < result_t)
            {
                vector3d inter_point = ray.R0 + ray.Rd * temp_t;
                if((inter_point - vector3d(0, Height, 0) - Center).l2norm() < Radius)
                {
                    result_t = temp_t;
                    result_n = vector3d(0, 1, 0);
                }
            }
        }
        //try round
        double a = ray.Rd.x * ray.Rd.x + ray.Rd.z * ray.Rd.z;
        double b = 2 * ray.Rd.x * (ray.R0.x - Center.x) + 2 * ray.Rd.z * (ray.R0.z - Center.z);
        double c = (ray.R0.x - Center.x) * (ray.R0.x - Center.x) + (ray.R0.z - Center.z) * (ray.R0.z - Center.z) - Radius * Radius;
        if(a != 0)
        {
            double delta = b*b - 4*a*c;
            if(delta >= 0)
            {
                //first
                temp_t = (- b - sqrt(delta))/(2*a);
                vector3d temp_inter = ray.R0 + ray.Rd * temp_t;
                if(temp_t > min_t && temp_t < result_t && temp_inter.y > Center.y && temp_inter.y < Center.y + Height)
                {
                    result_t = temp_t;
                    result_n = vector3d(temp_inter.x - Center.x, 0, temp_inter.z - Center.z);
                    result_n.normalize();
                }
                //second
                temp_t = (- b + sqrt(delta))/(2*a);
                temp_inter = ray.R0 + ray.Rd * temp_t;
                if(temp_t > min_t && temp_t < result_t && temp_inter.y > Center.y && temp_inter.y < Center.y + Height)
                {
                    result_t = temp_t;
                    result_n = vector3d(temp_inter.x - Center.x, 0, temp_inter.z - Center.z);
                    result_n.normalize();
                }
            }
        }

        if(result_t < 1e5)
        {
            t = result_t;
            n = result_n;
            return true;
        }
        return false;
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
        return false;
    }
    bool intersect_piece(Ray ray, const double min_t, double& t,vector3d& n, double ratio) const
    {
        vector3d piece_center = Center + vector3d(0, Height * ratio, 0);
        Plane piece(piece_center, vector3d(0, 1, 0), vector3d(0, 0, 0));
        double temp_t;
        vector3d temp_n;
        if(piece.intersect(ray, min_t, temp_t, temp_n))
        {
            vector3d inter_point = ray.R0 + ray.Rd * temp_t;
            if((inter_point - piece_center).l2norm() < Radius)
            {
                t = temp_t;
                n = vector3d(0, -1, 0);
                return true;
            }
        }
        return false;
    }
};

#endif