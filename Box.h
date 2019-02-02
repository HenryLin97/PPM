#ifndef BOX_H
#define BOX_H

#include "Ray.h"
#include "Object.h"

double sign(double x)
{
    if(x > 0)
        return 1;
    return -1;
}

class Box: public Object{
public:
    double Upper[3];//in order x, y, z
    double Lower[3];
    Box(vector3d up_corner, vector3d low_corner)
    {
        Upper[0] = up_corner.x;
        Upper[1] = up_corner.y;
        Upper[2] = up_corner.z;
        Lower[0] = low_corner.x;
        Lower[1] = low_corner.y;
        Lower[2] = low_corner.z;
        Kd = vector3d(0, 0, 0);
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
    }
    bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const;
    bool GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const;
    bool IsInside(vector3d point) const;
};

bool Box::intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
{
    double t_min_cache[3], t_max_cache[3];
    vector3d R0 = ray.R0;
    vector3d Rd = ray.Rd;
    //x axis
    if(Rd.x > 0)
    {
        t_min_cache[0] = -(R0.x - Lower[0])/Rd.x;
        t_max_cache[0] = -(R0.x - Upper[0])/Rd.x;
    }
    else if(Rd.x < 0)
    {
        t_max_cache[0] = -(R0.x - Lower[0])/Rd.x;
        t_min_cache[0] = -(R0.x - Upper[0])/Rd.x;
    }
    else
    {
        if(R0.x < Lower[0] || R0.x > Upper[0])
            return false;
        t_max_cache[0] = 1e10;
        t_min_cache[0] = -1e10;
    }
    //y axis
    if(Rd.y > 0)
    {
        t_min_cache[1] = -(R0.y - Lower[1])/Rd.y;
        t_max_cache[1] = -(R0.y - Upper[1])/Rd.y;
    }
    else if(Rd.y < 0)
    {
        t_max_cache[1] = -(R0.y - Lower[1])/Rd.y;
        t_min_cache[1] = -(R0.y - Upper[1])/Rd.y;
    }
    else
    {
        if(R0.y < Lower[1] || R0.y > Upper[1])
            return false;
        t_max_cache[1] = 1e10;
        t_min_cache[1] = -1e10;
    }
    //z axis
    if(Rd.z > 0)
    {
        t_min_cache[2] = -(R0.z - Lower[2])/Rd.z;
        t_max_cache[2] = -(R0.z - Upper[2])/Rd.z;
    }
    else if(Rd.z < 0)
    {
        t_max_cache[2] = -(R0.z - Lower[2])/Rd.z;
        t_min_cache[2] = -(R0.z - Upper[2])/Rd.z;
    }
    else
    {
        if(R0.z < Lower[2] || R0.z > Upper[2])
            return false;
        t_max_cache[2] = 1e10;
        t_min_cache[2] = -1e10;
    }
    //select max's min
    double t_max = t_max_cache[0];
    int t_max_index = 0;
    for(int i = 1; i < 3; i++)
    {
        if(t_max_cache[i] < t_max)
        {
            t_max = t_max_cache[i];
            t_max_index = i;
        }
    }
    //select min's max
    double t_min = t_min_cache[0];
    int t_min_index = 0;
    for(int i = 1; i < 3; i++)
    {
        if(t_min_cache[i] > t_min)
        {
            t_min = t_min_cache[i];
            t_min_index = i;
        }
    }
    if(t_max <= t_min)
    {
        return false;
    }
    if(t_min > min_t)//ray from outside
    {
        t = t_min;
        switch(t_min_index)
        {
            case 0: n = vector3d(-sign(Rd.x), 0, 0);break;
            case 1: n = vector3d(0, -sign(Rd.y), 0);break;
            case 2: n = vector3d(0, 0, -sign(Rd.z));break;
            default: break;
        }
        return true;
    }
    else if(t_max > min_t)//ray from inside
    {
        t = t_max;
        switch(t_max_index)
        {
            case 0: n = vector3d(-sign(Rd.x), 0, 0);break;
            case 1: n = vector3d(0, -sign(Rd.y), 0);break;
            case 2: n = vector3d(0, 0, -sign(Rd.z));break;
            default: break;
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool Box::GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const
{
    kd = Kd;
    ks = Ks;
    ka = Ka;
    return true;
}

bool Box::IsInside(vector3d point) const
{
    if(point.x < Lower[0])
        return false;
    if(point.x > Upper[0])
        return false;
    if(point.y < Lower[1])
        return false;
    if(point.y > Upper[1])
        return false;
    if(point.z < Lower[2])
        return false;
    if(point.z > Upper[2])
        return false;
    return true;
}

#endif