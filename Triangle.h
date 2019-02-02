#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Ray.h"
#include "Object.h"

class Triangle: public Object{
public:
    vector3d P[3];
    vector3d Normal;
    Triangle(vector3d p0, vector3d p1, vector3d p2, vector3d KD, OBJMT mt)
    {
        P[0] = p0;
        P[1] = p1;
        P[2] = p2;
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
        vector3d E1 = P[0] - P[1];
        vector3d E2 = P[0] - P[2];
        Normal = crossProd(E1, E2);
        Normal.normalize();
        MT = mt;
    }
    Triangle(vector3d normal, vector3d p0, vector3d p1, vector3d p2, vector3d KD, OBJMT mt)
    {
        P[0] = p0;
        P[1] = p1;
        P[2] = p2;
        Normal = normal;
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
        MT = mt;
    }
    bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
    {
        vector3d E1 = P[0] - P[1];
        vector3d E2 = P[0] - P[2];
        vector3d S = P[0] - ray.R0;
        //result = dot3d(1/det(ray.Rd, E1, E2), result);
        double temp = det(ray.Rd, E1, E2);
        vector3d result = {det(S, E1, E2)/temp, det(ray.Rd, S, E2)/temp, det(ray.Rd, E1, S)/temp};
        if(result.x - min_t > 0 && result.y >= 0 && result.z >= 0 && result.y + result.z <= 1)
        {
            t = result.x;
            n = Normal;
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
        //std::cout<<dot3d(point - P[0], Normal)<<"\n";
        if(dot3d(point - P[0], Normal) < 0)
        {
            return true;
        }
        return false;
    }
    void Scale(double ratio = 1)
    {
        P[0] = dot3d(ratio, P[0]);
        P[1] = dot3d(ratio, P[1]);
        P[2] = dot3d(ratio, P[2]);
    }
    void Translate(vector3d loc)
    {
        P[0] = P[0] + loc;
        P[1] = P[1] + loc;
        P[2] = P[2] + loc;
    }
    void Rotate(int axis, double theta, double Cos, double Sin)
    {
        P[0] = axis_rotate(axis, P[0], theta, Cos, Sin);
        P[1] = axis_rotate(axis, P[1], theta, Cos, Sin);
        P[2] = axis_rotate(axis, P[2], theta, Cos, Sin);
        Normal = axis_rotate(axis, Normal, theta, Cos, Sin);
    }
};

#endif