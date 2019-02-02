#ifndef OBJECT_H
#define OBJECT_H

#include "Ray.h"

enum OBJMT {MT_DIFF, MT_SPEC, MT_REFR};  // material types, used in radiance()

class Object{
public:
    vector3d Kd;
    vector3d Ks;
    vector3d Ka;
    vector3d ReflectCoe;
    vector3d RefractCoe;
    OBJMT MT;
    double RefractRatio;
    virtual bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const = 0;
    virtual bool GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const = 0;
    virtual bool IsInside(vector3d point) const = 0;//if the point is inside this object?
};


#endif