#ifndef PLANE_H
#define PLANE_H

#include"Ray.h"
#include"Object.h"
#include<CImg.h>

class Plane: public Object{
public:
    vector3d Point;
    vector3d Normal;
    double D;
    bool has_texture;
    //cv::Mat texture;
    cimg_library::CImg<unsigned char>* texture;
    vector3d text_down;//use to define the direction of texture
    vector3d text_right;
    double text_ratio;//make the texture larger                                                                                                                 
    Plane(vector3d point, vector3d normal, vector3d KD, OBJMT mt = MT_DIFF)
    {
        Point = point;
        Normal = normal;
        D = -normal.x * point.x - normal.y * point.y - normal.z * point.z;
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
        has_texture = false;
    }
    bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const;
    bool GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const;
    bool loadTexture(cimg_library::CImg<unsigned char>* text, vector3d down, vector3d right, double ratio);
    bool IsInside(vector3d point) const;
};

bool Plane::intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
{
    double temp_t = -(D + dot3d(ray.R0, Normal))/(dot3d(Normal, ray.Rd));
    if(temp_t - min_t > 0)
    {
        t = temp_t;
        n = Normal;
        return true;
    }
    else
    {
        return false;
    }
}

bool Plane::GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const
{
    if(has_texture == false)
    {
        ka = Ka;
        kd = Kd;
        ks = Ks;
        return true;
    }
    else
    {
        int rows = texture->width();
        int cols = texture->height();
        double x = dot3d(intersect_point, text_down)/text_ratio;
        double y = dot3d(intersect_point, text_right)/text_ratio;
        int i = (int)x % rows;
        int j = (int)y % cols;
        if(i < 0)
        {
            i += rows;
        }
        if(j < 0)
        {
            j += cols;
        }
        unsigned char b[3];
        for(int k = 0; k < 3; k++)
        {
            b[k] = *texture->data(i, j, 0, k);
        }
        kd = vector3d((double)(int)b[0]/255, (double)(int)b[1]/255, (double)(int)b[2]/255);
        ka = Ka;
        ks = vector3d(0, 0, 0);;
        return true;
    }
}

bool Plane::loadTexture(cimg_library::CImg<unsigned char>* text, vector3d down, vector3d right, double ratio = 1)
{
    has_texture = true;
    texture = text;
    text_down = down;
    text_right = right;
    text_ratio = ratio;
    return true;
}

bool Plane::IsInside(vector3d point) const
{
    vector3d direct = point - Point;
    double product = dot3d(direct, Normal);
    if(product > 0)
    {
        return false;
    }
    return true;
}

#endif