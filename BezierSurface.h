#ifndef BEZIERSURFACE_H
#define BEZIERSURFACE_H

#include"Object.h"
#include"MyMatrix.h"
#include"Bezier.h"
#include"Ray.h"
#include"Cylinder.h"
#include<ctime>

class BezierSurface: public Object{
public:
vector3d Center;
Bezier Curve;
double radius;
double Scale;
Cylinder* BoundCylinder;//this Cylinder is constructed in the small scene
Cylinder* SmallCylinder;//this is a small Cylinder used to init values
Cylinder* ThinCylinder;
    BezierSurface(vector3d center, Bezier curve, vector3d KD, OBJMT mt, double scale = 1)
    {
        Center = center;
        Curve = curve;
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
        Scale = scale;
        MT = mt;
        vector2d range = Curve.range_y();
        radius = Curve.max_x() * 1.2;
        BoundCylinder = new Cylinder(Center/Scale + vector3d(0, range.x, 0), range.y - range.x, radius, KD, mt);
        SmallCylinder = new Cylinder(Center/Scale + vector3d(0, range.x, 0), range.y - range.x, radius/1.6, KD, mt);
        ThinCylinder = new Cylinder(Center/Scale + vector3d(0, range.x, 0), range.y - range.x, radius/4.8, KD, mt);
        srand(time(0));
    }
    vector3d S(double u, double theta) const
    {
        vector2d Pu = Curve.val(u);
        return vector3d(Center.x/Scale + sin(theta) * Pu.x, Center.y/Scale + Pu.y, Center.z/Scale + cos(theta) * Pu.x);
    }
    vector3d F(const Ray ray, double t, double u, double theta) const
    {
        return ray.R0/Scale + ray.Rd * t - S(u, theta);
    }
    vector3d initVal(Ray ray, vector<vector3d>& init_stack) const
    {
        /*std::cout<<BoundCylinder->Radius<<"\n";
        std::cout<<BoundCylinder->Center;
        std::cout<<BoundCylinder->Height;*/
        Ray small_ray;
        small_ray.R0 = ray.R0/Scale;
        small_ray.Rd = ray.Rd;
        double temp_t;
        vector3d temp_n;
        vector3d inter_result;
        if(BoundCylinder->intersect(small_ray, 1e-4, temp_t, temp_n))
        {
            //propose 11 possible initial points
            int step_num = 10;
            double theta;
            init_stack.clear();
            for(int i = 0; i <= step_num; i++)
            {
                //(double)std::rand()/(double)RAND_MAX
                if(BoundCylinder->intersect_piece(small_ray, 1e-4, temp_t, temp_n, (double)i/(double)step_num))
                {
                    vector3d inter_point = small_ray.R0 + small_ray.Rd * temp_t - BoundCylinder->Center;
                    double z = inter_point.z;
                    double x = inter_point.x;
                    double r = sqrt(x*x + z*z);
                    if(r > 0)
                    {
                        if(x > 0)
                        {
                            theta = acos(z/r);
                            inter_result = vector3d(temp_t, 0.5, theta);
                        }
                        else
                        {
                            theta = acos(z/r);
                            inter_result = vector3d(temp_t, 0.5, -theta);
                        }
                    }
                    else
                    {
                        inter_result = vector3d(temp_t, 0.5, 0.5);
                    }
                    init_stack.push_back(inter_result);
                }
            }
            //try thin Cylinder
            if(ThinCylinder->intersect(small_ray, 1e-4, temp_t, temp_n))
            {
                if(temp_n.x == 0 && temp_n.z == 0)
                {
                    vector3d inter_point = small_ray.R0 + small_ray.Rd * temp_t - BoundCylinder->Center;
                    double z = inter_point.z;
                    double x = inter_point.x;
                    double r = sqrt(x*x + z*z);
                    if(r > 0)
                    {
                        if(x > 0)
                        {
                            theta = acos(z/r);
                            inter_result = vector3d(temp_t, 0.5, theta);
                        }
                        else
                        {
                            theta = acos(z/r);
                            inter_result = vector3d(temp_t, 0.5, -theta);
                        }
                    }
                    else
                    {
                        inter_result = vector3d(temp_t, 0.5, 0.5);
                    }
                }
                else if(temp_n.x > 0)
                {
                    theta = acos(temp_n.z);
                    inter_result = vector3d(temp_t, 0.5, theta);
                }
                else
                {
                    theta = acos(temp_n.z);
                    inter_result = vector3d(temp_t, 0.5, -theta);
                }
                init_stack.push_back(inter_result);
                inter_result.y = 0.2;
                init_stack.push_back(inter_result);
                inter_result.y = 0.9;
                init_stack.push_back(inter_result);
            }
            //end Thin
            //try small Cylinder
            if(SmallCylinder->intersect(small_ray, 1e-4, temp_t, temp_n))
            {
                if(temp_n.x == 0 && temp_n.z == 0)
                {
                    vector3d inter_point = small_ray.R0 + small_ray.Rd * temp_t - BoundCylinder->Center;
                    double z = inter_point.z;
                    double x = inter_point.x;
                    double r = sqrt(x*x + z*z);
                    if(r > 0)
                    {
                        if(x > 0)
                        {
                            theta = acos(z/r);
                            inter_result = vector3d(temp_t, 0.5, theta);
                        }
                        else
                        {
                            theta = acos(z/r);
                            inter_result = vector3d(temp_t, 0.5, -theta);
                        }
                    }
                    else
                    {
                        inter_result = vector3d(temp_t, 0.5, 0.5);
                    }
                }
                else if(temp_n.x > 0)
                {
                    theta = acos(temp_n.z);
                    inter_result = vector3d(temp_t, 0.5, theta);
                }
                else
                {
                    theta = acos(temp_n.z);
                    inter_result = vector3d(temp_t, 0.5, -theta);
                }
                init_stack.push_back(inter_result);
                inter_result.y = 0.2;
                init_stack.push_back(inter_result);
                inter_result.y = 0.9;
                init_stack.push_back(inter_result);
            }
            //end small
            BoundCylinder->intersect(small_ray, 1e-4, temp_t, temp_n);
            if(temp_n.x == 0 && temp_n.z == 0)
            {
                vector3d inter_point = small_ray.R0 + small_ray.Rd * temp_t - BoundCylinder->Center;
                double z = inter_point.z;
                double x = inter_point.x;
                double r = sqrt(x*x + z*z);
                if(r > 0)
                {
                    if(x > 0)
                    {
                        theta = acos(z/r);
                        inter_result = vector3d(temp_t, 0.5, theta);
                    }
                    else
                    {
                        theta = acos(z/r);
                        inter_result = vector3d(temp_t, 0.5, -theta);
                    }
                }
                else
                {
                    inter_result = vector3d(temp_t, 0.5, 0.5);
                }
            }
            else if(temp_n.x > 0)
            {
                theta = acos(temp_n.z);
                inter_result = vector3d(temp_t, 0.5, theta);
            }
            else
            {
                theta = acos(temp_n.z);
                inter_result = vector3d(temp_t, 0.5, -theta);
            }
            init_stack.push_back(inter_result);
            inter_result.y = 0.2;
            init_stack.push_back(inter_result);
            inter_result.y = 0.9;
            init_stack.push_back(inter_result);
            return inter_result;
        }
        //std::cout<<"here!\n";
        return vector3d(0, 0, 0);
    }
    bool intersect(const Ray ray, const double min_t, double& tt, vector3d& nn) const
    {
        double u = 0.5, theta = 0.5, t = 100;
        //init
        vector<vector3d> init_stack;
        vector3d temp1 = initVal(ray, init_stack);
        t = temp1.x * Scale; u = temp1.y; theta = temp1.z;
        if(t == 0)
        {
            return false;
        }
        //
        //we try the proposed init_values in the init_stack
        double result_t = 1e5;double result_u;double result_theta;
        //std::cout<<"init stack:"<<init_stack.size()<<"\n";
        if(init_stack.size() == 0)
        {
            std::cout<<"error!\n";
        }
        //std::cout<<init_stack.size()<<"\n";
        for(int j = 0; j < init_stack.size(); j++)
        {
            int Newton_num = 30;
            temp1 = init_stack[j];
            t = temp1.x; u = temp1.y; theta = temp1.z;
            for(int i = 0; i < Newton_num; i++)
            {
                vector3d current_F = F(ray, t, u, theta);
                double array1[] = {-current_F.x, -current_F.y, -current_F.z};
                MyMatrix<double> f(3, 1, array1);
                //std::cout<<f;
                vector2d partial_u = Curve.partial(u);
                vector2d p_u = Curve.val(u);
                double array2[] = {ray.Rd.x, -sin(theta) * partial_u.x, -cos(theta) * p_u.x, ray.Rd.y,-partial_u.y, 0, ray.Rd.z, -cos(theta) * partial_u.x, sin(theta) * p_u.x};
                MySquareMatrix<double> A(3, array2);
                //std::cout<<A;
                MyMatrix<double> delta_x(3, 1);
                delta_x = A.SolveEqu(f);
                t += delta_x.GetNumber(0, 0);
                u += delta_x.GetNumber(1, 0);
                theta += delta_x.GetNumber(2, 0);
                //std::cout<<"t:"<<t<<"\n";
                if(!(t < 1e6 && t > -1e6))
                {
                    break;
                }
            }
            if(F(ray, t, u, theta).l2norm() < 1e-2 && t * Scale > min_t && t < result_t && u <=1 && u >= 0)
            {
                result_t = t;
                result_u = u;
                result_theta = theta;
                if(fabs(theta - 3.14159265358979/2) < 0.2 && init_stack.size() <= 24)
                {
                    init_stack.push_back(vector3d(t-1, u, 3.14159265358979 - theta));
                }
                if(fabs(theta + 3.14159265358979/2) < 0.2 && init_stack.size() <= 24)
                {
                    init_stack.push_back(vector3d(t-1, u, -3.14159265358979 - theta));
                }
                if(fabs(theta - 3.14159265358979*1.5) < 0.2 && init_stack.size() <= 24)
                {
                    init_stack.push_back(vector3d(t-1, u, 3 * 3.14159265358979 - theta));
                }
            }
        }
        t = result_t; u = result_u; theta = result_theta;
        if(F(ray, t, u, theta).l2norm() > 1e-2 || t * Scale < min_t)
        {
            return false;
        }
        if(u > 1 || u < 0)
        {
            return false;
        }
        if(!(t > -1e5))
        {
            return false;
        }
        tt = t * Scale;
        vector2d temp = Curve.partial(u);
        vector2d temp2 = vector2d(-temp.y, temp.x);
        //nn = vector3d(Center.x/Scale + sin(theta) * temp2.x, Center.y/Scale + temp2.y, Center.z/Scale + cos(theta) * temp2.x);
        nn = vector3d(sin(theta) * temp2.x, temp2.y, cos(theta) * temp2.x);
        nn.normalize();
        //nn = vector3d(0, 1, 0);
        //Cautious! test code, remove later
        /*if(nn.dot(ray.Rd) > 0)
        {
            nn = vector3d(-nn.x, -nn.y, -nn.z);
        }*/
        //end Cautious!
        //std::cout<<"Intersect! t:"<<t<<"\n";
        return true;
    }
    bool GetCoefK(vector3d intersect_point, vector3d& kd, vector3d& ks, vector3d& ka) const
    {
        ka = Ka;
        kd = Kd;
        ks = Ks;
        //Cautious! testcode
        /*if(intersect_point.z < 80)
        {
            kd = vector3d(0.1, 0.1, 0.99);
        }*/
        //end Cautious
        return true;
    }
    bool IsInside(vector3d point) const
    {
        return false;
    }
};

#endif