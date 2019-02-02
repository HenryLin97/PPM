#ifndef BEZIER_H
#define BEZIER_H

#include<vector>
#include"Ray.h"
#include<math.h>

struct vector2d{
    double x;
    double y;
    vector2d(double xx, double yy)
    {
        x = xx; y = yy;
    }
};

int C(int n, int i)
{
    if(i < 0)
    {
        return 0;
    }
    int prod = 1;
    for(int k = 0; k < i; k++)
    {
        prod *= (n-k);
    }
    for(int k = 1; k < i+1; k++)
    {
        prod /= k;
    }
    return prod;
}

double power(double u, int i)
{
    double prod = 1;
    for(int k = 0; k < i; k++)
    {
        prod *= u;
    }
    return prod;
}

class Bezier{
public:
    int control_num;
    std::vector<vector2d> control_points;
    Bezier()
    {
        control_num = 0;
        control_points.clear();
    }
    void add_control(vector2d new_point)
    {
        control_points.push_back(new_point);
        control_num ++;
    }
    double B(int i, int n, double u) const
    {
        return C(n, i) * power(1-u, n-i) * power(u, i);
    }
    double Bp(int i, int n, double u) const
    {
        return n*(B(i-1, n-1, u) - B(i, n-1, u));
    }
    vector2d val(double u) const
    {
        vector2d result = vector2d(0, 0);
        int n = control_num - 1;
        for(int i = 0; i < control_num; i++)
        {
            vector2d P = control_points[i];
            double temp = B(i, n, u);
            result.x += (P.x * temp);
            result.y += (P.y * temp);
        }
        return result;
    }
    vector2d partial(double u) const
    {
        vector2d result = vector2d(0, 0);
        int n = control_num - 1;
        for(int i = 0; i < control_num; i++)
        {
            vector2d P = control_points[i];
            double temp = Bp(i, n, u);
            result.x += (P.x * temp);
            result.y += (P.y * temp);
        }
        return result;
    }
    double max_x() const
    {
        double result = -1e5;
        int num_step = 10;
        for(int i = 0; i < num_step; i++)
        {
            double temp = val((double)i/(double)num_step).x;
            if(temp > result)
            {
                result = temp;
            }
        }
        return result;
    }
    vector2d range_y() const
    {
        double temp1 = control_points[0].y;
        double temp2 = control_points[control_num - 1].y;
        if(temp1 > temp2)
        {
            return vector2d(temp2, temp1);
        }
        else
        {
            return vector2d(temp1, temp2);
        }
    }
};

#endif