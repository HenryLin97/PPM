#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include<fstream>
#include<iostream>
#include "Ray.h"
#include "Object.h"
#include "util.h"
#include "Sphere.h"
#include "Plane.h"
#include "Triangle.h"
#include "StlModel.h"
#include "Bezier.h"
#include "BezierSurface.h"
#include <CImg.h>

using namespace std;

int main()
{
    Bezier b1;
    b1.add_control(vector2d(0, 0));
    b1.add_control(vector2d(1, 1));
    cout<<b1.control_num<<endl;
    cout<<b1.B(0, 1, 0.5)<<endl;
    cout<<C(1, 0)<<" "<<C(1, 1)<<endl;
    cout<<power(0.5, 0)<<" "<<power(0.5, 1)<<endl;
    cout<<b1.val(0.5).x<<" "<<b1.val(0.5).y<<endl;
    BezierSurface* s11 = new BezierSurface(vector3d(0, 0, 0), b1, vector3d(0, 0, 0), MT_DIFF);
    Ray R;
    R.R0 = vector3d(-100, 0.5, 0);
    R.Rd = vector3d(1, 0, 0);
    double t;
    vector3d n;
    s11->intersect(R, 0, t, n);
    cout<<t<<endl;
    cout<<n.x<<" "<<n.y<<" "<<n.z<<"\n";
    return 0;
}