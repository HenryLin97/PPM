#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include<fstream>
#include<iostream>
#include "Ray.h"
#include "Object.h"
#include "util2.h"
#include "Sphere.h"
#include "Plane.h"
#include "Triangle.h"
#include "StlModel.h"
#include "Bezier.h"
#include "BezierSurface.h"
#include "Cylinder.h"
#include <CImg.h>

/*Compile command: 
g++ -o main2 main2.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
*/
using namespace cimg_library;
//using namespace std;
int Width = 1024;
int Height = 768;

void ConstructGlass(vector3d position, vector3d kd, OBJMT mt, double scale = 7)
{
    Bezier b1;
    b1.add_control(vector2d(0.1, 2.1));
    b1.add_control(vector2d(0.1, 1.1));
    b1.add_control(vector2d(0.1, 0.1));
    b1.add_control(vector2d(1, 0.1));
    b1.add_control(vector2d(1, 0.1));
    b1.add_control(vector2d(1.05, 0.05));
    b1.add_control(vector2d(1, 0));
    BezierSurface* s11 = new BezierSurface(position, b1, kd, mt, scale);
    Bezier b2;
    b2.add_control(vector2d(0.85, 5.1));
    b2.add_control(vector2d(1, 3.1));
    b2.add_control(vector2d(0.1, 3.1));
    b2.add_control(vector2d(0.1, 2.1));
    BezierSurface* s12 = new BezierSurface(position, b2, kd, mt, scale);
    Bezier b3;
    /*b3.add_control(vector2d(0.75, 5.1));
    b3.add_control(vector2d(0.7, 4.1));
    b3.add_control(vector2d(0.5, 3.1));*/
    b3.add_control(vector2d(0.1, 3.1));
    b3.add_control(vector2d(0.5, 3.1));
    b3.add_control(vector2d(0.7, 4.1));
    b3.add_control(vector2d(0.75, 5.1));
    BezierSurface* s13 = new BezierSurface(position, b3, kd, mt, scale);
    Bezier b5;
    /*b5.add_control(vector2d(0.85, 5.1));
    b5.add_control(vector2d(0.8, 5.2));*/
    b5.add_control(vector2d(0.75, 5.1));
    b5.add_control(vector2d(0.8, 5.2));
    b5.add_control(vector2d(0.85, 5.1));
    BezierSurface* s15 = new BezierSurface(position, b5, kd, mt, scale);
    ObjectQueue.push_back((Object*)s11);
    ObjectQueue.push_back((Object*)s12);
    ObjectQueue.push_back((Object*)s13);
    ObjectQueue.push_back((Object*)s15);
}

//begin
void ConstructScene()
{
    //Sphere* s1 = new Sphere(1e5, vector3d( 1e5+1,40.8,81.6), vector3d(.75,.25,.25),MT_DIFF);//left
    Plane* s1 = new Plane(vector3d(1, 40.8, 81.6), vector3d(1, 0, 0), vector3d(.95,.15,.15));
    //Sphere* s2 = new Sphere(1e5, vector3d(-1e5+99,40.8,81.6),vector3d(.25,.25,.75),MT_DIFF);//Right
	Plane* s2 = new Plane(vector3d(99,40.8,81.6), vector3d(-1, 0, 0), vector3d(.25,.25,.75));
	//cv::Mat texture = cv::imread("/home/lin/Pictures/timg.jpeg", CV_LOAD_IMAGE_COLOR);
	CImg<unsigned char>* texture = new CImg<unsigned char>("/home/lin/Pictures/index4.jpg");
    //s2->loadTexture(texture, vector3d(0, 1, 0), vector3d(0, 0, 1), 0.05);
    //Sphere* s3 = new Sphere(1e5, vector3d(50,40.8, 1e5),     vector3d(.75,.75,.75),MT_DIFF);//Back
	Plane* s3 = new Plane(vector3d(50,40.8,0), vector3d(0, 0, 1), vector3d(.75,.75,.75));
    //Sphere* s4 = new Sphere(1e5, vector3d(50,40.8,-1e5+170), vector3d(),           MT_DIFF);//Front
	Plane* s4 = new Plane(vector3d(50,40.8,170), vector3d(0, 0, -1), vector3d(0,0,0));
    //Sphere* s5 = new Sphere(1e5, vector3d(50, 1e5, 81.6),    vector3d(.75,.75,.75),MT_DIFF);//Bottomm
	Plane* s5 = new Plane(vector3d(40, 0, 71.6), vector3d(0, 1, 0), vector3d(.25,.75,.25));
    s5->loadTexture(texture, vector3d(1, 0, 0), vector3d(0, 0, 1), 0.1);
    //Sphere* s6 = new Sphere(1e5, vector3d(50,-1e5+81.6,81.6),vector3d(.75,.75,.75),MT_DIFF);//Top
	Plane* s6 = new Plane(vector3d(50, 81.6, 81.6), vector3d(0, -1, 0), vector3d(.75,.75,.75));
    Sphere* s7 = new Sphere(16.5,vector3d(27,16.5,47),       vector3d(1,1,1)*.999, MT_SPEC);//Mirror
    Sphere* s8 = new Sphere(8,vector3d(85,8,90),       vector3d(1,1,1)*.999, MT_REFR);//Glass
    Sphere* s9 = new Sphere(8.5, vector3d(50,8.5,60),        vector3d(1,1,1)*.999, MT_DIFF);//Middle
	Triangle* s10 = new Triangle(vector3d(73, 40, 88), vector3d(50, 50, 60), vector3d(27, 16.5, 47), vector3d(.25, .25, .75), MT_DIFF);
    Cylinder* c12 = new Cylinder(vector3d(70, 0, 80), 20, 10, vector3d(0.99, 0.99, 0.99), MT_DIFF);

    StlModel* s13 = new StlModel("dragon.obj.stl", vector3d(.99, .99, .1), MT_DIFF);
    //s13->Rotate(1, 1);
	s13->Scale(4);
	s13->Translate(vector3d(70, 0, 90));
    //s13->ConstructBoundingBox();
    s13->ConstructKdTree();
    ObjectQueue.push_back((Object*)s1);
    ObjectQueue.push_back((Object*)s2);
    ObjectQueue.push_back((Object*)s3);
    ObjectQueue.push_back((Object*)s4);
    ObjectQueue.push_back((Object*)s5);
    ObjectQueue.push_back((Object*)s6);
    ObjectQueue.push_back((Object*)s7);
    //ObjectQueue.push_back((Object*)s8);
    /*ObjectQueue.push_back((Object*)s9);
	ObjectQueue.push_back((Object*)s10);*/
	//ObjectQueue.push_back((Object*)s11);
    ObjectQueue.push_back((Object*)s13);
    //ConstructGlass(vector3d(50, 0, 120), vector3d(0.95, 0.95, 0.95), MT_REFR, 10);
}

int main(int argc, char *argv[]) {
	//cimg_library::CImg<float> image("/home/lin/Pictures/timg.jpeg");
    ConstructScene();
    srand(time(0));
	//cv::Mat texture = cv::imread("/home/lin/Pictures/timg.jpeg", CV_LOAD_IMAGE_COLOR);
	//texture.release();
    //hitpoints = NULL;
	// samps * 1000 photon paths will be traced
	int w=Width, h=Height;
    int samps = (argc==2) ? MAX(atoi(argv[1])/1000,1) : 1000;

	// trace eye rays and store measurement points
	Ray cam(vector3d(50,48,355.6), vector3d(0,-0.042612,-1).norm());
	vector3d cx=vector3d(w*.5135/h, 0, 0), cy = (cx%cam.Rd).norm()*.5135, *c=new vector3d[w*h], vw;
    int lost_count = 0;
    //Add depth here!
    Plane FocusPlane(vector3d(50,40.8,90), vector3d(0, 0, 1), vector3d(.75,.75,.75));
    int sample_num = 10;
    double focus_z = 169;
    //end depth
	for (int y=0; y<h; y++){
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y/(h-1));
		for (int x=0; x<w; x++) {
			pixel_index = x + y * w;
			vector3d d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5)+cam.Rd;
            double temp_t;
            vector3d temp_n;
            vector3d new_direction = vector3d(50, 48, focus_z) - (cam.R0 + d*140);
            new_direction.normalize();
            Ray ray1(cam.R0 + d * 140, new_direction);
            FocusPlane.intersect(ray1, 0.01, temp_t, temp_n);
            vector3d intersect_point = ray1.R0 + ray1.Rd * temp_t;
            //std::cout<<"d: "<<d.norm();
            double range = 5;
            for(int i = 0; i < sample_num; i++)
            {
                double xx = (random_num()-0.5)*range + 50;
                double yy = (random_num()-0.5)*range + 48;
                vector3d new_start = vector3d(xx, yy, focus_z);
                vector3d direction = (intersect_point - new_start).norm();
                //std::cout<<"direction: "<<direction;
			    trace(Ray(new_start, direction), 0, true, vector3d(), vector3d(1, 1, 1),0);
            }
            //debug code: remove later
            if(pixel_record[pixel_index] == 0)
            {
                std::cout<<"strange...hit point lost\n";
                trace(Ray(cam.R0 + d * 140, d.norm()), 0, true, vector3d(), vector3d(1, 1, 1),0, true);
                lost_count++;
            }
            //std::cout<<"here!\n";
		}
	}
    std::cout<<"total lost:"<<lost_count<<"\n";
	fprintf(stderr,"\n");
    //debug code; remove later
    /*for(int i = 0; i < 786432; i++)
    {
        if(pixel_record[i] == 0)
        {
            std::cout<<"strange...hit point lost\n";
        }
    }*/
	
	// build the hash table over the measurement points
	build_hash_grid(w,h); 
	
	// trace photon rays with multi-threading
	num_photon=samps; 
	vw=vector3d(1,1,1);
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0;i<num_photon;i++) {
		double p=100.*(i+1)/num_photon;
		fprintf(stderr,"\rPhotonPass %5.2f%%",p); 
		int m=1000*i; 
		Ray r; 
		vector3d f;
		for(int j=0;j<1000;j++){
            //std::cout<<j<<"\n";
			genp(&r,&f,m+j); 
			trace(r,0,0>1,f,vw,m+j);
		}
	}

	// density estimation
	for(auto item: hitpoints)
	{
		int i=item->pix;
		c[i]=c[i]+item->flux*(1.0/(PI*item->r2*num_photon*1000.0));
	}

	// save the image after tone mapping and gamma correction
	FILE* f = fopen("image10.ppm","w"); fprintf(f,"P3\n%d %d\n%d\n",w,h,255);
	for(int i = w * h - 1; i>=0; i--) 
	{
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	}
	fclose(f);

	return 0;
}