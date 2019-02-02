#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Ray.h"
#include "Object.h"
#define PI ((double)3.14159265358979)
#define ALPHA ((double)0.7) // the alpha parameter
#define MAX(x, y) ((x > y) ? x : y)
#define MIN(x, y) ((x < y) ? x : y)
#define MAX_DEPTH 8 //the max trace depth

std::vector<Object*> ObjectQueue; //this include all objects in the scene
vector3d point_light_loc(50,60,85);//this defines the location of the light
int light_type = 0;//set this number to change to different types of light
//0 for normal point light; 1 for circle light
double light_flux = 300;//set this to control the light or dark of the scene

int pixel_record[20786432];

//generate a random number in [0, 1]
double random_num()
{
	return (double)rand()/(double)(RAND_MAX);
}

struct ABBox{
    vector3d min;
    vector3d max;
    void fit(vector3d p)
    {
		min = element_min(p, min);
		max = element_max(p, max);
    }
    void reset()
    {
        min = vector3d(1e10, 1e10, 1e10);
        max = vector3d(-1e10, -1e10, -1e10);
    }
};

struct HPoint {
	vector3d f,pos,nrm,flux; 
	double r2; 
	unsigned int n; // n = N / ALPHA in the paper
	int pix;
};

unsigned int num_hash, pixel_index, num_photon;
double hash_s;
std::vector<HPoint*>* hash_grid;
std::vector<HPoint*> hitpoints;
ABBox hpbbox;

// spatial hash function
inline unsigned int my_hash(const int ix, const int iy, const int iz) {
	return (unsigned int)((ix*73856093)^(iy*19349663)^(iz*83492791))%num_hash;
}

void build_hash_grid(const int w, const int h) {
	// find the bounding box of all the measurement points
	hpbbox.reset();
	
	for(auto item: hitpoints)
	{
		hpbbox.fit(item->pos);
	}

	// heuristic for initial radius
	vector3d ssize = hpbbox.max - hpbbox.min;
	double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;
    //std::cout<<"irad:"<<irad<<"\n";
    /*char c;
    std::cin>>c;*/
	// determine hash table size
	// we now find the bounding box of all the measurement points inflated by the initial radius
	hpbbox.reset(); 
	int vphoton = 0;
	for(auto item: hitpoints)
	{
		item->r2=irad * irad; 
		item->n = 0;
		item->flux = vector3d(0, 0, 0);
		hpbbox.fit(item->pos-irad); 
		hpbbox.fit(item->pos+irad);
		vphoton++;
	}

	// make each grid cell two times larger than the initial radius
	hash_s=1.0/(irad*2.0); 
	num_hash = vphoton; 
    //std::cout<<"num_hash:"<<num_hash<<"\n";
    /*std::cin>>c;*/
	// build the hash table
	hash_grid=new std::vector<HPoint*>[num_hash];
	for (int i=0; i<num_hash;i++)
		hash_grid[i].clear();
	for(auto item: hitpoints)
	{ 
		HPoint *hp = item;
		vector3d BMin = ((hp->pos - irad) - hpbbox.min) * hash_s;
		vector3d BMax = ((hp->pos + irad) - hpbbox.min) * hash_s;
		for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
		{
			for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
			{
				for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
				{
					int hv=my_hash(ix,iy,iz);
                    //std::cout<<ix<<" "<<iy<<" "<<iz<<"\n";
                    //std::cout<<"hv:"<<hv<<"\n";
					hash_grid[hv].push_back(hp);
				}
			}
		}
	}
}

// tone mapping and gamma correction
int toInt(double x){
	return int(pow(1-exp(-x),1/2.2)*255+.5);
}

// find the closet interection
inline bool intersect(const Ray &ray,double &t,vector3d &normal, int &id){
    
	//determine the intersect point
    int ObjNum = ObjectQueue.size();
    //std::cout<<"ObjNum:"<<ObjNum<<"\n";
    double intersect_t = -1;
    vector3d intersect_normal;
    int intersect_index = -1;
    for(int i = 0; i < ObjNum; i++)
    {
        double temp_t = -1;
        vector3d temp_n;
        if(ObjectQueue[i]->intersect(ray, 1e-4, temp_t, temp_n))//there is an intersect point
        {
            if(temp_t < intersect_t || std::abs(intersect_t + 1) < 1e-5)
            {
                intersect_t = temp_t;
                intersect_normal = temp_n;
                normal = temp_n;
                intersect_index = i;
                t = temp_t;
                id = i;
            }
        }
    }
    if(intersect_index == -1)
    {
        return false;
    }
    else
    {
        //std::cout<<"true\n";
        return true;
    }
}

// generate a photon ray from the point light source
void genp(Ray* pr, vector3d* f, int i) {
    vector3d tmp(light_flux, light_flux, light_flux);
	*f = tmp*(PI*4.0);
	pr->R0 = point_light_loc;
	double r1 = random_num();
	double r2 = random_num();
	double theta = 2 * PI * r1;
	double Sin1 = sin(theta), Cos1 = cos(theta);
	if(light_type == 1)
	{
		r2 /= 10;
	}
	double Sin2 = 2 * sqrt((1 - r2)*r2);
	double Cos2 = 1 - 2*r2;
	pr->Rd = vector3d(Cos1 * Sin2, Cos2, Sin1 * Sin2);
}

void trace(const Ray &r,int dpt,bool m,const vector3d &fl,const vector3d &adj,int i, bool debug = false) 
{
	double t;
	int id;
    //std::cout<<"dpt:"<<dpt<<"\n";
	dpt++;//increase the depth
    vector3d n;//this is the normal vector of intersect point
	if(!intersect(r,t, n, id)||(dpt>=MAX_DEPTH))//if no intersect point or depth reaches MAX_DEPTH, stop tracing
	{
		return;
	}
	Object* obj = ObjectQueue[id]; 
	vector3d x=r.R0+r.Rd*t;//x is the intersect point
	if(debug)
	{
		std::cout<<x;
	}
	vector3d f, f1, f2;//f is the color vector, which is given by Kd
	obj->GetCoefK(x, f, f1, f2);
	vector3d adjust_normal = n;//this is the normal direction opposed to ray direction
	if(n.dot(r.Rd) >= 0)
	{
		adjust_normal = n*(-1);
	}
	double p=f.x>f.y&&f.x>f.z?f.x:f.y>f.z?f.y:f.z;

    /*std::cout<<"hit_points:"<<hitpoints<<"\n";
    char c;
    std::cin>>c;*/
	if (obj->MT == MT_DIFF) {

		double r1=2.*PI*random_num(),r2=random_num();
		double n_Sin = sqrt(r2);
		double n_Cos = sqrt(1 - r2);
		vector3d w = adjust_normal, u, v;//use 
		if(fabs(w.x) > 0.1)//w has x component
		{
			u = vector3d(0, 1, 0) % w;
			u.normalize();
		}
		else
		{
			u = vector3d(1, 0, 0) % w;
			u.normalize();
		}
		v = w % u; v.normalize();
		vector3d d = (u*cos(r1) + v*sin(r1))*n_Sin + w*n_Cos;
		d.normalize(); 

		if (m) {
			// eye ray
			// store the measurment point
			HPoint* hp=new HPoint; 
			hp->f=f.mul(adj); 
			hp->pos=x;
			hp->nrm=n; 
			hp->pix = pixel_index;
			hitpoints.push_back(hp);
			pixel_record[pixel_index]++;
		} 
		else 
		{
			vector3d hh = (x-hpbbox.min) * hash_s;
			int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
			{
				std::vector<HPoint*>& grid_list = hash_grid[my_hash(ix, iy, iz)];
                int count = 0;
				for(auto hit_point: grid_list)
				{
                    count++; 
					vector3d v = hit_point->pos - x;
					// check normals to be closer than 90 degree (avoids some edge brightning)
					//update the data in hitpoint
					if ((hit_point->nrm.dot(n) > 1e-3) && (v.dot(v) <= hit_point->r2) && hit_point->n < 4294967295) {
						// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
						double g = (hit_point->n*ALPHA+ALPHA) / (hit_point->n*ALPHA+1.0);
						hit_point->r2=hit_point->r2*g; 
						hit_point->n++;
						hit_point->flux=(hit_point->flux+hit_point->f.mul(fl)*(1./PI))*g;
					}
				}
			}
			if (random_num()<p) 
				trace(Ray(x,d),dpt,m,f.mul(fl)*(1./p),adj,i);
		}

	} else if (obj->MT == MT_SPEC) {
		// mirror
		trace(Ray(x, r.Rd-n*2.0*n.dot(r.Rd)), dpt, m, f.mul(fl), f.mul(adj),i);

	} else {
		// glass
		Ray lr(x,r.Rd-n*2.0*n.dot(r.Rd));//this is the spec ray
		bool into = (n.dot(adjust_normal)>-1e7);//use the direction of normal vector to judge into or out of
		double nc = 1.0, nt=1.5, nnt;
		if(into)
		{
			nnt = nc/nt;
		}
		else
		{
			nnt = nt/nc;
		}
		double ddn = r.Rd.dot(adjust_normal);
		double cos2t = 1-nnt*nnt*(1-ddn*ddn);
		// total internal reflection
		if (cos2t<0)
			return trace(lr,dpt,m,fl,adj,i, debug);
		
		double a=nt-nc, b=nt+nc, R0=a*a/(b*b);
		vector3d td;//refract direction
		if(into)
		{
			td = r.Rd * nnt - n*(ddn*nnt + sqrt(cos2t));
			td.normalize();
		}
		else
		{
			td = r.Rd * nnt + n*(ddn*nnt + sqrt(cos2t));
			td.normalize();
		}
		double c;
		if(into)
		{
			c = 1 + ddn;
		}
		else
		{
			c = 1 - td.dot(n);
		}
		double Re=R0+(1-R0)*c*c*c*c*c;
		Ray rr(x,td);//refract ray
		vector3d fa=f.mul(adj);
		if (m) {
			// eye ray (trace both rays)
			trace(lr,dpt,m,fl,fa*Re,i);
			trace(rr,dpt,m,fl,fa*(1.0-Re),i);
		} else {
			//random trace photon ray
			if(random_num() < Re)
			{
				trace(lr,dpt,m,fl,fa,i);
			}
			else
			{
				trace(rr,dpt,m,fl,fa,i);
			}
		}
	}
}
//end

#endif