#ifndef STLMODEL_H
#define STLMODEL_H

#include"Object.h"
#include"Ray.h"
#include"Triangle.h"
#include"Box.h"
#include"ObjectKdTree.h"
#include<fstream>

class StlModel: public Object{
public:
    std::vector<Triangle> TriangleList;
    int TriangleNum;
    bool has_Bounding_Box;
    bool has_kd_tree;
    ObjectKdTree* kd_tree;
    Box* Bounding_Box;
    StlModel(const char* filename, vector3d KD, OBJMT mt = MT_DIFF)
    {
        has_Bounding_Box = false;
        has_kd_tree = false;
        TriangleNum = LoadAsciiStl(filename, KD, mt);
        Kd = KD;
        Ka = vector3d(0, 0, 0);
        Ks = vector3d(0, 0, 0);
        ReflectCoe = vector3d(0, 0, 0);
        RefractCoe = vector3d(0, 0, 0);
        RefractRatio = 1;
        MT = mt;
        std::cout<<TriangleNum<<"\n";
    }
    bool intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
    {
        if(has_kd_tree)
        {
            double temp_t;
            vector3d temp_n;
            int c = kd_tree->intersect(ray, min_t, temp_t, temp_n);
            if(c == 0)
            {
                return false;
            }
            else
            {
                t = temp_t;
                n = temp_n;
                if(c %2 == 0)//ray into
                {
                    if(n.dot(ray.Rd) > 0)
                    {
                        n = vector3d(-n.x, -n.y, -n.z);
                    }
                }
                else//ray out
                {
                    if(n.dot(ray.Rd) < 0)
                    {
                        n = vector3d(-n.x, -n.y, -n.z);
                    }
                }
                return true;
            }
        }
        if(has_Bounding_Box)
        {
            double t1;
            vector3d n1;
            bool pre_test = Bounding_Box->intersect(ray, min_t, t1, n1);
            if(pre_test == false)
            {
                return false;
            }
        }
        double current_min = 1e5;
        bool flag = false;
        int count = 0;
        for(int i = 0; i < TriangleNum; i++)
        {
            double temp_t;
            vector3d temp_n;
            if(TriangleList[i].intersect(ray, min_t, temp_t, temp_n))
            {
                count++;
                if(temp_t < current_min)
                {
                    t = temp_t;
                    n = temp_n;
                    current_min = temp_t;
                    flag = true;
                    //std::cout<<"lower_z: "<<TriangleList[i].P[2].z<<"\n";
                }
            }
        }
        if(count %2 == 0)//ray into
        {
            if(n.dot(ray.Rd) > 0)
            {
                n = vector3d(-n.x, -n.y, -n.z);
            }
        }
        else//ray out
        {
            if(n.dot(ray.Rd) < 0)
            {
                n = vector3d(-n.x, -n.y, -n.z);
            }
        }
        return flag;
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
    void Scale(double ratio = 1)
    {
        for(int i = 0; i < TriangleNum; i++)
        {
            TriangleList[i].Scale(ratio);
        }
    }
    void Translate(vector3d loc)
    {
        for(int i = 0; i < TriangleNum; i++)
        {
            TriangleList[i].Translate(loc);
        }
    }
    void Rotate(int axis, double theta)
    {
        double Cos = cos(theta);
        double Sin = sin(theta);
        for(int i = 0; i < TriangleNum; i++)
        {
            TriangleList[i].Rotate(axis, theta, Cos, Sin);
        }
    }
    int LoadAsciiStl(const char* filename, vector3d KD, OBJMT mt);
    void ConstructBoundingBox()
    {
        vector3d upper_corner = TriangleList[0].P[0];
        vector3d lower_corner = upper_corner;
        for(int i = 0; i < TriangleNum; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                upper_corner = element_max(upper_corner, TriangleList[i].P[j]);
                lower_corner = element_min(lower_corner, TriangleList[i].P[j]);
            }
        }
        //std::cout<<"TriangleNum:"<<TriangleNum<<" lower_z:"<<lower_corner.z<<"\n";
        Bounding_Box = new Box(upper_corner, lower_corner);
        has_Bounding_Box = true;
    }
    void ConstructKdTree()
    {
        kd_tree = new ObjectKdTree(&TriangleList);
        has_kd_tree = true;
    }
};

int StlModel::LoadAsciiStl(const char* filename, vector3d KD, OBJMT mt)
{
    int i=0,j=0,cnt=0 ,pCnt=4;
    char a[100];
    char str[100];
    double x=0,y=0,z=0;
    vector3d tPoint;
    std::vector<vector3d> pointList;
 
    std::ifstream in;
    in.open(filename, std::ios::in);
	if (!in) 
	{ 
		return 0; 
	} 
	do 
	{ 
		i=0; 
		cnt=0; 
		in.getline(a,100, '\n'); 
		while(a[i]!='\0') 
		{ 
			if (!islower((int)a[i]) && !isupper((int)a[i]) && a[i]!=' ') 
			  break; 
			cnt++; 
			i++; 
		} 
		
		while(a[cnt]!='\0') 		
		{ 
			str[j]=a[cnt]; 
			cnt++; 
			j++; 
		} 
		str[j]='\0'; 
		j=0; 
		
		if (sscanf(str,"%lf%lf%lf",&x,&y,&z)==3) 
		{ 
			tPoint = {x, y, z};
            pointList.push_back(tPoint);
        }
        pCnt++;
    }while(!in.eof());

    int length = pointList.size();
    TriangleList.clear();
    for(int i = 0; i < length/4; i++)
    {
        Triangle temp(pointList[4*i], pointList[4*i + 1], pointList[4*i + 2], pointList[4*i + 3], KD, mt);
        TriangleList.push_back(temp);
    }
 
    return length/4;
}

#endif