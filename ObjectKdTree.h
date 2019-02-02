#ifndef OBJECTKDTREE_H
#define OBJECTKDTREE_H

#include"Ray.h"
#include"Triangle.h"
#include"Box.h"
#include<vector>

struct TreeNode{
    vector3d Min, Max;
    std::vector<Triangle>* Faces;
    TreeNode *Ls, *Rs;
    Box* BoundingBox;
    TreeNode(vector3d min, vector3d max)
    {
        Min = min;
        Max = max;
        BoundingBox = new Box(Max, Min);
    }
    bool inside(Triangle face)
    {
        vector3d v_min = element_min(face.P[0], face.P[1]);
        v_min = element_min(v_min, face.P[2]);
        vector3d v_max = element_max(face.P[0], face.P[1]);
        v_max = element_max(v_max, face.P[2]);
        if(v_min.x <= Max.x && v_min.y <= Max.y && v_min.z <= Max.z && v_max.x >= Min.x && v_max.y >= Min.y && v_max.z >= Min.z)
        {
            return true;
        }
        return false;
    }
};

class ObjectKdTree{
public:
TreeNode* root;
const int Max_Depth = 15;
const int Max_Faces = 10;
ObjectKdTree(std::vector<Triangle>* faces)
{
    vector3d upper_corner = (*faces)[0].P[0];
    vector3d lower_corner = upper_corner;
    for(int i = 0; i < faces->size(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            upper_corner = element_max(upper_corner, (*faces)[i].P[j]);
            lower_corner = element_min(lower_corner, (*faces)[i].P[j]);
        }
    }
    root = BuildTree(0, 0, faces, lower_corner, upper_corner);
}
TreeNode* BuildTree(int depth, int d, std::vector<Triangle>* faces, vector3d min, vector3d max)
{
    TreeNode *p = new TreeNode(min, max);
    vector3d maxL, minR;
    d %= 3;//this is the dimension we want to divide
    //now select control point of left son and right son
    if (d == 0) {
        maxL = vector3d((p->Min.x + p->Max.x) / 2, p->Max.y, p->Max.z);
        minR = vector3d((p->Min.x + p->Max.x) / 2, p->Min.y, p->Min.z);
    }
    else if (d == 1) {
        maxL = vector3d(p->Max.x, (p->Min.y + p->Max.y) / 2, p->Max.z);
        minR = vector3d(p->Min.x, (p->Min.y + p->Max.y) / 2, p->Min.z);
    }
    else {
        maxL = vector3d(p->Max.x, p->Max.y, (p->Min.z + p->Max.z) / 2);
        minR = vector3d(p->Min.x, p->Min.y, (p->Min.z + p->Max.z) / 2);
    }
    //select the faces we want
    p->Faces = new std::vector<Triangle>;
    p->Faces->clear();
    for(Triangle face: *faces)
    {
        if(p->inside(face))
        {
            p->Faces->push_back(face);
        }
    }
    int face_num = p->Faces->size();
    //construct son nodes
    if(depth < Max_Depth && face_num > Max_Faces)
    {
        p->Ls = BuildTree(depth+1, d+1, p->Faces, min, maxL);
        p->Rs = BuildTree(depth+1, d+1, p->Faces, minR, max);
    }
    else
    {
        p->Ls = NULL; p->Rs = NULL;
    }
    return p;
}
int intersect(const Ray ray, const double min_t, double& t, vector3d& n) const
{
    return intersect(ray, min_t, t, n, root);
}
int intersect(const Ray ray, const double min_t, double& t, vector3d& n, TreeNode* p) const
{
    if(p->Ls == NULL || p->Rs == NULL)//no child node is find
    {
        bool flag = false;
        int count = 0;
        double current_min = 1e5;
        for(auto face: *(p->Faces))
        {
            double temp_t;
            vector3d temp_n;
            if(face.intersect(ray, min_t, temp_t, temp_n))
            {
                if(p->BoundingBox->IsInside(ray.R0 + ray.Rd * temp_t))
                {
                    count++;
                }
                if(temp_t < current_min)
                {
                    t = temp_t;
                    n = temp_n;
                    current_min = temp_t;
                    flag = true;
                }
            }
        }
        return count;
    }
    else
    {
        double temp_t1, temp_t2;
        vector3d temp_n1, temp_n2;
        bool b1 = p->Ls->BoundingBox->intersect(ray, min_t, temp_t1, temp_n1);
        bool b2 = p->Rs->BoundingBox->intersect(ray, min_t, temp_t2, temp_n2);
        if(b1 && b2)//intersect two sons
        {
            int c1, c2;
            c1 = intersect(ray, min_t, temp_t1, temp_n1, p->Ls);
            c2 = intersect(ray, min_t, temp_t2, temp_n2, p->Rs);
            if(c1 > 0 && c2 > 0)
            {
                if(temp_t1 <= temp_t2)
                {
                    t = temp_t1;
                    n = temp_n1;
                }
                else
                {
                    t = temp_t2;
                    n = temp_n2;
                }
                return c1 + c2;
            }
            else if(c1 > 0)
            {
                t = temp_t1;
                n = temp_n1;
                return c1;
            }
            else if(c2 > 0)
            {
                t = temp_t2;
                n = temp_n2;
                return c2;
            }
            else
            {
                return 0;
            }
        }
        else if(b1)
        {
            return intersect(ray, min_t, t, n, p->Ls);
        }
        else if(b2)
        {
            return intersect(ray, min_t, t, n, p->Rs);
        }
        else
        {
            return 0;
        }
    }
}
};

#endif