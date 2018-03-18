#ifndef HITABLE_H
#define HITABLE_H
#include "ray.h"
struct hit_record {
	float t;
	vec3 p;
	vec3 normal;
};/*该结构体记录“撞点”处的信息：离光线起点的距离t、撞点的坐标向量p、撞点出的法向量normal。*/ 

class hitable {
public:/*hitable这个类表示能够被光线撞上的任何物体。比如，球体*/ 
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
	/*hit()在此被声明为虚函数，则hitable为抽象类。抽象类的子类中必须实现其虚函数*/
};
#endif // HITABLE_H
