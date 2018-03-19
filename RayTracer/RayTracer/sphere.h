#ifndef SPHERE_H
#define SPHERE_H

#include "hitable.h"

class sphere : public hitable {
	/*: public hitable表示sphere继承hitable。即：sphere为hitable的子类*/
public:
	sphere() {}
	sphere(vec3 cen, float r) : center(cen), radius(r) {}
	/*此处为使用初始化列表的构造函数来初始化成员变量*/
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	/*必须实现父类的虚函数。*/
	vec3 center;
	float radius;
};
bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	vec3 oc = r.orgin() - center;
	//得到向量A-C=oc
	float a = dot(r.direction(), r.direction());
	//B即为光线的方向向量。a = (B·B)=dot(B,B)
	float b = 2.0 * dot(oc, r.direction());
	// b=2*(B·(A-C))=2*dot(B,(A-C))
	float c = dot(oc, oc) - radius*radius;
	//c= (A-C)·(A-C) - R2=dot((A-C),(A-C)) - R2
	float discriminant = b*b - 4 * a*c;
	//判别式>0，有实根，光线撞上球；
	//判别式=0，有1个实根，光线撞上球，但是没有折射；
	//判别式<0，没有实根，光线没有撞上球；
	if (discriminant > 0) {
		float temp = (-b - sqrt(discriminant)) / (2.0*a);//小根
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
		temp = (-b + sqrt(discriminant)) / (2.0*a);//大根
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
	}
	/*判断小根和大根是否在范围内。首先判断小根是否在范围内。是：保存相关信息，然后直接返回；否，判断大根是否在范围内。也就是优先选小根，小根不行再考虑大根 */
	return false;
}
#endif // !SPHERE_H
#pragma once
