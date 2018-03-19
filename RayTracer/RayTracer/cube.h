
#ifndef CUBE_H
#define CUBE_H
#include "hitable.h"
class cube : public hitable{
	/*: public hitable表示cube继承hitable。即：cube为hitable的子类*/
public:
	cube() {}
	cube(vec3 position, vec3 r_x,vec3 r_y,vec3 r_z) : c_point(position), c_x(r_x),c_y(r_y),c_z(r_z) {}
	/*此处为使用初始化列表的构造函数来初始化成员变量*/
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	/*必须实现父类的虚函数。*/
	vec3 c_point, c_x, c_y, c_z;
};
bool cube::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	float x = c_x.length();
	float y = c_y.length();
	float z = c_z.length();
	float min_d = -1;
	vec3 N = vec3(0, 0, 0);

	vec3 oc = c_point - r.orgin();
	float denom;
	denom = dot(r.direction(), c_z);
	if (denom != 0) {
		float d = dot(oc, c_z) / denom;
		if (d > 0) {
			vec3 point = r.orgin() + d * r.direction();
			vec3 vex = point - c_point;
			float dis1 = dot(vex, c_x) / x;
			float dis2 = dot(vex, c_y) / y;
			if (dis1 > 0 && dis1 < x && dis2 > 0 && dis2 < y) {
				if (min_d != -1) {
					if (d < min_d) {
						min_d = d;
						N = unit_vector(-c_z);
					}
				}
				else {
					min_d = d;
					N = unit_vector(-c_z);
				}
			}
		}
	}

	oc = c_point + c_z - r.orgin();
	denom = dot(r.direction(), c_z);
	if (denom != 0) {
		float d = dot(oc, c_z) / denom;
		if (d > 0) {
			vec3 point = r.orgin() + d * r.direction();
			vec3 vex = point - c_point;
			float dis1 = dot(vex, c_x) / x;
			float dis2 = dot(vex, c_y) / y;
			if (dis1 > 0 && dis1 < x && dis2 > 0 && dis2 < y) {
				if (min_d != -1) {
					if (d < min_d) {
						min_d = d;
						N = unit_vector(c_z);
					}
				}
				else {
					min_d = d;
					N = unit_vector(c_z);
				}
			}
		}
	}

	oc = c_point - r.orgin();
	denom = dot(r.direction(), c_y);
	if (denom != 0) {
		float d = dot(oc, c_y) / denom;
		if (d > 0) {
			vec3 point = r.orgin() + d * r.direction();
			vec3 vex = point - c_point;
			float dis1 = dot(vex, c_x) / x;
			float dis2 = dot(vex, c_z) / z;
			if (dis1 > 0 && dis1 < x && dis2 > 0 && dis2 < z) {
				if (min_d != -1) {
					if (d < min_d) {
						min_d = d;
						N = unit_vector(-c_y);
					}
				}
				else {
					min_d = d;
					N = unit_vector(-c_y);
				}
			}
		}
	}

	oc = c_point + c_y - r.orgin();
	denom = dot(r.direction(), c_y);
	if (denom != 0) {
		float d = dot(oc, c_y) / denom;
		if (d > 0) {
			vec3 point = r.orgin() + d * r.direction();
			vec3 vex = point - c_point;
			float dis1 = dot(vex, c_x) / x;
			float dis2 = dot(vex, c_z) / z;
			if (dis1 > 0 && dis1 < x && dis2 > 0 && dis2 < z) {
				if (min_d != -1) {
					if (d < min_d) {
						min_d = d;
						N = unit_vector(c_y);
					}
				}
				else {
					min_d = d;
					N = unit_vector(c_y);
				}
			}
		}
	}
	oc = c_point - r.orgin();
	denom = dot(r.direction(), c_x);
	if (denom != 0) {
		float d = dot(oc, c_x) / denom;
		if (d > 0) {
			vec3 point = r.orgin() + d * r.direction();
			vec3 vex = point - c_point;
			float dis1 = dot(vex, c_y) / y;
			float dis2 = dot(vex, c_z) / z;
			if (dis1 > 0 && dis1 < y && dis2 > 0 && dis2 < z) {
				if (min_d != -1) {
					if (d < min_d) {
						min_d = d;
						N = unit_vector(-c_x);
					}
				}
				else {
					min_d = d;
					N = unit_vector(-c_x);
				}
			}
		}
	}
	oc = c_point + c_x - r.orgin();
	denom = dot(r.direction(), c_x);
	if (denom != 0) {
		float d = dot(oc, c_x) / denom;
		if (d > 0) {
			vec3 point = r.orgin() + d * r.direction();
			vec3 vex = point - c_point;
			float dis1 = dot(vex, c_y) / y;
			float dis2 = dot(vex, c_z) / z;
			if (dis1 > 0 && dis1 < y && dis2 > 0 && dis2 < z) {
				if (min_d != -1) {
					if (d < min_d) {
						min_d = d;
						N = unit_vector(c_x);
					}
				}
				else {
					min_d = d;
					N = unit_vector(c_x);
				}
			}
		}
	}
	if (min_d == -1) {
		return false;
	}
	else {
		rec.t = min_d;
		rec.p = r.point_at_parameter(rec.t);
		rec.normal = N;
		return true;

	}


}

#endif // !CUBE_H
#pragma once
