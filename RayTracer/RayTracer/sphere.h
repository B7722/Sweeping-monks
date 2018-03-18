#ifndef SPHERE_H
#define SPHERE_H

#include "hitable.h"

class sphere : public hitable {
	/*: public hitable��ʾsphere�̳�hitable������sphereΪhitable������*/
public:
	sphere() {}
	sphere(vec3 cen, float r) : center(cen), radius(r) {}
	/*�˴�Ϊʹ�ó�ʼ���б�Ĺ��캯������ʼ����Ա����*/
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	/*����ʵ�ָ�����麯����*/
	vec3 center;
	float radius;
};
bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	vec3 oc = r.orgin() - center;
	//�õ�����A-C=oc
	float a = dot(r.direction(), r.direction());
	//B��Ϊ���ߵķ���������a = (B��B)=dot(B,B)
	float b = 2.0 * dot(oc, r.direction());
	// b=2*(B��(A-C))=2*dot(B,(A-C))
	float c = dot(oc, oc) - radius*radius;
	//c= (A-C)��(A-C) - R2=dot((A-C),(A-C)) - R2
	float discriminant = b*b - 4 * a*c;
	//�б�ʽ>0����ʵ��������ײ����
	//�б�ʽ=0����1��ʵ��������ײ���򣬵���û�����䣻
	//�б�ʽ<0��û��ʵ��������û��ײ����
	if (discriminant > 0) {
		float temp = (-b - sqrt(discriminant)) / (2.0*a);//С��
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
		temp = (-b + sqrt(discriminant)) / (2.0*a);//���
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
	}
	/*�ж�С���ʹ���Ƿ��ڷ�Χ�ڡ������ж�С���Ƿ��ڷ�Χ�ڡ��ǣ����������Ϣ��Ȼ��ֱ�ӷ��أ����жϴ���Ƿ��ڷ�Χ�ڡ�Ҳ��������ѡС����С�������ٿ��Ǵ�� */
	return false;
}
#endif // !SPHERE_H
#pragma once
