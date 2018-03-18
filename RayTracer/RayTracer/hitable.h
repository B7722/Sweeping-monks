#ifndef HITABLE_H
#define HITABLE_H
#include "ray.h"
struct hit_record {
	float t;
	vec3 p;
	vec3 normal;
};/*�ýṹ���¼��ײ�㡱������Ϣ����������ľ���t��ײ�����������p��ײ����ķ�����normal��*/ 

class hitable {
public:/*hitable������ʾ�ܹ�������ײ�ϵ��κ����塣���磬����*/ 
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
	/*hit()�ڴ˱�����Ϊ�麯������hitableΪ�����ࡣ������������б���ʵ�����麯��*/
};
#endif // HITABLE_H
