#ifndef RAY_H  
#define RAY_H  
#include "vec3.h"  

class ray{
public:
	ray() {}
	ray(const vec3& a, const vec3& b) { A = a; B = b; }
	vec3 orgin() const { return A; }
	vec3 direction() const { return B; }
	vec3 point_at_parameter(float t) const { return A + t*B; }
	//��֪tʱ�����Ի�ù����ϸõ�����꣨������  

private:
	vec3 A;
	vec3 B;

	/*
	virtual ~ray();

	protected:

	*/
};

#endif // RAY_H  #pragma once
