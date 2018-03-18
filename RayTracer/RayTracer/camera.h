#ifndef CAMERA_H
#define CMAERA_H
#include "ray.h"

class camera {
public:
	camera() {
		lower_left_corner=vec3(-2.0, -1.0, -1.0);//屏幕的左下角
		horizontal = vec3(4.0, 0.0, 0.0);
		vertical = vec3(0.0, 2.0, 0.0);
		origin = vec3(0.0, 0.0, 0.0);
	}
	ray get_ray(float u,float v){return ray(origin, lower_left_corner + u*horizontal + v*vertical);
	/*产生光线，光线方向为光源到屏幕上的点的方向（因为此例的光源为（0，0，0），所以方向直接为屏幕点的坐标）*/
	}

	vec3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 origin;
};
#endif // !CAMERA_H
#pragma once
