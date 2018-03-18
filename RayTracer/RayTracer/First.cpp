#include<iostream>
#include"ray.h"
#include <fstream>
using namespace std;
bool hit_sphere(const vec3& center, float radius, const ray& r)

{

	vec3 oc = r.orgin() - center;
	//得到向量A-C=oc
	float a = dot(r.direction(), r.direction());
	//B即为光线的方向向量。a = (B·B)=dot(B,B)，
	float b = 2.0 * dot(oc, r.direction());
	// b=2*(B·(A-C))=2*dot(B,(A-C))，
	float c = dot(oc, oc) - radius*radius;
	//c= (A-C)·(A-C) - R2=dot((A-C),(A-C)) - R2
	float discriminant = b*b - 4 * (a*c);
	return (discriminant > 0);
	//判别式>0，有实根，光线撞上球，返回1；
	//判别式=0，有1个实根，光线撞上球，但是没有折射，返回0；
	//判别式<0，没有实根，光线没有撞上球，返回0；
}
vec3 color(const ray& r)
{

	if (hit_sphere(vec3(0, 0, -1), 0.5, r)) {
		return vec3(1, 0, 0);
	}
	//此处判断这条光线是否撞上球。
	//如果撞上，则执行return vec3(1, 0, 0)，将这个像素点设置为红色。
	//如果没有撞上，则不执行returnvec3(1, 0, 0)，程序向下运行，将这个像素点设置为背景色。
	vec3 unit_direction = unit_vector(r.direction());//单位方向（对方向向量进行标准化。） 
	float t = 0.5*(unit_direction.y() + 1.0);//标准化之后的y值在[-1,1]中y+1在[0,2]中0.5*(y+1)在[0,1]中 
	return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
	//t=0时，color=vec3(1,1,1)，乘以255后对应的rgb为（255，255，255）  
	//t=1时，color=vec3(0.5,0.7,1)，乘以255后对应的rgb为（127.5，178.5，255）  
	//如上两个颜色分别对应着“白色”和“浅蓝色”。  
	/*画面颜色=(1-t)*“白色”+ t*“浅蓝色”，即画面颜色为“白色”和“浅蓝色”（沿着y方向）的线性插值的结果。如果要换成x或者z方向，将上面的.y()改成.x()或者.z()即可。
	若要换其他颜色，设置对应的rgb值即可。rgb值和颜色的对应可参考word中“字体颜色设置”*/
}

int main() {
	int nx = 200;
	int ny = 100;

	ofstream outfile("Sphere_2D.txt", ios_base::out);
	outfile << "P3\n" << nx << " " << ny << "\n255\n";
	cout << "p3\n" << nx << " " << ny << "\n255\n";


	vec3 lower_left_corner(-2.0, -1.0, -1.0);//图像左下角坐标
	vec3 horizontal(4.0, 0.0, 0.0);//横向宽度
	vec3 vertical(0.0, 2.0, 0.0);//纵向高度
	vec3 origin(0.0, 0.0, 0.0);//原点位置


	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {
			float u = float(i) / float(nx);
			float v = float(j) / float(ny);
			ray r(origin, lower_left_corner + u*horizontal + v*vertical);//由画面中每个像素点在画面中的相对位置每个像素点对应的光线的方向向量从而确定画面的范围/大小。（完成第一步）
			vec3 col = color(r);//根据光线对每一个像素点上色。（完成第二步）
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);

			outfile << ir << " " << ig << " " << ib << "\n";
			std::cout << ir << " " << ig << " " << ib << "\n";
		}

	}
}

