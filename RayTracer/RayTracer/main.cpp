#include <iostream>  
#include <fstream>  
#include <limits>  
#include "sphere.h"  
#include "cube.h"
#include "hitablelist.h"  
#include "float.h" 
#include "camera.h"          

using namespace std;

vec3 color(const ray& r, hitable *world) {

	hit_record rec;
	if (world->hit(r, 0.0, (numeric_limits<float>::max)(), rec)) {
		return 0.5*vec3(rec.normal.x() + 1, rec.normal.y() + 1, rec.normal.z() + 1);
		/*有撞点：即为球体，将球体颜色设置为球在该点的单位法向量的色彩表映射值。另外numeric_limits<float>::max)()表示最大浮点数*/
	}
	else {
		vec3 unit_direction = unit_vector(r.direction());
		float t = 0.5*(unit_direction.y() + 1.0);
		return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue  
																	 // 无撞点：将颜色设置为背景色  
	}
}


int main() {
	int nx = 200;//图片宽
	int ny = 100;//图片长
	int ns = 100;//采样次数

	ofstream outfile("SeveralSpheres.txt", ios_base::out);
	outfile << "P3\n" << nx << " " << ny << "\n255\n";

	std::cout << "P3\n" << nx << " " << ny << "\n255\n";

	hitable *list[3];
	list[0] = new cube(vec3(0, -0.5, -2), vec3(1,0,0),vec3(0,1,0),vec3(0,0,1));
	list[1] = new sphere(vec3(0, -100.5, -1), 100);
	list[2] = new sphere(vec3(0, 0, -1), 0.5);
	hitable *world = new hitable_list(list, 3);
	/*将所有能够被撞击的物体信息保存在列表中*/

	camera cam;

	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {
			vec3 col(0,0,0);
			/*设置光线颜色*/
			for (int s = 0; s < ns; s++) {/*每个像素点区域采样ns次，此处ns=100*/

				float random = rand() % (100) / (float)(100);
				/*generate a random in range[0,1]。每个像素点的区域是以像素中心点为中心向外距离为1的范围。中心点位置+random相当于在这个像素点的区域内采样 */

				float u = float(i + random) / float(nx);
				float v = float(j + random) / float(ny);

				ray r = cam.get_ray(u, v);/*获得这个像素点区域随机采样点的颜色值。我们已经将颜色获得函数封装成一个叫做“camera”的类*/

				col += color(r, world);
				/*将这个像素点区域的所有ns个随机采样点的颜色值累加*/

			}

			col /= float(ns);
			/*将这个像素点区域的所有ns个随机采样点的颜色累加值除以ns获得其平均值，作为这个像素点区域最终的像素值。*/
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);

			outfile << ir << " " << ig << " " << ib << "\n";
			std::cout << ir << " " << ig << " " << ib << "\n";
		}
	}

}