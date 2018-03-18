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
		/*��ײ�㣺��Ϊ���壬��������ɫ����Ϊ���ڸõ�ĵ�λ��������ɫ�ʱ�ӳ��ֵ������numeric_limits<float>::max)()��ʾ��󸡵���*/
	}
	else {
		vec3 unit_direction = unit_vector(r.direction());
		float t = 0.5*(unit_direction.y() + 1.0);
		return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue  
																	 // ��ײ�㣺����ɫ����Ϊ����ɫ  
	}
}


int main() {
	int nx = 200;//ͼƬ��
	int ny = 100;//ͼƬ��
	int ns = 100;//��������

	ofstream outfile("SeveralSpheres.txt", ios_base::out);
	outfile << "P3\n" << nx << " " << ny << "\n255\n";

	std::cout << "P3\n" << nx << " " << ny << "\n255\n";

	hitable *list[3];
	list[0] = new cube(vec3(0, -0.5, -2), vec3(1,0,0),vec3(0,1,0),vec3(0,0,1));
	list[1] = new sphere(vec3(0, -100.5, -1), 100);
	list[2] = new sphere(vec3(0, 0, -1), 0.5);
	hitable *world = new hitable_list(list, 3);
	/*�������ܹ���ײ����������Ϣ�������б���*/

	camera cam;

	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {
			vec3 col(0,0,0);
			/*���ù�����ɫ*/
			for (int s = 0; s < ns; s++) {/*ÿ�����ص��������ns�Σ��˴�ns=100*/

				float random = rand() % (100) / (float)(100);
				/*generate a random in range[0,1]��ÿ�����ص�����������������ĵ�Ϊ�����������Ϊ1�ķ�Χ�����ĵ�λ��+random�൱����������ص�������ڲ��� */

				float u = float(i + random) / float(nx);
				float v = float(j + random) / float(ny);

				ray r = cam.get_ray(u, v);/*���������ص�����������������ɫֵ�������Ѿ�����ɫ��ú�����װ��һ��������camera������*/

				col += color(r, world);
				/*��������ص����������ns��������������ɫֵ�ۼ�*/

			}

			col /= float(ns);
			/*��������ص����������ns��������������ɫ�ۼ�ֵ����ns�����ƽ��ֵ����Ϊ������ص��������յ�����ֵ��*/
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);

			outfile << ir << " " << ig << " " << ib << "\n";
			std::cout << ir << " " << ig << " " << ib << "\n";
		}
	}

}