def irradiance (pmap,origin,direction,triangle_list,depth,RI) 
# pmap is a object of cPhotonMap; origin, direction, triangle_list are vectors 
	id=0 #id of intersected object 

	if (!intersect(origin, direction, triangle_list, t, id)) 
		return Vec() # if miss, return black
	
	const Triangle &obj = triangle_list[id]; # the hit object
 
	origin_new = origin + direction * t # t distance to intersection
	n = obj.norm
	nl = n.dot(direction) < 0 ? n : n * -1
	f = obj.surfaceColor


	if (++depth > 5) 
		return Vec() #R.R.
	
	if (obj.token == 'D')      
		Vec col
	    ###########################direct visualization of the photon map
		float color[3]
		pos[3] = {origin_new.x,origin_new.y,origin_new.z}
		normal[3] = {n.x, n.y, n.z}
		pmap.irradiance_estimate(color,pos,normal,0.1,100)
		col = Vec(color[0],color[1],color[2])
		return col
	
	elif (obj.token == 'S') ########################Ideal SPECULAR reflection 
		return f.mult(irradiance(pmap, origin_new, direction - n * 2 * n.dot(direction), triangle_list, depth, RI))
	elif (obj.token == 'T')
		reflRay(direction - n * 2 * n.dot(direction))# Ideal dielectric REFRACTION 
		bool into = n.dot(nl) > 0 # Ray from outside going in? 
		nc = 1
		nt = RI 
		nnt=into ? nc/nt : nt/nc, ddn=direction.dot(nl)

		if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) # Total internal reflection
			return f.mult(irradiance(pmap, origin_new,reflRay,triangle_list,depth, RI))

		tdir = (direction * nnt - n * ((into ? 1 : -1)*(ddn * nnt + sqrt(cos2t)))).norm();
		a = nt - nc
		b = nt + nc 
		R0 = a * a / (b * b) 
		c = 1 - (into ? -ddn : tdir.dot(n)) 
		Re = R0 + (1 - R0) * c * c * c * c * c
		Tr = 1 - Re 
		return obj.emissionColor + f.mult(irradiance(pmap, origin_new,reflRay,triangle_list,depth, RI)*Re+irradiance(pmap, origin_new,tdir,triangle_list,depth, RI)*Tr); 
	
	return Vec();
	 




def raytrace(pmap,pmap_caustic,origin,direction,triangle_list,light_list,depth,RI)
      
	id=0 # id of intersected object 
	
	if (!intersect(origin, direction, triangle_list, t, id)): # t is distance to intersection
		return Vec() # if miss, return black
	
	const Triangle &obj = triangle_list[id]        # the hit object
 
	origin_new = origin + direction * t
	n = obj.norm
	nl = n.dot(direction) < 0 ? n : n * -1 #法向量的正负
	f = obj.surfaceColor
   
	if ( ++depth > 20): 
		return obj.emissionColor # R.R. 
 
	if (obj.token == 'L'):
		return obj.emissionColor
	elif (obj.token == 'D' || obj.token == 'G'): # Ideal DIFFUSE reflection 
		col(0,0,0)


#****************   Direct visualization of Caustic ******************

	float color[3]
	pos[3] = {origin_new.x,origin_new.y,origin_new.z}
	normal[3] = {n.x, n.y, n.z}
	pmap_caustic.irradiance_estimate(color,pos,normal,0.1,100)
	col = Vec(color[0],color[1],color[2])

	
 #***************   Global illumination   ****************************
		nsamps = 200
		for i in(0 : nsamps-1):
	
			r1 = 2*PI*myrand()
			r2 = myrand()
			r2s = sqrt(r2)
     		w = nl
     		u = ((fabs(w.x) > .1 ? Vec(0,1,0) : Vec(1,0,0)) % w).norm()
     		v = w % u
     		d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1-r2)).norm()	
			col = col + irradiance(pmap,origin_new,d,triangle_list, 0, RI) * (1.0/nsamps)
		

#****************   Direct illumination  *****************************/
		for i in (1 : light_list.size()):
			factor = 1.0 / light_list[i].n_x / light_list[i].n_y
			for i in (i : j<light_list[i].n_x):
				for i in (1 : light_list[i].n_y):
					l_pos = light_list[i].pos - light_list[i].x_vec*0.5 + light_list[i].x_vec * (1./light_list[i].n_x*j) 
											  - light_list[i].y_vec*0.5 + light_list[i].y_vec * (1./light_list[i].n_y*k);
					d = (l_pos - origin_new)
					t_light = normalize(d)
					d = d.norm()
					id = 0
					if (!intersect(origin_new, d, triangle_list, t, id) || triangle_list[id].token=='L' || t>t_light):
						col = col + f.mult(light_list[i].color) * (d.dot(obj.norm)) * factor	
	
		
		if (obj.token == 'G'):
			col = col + f.mult(raytrace(pmap, pmap_caustic, origin_new,direction-n*2*n.dot(direction),triangle_list,light_list,depth,RI)); 
#**********************************************************************/		
		

		return col; 

	 elif (obj.token == 'S')  # Ideal SPECULAR reflection 
		return f.mult(raytrace(pmap, pmap_caustic, origin_new,direction - n * 2 * n.dot(direction),triangle_list,light_list,depth,RI)); 

   
	reflRay(direction - n * 2 * n.dot(direction)) # Ideal dielectric REFRACTION 
	into = n.dot(nl) > 0 # Ray from outside going in? 
	nc = 1
	nt = RI
	nnt = into ? nc/nt : nt/nc
	ddn = direction.dot(nl) 
	if ((cos2t=1 - nnt * nnt * (1 - ddn * ddn)) < 0) : # Total internal reflection
		return obj.emissionColor + f.mult(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI)); 
	

	tdir = (direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
	a = nt - nc
	b = nt + nc
	R0 = a * a / (b * b)
	c = 1 - ( into ? -ddn : tdir.dot(n)) 
	Re = R0 + (1 - R0) * c * c * c * c * c
	Tr = 1 - Re 
	return obj.emissionColor + f.mult(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI) * Re
                                     + raytrace(pmap, pmap_caustic, origin_new,tdir,triangle_list,light_list,depth,RI) * Tr)
