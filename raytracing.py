import numpy as np
import matplotlib.pyplot as plt

w = 400
h = 300

fov = 40.0
aspectratio = w/h
angle = tan(0.5*fov*PI/180.0)

AAFilter[AAKERNEL_SIZE][3] = {
        -0.52, 0.38, 0.128,
        0.41, 0.56, 0.119,
        0.27, 0.08, 0.294,
        -0.17, -0.29, 0.249,
        0.58, -0.55, 0.104,
        -0.31, -0.71, 0.106
    }

angle_x = 0
Rx[3][3] = {
        {1,0,0},
        {0,cos(angle_x*PI/180),-sin(angle_x*PI/180)},
        {0,sin(angle_x*PI/180),cos(angle_x*PI/180)}
    }
'''
    
    vector<Triangle> surfaces = parse(argv[1])
    vector<Light> lights

    Light light1, light2
    
    light1.pos = Vec(1.6,2.749,10.75)
    light1.color = Vec(0.5,0.5,0.5)
    light1.x_vec = Vec(1.2,0,0)
    light1.y_vec = Vec(0,0,1.2)
    light1.n_x = 9
    light1.n_y = 9
    
    lights.push_back(light1) 

    float total_area = 1.44
    float light_area[]={1.44}

    
    visualize_light(lights,surfaces)

    vector<Sphere> bound
    
    Sphere diamond(Vec(0,-0.75,10.75), 2.1)
    bound.push_back(diamond)
'''

    n_photons = 1000000
    n_photons_caustic = 1000000
    pmap_r(n_photons) = cPhotonMap()
    pmap_g(n_photons) = cPhotonMap()
    pmap_b(n_photons) = cPhotonMap()
    pmap_caustic_r(n_photons_caustic) = cPhotonMap()
    pmap_caustic_g(n_photons_caustic) = cPhotonMap() 
    pmap_caustic_b(n_photons_caustic) = cPhotonMap()
        
    
    for i in (0,lights.size()):
    
        while (pmap_r.stored_photons < (float)pmap_r.max_photons*light_area[i]/total_area) {
            offset_x = np.random.RandomState(0)
            offset_y = np.random.RandomState(0)
            Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
                                -lights[i].y_vec*0.5+lights[i].y_vec*offset_y
            Vec nl(0,-1,0)
            r1=2*PI*np.random.RandomState(0)
            r2=np.random.RandomState(0)
            r2s=sqrt(r2) 
            Vec w=nl, u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm(), 
            v=np.mod(w,u) 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            Vec color(15,15,15)
            photon_tracing(pmap_r,pos,d,surfaces,color,0,2.40)
        }

        print"finish generating red photon map"

        while (pmap_caustic_r.stored_photons<()pmap_caustic_r.max_photons*light_area[i]/total_area) {
            offset_x = np.random.RandomState(0)
            offset_y = np.random.RandomState(0)
            Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
                                -lights[i].y_vec*0.5+lights[i].y_vec*offset_y
            Vec nl(0,-1,0)
            r1=2*PI*np.random.RandomState(0), r2=np.random.RandomState(0), r2s=sqrt(r2) 
            Vec w=nl, u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm(), v=np.mod(w,u) 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            Vec color(1,1,1)
            photon_tracing_caustic(pmap_caustic_r,pos,d,surfaces,color,0,false, 2.40, bound)
        }

        print"finish generating caustic red photon map"

        while (pmap_g.stored_photons<(float)pmap_g.max_photons*light_area[i]/total_area) {
             offset_x = np.random.RandomState(0)
             offset_y = np.random.RandomState(0)
            Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
                                -lights[i].y_vec*0.5+lights[i].y_vec*offset_y
            Vec nl(0,-1,0)
             r1=2*PI*np.random.RandomState(0), r2=np.random.RandomState(0), r2s=sqrt(r2) 
            Vec w=nl, u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm(), v=np.mod(w,u) 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            Vec color(15,15,15)
            photon_tracing(pmap_g,pos,d,surfaces,color,0,2.43)
        }

        print"finish generating green photon map"

        while (pmap_caustic_g.stored_photons<(float)pmap_caustic_g.max_photons*light_area[i]/total_area) {
            offset_x = np.random.RandomState(0)
            offset_y = np.random.RandomState(0)
            Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
                                    -lights[i].y_vec*0.5+lights[i].y_vec*offset_y
            Vec nl(0,-1,0)
            r1=2*PI*np.random.RandomState(0), r2=np.random.RandomState(0), r2s=sqrt(r2) 
            w=nl
            u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm()
            v=np.mod(w,u) 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            Vec color(1,1,1)
            photon_tracing_caustic(pmap_caustic_g,pos,d,surfaces,color,0,false, 2.43, bound)
        }


        print"finish generating causitc green photon map"

        while (pmap_b.stored_photons<(float)pmap_b.max_photons*light_area[i]/total_area) {
             offset_x = np.random.RandomState(0)
             offset_y = np.random.RandomState(0)
            Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
                                    -lights[i].y_vec*0.5+lights[i].y_vec*offset_y
            Vec nl(0,-1,0)
            r1=2*PI*np.random.RandomState(0)
            r2=np.random.RandomState(0)
            r2s=sqrt(r2) 
            w=nl
            u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm() 
            v=np.mod(w,u) 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            Vec color(15,15,15)
            photon_tracing(pmap_b,pos,d,surfaces,color,0,2.46)
        }

        print"finish generating blue photon map"
        while (pmap_caustic_b.stored_photons<(float)pmap_caustic_b.max_photons*light_area[i]/total_area):
            offset_x = np.random.RandomState(0)
            offset_y = np.random.RandomState(0)
            Vec pos = lights[i].pos-lights[i].x_vec*0.5+lights[i].x_vec*offset_x 
                                    -lights[i].y_vec*0.5+lights[i].y_vec*offset_y
            Vec nl(0,-1,0)
            r1=2*PI*np.random.RandomState(0), r2=np.random.RandomState(0), r2s=sqrt(r2) 
            Vec w=nl, u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm(), v=np.mod(w,u) 
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            Vec color(1,1,1)
            photon_tracing_caustic(pmap_caustic_b,pos,d,surfaces,color,0,false, 2.46, bound)
        

        print"finish generating caustic blue photon map"

        
    
    pmap_r.balance()
    pmap_r.scale_photon_power(1.0/n_photons)   
    pmap_caustic_r.balance()
    pmap_caustic_r.scale_photon_power(1.0/n_photons_caustic)
    pmap_g.balance()
    pmap_g.scale_photon_power(1.0/n_photons)   
    pmap_caustic_g.balance()
    pmap_caustic_g.scale_photon_power(1.0/n_photons_caustic)
    pmap_b.balance()
    pmap_b.scale_photon_power(1.0/n_photons)       
    pmap_caustic_b.balance()
    pmap_caustic_b.scale_photon_power(1.0/n_photons_caustic)   

    visualize_light(lights,surfaces)
    
    

    Vec *c=new Vec[w*h]
    Vec output_color

def normalize(x):
    x /= np.linalg.norm(x)
    return x

def intersect_plane(O, D, P, N):
    # Return the distance from O to the intersection of the ray (O, D) with the 
    # plane (P, N), or +inf if there is no intersection.
    # O and P are 3D points, D and N (normal) are normalized vectors.
    denom = np.dot(D, N)
    if np.abs(denom) < 1e-6:
        return np.inf
    d = np.dot(P - O, N) / denom
    if d < 0:
        return np.inf
    return d

def intersect_sphere(O, D, S, R):
    # Return the distance from O to the intersection of the ray (O, D) with the 
    # sphere (S, R), or +inf if there is no intersection.
    # O and S are 3D points, D (direction) is a normalized vector, R is a scalar.
    a = np.dot(D, D)
    OS = O - S
    b = 2 * np.dot(D, OS)
    c = np.dot(OS, OS) - R * R
    disc = b * b - 4 * a * c
    if disc > 0:
        distSqrt = np.sqrt(disc)
        q = (-b - distSqrt) / 2.0 if b < 0 else (-b + distSqrt) / 2.0
        t0 = q / a
        t1 = c / q
        t0, t1 = min(t0, t1), max(t0, t1)
        if t1 >= 0:
            return t1 if t0 < 0 else t0
    return np.inf


def get_normal(obj, M):
    # Find normal.
    if obj['type'] == 'sphere':
        N = normalize(M - obj['position'])
    elif obj['type'] == 'plane':
        N = obj['normal']
    return N
    
def get_color(obj, M):
    color = obj['color']
    if not hasattr(color, '__len__'):
        color = color(M)
    return color

def judge_in_triangle(P0,P1,P2,P3):
    E12 = P2 - P1
    E13 = P3 - P1
    E10 = P0 - P1
    v11 = np.dot(E13,E13)
    v00 = np.dot(E12,E12)
    v20 = np.dot(E10,E12)
    v10 = np.dot(E13,E12)
    v21 = np.dot(E10,E13)
    u = (v11*v20-v10*v21)/(v00*v11-v10*v10)
    v = (v00*v21-v10*v20)/(v00*v11-v10*v10)
    if(u>=0 and v>=0 and u+v<=1):
        return True
    else:
        return False

def intersect_tetrahedron(O, D, P1, P2, P3, P4):
    E12 = P2 - P1
    E23 = P3 - P2
    E34 = P4 - P3
    E14 = P4 - P1
    N123 = normalize(np.cross(E12, E23))
    N234 = normalize(np.cross(E34, E23))
    N124 = normalize(np.cross(E14, E12))
    N134 = normalize(np.cross(E14, E34))
    min_d = np.inf
    vector = np.array([0, 0, 0])
    e = 0

    denom = np.dot(D, N123)
    if (denom != 0):
        d = np.dot(P1 - O, N123) / denom
        if d > e:
            point = O + d * D
            if judge_in_triangle(point, P1, P2, P3):
                if (d < min_d):
                    min_d = d
                    vector = N123

    denom = np.dot(D, N234)
    if (denom != 0):
        d = np.dot(P2 - O, N234) / denom
        if d > e:
            point = O + d * D
            if judge_in_triangle(point, P2, P3, P4):
                if (d < min_d):
                    min_d = d
                    vector = N234

    denom = np.dot(D, N124)
    if (denom != 0):
        d = np.dot(P1 - O, N124) / denom
        if d > e:
            point = O + d * D
            if judge_in_triangle(point, P1, P2, P4):
                if (d < min_d):
                    min_d = d
                    vector = N124

    denom = np.dot(D, N134)
    if (denom != 0):
        d = np.dot(P1 - O, N134) / denom
        if d > e:
            point = O + d * D
            if judge_in_triangle(point, P1, P3, P4):
                if (d < min_d):
                    min_d = d
                    vector = N134

    return dict(d=min_d, v=vector)

def intersect_cube(O, D, C, RX, RY, RZ):
    # Return the distance from O to the intersection of the ray (O, D) with the
    # cube (C, RX, RY, RZ), or +inf if there is no intersection.
    # O and C are 3D points, D is a normalized vector, RX, RY and RZ are ranged vectors.

    min_d = np.inf
    vector= np.array([0,0,0])
    X = np.linalg.norm(RX)
    Y = np.linalg.norm(RY)
    Z = np.linalg.norm(RZ)
    RX = normalize(RX)
    RY = normalize(RY)
    RZ = normalize(RZ)
    e = 0

    denom = np.dot(D,-RZ)
    if(denom != 0):
        d = np.dot(C - O, -RZ) / denom
        if d>e:
            point = O + d * D
            vex = point - C
            if 0 < np.dot(vex, RX)  < X and 0 < np.dot(vex, RY)  < Y:
                    if (d < min_d):
                        min_d = d
                        vector=-RZ


    denom = np.dot(D,RZ)
    if(denom != 0):
        d = np.dot(C + RZ - O, RZ) / denom
        if d>e:
            point = O + d * D
            vex = point - C-RZ
            if 0 < np.dot(vex, RX)< X and 0 < np.dot(vex, RY) < Y:
                if (d < min_d):
                    min_d = d
                    vector = RZ

    denom = np.dot(D, -RY)
    if(denom != 0):
        d = np.dot(C - O, -RY) / denom
        if d>e:
            point = O + d * D
            vex = point - C
            if 0 < np.dot(vex, RZ) < Z and 0 < np.dot(vex, RX) < X:
                if (d < min_d):
                    min_d = d
                    vector = -RY

    denom = np.dot(D, RY)
    if(denom != 0):
        d = np.dot(C + RY - O, RY) / denom
        if d>e:
            point = O + d * D
            vex = point - C-RY
            if 0 < np.dot(vex, RZ)  < Z and 0 < np.dot(vex, RX)  < X:
                if (d < min_d):
                    min_d = d
                    vector = RY

    denom = np.dot(D, -RX)
    if(denom != 0):
        d = np.dot(C - O, -RX) / denom
        if d>e:
            point = O + d * D
            vex = point - C
            if 0 < np.dot(vex, RZ) < Z and 0 < np.dot(vex, RY)  < Y:
                if (d < min_d):
                    min_d = d
                    vector = -RX

    denom = np.dot(D, RX)
    if(denom != 0):
        d = np.dot(C + RX - O, RX) / denom
        if d>e:
            point = O + d * D
            vex = point - C-RX
            if 0 < np.dot(vex, RZ)  < Z and 0 < np.dot(vex, RY) < Y:
                if (d < min_d):
                    min_d = d
                    vector = RX

    return dict(d = min_d,v = vector)

def intersect(O, D, obj):
    if obj['type'] == 'plane':
        return intersect_plane(O, D, obj['position'], obj['normal'])
    elif obj['type'] == 'sphere':
        return intersect_sphere(O, D, obj['position'], obj['radius'])
    elif obj['type'] == 'cube':
        return intersect_cube(O, D, obj['position'], obj['radius_X'], obj['radius_Y'], obj['radius_Z'])
    elif obj['type']=='tetrahedron':
        return intersect_tetrahedron(O,D,obj['p1'],obj['p2'],obj['p3'],obj['p4'])

def trace_ray(rayO, rayD):
        # Find first point of intersection with the scene.
        t = np.inf
        N = np.array([0, 0, 0])
        for i, obj in enumerate(scene):
            t_obj = intersect(rayO, rayD, obj)
            d = t_obj['d']
            n = t_obj['v']
            if d < t:
                t, obj_idx = d, i
                N = n
        # Return None if the ray does not intersect any object.
        if t == np.inf:
            return
        if np.array_equal(N, np.array([0, 0, 0])):
            return
        # Find the object.
        obj = scene[obj_idx]
        # Find the point of intersection on the object.
        M = rayO + rayD * t
        # Find properties of the object.
        color = get_color(obj, M)
        toL = normalize(L - M)
        toO = normalize(O - M)
        # Shadow: find if the point is shadowed or not.
        l_obj = [intersect(M + N * .0001, toL, obj_sh)
                 for k, obj_sh in enumerate(scene) if k != obj_idx]
        if l_obj:
            l = []
            for a, b in enumerate(l_obj):
                l.append(b['d'])
            if l and min(l) < np.inf:
                return
        # Start computing the color.
        col_ray = ambient
        # Lambert shading (diffuse).
        col_ray += obj.get('diffuse_c', diffuse_c) * max(np.dot(N, toL), 0) * color
        # Blinn-Phong shading (specular).
        col_ray += obj.get('specular_c', specular_c) * max(np.dot(N, normalize(toL + toO)),
                                                           0) ** specular_k * color_light
        return obj, M, N, col_ray

def add_cube(position, radius_X, radius_Y, radius_Z, color):
        return dict(type='cube', position=np.array(position), radius_X=np.array(radius_X), radius_Y=np.array(radius_Y),
                    radius_Z=np.array(radius_Z), color=np.array(color), reflection=.5)

def add_tetrahedron(p1, p2, p3, p4, color):
        return dict(type='tetrahedron', p1=np.array(p1), p2=np.array(p2), p3=np.array(p3), p4=np.array(p4),
                    color=np.array(color), reflection=.5)


def irradiance(pmap = cPhotonMap(),origin = [],direction = [],triangle_list = [],depth,RI) 
# pmap is a object of cPhotonMap origin, direction, triangle_list are vectors 
    id = 0 #id of intersected object 

    if (not intersect(origin, direction, triangle_list, t, id)) 
        return [0,0,0] # if miss, return black
    
    const Triangle &obj = triangle_list[id] ###################################### the hit object
 
    origin_new = origin + np.multiply(direction,t) # t distance to intersection
    n = obj.norm
    nl = n if n.np.dot(direction) < 0 else np.multiply(n,-1) # 内积小于零为正，否则为负
    sfco = obj.surfaceColor

    depth += 1
    if (depth > 5) 
        depth += 1
        return [0,0,0] # R.R.
    
    if (obj.token == 'D'):
        ########################### direct visualization of the photon map
        pos = [origin_new.x,origin_new.y,origin_new.z]
        normal = [n.x, n.y, n.z]
        pmap.irradiance_estimate(color,pos,normal,0.1,100)
        col = [color[0],color[1],color[2]]
        return col
    
    elif (obj.token == 'S') ######################## Ideal SPECULAR reflection 
        return sfco.np.multiply(irradiance(pmap, origin_new, direction - np.multiply(n,np.multiply(2,n.np.dot(direction))), triangle_list, depth, RI))

    elif (obj.token == 'T')
        reflRay(np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction)))))# Ideal dielectric REFRACTION 
        bool into = n.np.dot(nl) > 0 # Ray from outside going in? 
        nc = 1
        nt = RI 
        nnt = np.divide(nc,nt) if into else np.divide(nt,nc)
        ddn = direction.np.dot(nl)

        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) # Total internal reflection
            return f.np.multiply(irradiance(pmap, origin_new,reflRay,triangle_list,depth, RI))

        #tdir = (direction * nnt - n * ((1 if into else -1)*(np.multiply(ddn,nnt) + np.sqrt(cos2t)))).norm()
        #a = nt - nc
        #b = nt + nc 
        #R0 = a * a / (b * b) 
        #c = 1 - ( -ddn if into else tdir.np.dot(n)) 
        #Re = R0 + (1 - R0) * c * c * c * c * c
        #Tr = 1 - Re 

        #tdir = (np.multiply(direction,nnt) - n * ((1 if into else -1) * (np.multiply(ddn,nnt) + np.sqrt(cos2t)))).norm()
        tdir = np.substract(np.multiply(direction,nnt),np.multiply(n,np.multiply((1 if into else -1),np.add(np.multiply(ddn,nnt), np.sqrt(cos2t))))).norm()
        a = np.substract(nt,nc)
        b = np.add(nt,nc)
        R0 = np.divide(np.multiply(a,)a,np.multiply(b,b))
        c = np.substract(1, (~ddn if into else tdir.np.dot(n)) )
        Re = np.add(R0,np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(np.substract(1,R0), c),c) c)  c) c))
        Tr = np.substract(1,Re) 

        return obj.emissionColor + f.np.multiply(irradiance(pmap, origin_new,reflRay,triangle_list,depth, RI)*Re+irradiance(pmap, origin_new,tdir,triangle_list,depth, RI)*Tr)
    
    return [0,0,0]
     




def raytrace(pmap = cPhotonMap(),pmap_caustic = cPhotonMap(),origin = [],direction = [],triangle_list = [],light_list,depth,RI)
      
    id=0 # id of intersected object 
    
    if (~intersect(origin, direction, triangle_list, t, id)): # t is distance to intersection
        return [0,0,0] # if miss, return black
    
    obj = triangle_list[id]        # the hit object
 
    origin_new = origin + np.multiply(direction,t)
    n = obj.norm
    nl = n if n.np.dot(direction) < 0 else n * -1 #法向量的正负
    f = obj.surfaceColor
   
    depth += 1
    if (depth > 20): 
        depth += 1
        return obj.emissionColor # R.R. 
 
    if (obj.token == 'L'):
        return obj.emissionColor
    elif (obj.token == 'D' or obj.token == 'G'): # Ideal DIFFUSE reflection 
        [0,0,0]


#****************   Direct visualization of Caustic ******************

    
    pos = [origin_new.x,origin_new.y,origin_new.z]
    normal = [n.x, n.y, n.z]
    pmap_caustic.irradiance_estimate(color,pos,normal,0.1,100)
    col = [color[0],color[1],color[2]]

    
 #***************   Global illumination   ****************************
        nsamps = 200
        for i in(0 : nsamps):
    
            r1 = 2 * PI * np.random.RandomState(0)
            r2 = np.random.RandomState(0)
            r2s = np.sqrt(r2)
            w = nl
            u = np.mod(([0,1,0] if fabs(w.x) > .1 else [1,0,0]),w).norm()
            v = np.mod(w,u)
            d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * np.sqrt(1-r2)).norm()  
            col = col + irradiance(pmap,origin_new,d,triangle_list, 0, RI) * (1.0/nsamps)
        

#****************   Direct illumination  *****************************/
        for i in (1 : light_list.size()):
            factor = 1.0 / light_list[i].n_x / light_list[i].n_y
            for i in (i , j<light_list[i].n_x):
                for i in (1 , light_list[i].n_y):
                    l_pos = light_list[i].pos - light_list[i].x_vec*0.5 + light_list[i].x_vec * (1./light_list[i].n_x*j) 
                                              - light_list[i].y_vec*0.5 + light_list[i].y_vec * (1./light_list[i].n_y*k)
                    d = (l_pos - origin_new)
                    t_light = normalize(d)
                    d = d.norm()
                    id = 0
                    if (~intersect(origin_new, d, triangle_list, t, id) or triangle_list[id].token=='L' or t>t_light):
                        col = col + f.np.multiply(light_list[i].color) * (d.np.dot(obj.norm)) * factor  
    
        
        if (obj.token == 'G'):
            col = col + f.np.multiply(raytrace(pmap, pmap_caustic, origin_new,np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction)))),triangle_list,light_list,depth,RI))
#**********************************************************************/        
        

        return col 

     elif (obj.token == 'S')  # Ideal SPECULAR reflection 
        return f.np.multiply(raytrace(pmap, pmap_caustic, origin_new,direction - n * 2 * n.np.dot(direction),triangle_list,light_list,depth,RI))

   
    reflRay(np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction))))) # Ideal dielectric REFRACTION 
    into = n.np.dot(nl) > 0 # Ray from outside going in? 
    nc = 1
    nt = RI
    nnt = into ? nc/nt : nt/nc
    ddn = direction.np.dot(nl) 
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) : # Total internal reflection
        return obj.emissionColor + f.np.multiply(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI)) 
    

    tdir = np.substract(np.multiply(direction,nnt),np.multiply(n,np.multiply((1 if into else -1),np.add(np.multiply(ddn,nnt), np.sqrt(cos2t))))).norm()
    a = np.substract(nt,nc)
    b = np.add(nt,nc)
    R0 = np.divide(np.multiply(a,)a,np.multiply(b,b))
    c = np.substract(1, (~ddn if into else tdir.np.dot(n)) )
    Re = np.add(R0,np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(np.substract(1,R0), c),c) c)  c) c))
    Tr = np.substract(1,Re) 
    return obj.emissionColor + f.np.multiply(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI) * Re
                                     + raytrace(pmap, pmap_caustic, origin_new,tdir,triangle_list,light_list,depth,RI) * Tr)

#
# Photon mapping
#

import objlib

# define a photon object
class cPhoton:
    def __init__(self):
        self.pos = np.zeros(3)  # photon position
        self.power = np.zeros(3)   # photon power (uncompressed)
        self.plane  # splitting plane for kd-tree
        self.theta
        self.phi  # incoming direction

# define the nearest photons
class cNearestPhotons:
    def __int__(self):
        self.maxi
        self.found # how many photons have been found
        self.got_heap # to build a heap
        self.pos = np.zeros(3) # x,y,z
        self.dist2s = [] # x^2 + y^2 + z^2
        self.index = []
        self.photons = []

    def addPhoton(self, p):
        self.photons.append(p)


# define photon maps
class cPhotonMap(mxn_pho): # incoming parameter which means how many photons can be stored
    def __init__(self):
        #free(photons) ###############################

        self.stored_pho = 0 # how many photons are stored
        self.prev_scale = 1  
        self.mxn_photons = mxn_pho # user to set the maximum number of photon to be stored

        # 判断内存是否充足
        photons = cPhoton() # (cPhoton*) malloc(sys.getsizeof(cPhoton)*(mxn_Photons)+1)
        if(photons == 0)
            print "Out of memory initializing photon map"
            exit(-1)
        

        bbox_min[0] = bbox_min[1] = bbox_min[2] = 1e8
        bbox_max[0] = bbox_max[1] = bbox_max[2] = -1e8


        for i in (1,256)
            angle = i * (1.0 / 256.0) * PI # 将半球分成256份
            costheta[i] = cos(angle) 
            sintheta[i] = sin(angle)
            cosphi[i] = cos(2.0 * angle)
            sinphi[i] = sin(2.0 * angle)
    
    # define the direction of a photon
    def photon_dire(dire, pton):
        dire = []
        dire.append(0) = sintheta[pton.theta] * cosphi[pton.phi]
        dire.append(1) = sintheta[pton.theta] * sinphi[pton.phi]
        dire.append(2) = costheta[pton.theta]

    def irradiance_estimate(irrad = [],pos = [],normal = [],max_dist,num_ptons)
        irrad[0] = irrad[1] = irrad[2] = 0.0

        nrp = cNearestPhotons()
        nrp.dist2s = float(alloca(sys.getsizeof(float)*(num_ptons+1)))
        nrp.index = cPhoton(alloca(sys.getsizeof(cPhoton)*(num_ptons+1)))

        nrp.pos[0] = pos[0] #target
        nrp.pos[1] = pos[1]
        nrp.pos[2] = pos[2]
        nrp.maxi = num_ptons 
        nrp.found = 0
        nrp.got_heap = 0
        nrp.dist2s[0] = max_dist * max_dist

        locate_photons(nrp,1)
        
        if (nrp.found < 8):
            return
        pdir = []

        for i in (1, nrp.found+1):
            p = cPhoton()
            p = nrp.index[i]
            photon_dir(pdir,p)
            if ((pdir[0] * normal[0] + pdir[1] * normal[1] + pdir[2] * normal[2]) < 0.0):
                irrad[0] += p.power[0]
                irrad[1] += p.power[1]
                irrad[2] += p.power[2]
        
        tmp = (1.0/PI)/(nrp.dist2s[0])
        
        irrad[0] *= tmp
        irrad[1] *= tmp
        irrad[2] *= tmp

    
    # define the location of a photon
    def locate_photons(nrp = cNearestPhotons(),index): # index is a type of INT
        p = cPhoton.obj[index] # 需要一个cPhoton类型的
        #nrp = cNearestPhotons()

        # half_stored_photons has been define in cPhotonMap.balance
        # 从叶子节点开始构造树
        if (index < half_stored_photons): 
            dist1 = nrp.pos[p.plane] - p.pos[p.plane] # the distance between the target photon and the farest nearest photon
            
            if (dist1 > 0.0): # right-subtree
                locate_photons(nrp, 2 * index + 1) # 2i + 1
                if (dist1 * dist1 < nrp.dist2s[0]): #################################################################
                    locate_photons(nrp, 2 * index)
            else: # left-subtree (dist1 <= 0.0)
                locate_photons(nrp, 2 * index) # 2i
                if (dist1 * dist1 < nrp.dist2s[0]):
                    locate_photons(nrp, 2 * index + 1)
            

        dist1 = p.pos[0] - nrp.pos[0] # the distance between the target photon and the near photon||||x
        dist2s = dist1 * dist1 # to avoid the negative distance value |||| x^2
        dist1 = p.pos[1] - nrp.pos[1] # y
        dist2s += dist1 * dist1 # x^2 + y^2
        dist1 = p.pos[2] - nrp.pos[2] # z
        dist2s += dist1 * dist1 # x^2 + y^2 + z^2
        
        if (dist2s < nrp.dist2s[0]):
            # we found a photon, Insert it in the candidate list
            if (nrp.found < nrp.maxi):
                # heap is not full use array
                nrp.found += 1
                nrp.dist2s.append(dist2s) # = dist2s
                nrp.index.append(p)
            else:
                if (nrp.got_heap == 0): # need to build the heap
                    phot = cPhoton()
                    half_found = nrp.found >> 1 # Returns nrp.found with the bits shifted to the right by 1 place. 
                    for k in range(half_found,0,-1): # (int k=half_found k>=1k--)
                        parent = k 
                        phot = nrp.index[k] # set index-k photon as parent
                        dst2 = nrp.dist2s[k] # set the distance^2 of k photon as dst2
                        while (parent <= half_found): 
                            j = parent + parent 
                            if (j < nrp.found and nrp.dist2s[j] < nrp.dist2s[j+1]):
                                j += 1
                            if (dst2 >= nrp.dist2s[j]):
                                break
                            nrp.dist2s[parent] = nrp.dist2s[j]
                            nrp.index[parent] = nrp.index[j]
                            parent = j
                        
                        nrp.dist2s[parent] = dst2
                        nrp.index[parent] = phot
                    
                    nrp.got_heap = 1
                
                parent = 1 # for root
                j = 2 # start at index 2
                while (j <= nrp.found): 
                    if (j < nrp.found and nrp.dist2s[j] < nrp.dist2s[j+1]):
                        j += 1
                    if (dist2s > nrp.dist2s[j]):
                        break
                    nrp.dist2s[parent] = nrp.dist2s[j]
                    nrp.index[parent] = nrp.index[j]
                    parent = j
                    j += j
                
                nrp.index[parent] = p
                nrp.dist2s[parent] = dist2s
                
                nrp.dist2s[0] = nrp.dist2s[1]
    

    def store (power = [],pos = [], dire = []):
        if (stored_photons > mxn_photons):
            return
        stored_photons += 1
        node = cPhoton()

        for i in（0，3):
            node.pos[i] =  pos[i]
            if (node.pos[i] < bbox_min[i]):
                bbox_min[i] = node.pos[i]
            if (node.pos[i] > bbox_max[i]):
                bbox_max[i] = node.pos[i]
            node.power[i] = power[i]


        theta = int(acos(dire[2])*(256.0/PI))
        if (theta > 255):
            node.theta = 255
        else:
            node.theta = theta

        phi = int(atan2(dire[1],dire[0])*(256.0/(2.0*PI)))

        if (phi > 255):
            node.phi = 255
        elif: (phi < 0)
            node.phi = phi + 256 # (unsigned char)
        else:
            node.phi = phi # (unsigned char)


def scale_photon_power(scale):
        for (prev_scale,stored_photons + 1):
            photons[i].power[0] *= scale
            photons[i].power[1] *= scale
            photons[i].power[2] *= scale

        prev_scale = stored_photons

def balance():
        if (stored_photons > 1):
            pa1 = cPhoton() # Photon **pa1 = (Photon**) malloc(sizeof(Photon*)*(stored_photons+1))
            pa2 = cPhoton() # Photon **pa2 = (Photon**) malloc(sizeof(Photon*)*(stored_photons+1))

            for i in (0,stored_photons + 1):
                pa2[i] = photons[i]

            balance_segment(pa1, pa2, 1, 1, stored_photons)
            free(pa2)
            
            j = 1
            foo = 1
            foo_photon = photons[j]

            for i in (1, stored_photons + 1):
                d = pa1[j] - photons
                pa1[j] = 0
                if (d != foo):
                    photons[j] = photons[d]
                else:
                    photons[j] = foo_photon
                    if (i < stored_photons):
                        for foo in (,stored_photons + 1)
                            if (pa1[foo] != 0):
                                break
                        foo_photon = photons[foo]
                        j = foo
                    
                    continue
                
                j = d
            
            free(pa1)

        half_stored_photons = stored_photons / 2 - 1


    #define 
def swap(ph,a,b):
    ph2 = cPhoton()
    ph2 = ph[a]
    ph[a] = ph[b] 
    ph[b] = ph2


def median_split(p = cPhoton(),start, end, median, axis)
        left = start
        right = end

        while (right > left):
            v = p[right].pos[axis]
            i = left - 1
            j = right
            while true:
                i += 1
                while (p[i].pos[axis] < v):
                    pass
                j -= 1
                while (p[j].pos[axis] > v and j > left):
                    pass
                if (i >= j):
                    break
                swap(p,i,j)
            
            swap(p,i,right)
            if (i >= median):
                right = i - 1
            if (i <= median):
                left = i + 1
        
    

def balance_segment(pbal = cPhoton(), porg = cPhoton(), index, start, end):
        median = 1
        while ((4 * median) <= (end - start + 1)):
            median += median
        if ((3 * median) <= (end-start+1)):
            median += median
            median += start - 1
        else:
            median = end-median +1


        axis = 2
        if ((bbox_max[0]-bbox_min[0])>(bbox_max[1]-bbox_min[1]) and
            (bbox_max[0]-bbox_min[0])>(bbox_max[2]-bbox_min[2]))
            axis = 0
        elif ((bbox_max[1]-bbox_min[1])>(bbox_max[2]-bbox_min[2]))
            axis = 1

        median_split (porg, start, end, median, axis)
        pbal[index] = porg[median]
        pbal[index].plane = axis


        if (median > start):
            if (start<median-1):
                tmp = bbox_max[axis]
                bbox_max[axis] = pbal[index].pos[axis]
                balance_segment(pbal, porg, 2*index, start, median-1)
                bbox_max[axis] = tmp
            else:
                pbal[2*index] = porg[start]
            
        
        if (median<end):
            if (median+1<end):
                tmp = bbox_min[axis]
                bbox_min[axis] = pbal[index].pos[axis]
                balance_segment(pbal,porg, 2*index+1, median+1, end)
                bbox_min[axis] = tmp
            else:
                pbal[2*index+1] = porg[end]


def photon_tracing(pmap = cPhotonMap(),origin = [],direction = [], triangle_list = [],color=[],depth,RI)
    # t distance to intersection 
    id = 0 # id of intersected object 
    if (~intersect(origin, direction, triangle_list, t, id)):
        return  # if miss, return black
    
    obj = hit_object[id]
    #const Triangle &obj = triangle_list[id]      # the hit object
    depth += 1
    if (depth > 5):
        depth += 1
        return # R.R.
    origin_new = np.add(origin,np.multiply(direction, t))
    n = obj.norm
    nl = n if n.np.dot(direction) < 0 else np.multiply(n,-1)
    f = obj.surfaceColor
    direction.norm()
    
    
    if (f.x > f.y and f.x > f.z):
        p = f.x
    elif(f.y > f.z):
        p = f.y
    else:
        p = f.z

    
    if (obj.token == 'D' and obj.token == 'G'): #                 // Ideal DIFFUSE reflection 
        if (np.random.RandomState(0) > p):   # //absorb the photon
            color_tmp = [np.multiply(f.x,np.multiply(color.x,(1.0/(1-p)))),np.multiply(f.y,np.multiply(color.y,(1.0/(1-p)))),np.multiply(f.z,np.multiply(color.z,(1.0/(1-p))))]
            pos_tmp = [origin_new.x,origin_new.y,origin_new.z]
            dir_tmp = [direction.x,direction.y,direction.z]
            pmap.store(color_tmp,pos_tmp,dir_tmp)
            return
                        
        
        else:  # // trace another ray
            r1=2*PI*np.random.RandomState(0)
            r2=np.random.RandomState(0)
            r2s= np.sqrt(r2)

            w = nl 
            u = np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm()
            v = np.mod(w,u) 
            d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*np.sqrt(1-r2)).norm()
            photon_tracing(pmap, origin_new,d,triangle_list, f.np.multiply(color)*(1.0/p), depth, RI)

    elif (obj.token == 'S')    # // Ideal SPECULAR reflection 
        photon_tracing(pmap,origin_new,np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction)))),triangle_list,f.np.multiply(color)*(1.0/p),depth, RI) 
        
    elif (obj.token == 'T')   
       
        reflRay(np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction)))))#    // Ideal dielectric REFRACTION 
        bool into = n.np.dot(nl)>0              # Ray from outside going in? 
        nc=1
        nt=RI
        nnt=np.divide(nc,nt) if into else np.divide(nt,nc)
        ddn=direction.np.dot(nl)
        if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0):#// Total internal reflection

            photon_tracing(pmap,origin_new,reflRay,triangle_list, color, depth, RI) 
        
        tdir = np.substract(np.multiply(direction,nnt),np.multiply(n,np.multiply((1 if into else -1),np.add(np.multiply(ddn,nnt), np.sqrt(cos2t))))).norm()
        a = np.substract(nt,nc)
        b = np.add(nt,nc)
        R0 = np.divide(np.multiply(a,)a,np.multiply(b,b))
        c = np.substract(1, (~ddn if into else tdir.np.dot(n)) )
        Re = np.add(R0,np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(np.substract(1,R0), c),c) c)  c) c))
        Tr = np.substract(1,Re)
        P=np.add(.25,np.multiply(.5,Re))
        RP=np.divide(Re,P)
        TP=np.divide(Tr,np.substract(1,P))
        np.random.RandomState(0)<P #?   // Russian roulette 
        photon_tracing(pmap, origin_new,reflRay,triangle_list,np.multiply(color,RP),depth, RI):
        photon_tracing(pmap, origin_new,tdir,triangle_list,np.multiply(color,TP), depth, RI)



def photon_tracing_caustic(pmap = cPhotonMap(),origin=[],direction=[],triangle_list=[],color =[],depth, bool flag,RI,bound = []):
    #float t                               // distance to intersection 
    id=0 #                               // id of intersected object 
    
    if (~flag):
        for i in (0,i<bound.size()):
            if (~bound[i].intersect(origin, direction)):
                return


    
    if (~intersect(origin, direction, triangle_list, t, id)):
        return #  // if miss, return black

    depth += 1
    if (depth>5): 
        depth += 1
        return #//R.R.
    obj = triangle_list[id]#       // the hit object
    
    origin_new = origin + np.multiply(direction, t)
    n = obj.norm
    nl = n if n.np.dot(direction)<0 else n* -1
    f = obj.surfaceColor
    direction.norm()
    

    if (f.x > f.y and f.x > f.z):
        p = f.x
    elif(f.y > f.z):
        p = f.y
    else:
        p = f.z
    
    
    if (flag==false and obj.token != 'T'):
        return
    flag = true
    
    if (obj.token == 'D' or obj.token == 'G'):#                  // Ideal DIFFUSE reflection 
        if (depth>5 or np.random.RandomState(0)>p):#   //absorb the photon
            color_tmp = [f.x*color.x*(1.0/(1-p)),f.y*color.y*(1.0/(1-p)),f.z*color.z*(1.0/(1-p))]
            pos_tmp = [origin_new.x,origin_new.y,origin_new.z]
            dir_tmp = [direction.x,direction.y,direction.z]

            pmap.store(color_tmp,pos_tmp,dir_tmp)
            return
                        

        else: #  // trace another ray
            r1=2*PI*np.random.RandomState(0)
            r2=np.random.RandomState(0)
            r2s=np.sqrt(r2)
            
            w=nl
            u=np.mod(([0,1,0] if fabs(w.x)>.1 else [1,0,0]),w).norm()
            v=np.mod(w,u) 
            d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm()
            photon_tracing_caustic(pmap, origin_new,d,triangle_list, f.np.multiply(color)*(1.0/p), depth, flag,RI, bound)

    elif (obj.token == 'S'):#    // Ideal SPECULAR reflection 
        photon_tracing_caustic(pmap,origin_new,np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction)))),triangle_list,f.np.multiply(color)*(1.0/p),depth, flag,RI, bound)
        
    else if (obj.token == 'T'):  

        reflRay(np.substract(direction,np.multiply(n,np.multiply(2, n.np.dot(direction)))))#     // Ideal dielectric REFRACTION 
        bool into = n.np.dot(nl)>0 #                // Ray from outside going in? 
        nc=1
        nt=RI
        nnt=np.divide(nc,nt) if into else np.divide(nt,nc)
        ddn=direction.np.dot(nl) 
        if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0):#    {// Total internal reflection

            photon_tracing(pmap,origin_new,reflRay,triangle_list, color, depth, RI) 

        tdir = np.substract(np.multiply(direction,nnt),np.multiply(n,np.multiply((1 if into else -1),np.add(np.multiply(ddn,nnt), np.sqrt(cos2t))))).norm()
        a = np.substract(nt,nc)
        b = np.add(nt,nc)
        R0 = np.divide(np.multiply(a,)a,np.multiply(b,b))
        c = np.substract(1, (~ddn if into else tdir.np.dot(n)) )
        Re = np.add(R0,np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(np.substract(1,R0), c),c) c)  c) c))
        Tr = np.substract(1,Re)
        P=np.add(.25,np.multiply(.5,Re))
        RP=np.divide(Re,P)
        TP=np.divide(Tr,np.substract(1,P))
        np.random.RandomState(0)<P #?   // Russian roulette 
        photon_tracing_caustic(pmap, origin_new,reflRay,triangle_list,np.multiply(color,RP),depth, flag,RI, bound)
        photon_tracing_caustic(pmap, origin_new,tdir,triangle_list,np.multiply(color,TP), depth,flag,RI, bound)

    
# List of objects.
color_plane0 = 1. * np.ones(3)
color_plane1 = 0. * np.ones(3)
scene = [add_sphere([.75, .1, 1.], .6, [0., 0., 1.]),
         add_sphere([-.75, .1, 2.25], .6, [.5, .223, .5]),
         add_sphere([-2.75, .1, 3.5], .6, [1., .572, .184]),
         add_plane([0., -.5, 0.], [0., 1., 0.]),
    ]

# Light position and color.
L = np.array([5., 5., -10.])
color_light = np.ones(3)

# Default light and material parameters.
ambient = .05
diffuse_c = 1.
specular_c = 1.
specular_k = 50

depth_max = 5  # Maximum number of light reflections.
col = np.zeros(3)  # Current color.
O = np.array([0., 0.35, -1.])  # Camera.
Q = np.array([0., 0., 0.])  # Camera pointing to.
img = np.zeros((h, w, 3))

r = float(w) / h
# Screen coordinates: x0, y0, x1, y1.
S = (-1., -1. / r + .25, 1., 1. / r + .25)

# Loop through all pixels.
for i, x in enumerate(np.linspace(S[0], S[2], w)):
    if i % 10 == 0:
        print i / float(w) * 100, "%" # 进度
    for j, y in enumerate(np.linspace(S[1], S[3], h)):
        col[:] = 0
        Q[:2] = (x, y)
        D = normalize(Q - O)
        depth = 0
        rayO, rayD = O, D
        reflection = 1.
        # Loop through initial and secondary rays.
        while depth < depth_max:
            traced = trace_ray(rayO, rayD)
            if not traced:
                break
            obj, M, N, col_ray = traced
            # Reflection: create a new ray.
            rayO, rayD = M + N * .0001, normalize(rayD - 2 * np.dot(rayD, N) * N)
            depth += 1
            col += reflection * col_ray
            reflection *= obj.get('reflection', 1.)
        img[h - j - 1, i, :] = np.clip(col, 0, 1)

plt.imsave('fig.png', img)
