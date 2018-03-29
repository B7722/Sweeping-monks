import numpy as np
import matplotlib.pyplot as plt
import Data
objects_in_scene = Data.load('test.txt')

w = 800
h = 800


def normalize(x):
    x /= np.linalg.norm(x)
    return x

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
    # E12 = P2-P1
    # E23 = P3-P2
    # E03 = P3-P0
    # E01 = P1-P0
    # area123 = round(np.linalg.norm(np.cross(E12,E23))/2,3)
    # area023 = np.linalg.norm(np.cross(E03,E23))/2
    # area013 = np.linalg.norm(np.cross(E03,E01))/2
    # area012 = np.linalg.norm(np.cross(E01,E12))/2
    # area0 = round(area023+area012+area013,3)
    # print(area0,area123)
    # if(area123==area0):
    #     return True
    # else:
    #     return False

def intersect_tetrahedron(O,D,P1,P2,P3,P4):
    E12 = P2-P1
    E23 = P3-P2
    E34 = P4-P3
    E14 = P4-P1
    N123 = normalize(np.cross(E12,E23))
    N234 = normalize(np.cross(E34,E23))
    N124 = normalize(np.cross(E14,E12))
    N134 = normalize(np.cross(E14,E34))
    min_d = np.inf
    vector = np.array([0, 0, 0])
    e = 0

    denom = np.dot(D, N123)
    if (denom != 0):
        d = np.dot(P1 - O, N123) / denom
        if d > e:
            point = O + d * D
            if judge_in_triangle(point,P1,P2,P3):
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

    return dict(d = min_d,v = vector)

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



def intersect_plane(O, D, P, N):
    # Return the distance from O to the intersection of the ray (O, D) with the
    # plane (P, N), or +inf if there is no intersection.
    # O and P are 3D points, D and N (normal) are normalized vectors.
    denom = np.dot(D, N)
    if np.abs(denom) < 1e-6:
        d = np.inf
        return dict(d = d,v = np.array([0,0,0]))
    d = np.dot(P - O, N) / denom
    if d < 0:
        d = np.inf
        return dict(d = d,v = np.array([0,0,0]))
    return dict(d = d,v = N)


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
            if t0 < 0:
                t = t1
            else:
                t = t0
            # t= t1 if t0 < 0 else t0
            point = O + t * D
            v = point - S
            v = normalize(v)
            return dict(d = t,v = v)
    d = np.inf
    return dict(d = d,v = np.array([0,0,0]))


def intersect(O, D, obj):
    if obj['type'] == 'plane':
        return intersect_plane(O, D, obj['position'], obj['normal'])
    elif obj['type'] == 'sphere':
        return intersect_sphere(O, D, obj['position'], obj['radius'])
    elif obj['type'] == 'cube':
        return intersect_cube(O, D, obj['position'], obj['radius_X'], obj['radius_Y'], obj['radius_Z'])
    elif obj['type']=='tetrahedron':
        return intersect_tetrahedron(O,D,obj['p1'],obj['p2'],obj['p3'],obj['p4'])

# def get_normal(obj, M):
# 
#     # Find normal.
#     if obj['type'] == 'sphere':
#         N = normalize(M - obj['position'])
#     elif obj['type'] == 'plane':
#         N = obj['normal']
#     elif obj['type'] == 'cube':
#         # N = normalize(M - obj['position'])
#         C = obj['position']
#         T = obj['position']+obj['radius_X']+obj['radius_Y']+obj['radius_Z']
#         RX = obj['radius_X']
#         RY = obj['radius_Y']
#         RZ = obj['radius_Z']
#         if(np.dot(M-C,-RX)==0):
#             N = normalize(-RX)
#         elif(np.dot(M-C,-RY)==0):
#             N = normalize(-RY)
#         elif(np.dot(M-C,-RZ)==0):
#             N = normalize(-RZ)
#         elif (np.dot(M - T, RX) == 0):
#             N = normalize(RX)
#         elif (np.dot(M - T, RY) == 0):
#             N = normalize(RY)
#         elif (np.dot(M - T, RZ) == 0):
#             N = normalize(RZ)
#         else:
#             N = normalize(RZ)
#     return N


def get_color(obj, M):
    color = obj['color']
    if not hasattr(color, '__len__'):
        color = color(M)
    return color


def trace_ray(rayO, rayD):
    # Find first point of intersection with the scene.
    t = np.inf
    N = np.array([0,0,0])
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
    if np.array_equal(N,np.array([0,0,0])) :
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
    if l_obj :
        l = []
        for a,b in enumerate(l_obj):
            l.append(b['d'])
        if l and min(l) < np.inf:
            return
    # Start computing the color.
    col_ray = ambient
    # Lambert shading (diffuse).
    col_ray += obj.get('diffuse_c', diffuse_c) * max(np.dot(N, toL), 0) * color
    # Blinn-Phong shading (specular).
    col_ray += obj.get('specular_c', specular_c) * max(np.dot(N, normalize(toL + toO)), 0) ** specular_k * color_light
    return obj, M, N, col_ray


def add_sphere(position, radius, color):
    return dict(type='sphere', position=np.array(position),
                radius=np.array(radius), color=np.array(color), reflection=.5)


def add_plane(position, normal):
    return dict(type='plane', position=np.array(position),
                normal=np.array(normal),
                color=lambda M: (color_plane0
                                 if (int(M[0] * 2) % 2) == (int(M[2] * 2) % 2) else color_plane1),
                diffuse_c=.75, specular_c=.5, reflection=.25)


def add_cube(position, radius_X, radius_Y, radius_Z, color):
    return dict(type='cube', position=np.array(position), radius_X=np.array(radius_X), radius_Y=np.array(radius_Y),
                radius_Z=np.array(radius_Z), color=np.array(color), reflection=.5)

def add_tetrahedron(p1,p2,p3,p4,color):
    return dict(type='tetrahedron',p1 = np.array(p1),p2 = np.array(p2),p3 = np.array(p3),p4 = np.array(p4),color = np.array(color),reflection = .5)
# List of objects.
color_plane0 = 1. * np.ones(3)
color_plane1 = 0. * np.ones(3)
# scene = [add_tetrahedron([1,0.5,1],[1.75,.1,1],[1.25,.1,2],[1.5,1.1,1.5],[0., 0., 1.]),
#          add_cube([0, -0.5, 2], [.6, 0, 0], [0, 0.6, 0], [0, 0, .6], [0., 0., 1.]),
#          add_sphere([-.75, .1, 2.25], .6, [.5, .223, .5]),
#          add_sphere([-2.75, .1, 3.5], .6, [1., .572, .184]),
#          add_plane([0., -.5, 0.], [0., 1., 0.])
#          ]
scene = [add_plane([0., -.5, 0.], [0., 1., 0.])]
for i in range(len(objects_in_scene)):
    if (objects_in_scene[i].getName() == "Camera"):
        if(objects_in_scene[i].getPosition() !=[]):
            O = np.array(objects_in_scene[i].getPosition())
            Q = np.array(objects_in_scene[i].getOrientation())
    elif(objects_in_scene[i].getName() == "Tetrahedran"):
        point1 = objects_in_scene[i].getPoint1()
        point2 = objects_in_scene[i].getPoint2()
        point3 = objects_in_scene[i].getPoint3()
        point4 = objects_in_scene[i].getPoint4()
        color = objects_in_scene[i].getColor()
        scene.append(add_tetrahedron(point1,point2,point3,point4,color))
    elif(objects_in_scene[i].getName() =="Ball"):
        position = objects_in_scene[i].getPosition()
        r = objects_in_scene[i].getR()
        color = objects_in_scene[i].getColor()
        scene.append(add_sphere(position,r,color))
    elif(objects_in_scene[i].getName() =="Cuboid"):
        position = objects_in_scene[i].getPosition()
        x = objects_in_scene[i].getX()
        y = objects_in_scene[i].getY()
        z = objects_in_scene[i].getZ()
        color = objects_in_scene[i].getColor()
        scene.append(add_cube(position,x,y,z,color))
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
# O = np.array([-1., 0.35, -1.5])  # Camera.
# Q = np.array([0.5, 0., 0.])  # Camera pointing to.
img = np.zeros((h, w, 3))

r = float(w) / h
# Screen coordinates: x0, y0, x1, y1.
S = (-1., -1. / r + .25, 1., 1. / r + .25)

# Loop through all pixels.
for i, x in enumerate(np.linspace(S[0], S[2], w)):
    if i % 10 == 0:
        print
        i / float(w) * 100, "%"
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
