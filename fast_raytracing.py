"""
GPU Accelerated Raytracing

@author: Jiaming

"""

import numpy as np
import matplotlib.pyplot as plt
import math
from numba import cuda, float32, int8

w = 400
h = 300

def add_plane(scene, pos, color, normal):
    for i in range(scene.shape[0]):
        if scene[i,0,0] != 0:
            continue
        else:
            scene[i,0,0] = 1
            scene[i,1] = pos
            scene[i,2] = color
            scene[i,3] = normal
            break
        
def add_sphere(scene, pos, color, radius):
    for i in range(scene.shape[0]):
        if scene[i,0,0] != 0:
            continue
        else:
            scene[i,0,0] = 2
            scene[i,1] = pos
            scene[i,2] = color
            scene[i,3,0] = radius
            break
    


@cuda.jit(device = True)
def normalize(x):
    temp = 0
    for i in range(len(x)):
        temp += x[i]**2
    temp = math.sqrt(temp)
    for i in range(len(x)):
        x[i] /= temp
    return x

@cuda.jit(device = True)
def dot(A,B):
    temp = 0
    for i in range(len(A)):
        temp += A[i] * B[i]
    return temp

@cuda.jit(device = True)
def get_color(obj):
    r = obj[2,0]
    g = obj[2,1]
    b = obj[2,2]
    return r,g,b

@cuda.jit(device = True)
def get_normal(obj, M):
    N = cuda.local.array(shape=(3), dtype=float32)
    if obj[0,0] == 1:
        return obj[3,0], obj[3,1], obj[3,2]
    elif obj[0,0] == 2:
        N[0] = M[0] - obj[1,0]
        N[1] = M[1] - obj[1,1]
        N[2] = M[2] - obj[1,2]
        normalize(N)
        return N[0], N[1], N[2]
    
# intersect
@cuda.jit(device = True)
def intersect(O, D, obj):
    P = cuda.local.array(shape=(3), dtype=float32)
    N = cuda.local.array(shape=(3), dtype=float32)
    S = cuda.local.array(shape=(3), dtype=float32)
    R = cuda.local.array(shape=(1), dtype=float32)
    
    if obj[0,0] == 1:
        P[0] = obj[1,0]
        P[1] = obj[1,1]
        P[2] = obj[1,2]
        
        N[0] = obj[3,0]
        N[1] = obj[3,1]
        N[2] = obj[3,2]
        return intersect_plane(O, D, P, N)
    
    elif obj[0,0] == 2:
        S[0] = obj[1,0]
        S[1] = obj[1,1]
        S[2] = obj[1,2]
        
        R[0] = obj[3,0]
        return intersect_sphere(O, D, S, R)
    
@cuda.jit(device = True)
def intersect_plane(O, D, P, N):
    T = cuda.local.array(shape=(3), dtype=float32)
    denom = dot(D, N)
    if abs(denom) < 1e-6:
        return -1
    T[0] = P[0] - O[0]
    T[1] = P[1] - O[1]
    T[2] = P[2] - O[2]
    d = dot(T, N) / denom
    if d < 0:
        return -1
    return d

@cuda.jit(device = True)
def intersect_sphere(O, D, S, R):
    a = dot(D, D)
    OS = cuda.local.array(shape=(3), dtype=float32)
    OS[0] = O[0] - S[0]
    OS[1] = O[1] - S[1]
    OS[2] = O[2] - S[2]
    b = 2 * dot(D, OS)
    c = dot(OS, OS) - R[0] * R[0]
    disc = b * b - 4 * a * c
    if disc > 0:
        distSqrt = math.sqrt(disc)
        q = (-b - distSqrt) / 2.0 if b < 0 else (-b + distSqrt) / 2.0
        t0 = q / a
        t1 = c / q
        t0, t1 = min(t0, t1), max(t0, t1)
        if t1 >= 0:
            return t1 if t0 < 0 else t0
    return -1
    

@cuda.jit(device = True)
def trace_ray(rayO, rayD, obj, M, N, col_ray):
    toL = cuda.local.array(shape=(3), dtype=float32)
    toO = cuda.local.array(shape=(3), dtype=float32)
    LnO = cuda.local.array(shape=(3), dtype=float32)
    t = intersect(rayO, rayD, obj)
    for i in range(len(M)):
        M[i] = rayO[i] + rayO[i] * t
    # Find properties of the object.
    a1,a2,a3 = get_normal(obj, M)
    N[0] = a1
    N[1] = a2
    N[2] = a3
    color = get_color(obj)
    for i in range(len(toL)):
        toL[i] = L[i]
        toL[i] -= M[i]
    for i in range(len(toO)):
        toO[i] = O[i]
        toO[i] -= M[i]
    tmp = normalize(toL)
    toL = tmp
    tmp = normalize(toO)
    toO = tmp
    for i in range(len(toO)):
        LnO[i] = toL[i] +  toO[i]
    tmp = normalize(LnO)
    LnO = tmp

    # Lambert shading (diffuse).
    if obj[0,0] == 1:
        diff_c = .75
        spec_c = 0.5
    elif obj[0,0] == 2:
        diff_c = diffuse_c
        spec_c = specular_c
           
    for i in range(len(col_ray)):
        col_ray[i] = diff_c * max(dot(N, toL), 0) * color[i] + ambient
        col_ray[i] += spec_c * max(dot(N, LnO), 0) ** specular_k * color_light[i]
        
    return col_ray

@cuda.jit(device = True)
def test():
    if True:
        return 1,2,3,4
    else:
        return

@cuda.jit(device = True)
def load_object(obj, scene, idx):
    for a in range(len(scene[0])):
        for b in range(len(scene[0,1])):
            obj[a,b] = scene[idx,a,b]
    return obj


@cuda.jit
def fast_compute(img,xy,O,scene):
    col_d = cuda.shared.array(shape=(3), dtype=float32)
    Q_d = cuda.shared.array(shape=(3), dtype=float32)
    D_d = cuda.shared.array(shape=(3), dtype=float32)
    depth_d = cuda.shared.array(shape=(1), dtype=int8)
    reflection_d = cuda.shared.array(shape=(1), dtype=float32)
    rayO_d = cuda.shared.array(shape=(3), dtype=float32)
    rayD_d = cuda.shared.array(shape=(3), dtype=float32)
    obj = cuda.shared.array(shape=(numofattributes,numofdimensions), dtype=float32)
    obj_sh = cuda.shared.array(shape=(numofattributes,numofdimensions), dtype=float32)
    M = cuda.local.array(shape=(3), dtype=float32)
    N = cuda.local.array(shape=(3), dtype=float32)
    toL = cuda.local.array(shape=(3), dtype=float32)
    toO = cuda.local.array(shape=(3), dtype=float32)
    col_ray = cuda.local.array(shape=(3), dtype=float32)
    ray_sh = cuda.local.array(shape=(3), dtype=float32)


    i,j = cuda.grid(2)
    if i >= img.shape[0] and j >= img.shape[1]:
        return
    
    col_d[0] = 1.
    col_d[1] = 1.
    col_d[2] = 1.
    
    Q_d[0] = xy[i, j, 0]
    Q_d[1] = xy[i, j, 1]
    Q_d[2] = 0.
    
    D_d[0] = Q_d[0] - O[0]
    D_d[1] = Q_d[1] - O[1]
    D_d[2] = Q_d[2] - O[2]
    D_d = normalize(D_d) 
    
    
    rayO_d[0] = O[0]
    rayO_d[1] = O[1]
    rayO_d[2] = O[2]
    rayD_d = D_d
    reflection_d[0] = 1.
    cuda.syncthreads()
    

    #trace_ray(rayO_d, rayD_d, obj)
    t = -1
    idx = -1
    for i in range(len(scene)):
        load_object(obj, scene, i)
        dist = intersect(rayO_d, rayD_d, obj)
        if dist != -1:
            t = dist
            idx = i
    cuda.syncthreads()
    # Return None if the ray does not intersect any object.
    if t == -1:
        img[i,j,0] = 0
        img[i,j,1] = 0
        img[i,j,2] = 0
    else:
        load_object(obj, scene, idx)
        # Find the point of intersection on the object.
        toL = cuda.local.array(shape=(3), dtype=float32)
        toO = cuda.local.array(shape=(3), dtype=float32)
    
        t = intersect(rayO_d, rayD_d, obj)
        for i in range(len(M)):
            M[i] = rayO_d[i] + rayO_d[i] * t
            # Find properties of the object.
            a1,a2,a3 = get_normal(obj, M)
            N[0] = a1
            N[1] = a2
            N[2] = a3
            for i in range(len(toL)):
                toL[i] = L[i]
                toL[i] -= M[i]
            for i in range(len(toO)):
                toO[i] = O[i]
                toO[i] -= M[i]
            toL = normalize(toL)
            toO = normalize(toO)
        # Shadow: find if the point is shadowed or not.
        l = cuda.local.array(shape=(10000), dtype=float32)
        count = 0
        for k in range(len(scene)):
            if k != idx:
                load_object(obj_sh, scene, k)
                for a in range(len(ray_sh)):
                    ray_sh[a] = M[a] + N[a] * .0001
            l[count] = intersect(ray_sh, toL, obj_sh)
            count += 1
        for k in range(count):
            if l[k] == -1:
                l[k] = 10000
        minimum = -1
        if count > 0:
            minimum = l[0]
            for k in range(count):
                if minimum < l[k]:
                    minimum = l[k]
        if count > 0 and minimum < 10000:
            img[i,j,0] = 0
            img[i,j,1] = 0
            img[i,j,2] = 0
        else:
            while depth_d[0] < depth_max:
                col_ray = trace_ray(rayO_d, rayD_d, obj, M, N, col_ray)
                for i in range(len(rayO_d)):
                    rayO_d[i] += M[i] + N[i] * .0001
                tmp =  2 * dot(rayD_d, N)
                for i in range(len(rayD_d)):
                    rayD_d[i] -= (tmp * N[i])
                    rayD_d = normalize(rayD_d)
                depth_d[0] += 1
                for i in range(len(col_d)):
                    col_d[i] +=  col_ray[i]
            img[i,j,0] = col_d[0]
            img[i,j,1] = col_d[1]
            img[i,j,2] = col_d[2]
        
              
    cuda.syncthreads()
 

    
# List of objects. Each object is encoded into four attributes(Type, Position, Color, Radius/Normal), each attributes has three dimensions
maxnumofobjects = 10
numofattributes = 4
numofdimensions = 3
scene = np.zeros((maxnumofobjects,numofattributes,numofdimensions))


add_plane(scene, np.array([0., -.5, 0.]), np.array([1., 1., 1.]), np.array([0., 1., 0.]))
add_sphere(scene, np.array([.75, .1, 1.]), np.array([0., 0., 1.]) , .6)
add_sphere(scene, np.array([-.75, .1, 2.25]), np.array([.5, .223, .5]) , .6)
add_sphere(scene, np.array([-2.75, .1, 3.5]), np.array([1., .572, .184]) , .6)
#print scene

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


xy = np.zeros((h, w, 2))
for j, y in enumerate(np.linspace(S[1], S[3], h)):
    for i, x in enumerate(np.linspace(S[0], S[2], w)):
        xy[j,i,0] = x
        xy[j,i,1] = y

blockspergrid = (16,16)
threadsperblock_x = int(round(1.0*h/16 + 0.5))
threadsperblock_y = int(round(1.0*w/16 + 0.5))
threadsperblock = (threadsperblock_x, threadsperblock_y)
fast_compute[blockspergrid,threadsperblock](img,xy,O,scene)
#print img

plt.imsave('figGPU.png', img)

