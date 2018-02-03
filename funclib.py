#
# Functions Library
# General functions will be implemented in this file
#

import numpy as np
import objlib

def normalize(x):
    x /= np.linalg.norm(x)
    return x

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
#  scene, L, O, ambient, diffuse_c, specular_c, specular_k, color_light
def trace_ray(rayO, rayD, parameters):
    # Find first point of intersection with the scene.
    t = np.inf
    for i, obj in enumerate(parameters.scene.objs):
        t_obj = intersect(rayO, rayD, obj)
        if t_obj < t:
            t, obj_idx = t_obj, i
    # Return None if the ray does not intersect any object.
    if t == np.inf:
        return
    # Find the object.
    obj = parameters.scene.objs[obj_idx]
    # Find the point of intersection on the object.
    M = rayO + rayD * t
    # Find properties of the object.
    N = get_normal(obj, M)
    color = get_color(obj, M)
    toL = normalize(parameters.Light - M)
    toO = normalize(parameters.Camera - M)
    # Shadow: find if the point is shadowed or not.
    l = [intersect(M + N * .0001, toL, obj_sh) 
            for k, obj_sh in enumerate(parameters.scene.objs) if k != obj_idx]
    if l and min(l) < np.inf:
        return
    # Start computing the color.
    col_ray = parameters.ambient
    # Lambert shading (diffuse).
    col_ray += obj.get('diffuse_c', parameters.diffuse_c) * max(np.dot(N, toL), 0) * color
    # Blinn-Phong shading (specular).
    col_ray += obj.get('specular_c', parameters.specular_c) * max(np.dot(N, normalize(toL + toO)), 0) ** parameters.specular_k * parameters.color_light
    return obj, M, N, col_ray

def intersect(O, D, obj):
    if obj['type'] == 'plane':
        return objlib.intersect_plane(O, D, obj['position'], obj['normal'])
    elif obj['type'] == 'sphere':
        return objlib.intersect_sphere(O, D, obj['position'], obj['radius'])
