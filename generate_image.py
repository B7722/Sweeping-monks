#
# File that used to generate images, it was called main.py
# Parallel computing unfinished.
# Use gpu and multiple cpus for acceleration in the future.
# Tracing algorithm needs to be improved as well.
#

import funclib
import params

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
#from timeit import default_timer as timer


class Piece:
    def __init__(self, index, val):
        self.index = index
        self.val = val


def computepiece(inputdata):
    col = np.zeros(3) 
    Q = np.zeros(3)
    array = np.zeros((info.width, 3))
    for i, x in enumerate(np.linspace(info.Screen[0], info.Screen[2], info.width)):
        col[:] = 0 
        Q[:] = info.Camera_direction[:]
        Q[:2] = (x, inputdata[1])
        D = funclib.normalize(Q - info.Camera)
        depth = 0
        rayO, rayD = info.Camera, D
        reflection = 1.
        # Loop through initial and secondary rays.

        while depth < info.depth_max:
            traced = funclib.trace_ray(rayO, rayD, info)
            if not traced:
                break
            obj, M, N, col_ray = traced
            # Reflection: create a new ray.
            rayO, rayD = M + N * .0001, funclib.normalize(rayD - 2 * np.dot(rayD, N) * N)
            depth += 1
            col += reflection * col_ray
            reflection *= obj.get('reflection', 1.)
        array[i, :] = np.clip(col, 0, 1)
    return Piece(inputdata[0], array)
    
global info
info = params.Parameters()
img = np.zeros((info.height, info.width, 3))

#start = timer()
# Loop through all pixels.
inputdata = []
mapping = {};
for j, y in enumerate(np.linspace(info.Screen[1], info.Screen[3], info.height)):
    d = (j, y)
    inputdata.append(d)

p = multiprocessing.Pool(multiprocessing.cpu_count())
pieces = p.map(computepiece, inputdata)
#print(timer()-start)

for piece in pieces:
    img[info.height - piece.index - 1,:,:] = piece.val

plt.imsave('fig.png', img)

