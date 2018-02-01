#
# Main file, Beta version
# Parallel computing unfinished.
# Use gpu and multiple cpus for acceleration in the future.
# Tracing algorithm needs to be improved as well.
#

import funclib
import objlib
import parameters

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from timeit import default_timer as timer


class Data:
    def __init__(self, parameters, index, val):
        self.parameters = parameters
        self.index = index
        self.val = val


def computepixel(inputdata):
    col = np.zeros(3) 
    Q = np.zeros(3)
    array = np.zeros((inputdata.parameters.width, 3))
    for i, x in enumerate(np.linspace(inputdata.parameters.Screen[0], inputdata.parameters.Screen[2], inputdata.parameters.width)):
        col[:] = 0 
        Q[:] = inputdata.parameters.Camera_direction[:]
        Q[:2] = (x, inputdata.val)
        D = FunctionsLib.normalize(Q - inputdata.parameters.Camera)
        depth = 0
        rayO, rayD = inputdata.parameters.Camera, D
        reflection = 1.
        # Loop through initial and secondary rays.

        while depth < inputdata.parameters.depth_max:
            traced = FunctionsLib.trace_ray(rayO, rayD, inputdata.parameters)
            if not traced:
                break
            obj, M, N, col_ray = traced
            # Reflection: create a new ray.
            rayO, rayD = M + N * .0001, FunctionsLib.normalize(rayD - 2 * np.dot(rayD, N) * N)
            depth += 1
            col += reflection * col_ray
            reflection *= obj.get('reflection', 1.)
        array[i, :] = np.clip(col, 0, 1)
    return array
    
info = parameters.Parameters()
img = np.zeros((info.height, info.width, 3))

start = timer()
# Loop through all pixels.
inputdata = []
mapping = {};
for j, y in enumerate(np.linspace(info.Screen[1], info.Screen[3], info.height)):
    d = Data(info, j, y)
    temp = computepixel(d)
    mapping[j] = temp;
    #inputdata.append(d)

#p = multiprocessing.Pool(multiprocessing.cpu_count())
#data = p.map(computepixel, inputdata)

for j in range(info.height):
    temp = mapping[j]
    img[info.height - j - 1, :, :] = temp[:,:]

print(timer()-start)
plt.imsave('fig.png', img)
