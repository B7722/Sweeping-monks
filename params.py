#
# This file consists all the parameters needed to describe the 3D world.
# Defined default values only
# Need to add setters in the later version
#

import numpy as np
import objlib

class Parameters:
    def __init__(self):
        # Width & height of the output image
        self.width = 400
        self.height = 300

        # Parameters of the scene
        self.color_plane0 =  1. * np.ones(3)
        self.color_plane1 =  0. * np.ones(3)
        self.scene = [objlib.add_sphere([.75, .1, 1.], .6, [0., 0., 1.]),
                      objlib.add_sphere([-.75, .1, 2.25], .6, [.5, .223, .5]),
                      objlib.add_sphere([-2.75, .1, 3.5], .6, [1., .572, .184]),
                      objlib.add_plane([0., -.5, 0.], [0., 1., 0.],self.color_plane0, self.color_plane1),
                    ]

        # Light info
        self.Light = np.array([5., 5., -10.])
        self.color_light = np.ones(3)
        self.ambient = .05
        self.diffuse_c = 1.
        self.specular_c = 1.
        self.specular_k = 50
        self.depth_max = 5

        # Camera location & direction
        self.Camera = np.array([0., 0.35, -1.])
        self.Camera_direction = np.array([0., 0., 0.])

        # Screen coordinates: x0, y0, x1, y1.
        r = float(self.width) / self.height
        self.Screen = (-1., -1. / r + .25, 1., 1. / r + .25)
        
class Scene:
    def __init__(self):
        self.objs = []
    
    def add_obj(self, obj):
        self.objs.append(obj)
        

