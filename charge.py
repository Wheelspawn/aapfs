#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 17:31:33 2025

@author: nsage
"""

from vispy import scene

import typing
import numpy as np

import sys

from vispy import scene
from vispy.color import Color

CONST_ELECTRIC = 8.854 * 10**-12

class ChargedParticle(object):

    def __init__(self, pos: (float, float), charge: float):
        self.pos = np.array(pos)
        self.charge = charge

def computeElectricFieldAtR(charged_particle: list[ChargedParticle], r=(float,float)):
    
    ef = 1/(4*np.pi*CONST_ELECTRIC) * sum([p.charge/np.dot(np.abs(np.array(r)-p.pos),
                                                           np.abs(np.array(r)-p.pos)) for p in charged_particle])
    
    if (ef == np.inf):
        return np.nan
    else:
        return ef

p0 = ChargedParticle(pos=(0.0,0.0), charge=1.0)
p1 = ChargedParticle(pos=(-1.0,1.0), charge=-1.0)
p2 = ChargedParticle(pos=(0.0,1.0), charge=-1.0)
p3 = ChargedParticle(pos=(1.0,1.0), charge=1.0)

grid = np.array([[computeElectricFieldAtR([p0], (i,j))
                  for i in range(-2,3)] for j in range(2,-3,-1)],dtype=np.float32)

for r in grid:
    print(r)


canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)

# Set up a viewbox to display the cube with interactive arcball
view = canvas.central_widget.add_view()
view.bgcolor = '#efefef'
view.camera = 'turntable'
view.padding = 100

color = Color("#3f51b5")

cube = scene.visuals.Box(1, 1, 1, color=color, edge_color="black",
                         parent=view.scene)
if __name__ == '__main__' and sys.flags.interactive == 0:
    canvas.app.run()
    print('bluh')
