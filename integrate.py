#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:50:16 2025

@author: nsage
"""

import math
import typing
import numpy as np


class PhysicalMesh():
    def __init__(self, vertices = np.array([[0,0,0]]), pos=np.array([0,0,0])):
        self.vertices = vertices
        self.pos = pos

class Particle(PhysicalMesh):
    def __init__(self, vertices = np.array([[0,0,0]]), pos=np.array([0,0,0])):
        self.vertices = vertices
        self.pos = pos
        
class Circle(PhysicalMesh):
    def __init__(self, pos=np.array([0,0,0]), radius=1):
        self.pos = pos
        self.radius = radius
        self.vertices = np.array([[np.cos(math.radians(x*10))*self.radius*0.5,
                                   np.sin(math.radians(x*10))*self.radius*0.5,
                                   0.0] for x in range(0,36)]) + self.pos
    

def collides(a: Circle, b: Circle):
    return np.absolute(np.linalg.norm(a.pos-b.pos)) < a.radius + b.radius



def integrate():
    pass

a = Circle(pos=np.array([1,1,0]),radius=1)
b = Circle(pos=np.array([2,1,0]),radius=0.9)


