#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:50:16 2025

@author: nsage
"""

import math
import typing
import numpy as np
from multimethod import multimethod

class PhysicalMesh():
    def __init__(self, name, pos=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3)):
        self.name = name
        self.pos = pos
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
    
    def apply_force(self, force: np.array):
        pass

class Particle(PhysicalMesh):
    def __init__(self, name, pos=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3)):
        self.name = name
        self.pos = pos
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        
    def apply_force(self, lin_force: np.array, ang_force: np.array):
        self.lin_vel += lin_force / self.mass
        self.ang_vel += ang_force / self.mass
    
class Circle(PhysicalMesh):
    def __init__(self, name, pos=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3), radius=1.0, angle=0.0):
        self.name = name
        self.pos = pos
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        self.radius = radius
        self.angle = angle
        ''' self.vertices = np.array([[np.cos(math.radians(x*10))*self.radius*0.5,
                                       np.sin(math.radians(x*10))*self.radius*0.5,
                                       0.0] for x in range(0,36)]) + self.pos '''
    
    def apply_force(self, lin_force: np.array, ang_force: np.array):
        self.lin_vel += lin_force / self.mass
        self.ang_vel += ang_force / self.mass

def get_circle_force_components(m: PhysicalMesh, force: np.array, force_pos: np.array):
    # projection of force onto direction towards center of mass (linear component of the force)
    c_of_mass = m.pos - force_pos
    proj_f = (np.dot(c_of_mass, force) / (np.linalg.norm(c_of_mass) * np.linalg.norm(c_of_mass))) * c_of_mass
    
    # rejection (rotational component of the force)
    rej = force - proj_f
    
    return (proj_f,rej)

@multimethod
def collides(c1: Circle, c2: Circle):
    return np.absolute(np.linalg.norm(c1.pos-c2.pos)) < c1.radius + c2.radius

@multimethod
def collides(c: Circle, p: Particle):
    return np.absolute(np.linalg.norm(c.pos-p.pos)) < c.radius

def integrate(meshes: list[PhysicalMesh]):
    
    dt = 0.1
    e = 0.9
    
    for i in range(len(meshes)):
        meshes[i].pos += meshes[i].lin_vel * dt
        meshes[i].angle += np.linalg.norm(meshes[i].ang_vel) * dt
    
        for j in range(i+1, len(meshes)):
            if collides(meshes[i], meshes[j]):
                
                force_pos = meshes[j].pos + meshes[j].radius * (meshes[i].pos - meshes[j].pos)
                
                fc_1 = get_circle_force_components(meshes[i],meshes[i].lin_vel,force_pos)
                fc_2 = get_circle_force_components(meshes[j],meshes[j].lin_vel,force_pos)
                
                v_com = ((meshes[i].mass * meshes[i].lin_vel) + (meshes[j].mass * meshes[j].lin_vel))/(meshes[i].mass + meshes[j].mass)
                
                i_v2 = (1 + e) * v_com - e * meshes[i].lin_vel
                j_v2 = (1 + e) * v_com - e * meshes[j].lin_vel
                
                meshes[i].lin_vel = i_v2
                meshes[j].lin_vel = j_v2
                
                meshes[i].ang_vel = (1 + e) * v_com - e * meshes[i].ang_vel
                meshes[j].ang_vel = (1 + e) * v_com - e * meshes[j].ang_vel
                
                '''
                print(meshes[i].name, ", ",
                      meshes[j].name, "collide at: ",
                      meshes[i].pos, ", ",
                      meshes[j].pos)
                print("force_pos: ", force_pos)
                print("meshes[i].lin_vel: ", meshes[i].lin_vel)
                print("meshes[j].lin_vel: ", meshes[j].lin_vel)
                print(meshes[i].name, "force components: ", fc_1)
                print(meshes[j].name, "force components: ", fc_2)
                print("mesh i v2: ", i_v2)
                print("mesh j v2: ", j_v2)'''
        
        # print(meshes[0].name, " pos: ", meshes[0].pos, "   ", meshes[1].name, " pos: ", meshes[1].pos)

i = Circle(name="circle1",pos=np.array([0.5,0.5,0]),radius=0.5)
j = Circle(name="circle2",pos=np.array([1.25,2.5,0]),lin_vel=np.array([0.0,-0.15,0.0]),radius=0.5)

# force_a = np.array([0.0,0.5,0.0])
# force_a_pos = np.array([0.5,0.0,0.0])

for k in range(60):
    integrate([i,j])

