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

class TestMesh():
    def __init__(self):
        self.pos = pos = np.array([[0.0, 0.0, 0.0],
                                   [1.0, 0.0, 0.0],
                                   [1.0, 1.0, 0.0],
                                   [0.0, 1.0, 0.0]])
        
        self.edges = np.array([[0,1],
                               [1,2],
                               [2,3],
                               [3,0]])
        
        self.stiffness = np.array([[2.5],
                                   [2.5],
                                   [2.5],
                                   [2.5]])

class PhysicalMesh():
    def __init__(self, name, pos=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3), fixed=False):
        self.name = name
        self.pos = pos
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        self.fixed = fixed
    
    def apply_force(self, force: np.array):
        if self.fixed == True:
            return

class Particle(PhysicalMesh):
    def __init__(self, name, pos=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3), fixed=False):
        self.name = name
        self.pos = pos
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        self.fixed = fixed
        
    def apply_force(self, lin_force: np.array, ang_force: np.array):
        if self.fixed == True:
            return
        
        self.lin_vel += lin_force / self.mass
        self.ang_vel += ang_force / self.mass
    
class Circle(PhysicalMesh):
    def __init__(self, name, pos=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3), radius=1.0, angle=0.0, fixed = False):
        self.name = name
        self.pos = pos
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        self.radius = radius
        self.angle = angle
        self.fixed = fixed
        ''' self.vertices = np.array([[np.cos(math.radians(x*10))*self.radius*0.5,
                                       np.sin(math.radians(x*10))*self.radius*0.5,
                                       0.0] for x in range(0,36)]) + self.pos '''
    
    def apply_force(self, lin_force: np.array, ang_force: np.array):
        
        if self.fixed == True:
            return
        
        self.lin_vel += lin_force / self.mass
        self.ang_vel += ang_force / self.mass

def get_circle_force_components(m: PhysicalMesh, force: np.array, force_pos: np.array):
    # projection of force onto direction towards center of mass (linear component of the force)
    c_of_mass = m.pos - force_pos
    print("force: ", force)
    print("m.pos: ", m.pos)
    print("force_pos: ", force_pos)
    print("c_of_mass: ", c_of_mass)
    proj_f = (np.dot(c_of_mass, force) / (np.linalg.norm(c_of_mass) * np.linalg.norm(c_of_mass))) * c_of_mass + force_pos
    print("proj_f: ", proj_f)
    
    # rejection (rotational component of the force)
    rej = (force - proj_f) + force_pos
    print("rej: ", rej)
    
    return (proj_f,rej)

@multimethod
def closing(c1: Circle, c2: Circle):
    return np.absolute(np.linalg.norm(c1.pos-c2.pos)) - (c1.radius + c2.radius)

@multimethod
def collides(c1: Circle, c2: Circle):
    return closing(c1,c2) <= 0

@multimethod
def closing(c: Circle, p: Particle):
    return np.absolute(np.linalg.norm(c.pos-p.pos)) - c.radius
    
@multimethod
def collides(c: Circle, p: Particle):
    return closing(c,p) <= 0

def adjust_collision(m1: PhysicalMesh, m2: PhysicalMesh, dt, e=0.0001):
    c = closing(m1,m2)
    
    dt = 0.5
    
    # nothing to adjust; doesn't collide
    if c > 0:
        return
    
    # get previous position
    m1_p0 = m1.pos - m1.lin_vel*dt
    m2_p0 = m2.pos - m2.lin_vel*dt
    
    m1_p1 = m1.pos
    m2_p1 = m2.pos
    
    m1_p0_old = m1.pos - m1.lin_vel
    m2_p0_old = m2.pos - m2.lin_vel
    m1_p1_old = m1.pos
    m2_p1_old = m2.pos
    
    while (c >= e or c < 0):
        
        m1.pos = (m1_p1 - m1_p0)/2 + m1_p0
        m2.pos = (m2_p1 - m2_p0)/2 + m2_p0
        
        c = closing(m1,m2)
        if c > e:
            m1_p0 = m1.pos
            m2_p0 = m2.pos
        else:
            m1_p1 = m1.pos
            m2_p1 = m2.pos
    
    # the collision happens somewhere between t and t+1. We need to return this
    # number so we know how to apply the resulting forces properly in the integrator.
    return (np.linalg.norm(m1.pos)-np.linalg.norm(m1_p0_old))/(np.linalg.norm(m1_p1_old)-np.linalg.norm(m1_p0_old))

def integrate(meshes: list[PhysicalMesh], forces: list[np.array]):
    
    dt = 0.05
    e = 0.9
    
    for i in range(len(meshes)):
        i_pos_old = meshes[i].pos
        meshes[i].pos += meshes[i].lin_vel * dt
        # meshes[i].angle += np.linalg.norm(meshes[i].ang_vel) * dt
    
    for i in range(len(meshes)):
        for j in range(i+1, len(meshes)):
            
            if collides(meshes[i], meshes[j]):
                
                j_pos_old = meshes[j].pos
                
                # get the time interval of the collision
                col_dt = adjust_collision(meshes[i],meshes[j],dt)
                
                # contact point of meshes
                contact_point = meshes[i].pos + (meshes[j].pos - meshes[i].pos)/2
                print("contact_point: ", contact_point)
                print()
                
                # linear and rotational force vectors
                fc_1 = get_circle_force_components(meshes[j],meshes[i].lin_vel,contact_point)
                # fc_2 = get_circle_force_components(meshes[j],meshes[j].lin_vel,contact_point)
                
                forces.append( np.append( contact_point, meshes[i].pos ) )
                # forces.append( np.append( contact_point, meshes[j].pos ) )
                
                forces.append( np.append( contact_point, fc_1[0] ) )
                forces.append( np.append( contact_point, fc_1[1] ) )
                
                # forces.append( np.append( contact_point, fc_2[0] ) )
                # forces.append( np.append( contact_point, fc_2[1] ) )
                
                
                # velocity of center of mass
                v_com = ((meshes[i].mass * meshes[i].lin_vel) + (meshes[j].mass * meshes[j].lin_vel))/(meshes[i].mass + meshes[j].mass)
                
                # new velocities
                i_v2 = (1 + e) * v_com - e * meshes[i].lin_vel
                j_v2 = (1 + e) * v_com - e * meshes[j].lin_vel
                
                # new linear velocity
                meshes[i].lin_vel = i_v2
                meshes[j].lin_vel = j_v2
                
                # new angular velocity
                # meshes[i].ang_vel = (1 + e) * v_com - e * meshes[i].ang_vel
                # meshes[j].ang_vel = (1 + e) * v_com - e * meshes[j].ang_vel
                
                # position has been adjusted from meshes[i] to contact point over dt.
                # displace position by the velocity of the remaining period of time.
                meshes[i].pos += (1 - col_dt) * meshes[i].lin_vel
                meshes[j].pos += (1 - col_dt) * meshes[j].lin_vel
                
                

# i = Circle(name="circle1",pos=np.array([0.0, 0.0, 0.0]),lin_vel=np.array([0.0,0.2,0.0]),radius=1.0)
# j = Circle(name="circle2",pos=np.array([0.0, 1.7001, 0.0]),lin_vel=np.array([0.0,-0.1,0.0]),radius=1.0)

# print(closing(i,j))
# print(adjust_collision(i,j,1))
# print(closing(i,j))