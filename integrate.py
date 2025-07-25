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
from enum import Enum
import shapely
from shapely.geometry import LineString, Point
import copy

# which way the springy end of the spring points
# e.g. if the direction is TOP then if something is dropped onto the spring it activates the spring force and bounces
class SpringDirection(Enum):
    TOP = 1
    LEFT = 2
    BOTTOM = 3
    RIGHT = 4

class TestMesh():
    def __init__(self):
        self.verts = pos = np.array([[0.0, 0.0, 0.0],
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
    def __init__(self, name, verts=np.zeros(3), mass=1.0, lin_vel=np.zeros(3), ang_vel=np.zeros(3), fixed=False):
        self.name = name
        self.verts = verts
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        self.fixed = fixed
    
    def apply_force(self, force: np.array):
        if self.fixed == True:
            return
    
class Cube(PhysicalMesh):
    def __init__(self, name, verts=np.zeros((4,2)), mass=1.0, lin_vel=np.zeros(2), ang_vel=np.zeros(2), fixed=False):
        self.name = name
        self.verts = verts
        self.mass = mass
        self.lin_vel = lin_vel
        self.ang_vel = ang_vel
        self.fixed = fixed
    
    def apply_force(self, lin_force: np.array):
        if self.fixed == True:
            return
        
        self.lin_vel += lin_force / self.mass
        
    def apply_velocity(self, lin_vel: np.array):
        self.lin_vel += lin_vel
    
    def set_velocity(self, lin_vel: np.array):
        self.lin_vel = lin_vel    
    
    def displace(self, distance):
        self.verts += distance
    
    def center(self):
        # print([sum(self.verts[:,i])/self.verts.shape[0] for i in range(self.verts.shape[1])])
        return [sum(self.verts[:,i])/self.verts.shape[0] for i in range(self.verts.shape[1])]
    
    def height(self):
        # print("height: ", max(self.verts[:,1]) - min(self.verts[:,1]))
        return max(self.verts[:,1]) - min(self.verts[:,1])
    
    def width(self):
        # print("width: ", max(self.verts[:,0]) - min(self.verts[:,0]))
        return max(self.verts[:,0]) - min(self.verts[:,0])
    
    def bounding_box(self):
        return np.array([ [self.verts[:,0].min(),self.verts[:,1].min()],
                          [self.verts[:,0].max(),self.verts[:,1].min()],
                          [self.verts[:,0].max(),self.verts[:,1].max()],
                          [self.verts[:,0].min(),self.verts[:,1].max()] ] )

    
class FixedSpringCube(PhysicalMesh):
    def __init__(self, name, verts=np.zeros((4,2)), mass=1.0, lin_vel=np.zeros((4,2)), fixed=np.array([[0,0,0,0]]), k=1.0):
        self.name = name
        self.verts = verts
        self.orig_pos = copy.deepcopy(verts)
        self.mass = mass
        self.lin_vel = lin_vel
        self.fixed = fixed
        self.k = k
    
    def apply_force(self, force: np.array):
        self.lin_vel += (force / self.mass)
            
    def apply_spring_force(self):
        force = self.k * (self.orig_pos - self.verts)
        # print(force)
        self.lin_vel += (force / self.mass)[0]
    
    def apply_velocity(self, lin_vel: np.array):
        self.lin_vel += lin_vel
    
    def set_velocity(self, lin_vel: np.array):
        self.lin_vel = lin_vel
                
    def displace(self, distance):
        # print("verts: ", self.verts)
        # print("distance: ", distance)
        for i in range(len(self.fixed)):
            if self.fixed[i] != 1:
                # print(i)
                # print("diff: ", self.verts[i] + distance[i])
                self.verts[i] += distance
        # print()
            
    def center(self):
        # print([sum(self.verts[:,i])/self.verts.shape[0] for i in range(self.verts.shape[1])])
        return [sum(self.verts[:,i])/self.verts.shape[0] for i in range(self.verts.shape[1])]
    
    def height(self):
        # print("height: ", max(self.verts[:,1]) - min(self.verts[:,1]))
        return max(self.verts[:,1]) - min(self.verts[:,1])
    
    def width(self):
        # print("width: ", max(self.verts[:,0]) - min(self.verts[:,0]))
        return max(self.verts[:,0]) - min(self.verts[:,0])
    
    def bounding_box(self):
        return np.array([ [self.verts[:,0].min(),self.verts[:,1].min()],
                          [self.verts[:,0].max(),self.verts[:,1].min()],
                          [self.verts[:,0].max(),self.verts[:,1].max()],
                          [self.verts[:,0].min(),self.verts[:,1].max()] ] )

def intersects(m1: PhysicalMesh, m2: PhysicalMesh):
    for i in range(len(m1.verts)):
        for j in range(len(m2.verts)):
            
            # these two points in a and b make line segments.
            # by iterating each line segment in a and b and checking to see if the lines intersect,
            # we obtain a way to check if any arbitrary polygon intersects with another.
            # although this is in m*n time, it is still relatively efficient if the polygons
            # are small enough (triangles, cubes, etc).
            # We have also added the optimization trick of first checking the bounding sphere,
            # which eliminates this step if the meshes are not nearly close to each other.
            
            a_p1 = m1.verts[i]
            a_p2 = m1.verts[(i+1) % m1.verts.shape[0]]
            
            b_p1 = m2.verts[j]
            b_p2 = m2.verts[(j+1) % m2.verts.shape[0]]
            
            a = LineString([a_p1,a_p2])
            b = LineString([b_p1,b_p2])
            
            if shapely.intersects(a,b) == True:
                # print("True:  ", a_p1, a_p2, b_p1, b_p2)
                return True
            else:
                return False
                # print("False:  ", a_p1, a_p2, b_p1, b_p2)
            
    return False

def adjust_collision(m1: PhysicalMesh, m2: PhysicalMesh, dt, e=0.001):
    d = shapely.distance(shapely.Polygon(m1.verts),shapely.Polygon(m2.verts))
    i = shapely.intersection(shapely.Polygon(m1.verts),shapely.Polygon(m2.verts))
    # nothing to adjust; doesn't intersect
    if (d >= 0 and i.area == 0):
        return
    
    # get previous position
    m1_p0 = copy.deepcopy(m1.verts - m1.lin_vel*dt)
    m2_p0 = copy.deepcopy(m2.verts - m2.lin_vel*dt)
    
    m1_p1 = copy.deepcopy(m1.verts)
    m2_p1 = copy.deepcopy(m2.verts)
    
    # m1_p0_old = m1.verts - m1.lin_vel
    # m2_p0_old = m2.verts - m2.lin_vel
    # m1_p1_old = m1.verts
    # m2_p1_old = m2.verts
    
    while (d >= e or i.area > 0):
        
        # m1.verts = (m1_p1 + m1_p0)/2
        # m2.verts = (m2_p1 + m2_p0)/2
        
        m1.displace(-(m1_p1 - m1_p0)/2)
        m2.displace(-(m2_p1 - m2_p0)/2)
        
        m1_pm = copy.deepcopy(m1.verts)
        m2_pm = copy.deepcopy(m2.verts)
        
        print("m1: ", m1.name)
        print("m2: ", m2.name)
        
        print("m1_p0: ", m1_p0)
        print("m1_pm: ", m1_pm)
        print("m1_p1: ", m1_p1)
        
        print("")
        
        print("m2_p0: ", m2_p0)
        print("m2_pm: ", m2_pm)
        print("m2_p1: ", m2_p1)
        
        print("")
        
        d = shapely.distance(shapely.Polygon(m1.verts),shapely.Polygon(m2.verts))
        i = shapely.intersection(shapely.Polygon(m1.verts),shapely.Polygon(m2.verts))
        area = i.area
        
        print("d: ", d)
        print("area: ", area)
        
        print("\n")
        
        if np.abs(np.sum(m1_p0)) > 40:
            return
        
        if np.abs(np.sum(m1_pm)) > 40:
            return
        
        if np.abs(np.sum(m1_p1)) > 40:
            return
        
        if np.abs(np.sum(m1_p0)) > 40:
            return
        
        if np.abs(np.sum(m2_pm)) > 40:
            return
        
        if np.abs(np.sum(m2_p1)) > 40:
            return
        
        if d >= e:
            m1_p0 = copy.deepcopy(m1.verts)
            m2_p0 = copy.deepcopy(m2.verts)
        else:
            m1_p1 = copy.deepcopy(m1.verts)
            m2_p1 = copy.deepcopy(m2.verts)
    
    # the collision happens somewhere between t and t+1. We need to return this
    # number so we know how to apply the resulting forces properly in the integrator.
    # the ratio of the differences in magnitudes gives us this number.
    # return (np.linalg.norm(m1.verts)-np.linalg.norm(m1_p0_old))/(np.linalg.norm(m1_p1_old)-np.linalg.norm(m1_p0_old))
    
def integrate(meshes: list[PhysicalMesh]):
    
    dt = 0.04
    
    for a in range(len(meshes)):
        a_pos_old = copy.deepcopy(meshes[a].verts)
        if type(meshes[a]) == FixedSpringCube:
            meshes[a].apply_spring_force()
        meshes[a].displace(meshes[a].lin_vel * dt)
        # meshes[a].angle += np.linalg.norm(meshes[a].ang_vel) * dt
    
        for b in range(len(meshes)):
            if (a != b and intersects(meshes[a], meshes[b])):
                
                # print("intersects")
                
                # print("a: ", meshes[a].name)
                # print("b: ", meshes[b].name)
                
                # adjust_collision(meshes[a],meshes[b],dt)
                meshes[a].displace(-meshes[a].lin_vel * dt)
                meshes[b].displace(-meshes[b].lin_vel * dt)
                
                # meshes[a].displace(meshes[a].lin_vel * col)
                # meshes[b].displace(meshes[b].lin_vel * col)
                
                # from conservation of momentum
                # m_a * a_v1 + m_b * b_v1 =
                # m_a * a_v2 + m_b * b_v2
                # and
                # b_v2 - a_v2 = -(b_v1 - a_v1)
                # or 
                # b_v2 = a_v2 - (b_v1 - a_v1)
                # and
                # a_v2 = b_v2 + (b_v1 - a_v1)
                
                # m_a * a_v1 + m_b * b_v1 = m_a * (b_v2 + (b_v1 - a_v1)) + m_b * b_v2
                # m_a * a_v1 + m_b * b_v1 = m_a * b_v2 + m_a * (b_v1 - a_v1) + m_b * b_v2
                # m_a * a_v1 + m_b * b_v1 - m_a * (b_v1 - a_v1) = m_a * b_v2 + m_b * b_v2
                # m_a * a_v1 + m_b * b_v1 - m_a * (b_v1 - a_v1) = b_v2 * (m_a + m_b)
                # (m_a * a_v1 + m_b * b_v1 - m_a * (b_v1 - a_v1) ) / (m_a + m_b) = b_v2
                # (m_a * a_v1 + m_b * b_v1 - m_a * b_v1 + m_a * a_v1) / (m_a + m_b) = b_v2
                
                a_v1 = meshes[a].lin_vel
                b_v1 = meshes[b].lin_vel
                
                v_b2_minus_v_a2 = -(b_v1 - a_v1)
                
                # rearranging
                a_v2 = ( 2 * ( meshes[b].mass * b_v1 ) + a_v1 * ( meshes[a].mass - meshes[b].mass ) ) / (meshes[a].mass + meshes[b].mass)
                
                # b_v2 = (meshes[a].mass * a_v1 + meshes[b].mass * b_v1 - meshes[a].mass * a_v2) / meshes[b].mass
                b_v2 = a_v2 - (b_v1 - a_v1)
                
                meshes[a].set_velocity(a_v2)
                meshes[b].set_velocity(b_v2)
                
                # print(meshes[a].lin_vel)
                # print(meshes[b].lin_vel)
                # print()
                
                # print("a: ", meshes[a].verts)
                # print("b: ", meshes[b].verts)
                # print()
                

