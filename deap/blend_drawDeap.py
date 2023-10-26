#!/usr/bin/python
# -*- coding: utf-8 -*-
# works in Blender 3.3.12
import bpy
from bpy.ops import mesh
import math
from numpy import *
import numpy as np
import os, sys
import csv

rpsup = 8390
rpmt = 125 # hamamatsu R1408
rneck = 732
scale = 100

rpsup = rpsup/scale
rpmt = rpmt/scale
rneck = rneck/scale

PI = np.pi

# Function to create a cone at a specific location
def create_cone(location, orient):
    dirx = orient[0]
    diry = orient[1]
    dirz = orient[2]
    costheta = dirz
    theta = arccos(costheta)
    phi = arctan2(diry,dirx)
    
    mesh.primitive_cone_add(
        radius1=rpmt,
        radius2=rpmt,
        depth=1,
        end_fill_type='NGON',
        calc_uvs=True,
        enter_editmode=False,
        align='WORLD',
        location=location, 
        rotation = (0, theta, phi), ## Euler: rotate around y for theta, then z for phi
        scale=(1, 1, 1)
    )
    obj = bpy.context.active_object
    obj.name = 'pmtdisc'
    obj.data.name = 'pmtdisc'
    bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
    
def create_neck(location):    
    mesh.primitive_cone_add(
        radius1=rneck,
        radius2=rneck,
        depth=120, # 12801.6 mm in length
        end_fill_type='NGON',
        calc_uvs=True,
        enter_editmode=False,
        align='WORLD',
        location=location, 
        scale=(1, 1, 1)
    )
    obj = bpy.context.active_object
    obj.name = 'neck'
    obj.data.name = 'neck'
    bpy.ops.object.transform_apply(location=True, rotation=False, scale=True) 
    
    
def create_sphere(location):    
    bpy.ops.mesh.primitive_uv_sphere_add(
    radius=6000./scale,
    calc_uvs=True,
    enter_editmode=False,
    align='WORLD',
    location=(0.0, 0.0, 0.0), 
    scale=(1.0, 1.0, 1.0)
    )
    obj = bpy.context.active_object
    obj.name = 'sphere'
    obj.data.name = 'sphere'
    bpy.ops.object.transform_apply(location=True, rotation=False, scale=True)        

# Clear existing cones (optional, but useful if running the script multiple times)
for obj in bpy.data.objects:
    if obj.name.startswith('pmtdisc') or obj.name.startswith('neck') or obj.name.startswith('sphere'):
        bpy.data.objects.remove(obj)
    #if obj.name.startswith('neck'):
    #    bpy.data.objects.remove(obj)    

eventPosFit = []

positions = []
orientations = []
# Read positions from the CSV file and create cones
positions = []
csv_path = bpy.path.abspath(r"C:\deap\snoplus_positions.csv")  # assuming the CSV is in the same directory as the Blender file
with open(csv_path, 'r', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        x, y, z = map(float, row)
        positions.append((x, y, z))
        mag = sqrt(x*x+y*y+z*z)
        x0 = x/mag
        y0 = y/mag
        z0 = z/mag
        x = x0*rpsup; y = y0*rpsup; z = z0*rpsup;
        positions.append( (x, y, z) )
        orientations.append( (x0, y0, z0) )


# Create cones at the read positions
for i in range(7000):
    create_cone( positions[i], orientations[i] )

neck_boss_z = 6036.6
neck_location = (0,0, neck_boss_z/scale)

create_neck(neck_location)
create_sphere( (0.0,0.0,0.0) )
