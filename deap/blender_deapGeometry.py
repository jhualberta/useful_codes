#!/usr/bin/python
# -*- coding: utf-8 -*-
# works in Blender 3.3.12
# Jie Hu jhu9@ualberta.ca, 26 Oct, 2023
import bpy
from bpy.ops import mesh
import math
from numpy import *
import numpy as np
import os, sys

rav = 851 # mm
rpmt = 202./2 # Hamamatsu R5912
rneck = 173
lightguide = 1426

## here the length unit is scaled to 1 mm/100, to save some "out of world" problem
scale = 100 

rav = rav/scale
rpmt = rpmt/scale
rneck = rneck/scale
lightguide = lightguide/scale

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
        depth=lightguide,
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
        depth=1000/scale, ## assume length as 10000 mm
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
    
    
def create_avsphere(location):    
    bpy.ops.mesh.primitive_uv_sphere_add(
    radius=850./scale,
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
    if obj.name.startswith('pmtdisc') or obj.name.startswith('neck'):
        bpy.data.objects.remove(obj)
    #if obj.name.startswith('neck'):
    #    bpy.data.objects.remove(obj)    


## find this file in rat/data/DEAP-3600
file_pmtpos = 'C:\deap\Aksel_PMTpos_260.ratdb'
ff = open(file_pmtpos,'r')
pmtpos_offline = []
pmtposdirection = []
xpos = []
ypos = []
zpos = []

eventPosFit = []

if not os.path.isfile(file_pmtpos):
    print ('ERROR: need Aksel_PMTpos_260.ratdb to work')
    print ('INFO: inside TH2DPMTfancy.py, edit the variable pmtpos to the right place')
    sys.exit(1)

dead_pmt_id = [149,204]
# Read positions from the CSV file and create cones
positions = []
orientations = []
for l in ff:
    if l.startswith('x'): xpos = [float(val) for val in l.split('[')[1][:-4].split(',')]
    if l.startswith('y'): ypos = [float(val) for val in l.split('[')[1][:-4].split(',')]
    if l.startswith('z'): zpos = [float(val) for val in l.split('[')[1][:-4].split(',')]

    for i,(x,y,z) in enumerate(zip(xpos, ypos, zpos)):
        v = array([x,y,z])
        v_direction = v/sqrt(sum(v*v))
        v_scale = v_direction*rav
        pmtpos_offline.append(v_scale)
        pmtposdirection.append(v_direction)
        positions.append( (v_scale[0], v_scale[1], v_scale[2]) )
        orientations.append( (v_direction[0], v_direction[1], v_direction[2]) )

# Create cones at the read positions
for i in range(255):
    create_cone( positions[i], orientations[i] )
    
neck_location = (0,0,1000/scale)
create_neck(neck_location)
create_avsphere( (0,0,0) )