import os
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.patches import Circle
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from mpl_toolkits.mplot3d import art3d

rav = 851
rpmt = 202./2

ffpmtpos = './Aksel_PMTpos_260.ratdb'
ff = open(ffpmtpos,'r')
pmtpos_offline = []
pmtposdirection = []
xpos = []
ypos = []
zpos = []

eventPosFit = []

if not os.path.isfile(ffpmtpos):
    print ('ERROR: need Aksel_PMTpos_260.ratdb to work')
    print ('INFO: inside TH2DPMTfancy.py, edit the variable pmtpos to the right place')
    sys.exit(1)

dead_pmt_id = [149,204]
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

def rotation_matrix(d):
    """
    Calculates a rotation matrix given a position vector d.
    Here I first put a disc on (0,0,851) and its norm is (0,0,1)
    Then rotate the disc around y-axis by zenith-angle theta, and
    then around z-axis by azimuth-angle phi. Euler rotation matrices are used.
    You may have a more efficient way to plot.
    """
    d_x = d[0]; d_y = d[1]; d_z = d[2]; 
    theta = arccos(d_z)
    phi = arctan2(d_y, d_x) 
    # rotate around y and then around z
    rotateY = np.array( [[cos(theta), 0, sin(theta)],[0, 1, 0],[-sin(theta),0,cos(theta)]], dtype = np.float64)
    rotateZ = np.array( [[cos(phi), -sin(phi), 0],[sin(phi), cos(phi), 0], [0, 0, 1]], dtype = np.float64)
    
    M = matmul(rotateZ, rotateY)
    return M

def pathpatch_2d_to_3d(pathpatch, pmtpos):#, pmt204):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    z = rav#pmtpos[2]
    vectorPos = np.array(pmtpos) 
    vectorPos /= np.linalg.norm(vectorPos) #Make sure the vector is normalised

    path = pathpatch.get_path() #Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path) #Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D #Change the class
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    

    verts = path.vertices #Get the vertices in 2D

    M = rotation_matrix(vectorPos) #Get the rotation matrix
    #print M
    pathpatch._segment3d = np.array([np.dot(M, (x, y, z)) for x, y in verts])

def pathpatch_translate(pathpatch, delta):
    """
    Translates the 3D pathpatch by the amount delta.
    """
    pathpatch._segment3d += delta

ax = plt.axes(projection = '3d') #Create axes
setframe = 900
ax.set_xlim(-setframe, setframe)
ax.set_ylim(-setframe, setframe)
ax.set_zlim(-setframe, setframe)

#p = Circle((0,0), rpmt, facecolor = 'g') #Add a circle in the xy plane
#ax.add_patch(p)
#pathpatch_2d_to_3d(p, z = rav)#, pmt204)
#pathpatch_translate(p, (0.5, 0.5, 0))

#test_offline = pmtpos_offline[0:250] #[ pmtpos_offline[30], pmtpos_offline[149], pmtpos_offline[204] ]

pmtid = [i for i in range(255)]
#pmtid = [30, 149, 204]
ids = [str(i) for i in pmtid]

for pmtid in range(255):
    pmtpos = pmtpos_offline[pmtid]
    pmtx = pmtpos[0]; pmty = pmtpos[1]; pmtz = pmtpos[2];
    p = Circle((0,0), rpmt, facecolor = 'b', alpha = .3)
    if pmtid  == 149 or pmtid == 204:
        #print ("warning: dead pmt", pmtid)
        p = Circle((0,0), rpmt, facecolor = 'r', alpha = .6)
        label = str(pmtid)
        ax.text( pmtx, pmty, pmtz, label )
    #p = Circle((0,0), .2, facecolor = 'y', alpha = .2)
    ax.add_patch(p)
    pathpatch_2d_to_3d(p, pmtpos)
    textDirection = pmtposdirection[pmtid]
    label = str(pmtid)
    #ax.text( pmtx, pmty, pmtz, label )
    #ax.text( pmtx, pmty, pmtz, label )#, textDirection)
    #distance = sqrt( sum(pmtpos*pmtpos) )
    #pathpatch_translate(p, distance)

plt.show()
