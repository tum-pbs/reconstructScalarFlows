# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# Scene for generating a synthetic smoke plume simulation
#
# 1. adapt variable path 
# 2. example call "./manta ../scenes/simpleplume.py 0 100 1"
#

import os, sys, shutil, time, math, platform, datetime
import numpy as np
import _visualize as v
from manta import *
from enum import Enum

path = '/home/eckert/results/'
folderOut  = 'synthReal_%06d/'%int(sys.argv[1])
res = int(sys.argv[2])
version = int(sys.argv[3])

class drawKind(Enum):
	den3D = 1
	vel3D = 2
	den2D = 3
	vel2D = 4
	
# drawKind: 1: 3D den, 2: 3D vel, 3: 2D den, 4: 2D vel	
def saveVisGrid(grid, npy, output, dK, scale, negativeValues=False):
	if   dK==drawKind.den3D or dK==drawKind.den2D: copyGridToArrayReal(grid, npy)
	elif dK==drawKind.vel3D or dK==drawKind.vel2D: copyGridToArrayMAC(grid, npy)
	
	np.savez_compressed(path+folderOut+'tmp.npz', data=npy)
	if os.path.isfile(path+folderOut+output): os.remove(path+folderOut+output)
	os.rename(path+folderOut+'tmp.npz', path+folderOut+output)
	
	if   dK==drawKind.den3D: v.draw3DDensityGridNpy(path+folderOut+output, scale, True, negativeValues)
	elif dK==drawKind.vel3D: v.draw3DVelGridNpy(path+folderOut+output, scale, True)
	elif dK==drawKind.den2D: v.draw2DDensityNpy(path+folderOut+output, scale, True, negativeValues)
	elif dK==drawKind.vel2D: v.draw2DVelGridNpy(path+folderOut+output, scale, True)

def loadVel(grid, npy, filename):
	npy = np.load(filename)['data']
	copyArrayToGridMAC(npy, grid)

def loadDen(grid, npy, filename, scale=0, setBounds=False, bWidth=1):
	npy = np.load(filename)['data']
	copyArrayToGridReal(npy, grid)
	if scale!=0:  grid.multConst(scale)
	if setBounds: grid.setBound(0,bWidth)

if not os.path.exists(path+folderOut): os.makedirs(os.path.dirname(path+folderOut))

# solver params
factorY = 1.77
gs  = vec3(res,math.ceil(factorY*res),res) if math.ceil(factorY*res)%2==0 else vec3(res,math.ceil(factorY*res)+1,res)
s   = FluidSolver(name='main', gridSize = gs)

# prepare grids
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
density  = s.create(RealGrid)
pressure = s.create(RealGrid)

velNpy   = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 3], order='C')
denNpy   = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 1], order='C')

# noise field, tweak a bit for smoke source
noise = s.create(NoiseField, loadFromFile=False)
noise.posScale = vec3(120)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 3
noise.valOffset = 1.5#0.75
noise.timeAnim = 0.2

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))

flags.initDomain(boundaryWidth=1)
flags.fillGrid()
setOpenBound(flags, 1,'xXyYzZ',FlagOutflow|FlagEmpty) 
if version == 0 or version == 2:
	p0 = vec3(math.ceil(gs.x*0.44),0.,math.ceil(gs.z*0.38))
	p1 = vec3(math.ceil(gs.x*0.64),math.ceil(gs.x*0.068),math.ceil(gs.z*0.58))
elif version == 1:
	p0 = vec3(math.ceil(gs.x*0.42),0.,math.ceil(gs.z*0.4))
	p1 = vec3(math.ceil(gs.x*0.66),math.ceil(gs.x*0.1),math.ceil(gs.z*0.60))
srcInflow = ShapeDetails(s, 2, (p0+p1)/2., (p1-p0)/2., 0)
setInflowStructure(flags, srcInflow)
source = s.create(Box, p0=p0, p1=p1)

if (GUI):
	gui = Gui()
	#gui.show()

# restart simulation if already exists
startFrame=0
restart=False
for t in range(250, 0, -1):
	if os.path.isfile(path+folderOut+'velocity_%06d.npz'%t) and os.path.isfile(path+folderOut+'density_%06d.npz'%t):
		startFrame = t+1
		loadVel(vel, velNpy, path+folderOut+'velocity_%06d.npz'%t)
		loadDen(density, denNpy, path+folderOut+'density_%06d.npz' %t)
		restart=True 
		break
	
#main loop
for t in range(startFrame,250):
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=0.8)
	resetOutflow(flags=flags,real=density) 
	
	setWallBcs(flags=flags, vel=vel)   
	if version==2:
		addBuoyancy(density=density, vel=vel, gravity=vec3(0,-8e-5,0), flags=flags)
	else:	
		addBuoyancy(density=density, vel=vel, gravity=vec3(0,-4e-5,0), flags=flags)
	
	setVelInflow(vel, flags, vec3(0,0.006*res,0))
	
	solvePressure( flags=flags, vel=vel, pressure=pressure)
	
	densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=5, sigma=1.5)
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, strength=0.8)    
	density.setBound(0,1)
	
	saveVisGrid(density, denNpy, 'density_%06d.npz'%t,  drawKind.den3D, 2.5)
	saveVisGrid(vel,     velNpy, 'velocity_%06d.npz'%t, drawKind.vel3D, 0.2)
	
	mantaMsg('\nFrame %i %.1f' % (t, vel.getMaxAbs()))

	s.step()

