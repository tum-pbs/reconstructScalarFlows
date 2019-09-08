# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# render density to generate synthetic input images for reconstruction 
#
# 1. adapt variable path and pathCalib
# 2. determine calibration folder, e.g., calib20190813
# 3. example call "./manta ../scenes/reconstruct/renderVol.py calib20190813 0 100 1"
#

import os, sys, shutil, time, math, platform, datetime, glob
import numpy as np
import _visualize as v
from manta import *
from enum import Enum

class drawKind(Enum):
	den3D = 1
	vel3D = 2
	den2D = 3
	vel2D = 4
	
class reconKind(Enum):
	synth     = 1 # synthetic data + addSource,  real or orthographic cameras, no cutoff, true first density field
	synthReal = 2 # synthetic data + vel Inflow, real cameras, 				   cutoff
	real      = 3 # real data, real inflow       real cameras, 				   cutoff
	
# drawKind: 1: 3D den, 2: 3D vel, 3: 2D den, 4: 2D vel	
def saveVisGrid(grid, npy, output, dK, scaleVis, negativeValues=False):
	if   dK==drawKind.den3D or dK==drawKind.den2D: copyGridToArrayReal(grid, npy)
	elif dK==drawKind.vel3D or dK==drawKind.vel2D: copyGridToArrayMAC(grid, npy)
	
	np.savez_compressed(path+folderOut+'tmp.npz', data=npy)
	if os.path.isfile(path+folderOut+output): os.remove(path+folderOut+output)
	os.rename(path+folderOut+'tmp.npz', path+folderOut+output)
	
	if   dK==drawKind.den3D: v.draw3DDensityGridNpy(path+folderOut+output, scaleVis, rK == reconKind.synthReal, negativeValues)
	elif dK==drawKind.vel3D: v.draw3DVelGridNpy(path+folderOut+output, scaleVis, rK == reconKind.synthReal)
	elif dK==drawKind.den2D: v.draw2DDensityNpy(path+folderOut+output, scaleVis, rK == reconKind.synthReal, negativeValues)
	elif dK==drawKind.vel2D: v.draw2DVelGridNpy(path+folderOut+output, scaleVis, rK == reconKind.synthReal)

def loadVel(grid, npy, filename):
	npy = np.load(filename)['data']
	copyArrayToGridMAC(npy, grid)

def loadDen(grid, npy, filename, scale=0, setBounds=False, bWidth=1):
	npy = np.load(filename)['data']
	copyArrayToGridReal(npy, grid)
	if scale!=0:  grid.multConst(scale)
	if setBounds: grid.setBound(0,bWidth)

def loadImg(grid, npy, filename, scale=0, interpol=False, imgs=0):
	npy = np.load(filename)['data']
	copyArrayToGridReal(npy, grid)
	if interpol:   interpolateImgs(imgs,grid,scale)
	elif scale!=0: grid.multConst(scale)


###### read input arguments ######
if len(sys.argv)<4:
	print("three arguments required: calibFolder, captureFolder, and res")
	sys.exit()
	
path          = '/home/eckert/results/'
pathCalib     =  path+'/%s/%s_rays.txt'
folderOut     = 'synthReal_%06d/'%int(sys.argv[2])
calibFolder   = sys.argv[1]
res           = int(sys.argv[3])
version       = int(sys.argv[4])

pathCalib = pathCalib%(calibFolder,'%i')
rK = reconKind.synth

###### create grids ######
factorY = 1.77
gs = vec3(res,math.ceil(factorY*res),res) if math.ceil(factorY*res)%2==0 else vec3(res,math.ceil(factorY*res)+1,res)
s          = Solver(name='volume', gridSize = gs, dim=3) 
den        = s.create(RealGrid) 
flags      = s.create(FlagGrid)
flags.initDomain(boundaryWidth=1) 
flags.fillGrid()
setOpenBound(flags, 1, 'xXyYzZ', FlagOutflow|FlagEmpty) 

denNpy   = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 1], order='C')

###### setup inflow region ######
if version == 0 or version == 2:
	p0 = vec3(math.ceil(gs.x*0.44),0.,math.ceil(gs.z*0.38))
	p1 = vec3(math.ceil(gs.x*0.64),math.ceil(gs.x*0.068),math.ceil(gs.z*0.58))
elif version == 1:
	p0 = vec3(math.ceil(gs.x*0.42),0.,math.ceil(gs.z*0.4))
	p1 = vec3(math.ceil(gs.x*0.66),math.ceil(gs.x*0.1),math.ceil(gs.z*0.60))
	
srcInflow = ShapeDetails(s, 2, (p0+p1)/2., (p1-p0)/2., 0)
setInflowStructure(flags, srcInflow)#, 10, False, src0)
###### end setup inflow region ######

###### tomography part ######
## image parameters ##
angles     = [0,1,2,3,4] 
width      = 600
height     = math.ceil(width*1.77) if math.ceil(width*1.77)%2==0 else math.ceil(width*1.77)+1
stepSize   = 0.7
switchXY   = True#False

gsImgs  = vec3(width,height,len(angles))
sImgs   = Solver(name='images', gridSize = gsImgs, dim=3) if len(angles)>1 else Solver(name='main', gridSize = gsImgs, dim=2)
imgs    = sImgs.create(RealGrid)
imgsNpy = np.empty(shape=[len(angles), height, width, 1], order='C')


## create image class ##
orthographic = False and rK == reconKind.synth
if orthographic: angles = [0,1,2]
i = Image(sImgs,width,height,len(angles),pathCalib,switchXY,orthographic)
###### end tomography part ######

print(path + folderOut)
for file in glob.glob('density_*.npz'):
	t = int(file[-7:-4])
	if not os.path.isfile(path+folderOut+'imgs_%06d.npz'%t):
		denNpy = np.load(file)['data']
		copyArrayToGridReal(denNpy, den)
		i.render(imgs, den, stepSize, '', '', False, flags)
		saveVisGrid(imgs,      imgsNpy, 'imgs_%06d.npz'%t,drawKind.den2D, 3.0)
		