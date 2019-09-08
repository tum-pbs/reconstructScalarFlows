# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# Reconstruction of both density and velocity volume based on input images
#
# 0. make sure not to use more than ~4 threads (export OMP_NUM_THREADS=4)
# 1. adapt variable path, pathCalib, and captureFolderPath
# 2. example call "./manta ../scenes/reconstruct/reconDenVel.py calib20190813 0813_80_0085 100 3 8 5e-2 5e-4 5e-1 5e-2 1e-4 1e-3 0.8"
#	
import os, sys, shutil, time, math, platform, datetime
import numpy as np
import _visualize as v
import _writeJson as wJ
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
if len(sys.argv)<5:
	print("four arguments required: calibFolder, captureFolder, res, and reconKind")
	sys.exit()
	
	
path      = '/home/eckert/results/' 
calibFolder   = sys.argv[1]
captureFolder = sys.argv[2]
res           = int(sys.argv[3])# make sure this will result in an integer in the y-domain size

pathCalib = path+'/%s/%s_rays.txt'%(calibFolder,'%i')
captureFolderPath = 'input/%s/postprocessed/'%captureFolder


###### set parameters ######
if len(sys.argv)>5: scale = float(sys.argv[5])
else: 
	if   rK == reconKind.real:      scale = 6
	elif rK == reconKind.synthReal: scale = 1#2.5
	else:                           scale = 7.5

if len(sys.argv)>6: smoothDen = float(sys.argv[6])
else:
	if rK == reconKind.synthReal and int(captureFolder) == 3: smoothDen = 1e-2
	elif rK == reconKind.synthReal: smoothDen = 1e-4
	else: smoothDen = 2.5e-2
	
if len(sys.argv)>7: kinDen = float(sys.argv[7])
else:
	if rK == reconKind.synthReal and int(captureFolder) == 3: kinDen = 1e-3
	elif rK == reconKind.synthReal: kinDen = 1e-6
	else: kinDen = 5e-4
	
if len(sys.argv)>8:  smoothVel  = float(sys.argv[8])
else:
	smoothVel = 6e-1
	if rK == reconKind.synthReal:
		if   int(captureFolder)==2: smoothVel = 5e-1
		elif int(captureFolder)==3: smoothVel = 5e-1
		elif int(captureFolder)==4: smoothVel = 5e-1
		elif int(captureFolder)==6: smoothVel = 5e-1

if len(sys.argv)>9:  kinVel     = float(sys.argv[9])
else:
	kinVel = 5e-2
	if rK == reconKind.synthReal:
		if   int(captureFolder)==2: kinVel = 4e-2
		elif int(captureFolder)==3: kinVel = 5e-2
		elif int(captureFolder)==4: kinVel = 4e-2
		elif int(captureFolder)==6: kinVel = 4e-2

if len(sys.argv)>10:  smoothInfl = float(sys.argv[10])
else: smoothInfl = 1e-2

if len(sys.argv)>11: kinInfl    = float(sys.argv[11])
else: kinInfl = 1e-2

if len(sys.argv)>12: velInflowValue = float(sys.argv[12])
else:
	velInflowValue = res/73. if '_90_' in captureFolder else res/80.	
	if rK == reconKind.synthReal:
		velInflowValue = 0.006*res
		if   int(captureFolder)==2: velInflowValue = 0.95
		elif int(captureFolder)==3: velInflowValue = 0.98
		elif int(captureFolder)==4: velInflowValue = 0.9
		elif int(captureFolder)==6: velInflowValue = 1.05


###### setup reconstruction parameters ######
rK           = reconKind(int(sys.argv[4]))#reconKind.synthReal
startFrame   = 0 if rK == reconKind.real else 14
lastFrame    = startFrame+151
step         = 1
restartRecon = True
resFactor    = res/50.
ms = res/2 if rK != reconKind.synthReal else res

factorY = 1.77 if rK != reconKind.synth else 1.33
boxTest   = False and rK == reconKind.synth
###### end setup reconstruction parameters ######

###### create grids ######
gs = vec3(res,math.ceil(factorY*res),res) if math.ceil(factorY*res)%2==0 else vec3(res,math.ceil(factorY*res)+1,res)
s          = Solver(name='volume', gridSize = gs, dim=3) 
sO         = Solver(name='volumeHigh', gridSize = vec3(200,math.ceil(factorY*200),200), dim=3) 
vel        = s.create(MACGrid) 
velUpdate  = s.create(MACGrid) 
denPredict = s.create(RealGrid) 
den        = s.create(RealGrid) 
denTarget  = s.create(RealGrid) 
denHelp    = s.create(RealGrid) 
pressure   = s.create(RealGrid) 
src0       = s.create(RealGrid) if rK != reconKind.real else 0
flags      = s.create(FlagGrid)
den1O      = sO.create(RealGrid) if rK != reconKind.real else 0
flagsO     = sO.create(FlagGrid) if rK != reconKind.real else 0
bWidth=1 
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()
setOpenBound(flags, 1, 'xXyYzZ', FlagOutflow|FlagEmpty) 
###### end create grids ######

###### create numpy arrays ######
velNpy   = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 3], order='C')
denNpy   = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 1], order='C')
den1ONpy = np.empty(shape=[200, math.ceil(factorY*200), 200, 1], order='C') if rK != reconKind.real else 0
###### end create numpy arrays ######

###### setup path and filenames ######
folderO        = path+'synthReal/synthReal_%06d_7/'%int(captureFolder) if rK == reconKind.synth or rK == reconKind.synthReal else path+captureFolderPath
folderIn       = folderO + '%d/' % res if (res<200 and rK != reconKind.real) else folderO
prefix         = 'rDV'
suffix         = '%d_%s_%.1f_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e_%.1f' % (res, rK.name, scale, smoothDen, kinDen, smoothVel, kinVel, smoothInfl, kinInfl, velInflowValue)
folderOut      = prefix+'_%06d_%s/'%(int(captureFolder),suffix) if rK != reconKind.real else prefix+'_%s_%s/'%(captureFolder, suffix)
densityName    = 'density_%06d.npz'  
velocityName   = 'velocity_%06d.npz' 

if not os.path.exists(path+folderOut): os.makedirs(os.path.dirname(path+folderOut))
try:
	scriptname = 'reconDenVel.py'
	shutil.copy(os.path.abspath(os.path.dirname(sys.argv[0]))+'/'+scriptname, path+folderOut+scriptname)
except OSError as ecx:
	raise;
###### end setup path and filenames ######

###### setup inflow region ######
if rK == reconKind.real or rK == reconKind.synthReal:
	p0 = vec3(math.ceil(gs.x*0.44),0.,math.ceil(gs.z*0.38))
	p1 = vec3(math.ceil(gs.x*0.64),math.ceil(gs.x*0.068),math.ceil(gs.z*0.58))
	if rK == reconKind.synthReal and int(captureFolder) == 4:
		p0 = vec3(math.ceil(gs.x*0.42),0.,math.ceil(gs.z*0.4))
		p1 = vec3(math.ceil(gs.x*0.66),math.ceil(gs.x*0.1),math.ceil(gs.z*0.60))
else: 
	p0 = vec3(0,0,0)
	p1 = vec3(0,0,0)
	srcInflow = 0
	
if rK == reconKind.synthReal or rK == reconKind.real: 
	srcInflow = ShapeDetails(s, 2, (p0+p1)/2., (p1-p0)/2., 0)
	setInflowStructure(flags, srcInflow)#, 10, False, src0)
	#saveVisGrid(src0, denNpy, 'finalSrc_%06d.npz'%0,   drawKind.den3D, 0.5*scale)
	#sys.exit()
	if rK != reconKind.real: upsampleFlagGrid(flagsO, flags)
###### end setup inflow region ######

###### tomography part ######
## image parameters ##
angles     = [0,1,2,3,4] 
width      = res*6
height     = math.ceil(width*1.77) if math.ceil(width*1.77)%2==0 else math.ceil(width*1.77)+1
stepSize   = 0.7
minNumCams = 2
switchXY   = True#False
useDenTarget = False # use same denTarget instead of denP+denU where denU is recalculated in each bivariate OF PD iter
	
## tomography parameters ##
pdT = PDParams(s, 0.01, 100, 1, 10) #PDParams(parent,sigma, tau, theta, mxIter)
if rK == reconKind.real: wT = RegWeightsTomo(s, smoothDen*resFactor, kinDen*resFactor, smoothInfl*resFactor, kinInfl*resFactor) #RegWeightsTomo(parent,smooth,kinetic)
else:                    wT = RegWeightsTomo(s, smoothDen*resFactor, kinDen*resFactor, smoothInfl*resFactor, kinInfl*resFactor) #RegWeightsTomo(parent,smooth,kinetic)
#TomoParams(parent,t,path,threshVH,threshMask,stepSize,minNumCams,numPixels,numVoxels,angleWeight,pdParams,regWeights,shapeLimit);
tomoParams = TomoParams(s, startFrame, path+folderOut, 1e-9, 1e-4, stepSize, minNumCams, 0, 0, 1, pdT, wT, 0) 

pdT_trgt = PDParams(s, 0.01, 100, 1, 20) #PDParams(parent,sigma, tau, theta, mxIter)
if rK == reconKind.real: wT_trgt = RegWeightsTomo(s, 1e-2*smoothDen*resFactor, 1e-2*kinDen*resFactor, 0, 0) #RegWeightsTomo(parent,smooth,kinetic)
else:                    wT_trgt = RegWeightsTomo(s, 1e-1*smoothDen*resFactor, kinDen*resFactor, 0, 0) #RegWeightsTomo(parent,smooth,kinetic)
#TomoParams(parent,t,path,threshVH,threshMask,stepSize,minNumCams,numPixels,numVoxels,angleWeight,pdParams,regWeights,shapeLimit);
tomoParams_trgt = TomoParams(s, startFrame, path+folderOut, 1e-9, 1e-4, stepSize, minNumCams, 0, 0, 1, pdT_trgt, wT_trgt, 0) 

# source tomography parameters
if rK != reconKind.synth:
	# tomography for src
	if rK == reconKind.real: wT_firstDen = RegWeightsTomo(s, 1e-2*smoothDen*resFactor, 1e-2*kinDen*resFactor, 0, 0) #RegWeightsTomo(parent,smooth,kinetic)
	else:                    wT_firstDen = RegWeightsTomo(s, 0, 1e-6, 0, 0) #RegWeightsTomo(parent,smooth,kinetic)
	src_firstDen = ShapeDetails(s, 2, gs*vec3(0.5,0.05,0.5), gs*vec3(0.3, 0.1, 0.3), 3) #ShapeDetails(parent,shape,center,vec,radius)
	#TomoParams(parent,t,path,threshVH,threshMask,stepSize,numPixels,numVoxels,angleWeight,pdParams,regWeights,shapeLimit);
	tomoParams_firstDen = TomoParams(s, startFrame, path+folderOut, 1e-9, 1e-4, stepSize, minNumCams, 0, 0, 1, pdT, wT_firstDen, src_firstDen) 

## image FluidSolver to read in captured images ##
if rK == reconKind.real: 
	sImgsO   = Solver(name='imagesO', gridSize = vec3(1080, 1920, 5), dim=3) 
	imgsO    = sImgsO.create(RealGrid)
	imgsONpy = np.empty(shape=[5, 1920, 1080, 1], order='C')

# image FluidSolver for target images
gsImgs  = vec3(width,height,len(angles))
sImgs   = Solver(name='images', gridSize = gsImgs, dim=3) if len(angles)>1 else Solver(name='main', gridSize = gsImgs, dim=2)
imgs    = sImgs.create(RealGrid)
imgsNpy = np.empty(shape=[len(angles), height, width, 1], order='C')

## create image class ##
orthographic = False and rK == reconKind.synth
if orthographic: angles = [0,1,2]
i = Image(sImgs,width,height,len(angles),pathCalib,switchXY,orthographic)
###### end tomography part ######

###### optical flow part ######
pd = PDParams(s, 0.01, 10, 1, 10) #PDParams(parent, sigma, tau, theta, mxIter)
if rK == reconKind.real: w = RegWeightsOF(s, smoothVel*resFactor, kinVel*resFactor) #RegWeightsOF(parent,smooth,kinetic,kineticZ=-1,adaptiveKineticZ=false,sumZero=0)
else:                    w = RegWeightsOF(s, smoothVel*resFactor, kinVel*resFactor) #RegWeightsOF(parent,smooth,kinetic,kineticZ=-1,adaptiveKineticZ=false,sumZero=0)
#OFParams(parent,t,minSize,bivariate,path,s,strides,N,dim,heightOfSrc,pdParams,regWeights)
ofParams = OFParams(s, startFrame, ms, True, path+folderOut, gs, vec3(1, gs.x, gs.x*gs.y), gs.x*gs.y*gs.z, 3, srcInflow.getHeightOfSrc(), useDenTarget, velInflowValue, pd, w)
###### end optical flow part ######

###### viscosity approximation ######
visc     = 0.0000148 # air at 15dC = 1.48 * 10^-5   ##0.000001 # ca. "water" ##0.0 # off
timestep = 1./30.
alpha    = visc * timestep * float(res*res)
###### end viscosity approximation ######

###### prepare initial condition and setup target images ######
restartedRecon = False
if boxTest:
	source0 = s.create(Box, center=gs*vec3(0.5,0.5,0.5), size=gs*vec3(0.1))
	source0.applyToGrid(grid=src0, value=0.7)
	source1 = s.create(Box, center=gs*vec3(0.5,0.51,0.5), size=gs*vec3(0.1))
	source1.applyToGrid(grid=den, value=0.7)
else:
	# get current state (den and vel)
	if restartRecon:
		# find last successful reconstruction step, set this to startFrame
		for t in range(lastFrame-step, startFrame+step, -step):
			if os.path.isfile(path+folderOut+'velocity_%06d.npz'%t) and os.path.isfile(path+folderOut+'density_%06d.npz'%t):
				startFrame = t
				restartedRecon = True
				break	
		if restartedRecon:
			print('restartRecon %d'%(startFrame))
			sys.stdout.flush()
			loadVel(vel, velNpy, path+folderOut+'velocity_%06d.npz'%startFrame)
			loadDen(den, denNpy, path+folderOut+'density_%06d.npz' %startFrame)
	# or get first den field
	if not restartedRecon:
		if rK == reconKind.real or rK == reconKind.synthReal: 
			# target images for first den estimation
			if   rK == reconKind.real:      loadImg(imgsO, imgsONpy, folderO+'imgs_%06d.npz'%startFrame, scale, True, imgs)
			elif rK == reconKind.synthReal: loadImg(imgs,  imgsNpy,  folderO+'imgs_%06d.npz'%startFrame, scale)
				
			# first den estimation
			saveVisGrid(imgs, imgsNpy, 'imgsTarget_%06d.npz'%startFrame,   drawKind.den2D, 3.0)
			reconstructDensity(den, i, imgs, flags, tomoParams_firstDen, 0) # mask unnused
			
			# save estimated first den
			i.render(imgs, den, stepSize, '', '', True, flags)
			saveVisGrid(imgs, imgsNpy, 'imgsDiff_%06d.npz'%startFrame,     drawKind.den2D, 3.0, True)
			i.render(imgs, den, stepSize, '', '', False, flags)
			saveVisGrid(imgs, imgsNpy, 'imgsRendered_%06d.npz'%startFrame, drawKind.den2D, 3.0)
			saveVisGrid(den,  denNpy,  'density_%06d.npz'%startFrame,      drawKind.den3D, 3.0)
		elif rK == reconKind.synth: 
			loadDen(den, denNpy, folderIn+densityName%startFrame, scale, True, bWidth)
		
	# load target images (imgs)
	if rK == reconKind.real:        loadImg(imgsO, imgsONpy, folderO+'imgs_%06d.npz'%(startFrame+step), scale, True, imgs)
	elif rK == reconKind.synthReal: loadImg(imgs,  imgsNpy,  folderO+'imgs_%06d.npz'%(startFrame+step), scale)
	elif rK == reconKind.synth: 
		loadDen(den1O, den1ONpy, folderO+'densityH_%06d.npz'%(startFrame+step), scale, True, bWidth)
		i.render(imgs, den1O, stepSize, '', '', False, flagsO)	
	
	# setup source 
	if rK == reconKind.synth: loadDen(src0, denNpy, folderIn+densityName%0, scale)
###### end prepare initial condition and setup target images ######

wJ.writeJasonFile(path+folderOut,calibFolder,captureFolder,folderOut,rK,res,factorY,p0,p1,scale,tomoParams,tomoParams_firstDen,tomoParams_trgt,ofParams,angles,width,height,stepSize,minNumCams,orthographic,restartedRecon,startFrame)

for t in range(startFrame+step, lastFrame, step):
	start = time.time()
	ofParams.setT(t)
	tomoParams.setT(t)
	
	# apply inflow
	if rK == reconKind.synth: addSrc(den, src0)
	
	blurRealGrid(den, denHelp, 2*resFactor)
	reconstructDensity(denTarget, i, imgs, flags, tomoParams_trgt, denHelp) 
	predict(vel, denPredict, den, pressure, flags, tomoParams, ofParams, denTarget, 1e-3, alpha)
	
	# calculate update
	iter = reconstructDenVelMS(velUpdate, vel, i, imgs, denPredict, den, flags, ofParams, tomoParams, pressure, denTarget, denHelp) 
	
	# align predicted with update quantities
	align(vel, den, pressure, velUpdate, flags, tomoParams, ofParams, denTarget) # takes 4.6s for res=210
	
	################################################# takes 20s for res=210, depending on how many fields visualized
	# save and visualize grids
	saveVisGrid(den,       denNpy, 'density_%06d.npz'%t,     drawKind.den3D, 3.0)
	saveVisGrid(vel,       velNpy, 'velocity_%06d.npz'%t,    drawKind.vel3D, 0.2)
	saveVisGrid(imgs,      imgsNpy,'imgsTarget_%06d.npz'%t,  drawKind.den2D, 3.0)
	i.render(imgs, den, stepSize, '', '', False, flags)
	saveVisGrid(imgs,      imgsNpy, 'imgsRendered_%06d.npz'%t,drawKind.den2D, 3.0)
	#################################################

	# prepare for next time step
	if   rK == reconKind.real:      loadImg(imgsO, imgsONpy, folderO+'imgs_%06d.npz'%(t+step), scale, True, imgs)
	elif rK == reconKind.synthReal: loadImg(imgs,  imgsNpy,  folderO+'imgs_%06d.npz'%(t+step), scale)
	elif rK == reconKind.synth:
		if boxTest: break
		loadDen(den1O, den1ONpy, folderO+'densityH_%06d.npz'%(t+step), scale, True, bWidth)
		i.render(imgs, den1O, stepSize, '', '', False, flagsO)

	print('t=%03d, whole took %.4f s' %(t, time.time() - start))
	sys.stdout.flush()
	