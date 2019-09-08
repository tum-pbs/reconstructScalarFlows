# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# delete inflow area (lower 15 cells) from density, images, and velocity - save as npz and visualize as jpg
# calculate statistics (min, max, mean, std. deviation)
# calculate psnr of image differences
# 1. adapt variable path and pathCalib
# 2. if global statistics are desired and global mean values are available, adapt meanDenGlobal, meanVelGlobal, and meanImgsGlobal
# 3. example call "./manta ../scenes/reconstruct/postprocessRecons.py 0813_80_0085 7 0 1 1 1 1 1" for full postprocessing 
#              or "./manta ../scenes/reconstruct/postprocessRecons.py 0813_80_0085 7 0 0 0 0 1 1" for retrieving data statistics only
# global min, max, mean must be retrieved manually from the csv files
# global standard deviation must be calculated with math.sqrt(devDenGlobal/(cntDenGlobal)) in the end
#
import csv
import os, sys, shutil, time, math, glob
import numpy as np
import _visualize as v
from manta import *
from enum import Enum

class drawKind(Enum):
	den3D = 1
	vel3D = 2
	den2D = 3
	vel2D = 4

def saveVisGrid(grid, npy, output, dK, scale, negativeValues=False):
	if   dK==drawKind.den3D or dK==drawKind.den2D: copyGridToArrayReal(grid, npy)
	elif dK==drawKind.vel3D or dK==drawKind.vel2D: copyGridToArrayMAC(grid, npy)
	
	np.savez_compressed(path+folder+'tmp.npz', data=npy)
	if os.path.isfile(path+folder+output): os.remove(path+folder+output)
	os.rename(path+folder+'tmp.npz', path+folder+output)
	
	if   dK==drawKind.den3D: v.draw3DDensityGridNpy(path+folder+output, scale, negativeValues)
	elif dK==drawKind.vel3D: v.draw3DVelGridNpy(path+folder+output, scale)
	elif dK==drawKind.den2D: v.draw2DDensityNpy(path+folder+output, scale, negativeValues)
	elif dK==drawKind.vel2D: v.draw2DVelGridNpy(path+folder+output, scale)


path      = '/home/eckert/results/' 
meanDenGlobal = 0
meanVelGlobal = 0
meanImgsGlobal = 0

pathCalib = path+'/calib20190813/%s_rays.txt'%('%i')
folder = path + 'rDV_%s_100_real_8.0_5.0e-02_5.0e-04_5.0e-01_5.0e-02_1.0e-04_1.0e-03_0.8/'%sys.argv[1]

plus = int(sys.argv[2]) # number of cells to cut off above the original inflow area

overwrite = int(sys.argv[3]) # overwrite existing files? 
writeNpz = bool(int(sys.argv[4])) # write npz files to disc 
visualize = bool(int(sys.argv[5])) # visualize density and images 
visualizeVel = bool(int(sys.argv[6])) # visualize velocity
calcPSNR = bool(int(sys.argv[7])) # calculate psnr on images 
calcStats = bool(int(sys.argv[8])) # retrieve statistics for density, images, and velocity

os.chdir(folder)

res = 100
factorY = 1.77
gs = vec3(res,math.ceil(factorY*res),res) if math.ceil(factorY*res)%2==0 else vec3(res,math.ceil(factorY*res)+1,res)
s  = Solver(name='volume', gridSize = gs, dim=3) 
den   = s.create(RealGrid) 
vel   = s.create(MACGrid) 
denNpy  = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 1], order='C')
velNpy   = np.empty(shape=[int(gs.z), int(gs.y), int(gs.x), 3], order='C')

angles  = [0,1,2,3,4] 
width   = res*6
height  = math.ceil(width*1.77) if math.ceil(width*1.77)%2==0 else math.ceil(width*1.77)+1
gsImgs  = vec3(width,height,len(angles))
sImgs   = Solver(name='images', gridSize = gsImgs, dim=3) if len(angles)>1 else Solver(name='main', gridSize = gsImgs, dim=2)
imgs0   = sImgs.create(RealGrid)
imgs1   = sImgs.create(RealGrid)
imgs2   = sImgs.create(RealGrid)
imgsNpy = np.empty(shape=[len(angles), height, width, 1], order='C')

switchXY = True#False
i = Image(sImgs,width,height,len(angles),pathCalib,switchXY,False)

sImgsO   = Solver(name='imagesO', gridSize = vec3(width, height, 5), dim=3) 
imgsO    = sImgsO.create(RealGrid)
imgsONpy = np.empty(shape=[5, height, width, 1], order='C') 

p0 = vec3(math.ceil(gs.x*0.44),0.,math.ceil(gs.z*0.38))
p1 = vec3(math.ceil(gs.x*0.64),math.ceil(gs.x*0.068),math.ceil(gs.z*0.58))

velStats = [0,0,0,0,0,0]
imgsStats = [0,0,0,0,0,0]
denStats = [0,0,0,0,0,0]

untilT = 150
N = gs.x*gs.y*gs.z

devDenGlobal = 0
devVelGlobal = 0
devImgsGlobal = 0
cntDenGlobal = 0
cntVelGlobal = 0
cntImgsGlobal = 0

for file in glob.glob('*.npz'):
	t = int(file[-7:-4])
	if t>untilT: continue
	if 'velocity_0' in file or 'velocity_noInflow_0' in file:
		velFileName = folder+'velocity_noInflow_%06d.npz'%(t)
		velLoaded = False
		if writeNpz and (overwrite or not os.path.isfile(velFileName)):
			velNpy = np.load(file)['data']
			copyArrayToGridMAC(velNpy, vel)
			deleteInflowVel(vel, p0, p1, plus)
			copyGridToArrayMAC(vel, velNpy)
			np.savez_compressed(velFileName, data=velNpy)
			velLoaded = True
		if visualizeVel and (overwrite or not os.path.isfile(velFileName[:-3]+'jpg')): v.draw3DVelGridNpy(velFileName, 0.2, False)
		if calcStats:
			if not velLoaded:
				velNpy = np.load(velFileName)['data']
				copyArrayToGridMAC(velNpy, vel)
			max = vel.getMax()
			min = vel.getMin()
			nCurrent = vel.getCellCount()
			meanUnscaled = vel.getMeanUnscaled()
			if velStats[0]< max: velStats[0]=max
			if velStats[1]> min: velStats[1]=min
			velStats[2]=velStats[2]+meanUnscaled
			velStats[3]=velStats[3]+nCurrent
			velStats[4]=velStats[4]+N
			if t==untilT:
				meanSum = velStats[2]/velStats[4]
				stdDevSum = 0
				stdDevSumGlobal = 0
				for t2 in range(1,untilT+1):
					velNpy = np.load(folder+'velocity_noInflow_%06d.npz'%t2)['data']
					copyArrayToGridMAC(velNpy, vel)
					stdDevSum = stdDevSum + vel.getStdDeviationUnscaled(meanSum)
					if meanVelGlobal: 
						devVelGlobal = devVelGlobal + vel.getStdDeviationUnscaled(meanVelGlobal)
						cntVelGlobal = cntVelGlobal + N
				stdDevSum = math.sqrt(stdDevSum/(velStats[4]))
				with open(path +'velStatsSingle.csv', 'a', newline='') as csvFile:
					writer = csv.writer(csvFile)
					writer.writerow([sys.argv[1], velStats[0], velStats[1], meanSum, stdDevSum, velStats[3]/(untilT*100*178)])
				csvFile.close()
	elif 'density_0' in file or 'density_noInflow_0' in file:
		denFileName = folder+'density_noInflow_%06d.npz'%(t)
		denLoaded = False
		if writeNpz and (overwrite or not os.path.isfile(denFileName)):
			denNpy = np.load(file)['data']
			copyArrayToGridReal(denNpy, den)
			deleteInflowDen(den, p0, p1, plus)
			copyGridToArrayReal(den, denNpy)
			np.savez_compressed(denFileName, data=denNpy)
			denLoaded = True
		if visualize and (overwrite or not os.path.isfile(denFileName[:-3]+'jpg')): v.draw3DDensityGridNpy(denFileName, 3.0, False, False)
		
		# handle images! also from input
		imgFileName = folder+'imgsRendered_noInflow_%06d.npz'%(t)
		imgInputFileName = folder+'imgsTarget_noInflow_%06d.npz'%(t)
		imgs1Loaded = False 
		imgs2Loaded = False
		if writeNpz and (overwrite or not os.path.isfile(imgFileName) or not os.path.isfile(imgInputFileName)):
			if overwrite or not os.path.isfile(imgFileName):
				denNpy = np.load(denFileName)['data']
				copyArrayToGridReal(denNpy, den)
				i.render(imgs1, den, 0.7, '', '', False, 0)
				copyGridToArrayReal(imgs1, imgsNpy)
				np.savez_compressed(imgFileName, data=imgsNpy)
			else:
				imgsNpy = np.load(imgFileName)['data']
				copyArrayToGridReal(imgsNpy, imgs1)
			imgs1Loaded = True
		
			if os.path.isfile(folder+'imgsRendered_%06d.npz'%t) and (overwrite or not os.path.isfile(imgInputFileName)):
				imgsNpy = np.load(folder+'imgsRendered_%06d.npz'%t)['data']
				copyArrayToGridReal(imgsNpy, imgs0)
				imgsONpy = np.load(folder+'imgsTarget_%06d.npz'%(t))['data']
				copyArrayToGridReal(imgsONpy, imgs2) 
					
				deleteInflowImg(imgs0, imgs1, imgs2)			
				copyGridToArrayReal(imgs2, imgsNpy)
				np.savez_compressed(imgInputFileName, data=imgsNpy)
				imgs2Loaded = True
		if visualize and (overwrite or not os.path.isfile(imgFileName[:-3]+'jpg')): v.draw2DDensityNpy(imgFileName, 3.0, False, False)
		if visualize and os.path.isfile(imgInputFileName) and (overwrite or not os.path.isfile(imgInputFileName[:-3]+'png')): v.draw2DDensityNpy(imgInputFileName, 3.0, False, False)
		if calcStats:
			if not imgs1Loaded:
				imgsNpy = np.load(imgFileName)['data']
				copyArrayToGridReal(imgsNpy, imgs1)
			max = imgs1.getMax()
			min = imgs1.getMin()
			nCurrent = imgs1.getCellCount()
			meanUnscaled = imgs1.getMeanUnscaled()
			if imgsStats[0]< max: imgsStats[0]=max
			if imgsStats[1]> min: imgsStats[1]=min
			imgsStats[2]=imgsStats[2]+meanUnscaled
			imgsStats[3]=imgsStats[3]+nCurrent
			imgsStats[4]=imgsStats[4]+N
			if not denLoaded:
				denNpy = np.load(denFileName)['data']
				copyArrayToGridReal(denNpy, den)
			max = den.getMax()
			min = den.getMin()
			nCurrent = den.getCellCount()
			meanUnscaled = den.getMeanUnscaled()
			if denStats[0]< max: denStats[0]=max
			if denStats[1]> min: denStats[1]=min
			denStats[2]=denStats[2]+meanUnscaled
			denStats[3]=denStats[3]+nCurrent
			denStats[4]=denStats[4]+N
			if t==untilT:
				meanSum = imgsStats[2]/imgsStats[4]
				stdDevSum = 0
				for t2 in range(1,untilT+1):
					imgsNpy = np.load(folder+'imgsRendered_noInflow_%06d.npz'%t2)['data']
					copyArrayToGridReal(imgsNpy, imgs1)
					stdDevSum = stdDevSum + imgs1.getStdDeviationUnscaled(meanSum)
					if meanImgsGlobal: 
						devImgsGlobal = devImgsGlobal + imgs1.getStdDeviationUnscaled(meanImgsGlobal)
						cntImgsGlobal = cntImgsGlobal + N
				stdDevSum = math.sqrt(stdDevSum/(imgsStats[4]))
				with open(path + 'imgsStatsSingle.csv', 'a', newline='') as csvFile:
					writer = csv.writer(csvFile)
					writer.writerow([sys.argv[1], imgsStats[0], imgsStats[1], meanSum, stdDevSum, imgsStats[3]/(untilT*100*178)])
				csvFile.close()
				
				meanSum = denStats[2]/denStats[4]
				stdDevSum = 0
				for t2 in range(1,untilT+1):
					denNpy = np.load(folder+'density_noInflow_%06d.npz'%t2)['data']
					copyArrayToGridReal(denNpy, den)
					stdDevSum = stdDevSum + den.getStdDeviationUnscaled(meanSum)
					if meanDenGlobal: 
						devDenGlobal = devDenGlobal + den.getStdDeviationUnscaled(meanDenGlobal)
						cntDenGlobal = cntDenGlobal + N
				stdDevSum = math.sqrt(stdDevSum/(denStats[4]))				
				with open(path +'denStatsSingle.csv', 'a', newline='') as csvFile:
					writer = csv.writer(csvFile)
					writer.writerow([sys.argv[1], denStats[0], denStats[1], meanSum, stdDevSum, denStats[3]/(untilT*100*178)])
				csvFile.close()
		if calcPSNR:
			if not imgs1Loaded:
				imgsNpy = np.load(imgFileName)['data']
				copyArrayToGridReal(imgsNpy, imgs1)
			if not imgs2Loaded:
				imgsONpy = np.load(imgInputFileName)['data']
				copyArrayToGridReal(imgsONpy, imgs2)
			imgs1.sub(imgs2)
			with open(path +'psnr.csv', 'a', newline='') as csvFile:
				writer = csv.writer(csvFile)
				writer.writerow([sys.argv[1], psnr_den(imgs1, imgs2.getMaxAbs())])
			csvFile.close()
				
if meanVelGlobal: 
	devVelobal# = math.sqrt(devDenGlobal/(cntDenGlobal))				
	with open(path +'velStatsAll.csv', 'a', newline='') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow([sys.argv[1], devVelobal])
	csvFile.close()	
if meanImgsGlobal: 
	devImgsGlobal# = math.sqrt(devDenGlobal/(cntDenGlobal))				
	with open(path +'imgsStatsAll.csv', 'a', newline='') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow([sys.argv[1], devImgsGlobal])
	csvFile.close()	
if meanDenGlobal: 
	devDenGlobal# = math.sqrt(devDenGlobal/(cntDenGlobal))				
	with open(path +'denStatsAll.csv', 'a', newline='') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow([sys.argv[1], devDenGlobal])
	csvFile.close()	