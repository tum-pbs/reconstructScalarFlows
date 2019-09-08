# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# plot density and velocity fields 
#

import math, time, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import MaxNLocator

import uniio

# some global variables to reuse axes for all figures
didSetupV   = False
didSetupD   = False
didSetupI   = False

fsGlobal = (360*3.4/96, 720/96)

fV, axGlobalV = plt.subplots(1, 3, figsize=fsGlobal, dpi=96, facecolor='b')
fD, axGlobalD = plt.subplots(1, 3, figsize=fsGlobal, dpi=96)

XYVGlobal = 0
XYDGlobal = 0
XYIGlobal = 0
fI = 0
axGlobalI = 0

def setupV(gs):
	global didSetupV
	didSetupV = True
	
	global XYVGlobal
	
	x = np.linspace(0, gs[0], gs[0]) 
	y = np.linspace(0, gs[1], gs[1]) 
	ySquare = np.linspace(0, gs[0], gs[0]) 
	XYVGlobal = [np.meshgrid(x, y), np.meshgrid(x, ySquare), np.meshgrid(x, y)]
	
	for c in range(0,3):
		square = (c==2)
		gsY = gs[0] if square else gs[1]
		
		pos = [0.029+0.33*c,0.,0.3,1]
		axGlobalV[c].set_position(pos,'both')
		
		axGlobalV[c].get_yaxis().set_major_locator(MaxNLocator(integer=True))
		axGlobalV[c].get_xaxis().set_major_locator(MaxNLocator(integer=True))
		axGlobalV[c].set_aspect('equal')
		axGlobalV[c].xaxis.set_ticks(np.arange(0, gs[0]+1, math.ceil((gs[0]-1)/5)))
		if square: axGlobalV[c].yaxis.set_ticks(np.arange(0, gsY+2, math.ceil((gsY-1)/5)))
		else: axGlobalV[c].yaxis.set_ticks(np.arange(0, gsY+2, math.ceil((gsY-1)/3))+1)
		axGlobalV[c].set_xlim([0, gs[0]])
		axGlobalV[c].set_ylim([0, gsY+1]) 
		axGlobalV[c].spines['bottom'].set_color('#dddddd')
		axGlobalV[c].spines['top'].set_color('#dddddd') 
		axGlobalV[c].spines['right'].set_color('#dddddd')
		axGlobalV[c].spines['left'].set_color('#dddddd')
		axGlobalV[c].tick_params(axis='x', colors='#dddddd')
		axGlobalV[c].tick_params(axis='y', colors='#dddddd')
		axGlobalV[c].yaxis.label.set_color('#dddddd')
		axGlobalV[c].xaxis.label.set_color('#dddddd')
		axGlobalV[c].set_axis_bgcolor('black')
		
def setupD(gs):
	global didSetupD
	didSetupD = True
	
	global XYDGlobal
	
	x = np.linspace(0, gs[0]+1, gs[0]+1) 
	y = np.linspace(0, gs[1]+1, gs[1]+1) 
	ySquare = np.linspace(0, gs[0]+1, gs[0]+1) 
	XYDGlobal = [np.meshgrid(x, y), np.meshgrid(x, ySquare), np.meshgrid(x, y)]
	
	for c in range(0,3):
		square = (c==2)
		gsY = gs[0] if square else gs[1]
		
		pos = [0.029+0.33*c,0.,0.3,1]
		axGlobalD[c].set_position(pos,'both')
		axGlobalD[c].get_yaxis().set_major_locator(MaxNLocator(integer=True))
		axGlobalD[c].get_xaxis().set_major_locator(MaxNLocator(integer=True))
		axGlobalD[c].set_aspect('equal')
		axGlobalD[c].xaxis.set_ticks(np.arange(0, gs[0]+2, math.ceil((gs[0])/5)))
		if square: axGlobalD[c].yaxis.set_ticks(np.arange(0, gsY+2, math.ceil((gsY)/5)))
		else: axGlobalD[c].yaxis.set_ticks(np.arange(0, gsY+2, math.ceil((gsY-1)/3))+1)
		axGlobalD[c].set_xlim([0, gs[0]+1])
		axGlobalD[c].set_ylim([0, gsY+1]) 
		axGlobalD[c].spines['bottom'].set_color('#dddddd')
		axGlobalD[c].spines['top'].set_color('#dddddd') 
		axGlobalD[c].spines['right'].set_color('#dddddd')
		axGlobalD[c].spines['left'].set_color('#dddddd')
		axGlobalD[c].tick_params(axis='x', colors='#dddddd')
		axGlobalD[c].tick_params(axis='y', colors='#dddddd')
		axGlobalD[c].yaxis.label.set_color('#dddddd')
		axGlobalD[c].xaxis.label.set_color('#dddddd')
		
def setupI(gs, numImgs):
	global didSetupI
	didSetupI = True
	
	global XYIGlobal
	global fI 
	global axGlobalI
	
	# images
	fI, axGlobalI = plt.subplots(1, numImgs, figsize=(540*numImgs/96, 960/96), dpi=96)
	
	x = np.linspace(0, gs[0]+1, gs[0]+1) 
	y = np.linspace(0, gs[1]+1, gs[1]+1) 
	XYIGlobal = []
	
	for c in range(0,numImgs):
		XYIGlobal.append(np.meshgrid(x, y))
	
		pos = [0.017+c/numImgs,0.,1/numImgs-0.025,1]
		axGlobalI[c].set_position(pos,'both')
		axGlobalI[c].get_yaxis().set_major_locator(MaxNLocator(integer=True))
		axGlobalI[c].get_xaxis().set_major_locator(MaxNLocator(integer=True))
		axGlobalI[c].set_aspect('equal')
		axGlobalI[c].xaxis.set_ticks(np.arange(0, gs[0]+2, math.ceil((gs[0])/5)))
		axGlobalI[c].yaxis.set_ticks(np.arange(0, gs[1]+2, math.ceil((gs[1])/3)))
		axGlobalI[c].set_xlim([0, gs[0]+1])
		axGlobalI[c].set_ylim([0, gs[1]+1]) 
		
def adjustSubplots(image, c):
	if c==0:
		c2=c
	elif c==1: 
		c2=2
	elif c==2:
		c2=1
		image = np.transpose(image, (1, 0, 2))

	return image, c2

	
def getVelCmps(image, c2):
	if c2==0 :
		U = image[:,:,0]
		V = image[:,:,1]
	elif c2==2:
		U = image[:,:,0]
		V = image[:,:,2]
	else: 
		U = image[:,:,2]
		V = image[:,:,1]

	return U, V
	
def draw2DDensityNpy(filenameIn, scale, synthReal, negativeValues=False, dark=True, useColorMap=True):
	useMap = 'viridis' if synthReal else plt.cm.gnuplot2
	if not useColorMap: scale = 1
	colorMapToUse = useMap if useColorMap else 'gray'
	if dark: plt.style.use('dark_background')
	
	content = np.load(filenameIn)['data']
	
	shapeAdj = [content.shape[2], content.shape[1]]
	if not didSetupI: setupI(shapeAdj,content.shape[0])

	vmin = -scale if negativeValues else 0
	
	plots = []
	for c in range(content.shape[0]):
		image = content[c,:,:,0]

		plots.append(axGlobalI[c].pcolormesh(XYIGlobal[c][0], XYIGlobal[c][1],  image[:,:], cmap = colorMapToUse, vmin = vmin, vmax = scale))
	
	fI.savefig(filenameIn[:-3]+'jpg')
	for p in plots: p.remove()

def draw3DDensityGridNpy(filenameIn, scale, synthReal, negativeValues=False):
	plt.style.use('dark_background')
	useMap = 'viridis' if synthReal else plt.cm.gnuplot2
	
	content = np.load(filenameIn)['data']

	if not didSetupD: setupD(content.shape)

	vmin = -scale if negativeValues else 0
	
	plots = []
	for c in range(3):
		image = np.average(content, axis=c)
		
		image, c2 = adjustSubplots(image, c)
		
		plots.append(axGlobalD[c2].pcolormesh(XYDGlobal[c][0], XYDGlobal[c][1], image[:,:,0], cmap = useMap, vmin = vmin, vmax = scale))

	fD.savefig(filenameIn[:-3]+'jpg')
	for p in plots: p.remove()

def draw3DVelGridNpy(filenameIn, scale, synthReal):
	plt.style.use('dark_background')
	useMap = 'cool_r' if synthReal else 'spring'
	
	content = np.load(filenameIn)['data']

	if not didSetupV: setupV(content.shape)

	plots = []
	for c in range(3):
		image = np.average(content, axis=c)
		
		image, c2 = adjustSubplots(image, c)
		
		U, V = getVelCmps(image, c2)
		plots.append(axGlobalV[c2].quiver(XYVGlobal[c][0], XYVGlobal[c][1], U, V, np.arctan2(np.abs(V),np.abs(U)), cmap=useMap, headwidth=0, headlength=0, headaxislength=0, angles='xy', scale_units='xy', scale=scale, width=0.004))

	fV.savefig(filenameIn[:-3]+'jpg')
	for p in plots: p.remove()
	
	