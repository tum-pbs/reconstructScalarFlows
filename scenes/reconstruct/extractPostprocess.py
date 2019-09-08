# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# extract and post-process single frames from captured images
# videos in path must have naming cam1.mp4 - cam4.mp4
#
# 1. adapt variable path
# 2. example call "./manta ../scenes/reconstruct/extractPostprocess.py 5 1 1 0"
#
import sys,cv2,os,numpy,re,math
import numpy as np

import sys
import _visualize as v
from manta import *

if len(sys.argv) < 4:
    print("USAGE: extractPostprocess.py PATH CAMCOUNT [STARTFRAME] [FRAMECOUNT]\n")
    print("CAMCOUNT        Number of cameras / videos")
    print("PREPIMAGES      Assemble, save, and visualize images in numpy arrays.\n")
    print("ASSEMBLEUNPROC  Do prepImages for unprocessed, raw input images as well.\n")
    print("ONLYPREP        If false, images are extracted, denoised, background-subtractd, or thresholded.")
    sys.exit(1)

path           = '/home/eckert/results/input/0813_80_0085/' #Path to the directory containing the video files
camCount       = int(sys.argv[1])
prepImages     = bool(int(sys.argv[2])) #True	
assembleUnproc = bool(int(sys.argv[3])) #True	
onlyPrep       = bool(int(sys.argv[4])) #False  	

startframe = 0
numOfFrames  = 162
averBrig = 0
maxBrig = 0
perfectDenoise = False

def extractFramesVideo(vName,path):
	if not os.path.isfile(path+"frame%04d.png" % (numOfFrames-1)):
		vidcap = cv2.VideoCapture(vName)
		success,image = vidcap.read()
		framenum = 0
		count = 0
		success = True
		while success:
			success,image = vidcap.read()
			framenum += 1
			if (framenum < startframe):
				continue
			if (numOfFrames > 0 and count >= numOfFrames):
				break

			if success:
				if not os.path.isfile(path+"frame%04d.png" % count):
					image = cv2.rotate(image,cv2.ROTATE_90_COUNTERCLOCKWISE)
					image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
					cv2.imwrite(path+"frame%04d.png" % count, image) 
				count += 1

def denoise(f,folderIn,folderOut):
	if not os.path.isfile(folderOut+f):
		image_gray = cv2.imread(folderIn+f, cv2.IMREAD_GRAYSCALE)
		filename = folderIn+f[:-8]+"%04d.png"
		currFrame = int(f[-8:-4])
		denoiseStrength = 3
		if perfectDenoise and int(f[-8:-4])>1 and os.path.isfile(filename%(currFrame+1)) and os.path.isfile(filename%(currFrame+2)) :
			images = [cv2.imread(filename%(currFrame-2), cv2.IMREAD_GRAYSCALE), cv2.imread(filename%(currFrame-1), cv2.IMREAD_GRAYSCALE), image_gray, cv2.imread(filename%(currFrame+1), cv2.IMREAD_GRAYSCALE), cv2.imread(filename%(currFrame+2), cv2.IMREAD_GRAYSCALE)]
			image_gray = cv2.fastNlMeansDenoisingMulti(images,math.floor(len(images)/2),len(images),None,denoiseStrength,7,21)
		else:
			image_gray = cv2.fastNlMeansDenoising(image_gray,None,denoiseStrength,7,21)
		cv2.imwrite(folderOut+f, image_gray)

def separateBackground(img1,imgCor,folderIn,folderOut):
	img2 = cv2.imread(imgCor, cv2.IMREAD_GRAYSCALE)
	return cv2.subtract(img1, img2)
	
for c in range(1,camCount+1):
	imgFolder = path+"cam%d/"%(c)
	denoisedFolder = path+"denoised%d/"%(c)
	postprocFolder = path+"postprocessed%d/"%(c)
	
	if not onlyPrep: 
		print("Extract frames for ", c)
		if not os.path.exists(imgFolder): os.makedirs(os.path.dirname(imgFolder))
		extractFramesVideo(path+"cam%d.mp4"%(c),imgFolder)
		print("Done")
	
	# read config file 
	if not onlyPrep:
		print("Denoise for ", c)
		if not os.path.exists(denoisedFolder): os.makedirs(os.path.dirname(denoisedFolder))
		for g in os.listdir(imgFolder):
			if g.endswith(".png"):
				denoise(g,imgFolder,denoisedFolder)
		print("Done")
	if not onlyPrep:
		print("Separate background and threshold for ", c)
		if not os.path.exists(postprocFolder): os.makedirs(os.path.dirname(postprocFolder))
		for f in os.listdir(denoisedFolder):
			if os.path.exists(denoisedFolder+"/frame0000.png"):
				if f.endswith(".png") and not os.path.isfile(postprocFolder+f): 
					img1 = cv2.imread(denoisedFolder+f, cv2.IMREAD_GRAYSCALE)
					img1 = separateBackground(img1,denoisedFolder+"/frame0000.png",denoisedFolder,postprocFolder)
					ret,img1 = cv2.threshold(img1,8,255,cv2.THRESH_TOZERO)
					cv2.imwrite(postprocFolder+f, img1)
		print("Done")

		
# now prepare images for reading in
if prepImages:
	print("Prep Images Start")
	scale = 1
	frameOffset = 11

	width = 1080
	height = 1920
	angles = 5
	
	folderIn 		= path+'cam'+"%d\\"
	folderOut 		= path+'cam'+"\\"
	if not os.path.exists(folderOut): os.makedirs(os.path.dirname(folderOut))
	
	folderInProc 	= path+'postprocessed'+"%d\\"
	folderOutProc	= path+'postprocessed'+"\\"
	if not os.path.exists(folderOutProc): os.makedirs(os.path.dirname(folderOutProc))

	imgName         = 'imgsUnproc_%06d.npz'
	imgNamePng      = 'imgsUnproc_%06d.jpg' 
	
	imgNameProc     = 'imgs_%06d.npz'
	imgNamePngProc  = 'imgs_%06d.jpg'

	imgsArray = np.empty(shape=[angles, height, width, 1], order='C')
		
	if assembleUnproc:
		print("Prep Images Unproc")
		for t in range(0,numOfFrames):
			if not os.path.isfile(folderIn%(1)+"frame%04d.png"%t): break
			if not os.path.isfile(folderOut+imgName%t): 
				for i in range(angles):
					img = cv2.imread(folderIn%(i+1)+"frame%04d.png"%t, cv2.IMREAD_GRAYSCALE)
					imgsArray[i, :, :, 0] = cv2.flip(img, 0)
				imgsArray[:,:,:,:] = scale*(1./255.)*imgsArray[:,:,:,:]
				np.savez_compressed(folderOut+imgName%t, data=imgsArray)
				
			if not os.path.isfile(folderOut+imgNamePng%t):
				v.draw2DDensityNpy(folderOut+imgName%t, 0.8, False, False, True, False) 
			
	print("Prep Images Proc")		
	for t in range(frameOffset,numOfFrames):
		if not os.path.isfile(folderInProc%(1)+"frame%04d.png"%t): break
		if not os.path.isfile(folderOutProc+imgNameProc%(t-frameOffset)): 
			for i in range(angles):
				img = cv2.imread(folderInProc%(i+1)+"frame%04d.png"%t, cv2.IMREAD_GRAYSCALE)
				imgsArray[i, :, :, 0] = cv2.flip(img, 0)
			imgsArray[:,:,:,:] = scale*(1./255.)*imgsArray[:,:,:,:]
			np.savez_compressed(folderOutProc+imgNameProc%(t-frameOffset), data=imgsArray)
			
		if not os.path.isfile(folderOutProc+imgNamePngProc%(t-frameOffset)):
			v.draw2DDensityNpy(folderOutProc+imgNameProc%(t-frameOffset), 0.8, False, False, True, True) 
	print("Prep Images End")