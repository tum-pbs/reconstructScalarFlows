# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# @author: Marie-Lena Eckert http://marielenaeckert.com/
#
#
# write json file with summary of reconstruction parameters
#
    
import os, sys, shutil, time, math, platform, datetime
from decimal import Decimal
from enum import Enum

class reconKind(Enum):
	synth     = 1 # synthetic data + addSource,  real or orthographic cameras, no cutoff, true first density field
	synthReal = 2 # synthetic data + vel Inflow, real cameras, 				   cutoff
	real      = 3 # real data, real inflow       real cameras, 				   cutoff

def writeJasonFile(path,calibFolder,captureFolder,folderOut,rK,res,factorY,p0,p1,scale,tomoParams,tomoParams_firstDen,tomoParams_trgt,ofParams,angles,width,height,stepSize,minNumCams,orthographic,restartedRecon,startFrame):
    now = datetime.datetime.now()

    if restartedRecon: file = open(path+'description_%06d.json'%startFrame,'w') 
    else: file = open(path+'description.json','w') 
    file.write('{\n')
    file.write('\t"grids": [\n')
    file.write('\t\t"density",\n')
    file.write('\t\t"density_noInflow",\n')
    file.write('\t\t"velocity",\n')
    file.write('\t\t"velocity_noInflow",\n')
    file.write('\t\t"imgsRendered",\n')
    file.write('\t\t"imgsRendered_noInflow",\n')
    file.write('\t\t"imgsTarget",\n')
    file.write('\t\t"imgsTarget_noInflow",\n')
    file.write('\t],\n')
    file.write('\t"creation_date": "%s",\n'%now.strftime("%Y-%m-%d %H:%M:%S")) 
    file.write('\t"calibration": "%s",\n'%calibFolder)
    file.write('\t"capture": "%s",\n'%captureFolder)
    file.write('\t"resolution": %d,\n'%res)
    file.write('\t"y-factor": %.2f,\n'%factorY)
    file.write('\t"dimension": %d,\n'%3)
    file.write('\t"inflow_start": "(%d,%d,%d)",\n'%(p0.x,p0.y,p0.z))
    file.write('\t"inflow_end": "(%d,%d,%d)",\n'%(p1.x,p1.y,p1.z))
    if rK == reconKind.synth:   file.write('\t"reconstruction_type": "synthetic",\n')
    elif rK == reconKind.synthReal: file.write('\t"reconstruction_type": "synthetic with real settings",\n')
    elif rK == reconKind.real: file.write('\t"reconstruction_type": "real",\n')
    else: file.write('\t"reconstruction_type": "real",\n')
    file.write('\t"output_folder": "%s",\n'%folderOut)
    file.write('\t"upscale_inputData": %.1f,\n'%scale)
    file.write('\t"of_smoothness_weight": %.1e,\n'%Decimal(ofParams.getSmooth()))
    file.write('\t"of_kinetic_weight": %.1e,\n'%Decimal(ofParams.getKinetic()))
    file.write('\t"vel_inflow": %.1e,\n'%Decimal(ofParams.getInflowValue()))
    file.write('\t"initial_tomography_smoothness_weight": %.1e,\n'%Decimal(tomoParams_firstDen.getSmooth()))
    file.write('\t"initial_tomography_kinetic_weight": %.1e,\n'%Decimal(tomoParams_firstDen.getKinetic()))
    file.write('\t"target_tomography_smoothness_weight": %.1e,\n'%Decimal(tomoParams_trgt.getSmooth()))
    file.write('\t"target_tomography_kinetic_weight": %.1e,\n'%Decimal(tomoParams_trgt.getKinetic()))
    file.write('\t"tomography_smoothness_weight": %.1e,\n'%Decimal(tomoParams.getSmooth()))
    file.write('\t"tomography_kinetic_weight": %.1e,\n'%Decimal(tomoParams.getKinetic()))
    file.write('\t"tomography_smoothness_weight_inflow": %.1e,\n'%Decimal(tomoParams.getSmoothInflow()))
    file.write('\t"tomography_kinetic_weight_inflow": %.1e,\n'%Decimal(tomoParams.getKineticInflow()))
    file.write('\t"tomography_thresh_visual_hull": %.1e,\n'%Decimal(tomoParams.getThreshVH()))
    file.write('\t"tomography_thresh_mask": %.1e,\n'%Decimal(tomoParams.getThreshMask()))
    file.write('\t"tomography_num_angles": %d,\n'%len(angles))
    file.write('\t"tomography_img_width": %d,\n'%width)
    file.write('\t"tomography_img_height": %d,\n'%height)
    file.write('\t"tomography_ray_stepsize": %.1f,\n'%stepSize)
    file.write('\t"tomography_minimal_number_seeing_cameras": %d,\n'%minNumCams)
    file.write('\t"tomography_orthographic": %s\n'%str(orthographic).lower()) 

    file.write('}\n')
    file.close() 
    