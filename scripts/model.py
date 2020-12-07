#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 16:14:30 2018

@author: ad
"""

'''
this file is the main file to run the full model, calls other scripts for information
'''

import os
import numpy as np
from math import radians, degrees, sin, asin, sinh, cos, tan, exp, sqrt, tanh, atanh, cosh, atan, atan2, pi
import random

import params
import bathy
import shoreline
import initialise
import output

#-------------------------------------- Model Initialisation ---------------------------------------

#coastline coordinates
coastx, coasty, beachwidths = gridgen(coastline_x, coastline_y, widthsx, widthsy)
yinterpolated.append(coasty)

#beachwidths[-1] = coasty[-1]

zerop = list(coastx).index(0)

#upper and lower plotting bounds
baseline_y = [-10000 for i in range(len(coastx))]
upper_line = [10000 for i in range(len(coastx))]

#initialise Gaussian tracking grid
gaussy = np.array([float(0) for i in coasty])
gaussx = coastx[:]

A = xshoreprofile(d50arr) #cross shore profile shape factor
d50 = np.mean(d50arr)

#Data Importing

data_dict, dateList = dataimport(file_name)
annualumg, annualumh, sedimnet, longshore_in, shore_change = annualvol(start_year, dateList)

#Storm Generation

stormdates, stormvols, yearvol = genstorm(river_range, storm_frequency, storm_length, width)
storm_index = 0 #tracks occurrence of storm events
input_hours = 0

iteration_count = int(analysis_p/delta_t) #number of one hour steps in simulation
iteration = 0 #tracks iteation count

#Storage Lists and Arrays for output

startd = dateList[0]
coastline_new = []
previous = 1990
analysisyears.append(previous)
print('1990')
Qumgeni = {}
Qumhlanga = {}

#annualumg = {}
#annualumh = {}
longshore_sup = {}
shore_change = {}
yinterpolated = []
river_input = 0

#------------------------------------- Iterative Model Loop  -----------------------------------------

for dataset in data_dict:
    
    if dataset == dt.datetime(dataset.year, 1, 1, 0, 0):
        plotindices.append(dateList.index(dataset))
    
    if dataset == simenddate:
        break
    
    if dataset.year != previous:
        print(dataset.year)
        analysisyears.append(dataset.year)
        yearindex += 1
    
    Hs_dict = data_dict[dataset][0] #wave data import from dictionary
    Tp_dict = data_dict[dataset][1]
    D_dict = data_dict[dataset][2]
    
    #--------------- Check storm occurrence ----------------------#
    if dataset == stormdates[storm_index]:
        gaussian, input_hours, storm_index = stormsed(dataset, storm_index) #storm sediment nourishment
        gaussy += gaussian #CHECK, DELETE FOR MODEL
        #break
        
    if input_hours == 0: #check to ensure storm length is applied correctly
        gaussian = np.array([0 for i in coastx])
    else:
        input_hours -= 1

    gaussy += gaussian #update gaussian plume with new inflow if necessary
    
#    if np.any(gaussy != 0): #gaussian input from storms, require decoupling
#    
#        #-------- Diffuse coastline without plume -----------#
#        coastyadj = coasty - gaussy #remove gaussian plume from coastline
#        normal_seg, normal_node = slopes(coastx, coastyadj) #calculate segment orientations
#        wave_direc = normal_seg - radians(D_dict) #relative wave angles for each segment
#        ktrans, thetab, h, unfinished_kr = transform(wave_direc, Hs_dict, Tp_dict) #transform deepwater waves to shallow water
#    
#        Q = lst(Hs_dict, thetab, h, ktrans) #lst rates, coastline without gaussian
#        Q = boundary(Q, D_dict) #applies boundary conditions to ends of domain
#        
#        #------- Diffuse and advect plume  ------------------#
#        
#        normal_seg_g, normal_node_g = slopes(gaussx, gaussy) #calculate segment orientations
#        wave_direc_g = normal_seg_g - radians(D_dict) #relative wave angles for each segment
#        ktrans_g, thetab_g, h_g, unfinished_kr_g = transform(wave_direc_g, Hs_dict, Tp_dict) #transform deepwater waves to shallow water
#
#        Q_g = lst(Hs_dict, thetab_g, h_g, ktrans_g, wave_direc_g) #lst rates, coastline + gaussian
#        Q_g = boundary(Q_g, D_dict) #applies boundary conditions to ends of domain
#        
#        changexg, changeyg = advect(gaussy, Q_g, ktrans_g, h_g, D_dict, Hs_dict, Tp_dict, normal_node_g) #advect gauss plume
#        gaussy = gaussupdate(gaussx, gaussy, changexg, changeyg)
#        #gaussy += gausschangey

    #else: #----------- Diffuse coastline, no plume involved -----------#
        
    normal_seg, normal_node = slopes(coastx, coasty) #calculate segment orientations
    wave_direc = normal_seg - radians(D_dict) #relative wave angles for each segment
    ktrans, thetab, h, unfinished_kr = transform(wave_direc, Hs_dict, Tp_dict) #transform deepwater waves to shallow water
        
    Q = lst(Hs_dict, thetab, h, ktrans, wave_direc) #calculate longshore sediment transport rates
    Q = boundary(Q, D_dict) #applies boundary conditions to ends of domain
    
    #--------- Reconstruct coastline with plumes --------#
    delta_y = change(Q, normal_node) #calculates coastline changes along node normals
    coastynew = updatecoast(delta_y, coasty) #updates coastline using all changes
    #coastynew = coastynew2 + gaussy #add advected plume back in. (0 array if no plume)
    
    #-------------- Adjust transport rates if necessary  -------------#
    if np.any(coastynew <= beachwidths) == True:
        Qadj = correct(Q, coasty) #adust Q rates
        deltay2 = change(Qadj, normal_node) #insert function to calculate change
        coastynew2 = updatecoast(deltay2, coasty) #update coast again
        coasty = coastynew2
    else:
        coasty = coastynew
        
    #--------------- Append lists to storage list  ------------------#
    yinterpolated.append(coasty) #appends to list of coastlines
    #gaussplumes.append(gaussy)
    
    iteration += 1
    previous = dataset.year
    
print('')
print('Simulation complete. Unstable steps: ' + str(unfinished_kr))
print('')

plotindices[-1] = plotindices[-1] - 1

plotcoast(coastx, yinterpolated, beachwidths, baseline_y, upper_line, plotindices, 5)

#plotgauss(gaussplumes, coastx, baseline_y, upper_line, 360)

umgin, umhout, sedimnet, longshore_in, longshore_out = volcalcs(annualumg, annualumh, longshore_sup, shore_change)

meanvol, stdvol = annualmean(sedimnet)

print(meanvol, stdvol)
