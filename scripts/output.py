#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 23:24:28 2018

@author: ad
"""
import numpy as np
#----------------------------------  Volume Calculations  ------------------------------#

def volcalcs(annualumg, annualumh, longshore_sup, shore_change):
    '''
    Calculates annual volumes and shoreline changes
    '''
    
    longshore_out = {}
    umgkeys = list(annualumg.keys())
    sedimnet = {}

    for year in range(len(analysisyears)):
        umgin = np.sum(np.array(annualumg[umgkeys[year]]))
        umhout = np.sum(np.array(annualumh[umgkeys[year]]))
        sedimnet[analysisyears[year]] = umhout - umgin
        longshore_in[analysisyears[year]] = umgin
        longshore_out[analysisyears[year]] = umhout
        
    return umgin, umhout, sedimnet, longshore_in, longshore_out
    
#----------------------------------  Plotting Function  --------------------------#

def plotcoast(xcoord, ylist, beachwidths, baselow, baseup, indices, freq):
    '''
    Plots coastline change at given intervals using freq
    '''
    
    import matplotlib.pyplot as itera

    left_plot = 0
    right_plot = len(ylist[0])

    itera.figure('Iterated Coastline', figsize=(24,4))
    itera.clf()
    #itera.figure(1)
    itera.xlabel('Distance Alongshore $(m)$', fontsize = 13)
    itera.ylabel('Distance from Baseline $(m)$', fontsize = 13)
    #itera.axis('equal')
    itera.ylim(0, 600) #max(coastline_y)*1.5)
    itera.xlim(0, max(xcoord)) 
    itera.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #itera.xticks(coastx)
    itera.fill_between(xcoord, beachwidths, baselow, color = 'green', alpha = 0.7) #return
    itera.fill_between(xcoord, ylist[0], beachwidths, color = 'yellow')
    itera.plot(xcoord[left_plot:right_plot], baseup[left_plot:right_plot], '-')
    itera.grid(alpha = 0.5)
    colors = ['b', 'g', 'r', 'm', 'y', 'k', 'w']
    itera.plot(xcoord[left_plot:right_plot], baselow[left_plot:right_plot], '-') #constant
    itera.fill_between(xcoord, baseup, ylist[0], color='cyan', alpha = 0.8) #return

    for newl in range(0, len(indices), freq): #, 8760): #, 30*24): #len(y_val_dict), 87600):
        graph = indices[newl]
        if graph == 0:
            itera.plot(xcoord[left_plot:right_plot], ylist[0][left_plot:right_plot], 'k-', label = 'Initial Coastline',linewidth=2) #marker = 'o'
            plot1 = 1
        elif graph != 0 and graph != len(ylist)-1 and plot1 == 1:
            itera.plot(xcoord[left_plot:right_plot], ylist[graph][left_plot:right_plot], 'r', linewidth=1.2, label='Iterated Coastline')
            plot1 = 0
        elif graph != 0 and graph != len(ylist) - 1 and plot1 != 1:
            itera.plot(xcoord[left_plot:right_plot], ylist[graph][left_plot:right_plot], 'r', linewidth=1.2)
        elif graph == len(ylist)-1:
            itera.plot(xcoord[left_plot:right_plot], ylist[graph][left_plot:right_plot], 'b', linewidth=2, label='Final Coastline')
            
    itera.legend(loc= 'lower left')
    
def plotgauss(gaussplumes, coastx, baselow, baseup, freq):
    '''
    plots evolution of gaussian plumes
    '''
    
    import matplotlib.pyplot as itera
    
    #itera.figure(2)
    itera.figure('Nourishment Advection', figsize=(20,4))
    itera.xlabel('Distance Alongshore (m)')
    itera.ylabel('Coastal Outcrop (m)')
    itera.ylim(0, 2) 
    itera.xlim(0, max(coastx)) 
    itera.title('Nourishment Advection')
    itera.fill_between(coastx, baseup, gaussplumes[0], color = 'cyan')
    itera.fill_between(coastx, gaussplumes[0], baselow, color = 'yellow')
    
    for line in range(0, len(gaussplumes), freq):
        itera.plot(coastx, gaussplumes[line], linewidth = 1.5, color = 'k')
        
def annualmean(array):
    '''
    Returns mean and st. deviation for net sediment into model domain over period
    '''
    
    vollist = np.array([array[i] for i in array])
        
    meanvol = np.mean(vollist)
    stdvol = np.std(vollist)
    
    return meanvol, stdvol

