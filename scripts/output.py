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
    itera.ylim(0, 2) #max(coastline_y)*1.5)
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

#plotdat = {}
#
#inputvol = np.arange(0, 600001, 100000)
#
#plotdat[0] = [-487653.6, -462197.2, -472541.9, -488409.2, -482426.5, -506907, -448531.1, -497182.8, -499022.2, -441875.9]
#plotdat[1] = [-381939]
#plotdat[2] = [-282007]
#plotdat[3] = [-180935, -154559.8, -178870.7, -173191.9, -159326, -167676.5, -180307.6, -180879.5, -196999.4]
#plotdat[4] = [-82143]
#plotdat[5] = [17791]
#plotdat[6] = [119064.6, 128014.8, 106729.1, 116447.5, 100108.4, 137048.3, 132324.9, 118093.3, 127761.3, 119256.9]
#
#means = np.array([np.mean(plotdat[i]) for i in range(7)])
#
#yt = np.arange(-600000, 200001, 100000)
#
#bestf = 0.9989*inputvol - 479834
#x = np.arange(-200000, 800001, 100000)
#y = np.array([0 for i in x])
    
#xer = [inputvol[0], inputvol[3], inputvol[6]]
#yer = [means[0], means[3], means[6]]
#
##yelo = np.array([-506907, -196999.4, 100108.4])
#yelo = np.array([28232.25999999995, 22249.800000000017, 20376.509999999995])
#yehi = np.array([36798.84, 20189.8 , 16563.39])
##yehi = np.array([-441875.9, -154559.8, 137048.3])
#
#
#yv = [yelo, yehi]
#
#import matplotlib.pyplot as plt
#
#plt.figure('compplot', figsize=(10,8))
#plt.grid(alpha = 0.3)
#plt.xlabel('Annual Sediment Supply $(m^3/year)$', fontsize = 15)
#plt.ylabel('Net Sediment Change $(m^3)$', fontsize = 15)
#plt.yticks(yt)
#plt.xticks(inputvol)
#plt.xlim(-50000, 650000)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
##plt.scatter(np.array([inputvol[0] for i in plotdat[0]]), plotdat[0])
##plt.scatter(np.array([inputvol[3] for i in plotdat[3]]), plotdat[3])
##plt.scatter(np.array([inputvol[-1] for i in plotdat[0]]), plotdat[6])
#plt.scatter(inputvol, means, color = 'k', s = 30)
#plt.scatter(418333.3333333333, -61960.83333333337, color = 'k', s=95, marker='^', label='Corbella & Stretch (2012) - 1')
#plt.scatter(461330.9641, -19010.5, color = 'k', s = 95, marker = 's', label='Corbella & Stretch (2012) - 2')
#plt.scatter(429924.7899, -50382.12739, color = 'k', s = 95, marker='P', label='CSCM (2017)')
#plt.plot(x, y, 'k-', linewidth = 0.5)
#plt.plot(inputvol, bestf, 'k--', linewidth=0.5, label=r'$y = 0.9989x + 479834$')
#plt.errorbar(xer, yer, yerr = yv, fmt = 'o', color = 'k', linewidth = 1)
#plt.legend(loc = 'lower right')


