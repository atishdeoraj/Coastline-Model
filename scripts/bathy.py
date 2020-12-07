#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 12:43:42 2018

@author: ad
"""
import numpy as np
from math import radians, degrees, sin, asin, sinh, cos, tan, exp, sqrt, tanh, atanh, cosh, atan, atan2, pi

#FUNCTIONS ----------------------------------------------------------------------------------

#Cross-shore Profile

def xshoreprofile(diamaters):
    '''
    calculates cross-shore profile constant A using an array of sediment diamaters
    '''
    d50 = np.mean(np.array(diamaters))

    if d50 < 0.4:
        profilesp = 0.41*(d50**(0.94))
    
    elif d50 >= 0.41 and d50 < 10:
        profilesp = 0.23*(d50**(0.32)) #function is depth(h) = Ay^(2/3) where y is distance away from shore

    return profilesp

#-------------------------------------------------------------------------------------
    
#Iterative Loop Calculations
    
#-----------------------  Coastline Orientation Calculation  ----------------------#
    
def slopes(x, y):
    '''
    Calculates individual grid cell orientation in radians
    '''
    
    coast_orien_rad = radians(coast_orien)
    
    ydistance = np.array(y[1:]) - np.array(y[:-1])
    xsegments = np.array(x[1:]) - np.array(x[:-1])
    
    slopess = 2*pi - np.arctan2(ydistance, xsegments) #cell inclination along shoreline
    normal_slope = slopess + coast_orien_rad + (pi/2) - 2*pi #cell normal seaward
    
    normals = ((normal_slope[1:] + normal_slope[:-1])/2) - pi/2 - coast_orien_rad #node normal
    normals = np.insert(normals, 0, 0)
    normals = np.round(np.insert(normals, len(normals), 0), 3)
    
    return normal_slope, normals
    
#----------------------------   Wave Transformation   --------------------------#
    
def transform(wavedirecs, Hs_dict, Tp_dict):
    '''
    Transforms deepwater waves to shallow or breaking conditions for shoaling and refraction,
    returns transformation coefficients and new wave angles
    '''
    unfinished_kr = 0
    
    Lo = (g*Tp_dict**2)/(2*pi) #deepwater wavelength
    
    steepness = Hs_dict/Lo #deepwater steepness
    
    if steepness <= 0.142:
        Hs_dict = Hs_dict
    else:
        Hs_dict = 0.142*Lo

    Ko = (4*(pi**2))/(g*(Tp_dict**2)) #deepwater wave number, constant
    c_deep = ((g*Tp_dict)/(2*pi)) #constant for all segments, only dependent on period

    h = np.array([4 for i in wavedirecs]) #initial assumption of 4m depth of breaking

    E_deep = np.around((1/8)*rho_water*g*(Hs_dict**2)*0.5*c_deep*abs(np.cos(wavedirecs)), 4)
    
    loop_count = 0

    while True:
        
        Lb_init = Tp_dict*np.sqrt(g*h*(Ko + (1 + 0.6522*Ko + 0.4622*(Ko**2) + 0.0864*(Ko**4) + 0.0675*(Ko**5))**-1)**-1)

        k_b = 2*pi/Lb_init
        n = 0.5*(1 + (2*k_b*h)/np.sinh(2*k_b*h))
        
        c_phase = ((g*Tp_dict)/(2*pi))*np.tanh(k_b*h) #transitional water conditions
        c_groupbv = n*c_phase
        
        theta_b = np.arcsin((c_phase*np.sin(wavedirecs))/c_deep)

        Ksh = np.sqrt((c_deep*0.5)/(c_phase*n))
        Kr = np.sqrt(abs(np.cos(wavedirecs))/abs(np.cos(theta_b)))
        K_trans = Ksh*Kr

        Hb = K_trans*Hs_dict
        E_shallow = np.around((1/8)*rho_water*g*(Hb**2)*n*c_phase*np.cos(abs(theta_b)), 4)
        breaker_ind = Hb/h
    
        if np.around(np.mean(breaker_ind), decimals=4) != 0.78:
            loop_count += 1
            h = Hb/0.78
            if loop_count > 100:
                unfinished_kr += 1
                break
        elif np.around(np.mean(breaker_ind), decimals=4) == 0.78:
            break
        
    anglediff = wavedirecs - theta_b
    
    np.place(theta_b, abs(anglediff)> pi/2, 0) #wave transformation limit checks
    #np.place(theta_b, abs(anglediff) > pi/4, 0) 
     
    return K_trans, theta_b, h, unfinished_kr
    
    
#---------------------------  Longshore Sediment Transport  -------------------------#
    
def lst(hsig, angles, h, ktrans, waved):
    '''
    Computes longshore sediment transport rates
    '''
    
    Hbs = hsig*ktrans
    tanbeta = (A**1.5)/(np.sqrt(h)) #use for kamphius
    sintheta = abs(np.sin(2*angles))**0.6 #if wave transformation is used, kamphius
    signs = np.sign(waved) 

    Q = signs*(const*(Hbs**2)*(Tp_dict**1.5)*(tanbeta**0.75)*((d50/1000)**-0.25)*(sintheta))/denom #theta_dict[ran_slope1]

    return Q
    
#-----------------------------  Boundary Conditions  ----------------------------------#
    
def boundary(Q, D_dict):
    '''
    Applies boundary conditions to each end of the domain
    '''
    
    inrate = 0/(365.25*24*3600)
    
    if D_dict <= coast_orien+90:
        
        Q = np.insert(Q, 0, Q[0]) #umhlanga, does not change
        Q = np.insert(Q, len(Q), -inrate) #umgeni, variable
        
        #longshore_input += Qdd[-1]*byp_dd*3600*delta_t

        annualumg['sedim' + str(dataset.year)].append(Q[-1]*delta_t*3600)
        annualumh['sedim' + str(dataset.year)].append(Q[zerop-1]*delta_t*3600)

    elif D_dict > coast_orien+90:
        
        Q = np.insert(Q, 0, Q[0]) #umhlanga, does not change
        Q = np.insert(Q, len(Q), -inrate) #umgeni, variable
        
        #longshore_input += Qud[0]*byp_ud*3600*delta_t
        
        annualumg['sedim' + str(dataset.year)].append(Q[-1]*delta_t*3600)
        annualumh['sedim' + str(dataset.year)].append(Q[zerop-1]*delta_t*3600)
        
    return Q

#------------------------------  Storm Sediment Supply  ---------------------------------#
    
def stormsed(date, storm_index):
    '''
    Generates storm sediment distribution curve for model domain
    '''

    print("Storm Event" + ' - ' + str(date))
    
    hourly_input = (stormvols[date]/storm_length)*delta_t
    dy = (hourly_input)/((closure_depth + berm_height)) #area under curve
    mu = coastx[-1] - width/2 #coastlinex_grid[(len(coastlinex_grid) - len(coastxgroyne))-1] #coastline_x[-10]
    stdev = width/8 #in metres
    prob2 = (1/(np.sqrt(2*pi)*stdev))*np.exp((-(coastx-mu)**2)/(2*(stdev**2)))
    ydistrib = prob2*dy
    ydistrib = np.round(ydistrib, 3)    
    
    input_hours = storm_length
    storm_index += 1
        
    if storm_index > len(stormdates)-1:
        storm_index = 0

    return ydistrib, input_hours, storm_index

#-------------------------------------------------------------------
#import shoreline
#
#x = coastx
#y = coasty
#
#ys = []
#ys.append(y)
#
#for i in range(0, 60):
#    
#    hs = 1.5
#    d = 150
#    tp = 10
#    
#    normal_seg_g2, normal_node_g2 = slopes(x, y) #calculate segment orientations
#    wave_direc_g2 = normal_seg_g2 - radians(d) #relative wave angles for each segment
#    ktrans_g2, thetab_g2, h_g2, unfinished_kr_g2 = transform(wave_direc_g2, hs, tp) #transform deepwater waves to shallow water
#
#    Q_g2 = lst(hs, thetab_g2, h_g2, ktrans_g2, wave_direc_g2) #lst rates, coastline + gaussian
#    Q_g2 = boundary(Q_g2, d) #applies boundary conditions to ends of domain
#      
#    delty = change(Q_g2, normal_node_g2)
#    ynew = updatecoast(delty, delty, y)
#    ynew2 = advect(x, ynew, Q_g2, ktrans_g2, h_g2, d, hs, tp, normal_node_g2)#advect gauss plume
#    np.place(ynew2, ynew2<0, 0)
#    #y=ynew
##    if np.any(ynew < beachwidths):
##        print('correction')
##        Qad = correct(Q_g2, y)
##        ynew2 = advect(x, y, Qad, ktrans_g2, h_g2, d, hs, tp, normal_node_g2) #np.round(diffuse(x, y, Qad, normal_node_g2), 3)
##        volcheck(Qad, y, ynew2)
##        y = ynew2
##    else:
#    volcheck(Q_g2, y, ynew2)
#    y = ynew2
#        
#    #y = advect(x, y, Q_g2, ktrans_g2, h_g2, d, hs, tp, normal_node_g2)
#    ys.append(y)
#
#import matplotlib.pyplot as plt
#
#plt.figure('Nourishment Advection', figsize=(20,4))
#for i in range(0, len(ys), 20):
#    plt.plot(x, ys[i])