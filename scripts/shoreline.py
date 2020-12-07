#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 17:55:39 2018

@author: ad
"""

#------------------------------  Coastline Change  --------------------------------------#
        
def change(Q, normal_node):
    '''
    Calculates coastline changes along grid nodes
    '''
    
    Qdelta = Q[1:]-Q[:-1]
    
    deltayq = np.nan_to_num(((-delta_t*3600)/(closure_depth + berm_height))*(Qdelta/stagger))
    
    return deltayq

def updatecoast(ycomp, coasty):
    '''
    returns updates coastline
    '''
    
    coastline_newy = coasty + ycomp
    
    return np.round(coastline_newy, 3)

#--------------------------  Seawall Boundary Code  ---------------------------#

def correct(Q, coasty):
    '''
    Corrects transport rates where the hard boundary has been violated
    '''
    
    Qdup = Q[:]
    coastydup = coasty[:]
        
    ysbeg = 0
    ysend = len(coastydup) - 1
    
    Bq = (delta_t*3600)/((closure_depth + berm_height)*stagger)
    
    i = ysbeg
    k = 0
    breaker = 0
    initial = 0
    
    #CODE-----------------
    
    while True:
        
        if Qdup[i] > 0 or initial == 0 or initial == 1: #if Q is positive, i.e. iterate from left to right
            
            while Qdup[i+1] >= 0: #10
                yip = coastydup[i] - Bq*(Qdup[i+1] - Qdup[i])
                ysi = beachwidths[i]
                if yip < ysi:
                    diff = ysi - yip
                    Qdup[i+1] = Qdup[i+1] - diff/Bq
                i += 1
                if i == ysend + 1:
                    breaker = 1
                    break #breaks out of inner loop
                    
            if i == ysend + 1 or breaker == 1:
                break
            
            k = i
            i += 1
            if i == ysend + 1 or breaker == 1:
                yip = coastydup[i-1] - Bq*(Qdup[i] - Qdup[i-1])
                breaker = 1 
                break #breaks out of main loop
                
            if i == ysend:
                i += 1
                for points in range(i-1, k, -1): 
                    yip = coastydup[points] - Bq*(Qdup[points+1] - Qdup[points])
                    ysi = beachwidths[points]
                    if yip <= ysi and points >= ysbeg:
                        diff = ysi - yip
                        Qdup[points] = Qdup[points] + diff/Bq
                
            initial = 1

        else:
            k = ysbeg - 1
            if ysbeg == 1:
                k = 1
            
        if i == ysend + 1:
            break
        
        while Q[i+1] < 0: #20
            i += 1
            if i == ysend:
                if Qdup[i+1] <=0:
                    yip = coastydup[i] - Bq*(Qdup[i+1] - Qdup[i])
                    ysi = beachwidths[i]
                    if yip < ysi:
                        diff = ysi - yip
                        Qdup[i] = Qdup[i] + diff/Bq
                    for points in range(i-1, k, -1): 
                        yip = coastydup[points] - Bq*(Qdup[points+1] - Qdup[points])
                        ysi = beachwidths[points]
                        if yip < ysi and points >= ysbeg:
                            diff = ysi - yip
                            Qdup[points] = Qdup[points] + diff/Bq
                    i += 1
                    if i >= ysend +1:
                        breaker = 1
                        break #breaks out of internal loop
                    
        if breaker == 1:
            break #breaks out of main loop
                        
        yip = coastydup[i] - Bq*(Qdup[i+1] - Qdup[i])
        ysi = beachwidths[i]
    
        if yip < ysi:
            diff = ysi - yip
            qdiff = Qdup[i+1] - Qdup[i]
            Qdup[i] = Qdup[i] - (diff/Bq)*(Qdup[i]/qdiff)
            Qdup[i+1] = Qdup[i+1] - (diff/Bq)*(Qdup[i+1]/qdiff)
    
        for points in range(i-1, k, -1): 
            yip = coastydup[points] - Bq*(Qdup[points+1] - Qdup[points])
            ysi = beachwidths[points]
            if yip < ysi and points >= ysbeg:
                diff = ysi - yip
                Qdup[points] = Qdup[points] + diff/Bq
        i += 1
        if i >= ysend + 1:
            break #breaks out of main loop

    return Qdup

#--------------------------------- Advection and Diffusion of Sand Supply -------------------------------#

def diffuse(gaussx, gaussy, Qg, normal):
    '''
    Diffuses gaussian plume
    '''
    
    Qdelta = Qg[1:] - Qg[:-1]
    delta = Qdelta/((closure_depth + berm_height)*stagger)
    
    deltaf = -(delta_t*3600)*delta
    
    gaussnewy = gaussy + deltaf

    return gaussnewy

#-----------------------------------------------------#
    
def advect(gaussx, gaussy, Qg, ktransg, hg, D, Hs, Tp, normal):
    '''
    Advects sediment inputs by decoupling processes and summing coastline coordinates
    '''
    dydx = (gaussy[2:] - gaussy[:-2])/(stagger*2)
    dydx = np.insert(dydx, 0, (gaussy[1] - gaussy[0])/(2*stagger))
    dydx = np.insert(dydx, len(dydx), (gaussy[-1] - gaussy[-2])/(2*stagger))
    
    tanbeta = (A**1.5)/(np.sqrt(np.mean(hg))) #use for kamphius
    angle = np.radians((coast_orien+90) - D)
    sintheta = abs(np.sin(2*angle))**0.6 #if wave transformation is used, kamphius
    signs = np.sign(np.sin(2*angle))
    Sx0 = signs*(const*(np.mean(Hs*ktransg)**2)*(Tp**1.5)*(tanbeta**0.75)*((d50/1000)**-0.25)*(sintheta))/denom
    
    B =  10 #max(gaussy)
    
    celerity = (Sx0/((closure_depth + berm_height)*B))*dydx
    
    deltatot = -(delta_t*3600)*celerity

    gaussnewy = gaussy + deltatot

    return gaussnewy

        
