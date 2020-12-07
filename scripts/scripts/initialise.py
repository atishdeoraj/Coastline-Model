#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 12:49:06 2018

@author: ad
"""

import numpy as np
import os
import random
import datetime as dt

#Coastline Grid Generation

#Generate coastline grid sequence at defined stagger using interpolation

def gridgen(xcoord, ycoord, xwidth, ywidth):
    '''
    Generates equally spaced grids for coastline and erosion limits
    '''
    coastx = np.arange(min(xcoord), max(xcoord)+stagger, stagger) #these values remain constant throughout iterations by constant interpolation
    coasty = np.interp(coastx, xcoord, ycoord) #interpolate original grid into more detailed grid using initial y coords

    ywidths = np.interp(coastx, xwidth, ywidth)
    
    beachwidth = coasty - ywidths
    
    return coastx, coasty, beachwidth

#----------------------------------------------------------------------------------------------

#Data Import 
    
def dataimport(filename):
    '''
    Imports data from file and generates a time series at defined time step
    '''
    
    waveHeight = []
    wavePeriod = []
    waveDirection = []
    
    f = open(os.path.expanduser("~ad/MSc/data/Height/List_" + str(file_name) + '.dat'))
    for point in f:
        waveHeight.append(float(point))

    gx = open(os.path.expanduser("~ad/MSc/data/Period/List_" + str(file_name)+ '.dat'))
    for element in gx:
        wavePeriod.append(float(element))
    
    hx = open(os.path.expanduser("~ad/MSc/data/Direction/List_" + str(file_name)+ '.dat'))
    for element in hx:
        waveDirection.append(float(element))
    
    print('')
    if len(waveHeight) == len(waveDirection) == len(wavePeriod):
        print('Imported data okay.')
    else:
        print('Error importing wave data.')
    print('')

    start_year = 1990                   #   enter value here, rest is based on this
    end_year = start_year + 100         #   produces essentially a 101 year count, counting the first year 
    
    #creates a list of all time steps between start and end date using a calendar

    import datetime as dt

    steps = int(24/delta_t)

    sDate = dt.datetime(start_year,1,1)
    eDate = dt.datetime(end_year+1,1,1)
    
    gap = (eDate-sDate + dt.timedelta(seconds=delta_t*3600)).days
    dateList = [sDate+dt.timedelta(seconds=delta_t*3600*i) for i in range(steps*gap)]
            
    dateList2 = [sDate+dt.timedelta(seconds=6*3600*i) for i in range(4*gap)]

    #links Hs, Tp and D data to each timestep in datetime dictionary

    data_dict1 = {} #stores in the form: time step: Hs, Tp, D

    for i in range(len(dateList2)): #accounts for timestep thats smaller than 6 hours
    
        Hs = waveHeight[i]
        Tp = wavePeriod[i]
        D = waveDirection[i]

        data_dict1[dateList2[i]] = Hs, Tp, D

    keys1 = sorted(data_dict1.keys())

    data_dict = {} #creates a dictionary that fills in steps between 6 hour steps

    sub_steps = int(data_step/delta_t)
    print('Generating data for dictionary.')
    print('')

    for date in range(0, len(dateList), sub_steps): #last timestep is len(dateList)-substeps
    
        keyc = 0
    
        key_prev = dateList[date] #0

        Hs_time = data_dict1[key_prev][0]
        Tp_time= data_dict1[key_prev][1]
        D_time = data_dict1[key_prev][2]
    #print(key_prev, str(Hs_time) + ', ' + str(Tp_time) + ', ' + str(D_time)) 

        if date == len(dateList)-sub_steps:
        
            Hs_next = Hs_time
            Tp_next = Tp_time
            D_next = D_time
         
        else:
            key_next = dateList[date+sub_steps]

            Hs_next = data_dict1[key_next][0]
            Tp_next = data_dict1[key_next][1]

        temp_Hs = []
        temp_Tp = []
        temp_D = []
    
        #obtain interpolated values as inputs for each smaller timestep
        for step in range(sub_steps):
            Hs_step = Hs_time + step*((Hs_next - Hs_time)*(delta_t/data_step))
            temp_Hs.append(Hs_step)
            Tp_step = Tp_time + step*((Tp_next - Tp_time)*(delta_t/data_step))
            temp_Tp.append(Tp_step)
            temp_D.append(D_time)
        
        if date != len(dateList)-sub_steps:
            for key in range(sub_steps):
                data_dict[dateList[date+keyc]] = temp_Hs[key], temp_Tp[key], temp_D[key]
                keyc += 1
        
        else:
            for key in range(sub_steps):
                data_dict[dateList[date+keyc]] = temp_Hs[key], temp_Tp[key], temp_D[key]
                keyc += 1
                data_dict[dateList[len(dateList)-1]] = temp_Hs[0], temp_Tp[0], temp_D[0]

    waveHeight = np.array(waveHeight)
    wavePeriod = np.array(wavePeriod)
    waveDirection = np.array(waveDirection)

    Hsmean = np.mean(waveHeight)
    Tpmean = np.mean(wavePeriod)
    Dmean = np.mean(waveDirection)

    #lists that may be cleared
    waveDirection = []
    wavePeriod = []
    waveHeight = []
    dateList2 = []
    data_dict1 = {}
    
    return data_dict, dateList

#-----------------------------------  Annual Volume Lists ------------------------------------#

def annualvol(startyear, dateList):
    '''
    Generate annual volume dictionaries for storage and output
    '''
    simyears = []

    if analysis_p <= 8760:
        simyears.append(start_year)
    elif analysis_p > 8760:
        for year in range(start_year, dateList[analysis_p].year+1):
            simyears.append(year)
        
    annualumg = {}
    annualumh = {}
    sedimnet = {}
    longshore_in = {}
    shore_change = {}

    for year in simyears: 
        annualumg['sedim' + str(year)] = []
        annualumh['sedim' + str(year)] = []
        sedimnet[str(year)] = []
        longshore_in[str(year)] = []
        shore_change[str(year)] = []
    
    return annualumg, annualumh, sedimnet, longshore_in, shore_change

#------------------------------------- Storm Series Generation  ------------------------------------------------

def genstorm(volrange, freq, length, width):
    '''
    Generates a storm sequence for the specified conditions in params
    '''
    
    simyears = []

    if analysis_p <= 8760:
        simyears.append(start_year)
    elif analysis_p > 8760:
        for year in range(start_year, dateList[analysis_p].year+1):
            simyears.append(year)

    #generate sequence of storm events for each year

    seasrain = np.array([188.8, 265.9, 130.6, 81.4]) #rainfall values for spring, summer, autumn, winter
    rainprob = seasrain/np.sum(seasrain)

    if storm_frequency == 1:
        stormprob = np.array([0, 1, 0, 0])
        rainprob = stormprob
    else:
        stormprob = np.round(rainprob*storm_frequency, 0) #calculate distribution of storm for seasons

    if np.sum(stormprob) != storm_frequency:
        stormprob = np.round(rainprob*storm_frequency, 1)
        remainders = np.round(np.sum(stormprob%1), 0)
        stormprob[1] += remainders
    
    elif np.sum(stormprob) == storm_frequency:
        stormprob = stormprob

    stormdates = []
    stormvolumes = {}
    yearvol = []

    for year in simyears: #range(dateList[0].year, dateList[analysis_p].year+1):
    
        annualvol = random.randint(volrange[0], volrange[1]) #select random annual supply volume in specified range
        yearvol.append(annualvol)
        seasvol = (annualvol*stormprob)/np.sum(stormprob) #calculates sediment release for individual storms in given seasons
    
        for i in range(int(stormprob[0])): #spring
            date = dateList[random.randint(14569, 16776)]
            stormdates.append(dt.datetime(year, date.month, date.day, date.hour))
            stormvolumes[dt.datetime(year, date.month, date.day, date.hour)] = seasvol[0]/stormprob[0]
        for i in range(int(stormprob[1])): #summer
            date = dateList[random.randint(8016, 10176)]
            stormdates.append(dt.datetime(year, date.month, date.day, date.hour))
            stormvolumes[dt.datetime(year, date.month, date.day, date.hour)] = seasvol[1]/stormprob[1]
        for i in range(int(stormprob[2])): #autumn
            date = dateList[random.randint(10176, 12384)]
            stormdates.append(dt.datetime(year, date.month, date.day, date.hour))
            stormvolumes[dt.datetime(year, date.month, date.day, date.hour)] = seasvol[2]/stormprob[2]
        for i in range(int(stormprob[3])): #winter
            date = dateList[random.randint(12384, 14569)]
            stormdates.append(dt.datetime(year, date.month, date.day, date.hour))
            stormvolumes[dt.datetime(year, date.month, date.day, date.hour)] = seasvol[3]/stormprob[3]

    stormdates.sort()
    
    return stormdates, stormvolumes, yearvol









