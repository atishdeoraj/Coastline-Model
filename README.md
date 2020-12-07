# Coastline-Model
This repository contains the files required to run the coastline model developed by Atish Deoraj. 

*This software is free and and may be redistributed or modified to any extent by the user. This software has been developed for research purposes and
the developers take no responsibility should the code not perform in the desired way.

File Breakdown:

model.py        -   This is the main file that calls on functions from the partner files when running the model.
params.py       -   Contains all input information required such as the wave sequence, coastlie coordinates and sediment information.
initialise.py   -   This file initialises the simulation by generating the computational grid and interpolating wave sequences
bathy.py        -   Contains functions such as the wave transformation
shoreline.py    -   Contains functions for shoreline change
output.py       -   Contains code for generating output plots


