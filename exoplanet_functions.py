# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 10:56:58 2024

@author: Evelyn Macdonald
"""

import numpy as np 
import matplotlib.pyplot as plt
from exoplasim import gcmt
import glob


#%%
def map2d(dataset,var,time=None,symmetric=None, newfig=True, **kwargs):
    '''
    Plot a pcolormesh of a variable.
    
    dataset: a year of exoplasim data.
    var: str
        climate variable with 2 spatial dimensions that you want to plot. Examples: 'ts' (surface temperature), 'pr' (precipitation), 'clt' (cloud cover),
        'prw' (water vapour), 'sit' (sea ice thickness), 'snd' (snow depth)
    time: None or int
        if None, average all the months of the year. If int, take that month only (e.g. 0 for January)    
    symmetric: float (optional)
        midpoint of colorbar
    **kwargs: additional arguments to pass to pcolormesh (optional)
    '''    
    
    lat = dataset.variables['lat']
    lon = dataset.variables['lon']
    
    if var == 'albt':     #top of atmosphere albedo
       albedo_toa = -(dataset.variables['rsut'])/((dataset.variables['rst']) - (dataset.variables['rsut']))
       variable = gcmt.make2d(albedo_toa,lat=lat,lon=lon,time=time)
       kwargs['cmap'] = plt.cm.Grays_r
       kwargs['vmax']=1
       kwargs['vmin']=0
       
    elif var == 'albs':    #surface albedo  
       albedo_surf = -(dataset.variables['ssru'])/((dataset.variables['rss']) - (dataset.variables['ssru'])) 
       variable = gcmt.make2d(albedo_surf,lat=lat,lon=lon,time=time)
       kwargs['cmap'] = plt.cm.Grays_r
       kwargs['vmax']=1
       kwargs['vmin']=0

    else:
        variable = gcmt.make2d(variable=dataset.variables[var],lat=lat,lon=lon,time=time)
    
    if var == 'sit':
        kwargs['vmax']=9
        kwargs['vmin']=0
        kwargs['cmap']=plt.cm.Blues_r
        
        
    if var == 'ts':
        symmetric=273   #plot temperature to be symmetric around 273K so you can see where it's frozen
        kwargs['cmap'] = plt.cm.seismic    #red-blue colormap
    if symmetric:
        if 'vmax' not in kwargs:
            vmax = symmetric + np.nanmax(abs(variable-symmetric))
        else:
            vmax = kwargs['vmax']
        vmin = 2*symmetric - vmax
        kwargs['vmax']=vmax
        kwargs['vmin']=vmin
    
    if newfig == True:
        plt.figure()
    pc = plt.pcolormesh(lon,lat, variable, shading='gouraud', **kwargs)
    plt.xlabel('Longitude',fontsize=16)
    plt.ylabel('Latitude',fontsize=16)
    if time == None:
        plt.title(var+' average',fontsize=16)
    else:
        plt.title(var+' month '+str(time),fontsize=16)
    plt.colorbar()
    return pc

#%%
def globalaverage(dataset,var,**kwargs):
    #average the variable over the whole planet
    lat = dataset.variables['lat']
    lon = dataset.variables['lon']
    if var == 'albt':     #top of atmosphere albedo
       variable = -(dataset.variables['rsut'])/((dataset.variables['rst']) - (dataset.variables['rsut']))   #what shape is this?
       v = [gcmt.spatialmath(variable,lat,lon,time=m) for m in range(len(dataset.variables['ts']))]
       
    elif var == 'albs':    #surface albedo  
       variable = -(dataset.variables['ssru'])/((dataset.variables['rss']) - (dataset.variables['ssru']))     
       v = [gcmt.spatialmath(variable,lat,lon,time=m) for m in range(len(dataset.variables['ts']))]

    else:
        v = [gcmt.spatialmath(dataset.variables[var],lat,lon,time=m) for m in range(len(dataset.variables[var]))]
    # y = gcmt.spatialmath(dataset.variables[var],lat,lon)
    
    plt.plot(np.arange(0,12,1),v,**kwargs)
    plt.xlabel('Month')
    plt.ylabel(str(var))
   

#%%
def profile(dataset,var,latindex=16,lonindex=32,psurf=1,month='avg',**kwargs):
    #Vertical profile of a variable at the desired latitude and longitude indices.
    #dataset: a year of exoplasim data
    #var: name of an output variable with 3 spatial dimensions, e.g. 'ta' (air temperature) or 'hus' (specific humidity)
    #month: either 'avg' for a yearly average, or a number from 0-11 for a specific month.

    variable = dataset.variables[var]  
    if month == 'avg':
        variable = np.mean(variable,axis=0)  #average variable in time dimension
    else:
        variable = variable[month]
    prof = variable[:,latindex,lonindex]  #value of the variable at the desired latitude and longitude for each vertical level
    pressure = psurf*dataset.variables['lev']
    
    plt.plot(prof,pressure)
    plt.gca().invert_yaxis()
    plt.yscale('log')
    plt.ylabel('Pressure (bar)')
    plt.xlabel(var)
    return(prof)

#%%
def zonal_mean(dataset, var, month='avg'):
    #month should either be 'avg' for a whole year average, or a number for a specific month
    lat = dataset.variables['lat']
    lon = dataset.variables['lon']
    variable = dataset.variables[var]
    if month == 'avg':   #average over the whole year
        mean = np.mean(gcmt.lonmean(variable, lon), axis=0)   #zonal mean averaged over the whole year
    else:
        mean = gcmt.lonmean(variable[month], lon)    #take that specific month only
    plt.plot(lat,mean)
    plt.xlabel('Latitude')
    plt.ylabel(var)

#%%
'''
Some examples using the functions defined above
'''
   
#load a data file 
data = gcmt.load('datafiles/earth.npz')
 
#plot the monthly average temperature 
globalaverage(data,'ts')

#average of all the months 
map2d(data,'ts')   

#month 11 
map2d(data,'ts',11)

# average air temperature as a function of vertical level
plt.figure()
profile(data,'ta')

plt.figure()
for month in range(12):
    zonal_mean(data,'ts',month)



