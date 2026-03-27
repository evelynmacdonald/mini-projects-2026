# ---------------------------------------------------------------
# Script to convert netCDF climate files into PSG GCM files
# netCDF files in LMD climate model
# Villanueva, Fauchez - NASA Goddard Space Flight Center
# February 2021
# Adapted to convert files from ExoPlaSim 
# Evelyn Macdonald
# August 2025
# ---------------------------------------------------------------
import sys
import numpy as np



#%%
def load_config(configpath):   #open the exoplasim configuration file and extract relevant information
	with open(configpath,"r") as cfgf:
		cfg = cfgf.read().split("\n")
	flux=float(cfg[1])
	if cfg[2] == None:
		startemp = 5777
	else:
		startemp = float(cfg[2])
	gases = cfg[4].split("&")
	for gas in gases:
		species = gas.split("|") 
		amt = float(species[1])
		if species[0] == 'pN2':
			pN2 = amt
		elif species[0] == 'pCO2':
			pCO2 = amt
	gascon = float(cfg[5])
	rotationperiod = float(cfg[9])
	gravity = float(cfg[13])
	radius = float(cfg[14])
    
    
#%%
def convert(filename,path='datafiles/',configpath='datafiles/',templatepath='psg_template.txt',itime=0,surfaces=None):
    savename = filename
    if surfaces != None:
        savename += '_surf'
    fileout = 'psg_configs/'+savename+'.dat'  #output file we're going to save
    templ = {}
    with open(templatepath,"r") as f:     #load template psg config file
        templ = dict([line.strip('<\n').split('>') for line in f])
    with open(configpath+filename+'.cfg',"r") as cfgf:     #load exoplasim config file
        cfg = cfgf.read().split("\n")
    
    #extract 3D data from .npz file
    ptotal= float(cfg[6])  #prescribed surface pressure
    file = np.load(path+filename+'.npz')
    print(file)

    lat = file.f.lat
    lon = file.f.lon
    ps = file.f.ps[0]/1000
    print(ps.shape)
    tsurf = file.f.ts[itime]
    print(tsurf.shape)
    temp = file.f.ta[itime,::-1]   #3d variables need to be inverted in the pressure dimension
    u = file.f.ua[itime,::-1]
    v = file.f.va[itime,::-1]
    press3D = (file.f.flpr[0,::-1])/100000
    h2o = file.f.hus[itime,::-1]*press3D*28/18   
    h2o[h2o==0] = 1e-30   #fill zeros with a very small number so the log doesn't make NaNs and crash everything

    #extract mixing ratios 
    gaseslist = ['N2','CO2','H2O']
    gases = {}
    units = {}
    gases['N2'] = 98
    units['N2'] = 'pct'
    gases['CO2'] = 1
    units['CO2'] = 'pct'
    gases['H2O'] = 1
    units['H2O'] = 'scl'

    #edit psg config file
    templ['GENERATOR-RESOLUTION'] = 100
    if cfg[2] != 'None':
        templ['OBJECT-STAR-TEMPERATURE'] = int(cfg[2])
        templ['OBJECT-STAR-RADIUS'] = 0.7   #assuming 4500K K-dwarf
    else:
        templ['OBJECT-STAR-TEMPERATURE'] = 5777   #Sun is the default star
    templ['ATMOSPHERE-PRESSURE'] = ptotal
    templ['ATMOSPHERE-WEIGHT'] = 8314/float(cfg[5])
    templ['ATMOSPHERE-LAYERS'] = len(file.f.lev)
    templ['ATMOSPHERE-NGAS'] = str(len(gases))
    templ['ATMOSPHERE-GAS'] = 'N2,CO2,H2O'
    templ['ATMOSPHERE-ABUN'] = ','.join(str(gases[gas]) for gas in gases)
    templ['ATMOSPHERE-UNIT'] = ','.join(str(units[gas]) for gas in gases)
    templ['ATMOSPHERE-LAYERS-MOLECULES'] = ''   
    templ['ATMOSPHERE-NAERO'] = 0

    if surfaces == None:
        templ['SURFACE-NSURF'] = 0
    else:
        templ['SURFACE-NSURF'] = len(surfaces)
        templ['SURFACE-SURF'] = 'Snow,Ocean'
        templ['SURFACE-TYPE'] = 'Albedo_GSFC,Albedo_GSFC'
        templ['SURFACE-ABUN'] = '1,1'
        templ['SURFACE-UNIT'] = 'scl,scl'


        snow = np.zeros_like(file.f.lsm[itime])  #ice and snow are the same for now
        snow[np.where(file.f.sic[itime]>0)] = 1   
        ocean = np.ones_like(file.f.lsm[itime])   #ocean everywhere except where there's snow
        ocean[np.where(snow==1)] = 0
    
    #GCM parameters
    vars = '<ATMOSPHERE-GCM-PARAMETERS>'
    vars = vars + str("%d,%d,%d,%.1f,%.1f,%.2f,%.2f" %(len(lon),len(lat),press3D.shape[0],lon[0],lat[0],lon[1]-lon[0],lat[1]-lat[0]))
    vars = vars + ',Winds,Tsurf,Psurf,Temperature,Pressure,H2O'
    if surfaces != None:
        vars = vars + ',Surf_Snow,Surf_Ocean'
    print(vars)
    templ['ATMOSPHERE-GCM-PARAMETERS'] = ''


	
    with open(fileout,'w') as fw:
        for i in templ: fw.write('<'+i+'>'+str(templ[i])+'\n')      
        fw.write(vars+'\n')
    # with open(fileout,'ab') as fb:
    with open(fileout,'ab') as fb:
        if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
        else: bc=fb.write('<BINARY>')
        fb.write(np.asarray(u,order='C',dtype=np.float32))
        fb.write(np.asarray(v,order='C',dtype=np.float32))
        fb.write(np.asarray(tsurf,order='C',dtype=np.float32))
        fb.write(np.log10(np.asarray(ps,order='C',dtype=np.float32)))
        fb.write(np.asarray(temp,order='C',dtype=np.float32))
        fb.write(np.log10(np.asarray(press3D,order='C',dtype=np.float32)))
        fb.write(np.log10(np.asarray(h2o,order='C',dtype=np.float32)))
        
        if surfaces != None:   #currently the options are None or snow and ocean
            fb.write(np.asarray(snow,order='C',dtype=np.float32))
            fb.write(np.asarray(ocean,order='C',dtype=np.float32))

        if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
        else: bc=fb.write('</BINARY>')
        fb.close()
    return(templ)
    
data = convert('earth',surfaces=['snow','ocean'])      


