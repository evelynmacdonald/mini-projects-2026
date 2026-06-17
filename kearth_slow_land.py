import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth_slow_land',modelname='kearth_slow_land')
kearth.configure(desertplanet=True,startemp=4500,rotationperiod=10)
kearth.exportcfg()
kearth.runtobalance()



