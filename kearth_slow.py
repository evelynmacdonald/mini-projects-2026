import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth_slow',modelname='kearth_slow')
kearth.configure(aquaplanet=True,startemp=4500,rotationperiod=10)
kearth.exportcfg()
kearth.runtobalance()



