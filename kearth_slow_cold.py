import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth_slow_cold',modelname='kearth_slow_cold')
kearth.configure(aquaplanet=True,startemp=4500,rotationperiod=10,flux=800)
kearth.exportcfg()
kearth.runtobalance()



