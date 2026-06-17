import exoplasim as exo
kearth = exo.TLaquaplanet(ncpus=8,workdir='kearth_tl',modelname='kearth_tl')
kearth.configure(startemp=4500,rotationperiod=120)
kearth.exportcfg()
kearth.runtobalance()



