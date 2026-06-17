import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth_land',modelname='kearth_land')
kearth.configure(desertplanet=True,startemp=4500)
kearth.exportcfg()
kearth.runtobalance()



