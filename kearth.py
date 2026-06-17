import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth',modelname='kearth')
kearth.configure(aquaplanet=True,startemp=4500)
kearth.exportcfg()
kearth.runtobalance()



