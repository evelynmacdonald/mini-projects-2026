import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth_cold',modelname='kearth_cold')
kearth.configure(aquaplanet=True,startemp=4500,flux=800)
kearth.exportcfg()
kearth.runtobalance()



