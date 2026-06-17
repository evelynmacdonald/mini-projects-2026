import exoplasim as exo
earth = exo.Earthlike(ncpus=8,workdir='earth',modelname='earth')
earth.configure(aquaplanet=True)
earth.exportcfg()
earth.runtobalance()



