import exoplasim as exo
kearth = exo.Earthlike(ncpus=8,workdir='kearth_continent_run',modelname='kearth_continent')
kearth.configure(landmap='SC0.30T21_surf_0172.sra',startemp=4500)
kearth.exportcfg()
kearth.runtobalance()
kearth.finalize('kearth_continent',allyears=False,keeprestarts=False)



