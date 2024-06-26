import sys, os
sys.path.insert(0,os.path.join('../'))
import gwcatpy as gwcatpy
import json
import ciecplib

# print(ciecplib.get("https://ldas-jobs.ligo.caltech.edu/~duncan.macleod/hello.html"))

verbose=True
forceupdate=True
devMode=False
useLocal=False
update=True

dataDir='data/'
fileIn=os.path.join(dataDir,'gwosc_gracedb_empty.json')
if not os.path.isfile(fileIn):
    fileIn='gwosc_gracedb_empty.json'
    print('Reading from empty file: {}'.format(fileIn))
else:
    print('Reading from file: {}'.format(fileIn))

if devMode:
    devStr='dev'
    sess=ciecplib.Session("LIGO")
else:
    devStr='public'
    sess=None

gc=gwcatpy.GWCat(fileIn,dataDir=dataDir,verbose=verbose,mode=devStr)

if update:
    knownEvents=gc.getTimestamps()
    
    # print('\n\n*****\nReading GWTC...\n*****\n\n')
    gwtcdata=gwcatpy.gwosc.getGWTC(useLocal=useLocal,export=True,dirOut=dataDir,verbose=verbose,devMode=devMode,catalog='O4_Discovery_Papers',sess=sess)
    print('\n\n*****\nImporting O4_Discovery_Papers...\n*****\n\n')
    gc.importGWTC(gwtcdata,verbose=verbose, devMode=devMode)
    
    # gwtcm1data=gwcatpy.gwosc.getGWTC(useLocal=useLocal,export=True,dirOut=dataDir,verbose=verbose,devMode=devMode,catalog='GWTC-1-marginal',sess=sess)
    # print('\n\n*****\nImporting GWTC...\n*****\n\n')
    # gc.importGWTC(gwtcm1data,verbose=verbose, devMode=devMode)
    # 
    # gwtcm1data=gwcatpy.gwosc.getGWTC(useLocal=useLocal,export=True,dirOut=dataDir,verbose=verbose,devMode=devMode,catalog='GWTC-1-marginal',sess=sess)
    # print('\n\n*****\nImporting GWTC...\n*****\n\n')
    # gc.importGWTC(gwtcm1data,verbose=verbose, devMode=devMode)
    
    # print('\n\n*****\nImporting O3 Discovery Papers...\n*****\n\n')
    # o3discdata=gwcatpy.gwosc.getGWTC(useLocal=useLocal,export=True,dirOut=dataDir,verbose=verbose,catalog='O3_Discovery_Papers',sess=sess,devMode=devMode)
    # # o3discdata=gwcatpy.gwosc.getGWTC(useLocal=useLocal,export=True,dirOut=dataDir,verbose=verbose,catalog='O3_Discovery_Papers',sess=sess,url='data/local-mirror/O3_Discovery_Papers.json')
    # print('\n\n*****\nImporting O3 Discovery Papers...\n*****\n\n')
    # gc.importGWTC(o3discdata,verbose=verbose, devMode=devMode)
    
    # print('\n\n*****\nReading GraceDB...\n*****\n\n')
    # # gdb=json.load(open(os.path.join(dataDir,'gracedb.json')))
    # gdb=gwcatpy.gracedb.getSuperevents(export=True,dirOut=dataDir,verbose=verbose,knownEvents=knownEvents,forceUpdate=forceupdate)
    # 
    # 
    # print('\n\n*****\nimporting GraceDB...\n*****\n\n')
    # gc.importGraceDB(gdb,verbose=verbose,forceUpdate=forceupdate)
    
    # print('\n\n*****\nmatching GraceDB entries...\n*****\n\n')
    # gc.matchGraceDB(verbose=verbose)
    # 
    # print('\n\n*****\nremoving unnecessary GraceDB candidates\n*****\n\n')
    # gc.removeCandidates(verbose=verbose)
    
    # print('\n\n*****\nAdding manual references...\n*****\n\n')
    # gc.addRefs(verbose=verbose)

print('\n\n*****\nUpdating data from H5\n*****\n\n')
gc.updateH5(verbose=verbose,forceUpdate=forceupdate)

print('\n\n*****\nsetting precision...\n*****\n\n')
gc.setPrecision(extraprec=1,verbose=verbose)

print('\n\n*****\nUpdating maps\n*****\n\n')
gc.updateMaps(verbose=verbose,forceUpdate=forceupdate)

# print('\n\n*****\nPlotting maps\n*****\n\n')
# gc.plotMapPngs(verbose=verbose)
# 
# print('\n\n*****\nUpdating gravoscope\n*****\n\n')
# tilesurl='https://ligo.gravity.cf.ac.uk/~chris.north/LVC/gwcatpydev/'
# gc.makeGravoscopeTiles(verbose=verbose,maxres=6,tilesurl=tilesurl)

gc.exportJson(os.path.join(dataDir,'gwosc_gracedb_O4.json'))

gwcatpy.json2jsonp(os.path.join(dataDir,'gwosc_gracedb_O4.json'),os.path.join(dataDir,'gwosc_gracedb_O4.jsonp'))