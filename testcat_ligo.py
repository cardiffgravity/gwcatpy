import sys, os
sys.path.insert(0,os.path.join('../'))
import gwcatpydev as gwcatpy
import json
# import ciecplib

# print(ciecplib.get("https://ldas-jobs.ligo.caltech.edu/~duncan.macleod/hello.html"))

verbose=True
forceupdate=False
devMode=True
useLocal=False
update=False

dataDir='data/'
fileIn=os.path.join(dataDir,'gwosc_gracedb_GWTC2.json')
if not os.path.isfile(fileIn):
    fileIn='gwosc_gracedb_empty.json'
    print('Reading from empty file: {}'.format(fileIn))
else:
    print('Reading from file: {}'.format(fileIn))

gc=gwcatpy.GWCat(fileIn,dataDir=dataDir,verbose=verbose,mode='dev')

if update:
    print('\n\n*****\nReading GWTC...\n*****\n\n')
    gwtcdata=gwcatpy.gwosc.getGWTC(useLocal=useLocal,export=True,dirOut=dataDir,verbose=verbose,devMode=devMode)
    knownEvents=gc.getTimestamps()

    print('\n\n*****\nReading GraceDB...\n*****\n\n')
    # gdb=json.load(open(os.path.join(dataDir,'gracedb.json')))
    gdb=gwcatpy.gracedb.getSuperevents(export=True,dirOut=dataDir,verbose=verbose,knownEvents=knownEvents,forceUpdate=forceupdate)
                
    print('\n\n*****\nImporting GWTC...\n*****\n\n')
    gc.importGWTC(gwtcdata,verbose=verbose, devMode=devMode)
    # print('\n\n*****\nimporting GraceDB...\n*****\n\n')
    # gc.importGraceDB(gdb,verbose=verbose,forceUpdate=forceupdate)
    # print('\n\n*****\nmatching GraceDB entries...\n*****\n\n')
    gc.matchGraceDB(verbose=verbose)
    print('\n\n*****\nremoving unnecessary GraceDB candidates\n*****\n\n')
    gc.removeCandidates(verbose=verbose)
    print('\n\n*****\nAdding manual references...\n*****\n\n')
    gc.addRefs(verbose=verbose)

print('\n\n*****\nUpdating data from H5\n*****\n\n')
gc.updateH5(verbose=verbose,forceUpdate=forceupdate)

print('\n\n*****\nsetting precision...\n*****\n\n')
gc.setPrecision(extraprec=1,verbose=verbose)

print('\n\n*****\nUpdating maps\n*****\n\n')
gc.updateMaps(verbose=verbose,forceUpdate=forceupdate)

# print('\n\n*****\nPlotting maps\n*****\n\n')
# gc.plotMapPngs(verbose=verbose)
# 
print('\n\n*****\nUpdating gravoscope\n*****\n\n')
tilesurl='https://ligo.gravity.cf.ac.uk/~chris.north/LVC/gwcatpydev/'
gc.makeGravoscopeTiles(verbose=verbose,maxres=6,tilesurl=tilesurl)

gc.exportJson(os.path.join(dataDir,'gwosc_gracedb_GWTC2.json'))

gwcatpy.json2jsonp(os.path.join(dataDir,'gwosc_gracedb_GWTC2.json'),os.path.join(dataDir,'gwosc_gracedb_GWTC2.jsonp'))