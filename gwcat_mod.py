import json
import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.io import fits
import healpy as hp
import h5py
import os
import requests
from . import gwcat_gwosc as gwosc
from . import gwcat_gracedb as gracedb
from . import plotloc
from astropy.table import Table
import astropy_healpix as ah
from pycbc.waveform import get_td_waveform
import ciecplib
# import gracedb
# import gwosc

def json2jsonp(fileIn,fileOut=None,verbose=False):
    """Read Json file and convert to Jsonp
    Inputs:
        * fileIn [string]: json filename
        * fileOut [string, optional]: jsonp filename. Default = input with json->jsonp
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * None
    """
    if fileOut==None:
        fileOut=fileIn.replace('json','jsonp')
    if verbose:print('reading JSON from {}'.format(fileIn))
    if verbose:print('writing JSONP to {}'.format(fileOut))
    fIn=open(fileIn,'r')
    fOut=open(fileOut,'w')
    lines=fIn.readlines()
    lines[0]='catdata('+lines[0]
    lines[-1]=lines[-1]+');'
    for l in lines:
        fOut.write(l)
    fOut.close()
    return

def compareElement(el1,el2,verbose=False):
    """Convert event from pandas dataframe format to json-style object
    Inputs:
        * el1 [object]: element 1 to compare
        * el1 [object]: element 2 to compare
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * Boolean or None:
            * True: Elements are the same
            * False: Elements are not the same
            * None: Elements can't be compared
    """

    if type(el1)!=type(el2):
        # types are different
        if verbose: print('inconsistent types [{},{}]'.format(type(el1),type(el2)))
        return(False)
    elif type(el1)==dict and type(el2)==dict:
        for el in el1:
            if el in el2:
                # if verbose: print('Comparing elements {}'.format(el))
                return(compareElement(el1[el],el2[el],verbose=verbose))
            else:
                # if verbose: print('element {} not in dict 2'.format(el))
                return(False)
        for el in el2:
            if el not in el1:
                # if verbose: print('element {} not in dict 1'.format(el))
                return(False)
    else:
        try:
            # if verbose: print('Comparing elements {}'.format(el1))
            if el1==el2:
                # if verbose: print('{}=={}'.format(el1,el2))
                return(True)
            else:
                # if verbose: print('{}!={}'.format(el1,el2))
                return(False)
        except:
            print('ERROR: unable to compare [{} , {}] / [{},{}]'.format(el1,el2),type(el1),type(el2))
            return()

    return()

def dataframe2jsonEvent(evIn,params,verbose=False):
    """Convert event from pandas dataframe format to json-style object
    Inputs:
        * dataIn [object]:pandas-style event record
        * params [dictionary]: units information (merge with datadict)
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * Json-style event object:
    """
    evOut={}
    for p in range(len(params)):
        param=params[p]
        if verbose: print('Converting value {}'.format(param))
        if param in evIn:
            try:
                pv=evIn[param+'_valtype']
            except:
                if verbose: print('ERROR: No valtype found for {}'.format(param))
                pv=''
                pass
            try:
                np.isnan(pv)
                pass
            except:
                evOut[param]={}
            if pv=='value':
                evOut[param]['best']=evIn[param]
            elif pv=='bestfit':
                evOut[param]['best']=evIn[param]
                try:
                    evOut[param]['err']=[evIn[param+'_errp'],evIn[param+'_errm']]
                except:
                    if verbose: print('WARNING: no err found form {}'.format(param))
            elif pv=='range':
                try:
                    evOut[param]['lim']=[evIn[param+'_errp'],evIn[param+'_errm']]
                except:
                    if verbose: print('WARNING: no lim found form {}'.format(param))
            elif pv=='upper':
                evOut[param]['upper']=evIn[param]
            elif pv=='lower':
                evOut[param]['lower']=evIn[param]
            else:
                if verbose: print('ERROR: unknown type for {} [{}]'.format(param,pv))
    return(evOut)

def getManual(loc='',verbose=True,export=False,
    dirOut=None,fileOut=None,indent=2):
    """Get data from manually constructed data-file
    Inputs:
        * loc [string, optional]: filename to load. Default='data/manual-events.json'
        * verbose [boolean, optional]: set for verbose output. Default=False
        * export [boolean, optional]: set to for export data to new file. Default=False
        * dirOut [string, optional]: directory to export data to (if export=-True). Default=../../data/
        * fileOut [string, optional]: file to export data to (if export=-True). Default=manual-events.json
        * indent [integer, optional]: json indent in exported file (if export=True). Default=2
    Outputs:
        * Json-style object containing data in datafile
    """
    if loc=='':
        loc='data/manual-events.json'
    if verbose: print('Retrieving Manual data from {}'.format(loc))
    manread=json.load(open(loc,'r'))
    mandata={'meta':{'retrieved':Time.now().isot,'src':loc}}
    for s in manread:
        mandata[s]=manread[s]

    if verbose: print('Retrieved data for {} events'.format(len(mandata['data'])))

    if export:
        if dirOut==None:
            dirOut='../../data/'
        if fileOut==None:
            fileOut='manual-events.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(mandata,fOut,indent=indent)
        fOut.close()
    return mandata

def setPrec(v,prec):
    """Set precision of variable
    Inputs:
        * v [float]: value to set precision of
        * prec [integer]: precision to use
    Outputs:
        * [float]: value at set precision
    """
    if prec:
        precstr='{:.'+'{}'.format(prec)+'g}'
    else:precstr='{}'
    return float(precstr.format(v))

class GWCat(object):
    """Cataloge object containing useful methods for editing and accessing
    """
    def __init__(self,fileIn='../data/events.json',statusFile='status.json',
        dataDir='data/',baseurl='https://data.cardiffgravity.org/gwcat-data/',dataurl='https://ligo.gravity.cf.ac.uk/~chris.north/gwcat-data/',verbose=False,mode='public'):
        """Initialise catalogue from input file and set basic parameters
        Inputs:
            * fileIn [string, optional]: filename to load from. Default=../data/events.json
            * statusFile [string, optional]: filename to use for data status. Default=status.json
            * dataDir [string, optional]: Directory to store data in. Default=data/
            * baseurl [string, optional]: URL to use for absolute URLs. Default=https://data.cardiffgravity.org/gwcat-data/
            * dataurl [string, optional]: URL to use for absolute URLs for data files. Default='https://ligo.gravity.cf.ac.uk/~chris.north/gwcat-data/
            * verbose [boolean, optional]: set for verbose output. Default=False
            * mode [string, optional]: mode. Default=public
        Outputs:
            * None
        """
        self.dataDir=dataDir
        if not os.path.isdir(dataDir):
            # make new directory
            os.mkdir(dataDir)
        if baseurl[-1]!='/':baseurl=baseurl+'/'
        self.baseurl=baseurl
        self.dataurl=dataurl
        self.statusFile=os.path.join(dataDir,statusFile)

        self.mode=mode
        if mode=='dev':
            self.devMode=True
            self.sess=ciecplib.Session("LIGO")
        else:
            self.devMode=False
            self.sess=None
        self.getStatus()
        eventsIn=json.load(open(fileIn))
        self.data=eventsIn['data']
        self.datadict=eventsIn['datadict']
        self.cols=list(self.datadict.keys())
        self.links=eventsIn['links']
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        if 'meta' in eventsIn:
            self.meta=eventsIn['meta']
            self.meta['created']=Time.now().isot
        else:
            self.meta={'created':Time.now().isot}
        return

    def getEvent(self,ev):
        """Get event from database by event name
        Inputs:
            * ev [string]: event
        Outputs:
            * [object] json-style event object
            * (None if no event not in database)
        """
        if ev in self.data:
            return(self.data[ev])
        else:
            return
    def getParameter(self,ev,param):
        """Get parameter for event from database
        Inputs:
            * ev [string]: event name
            * param [string]: parameter name
        Outputs:
            * [object] json-style parameter object
            * (None if no event not in database)
        """
        pOut=None
        if ev in self.data:
            if param in self.data[ev]:
                pOut=self.data[ev][param]
        return(pOut)
    def getTimestamps(self):
        """Get timestamps of all events in database
        Inputs:
            * None
        Outputs:
            * [object] list of timestamps (of creation) in database, listed by event name
        """
        evTimes={}
        for ev in self.data:
            try:
                evTimes[ev]=Time(self.data[ev]['meta']['created_date'])
            except:
                # print('no created_date for {}'.format(ev))
                pass
        return(evTimes)

    def getStatus(self):
        """Load status from status file (set by self.statusFile). Initialise file if not present. Store in self.status.
        Inputs:
            * None
        Outputs:
            * None
        """
        try:
            self.status=json.load(open(self.statusFile))
        except:
            print("Initialising status")
            self.status={}
            json.dump(self.status,open(self.statusFile,'w'))
        return

    def saveStatus(self):
        """Save status (in self.status) for status file (set by self.statusFile).
        Inputs:
            * None
        Outputs:
            * None
        """
        json.dump(self.status,open(self.statusFile,'w'),indent=4)
        return

    def updateStatus(self,ev,desc='',statusIn=None,verbose=False):
        """Update status for event and save to status file
        Inputs:
            * ev [string]: event name
            * desc [string, optional]: Description to use in verbose output
            * statusIn [opject, optional]: status values to add for event. Default: load from metadata of event
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        # check if it's in status
        if type(ev)!=str:
            print('ERROR: ev must be string object:',type(ev))
            return
        if statusIn==None:
            # get status from event data
            try:
                event=self.getEvent(ev)
                status=event['meta']
            except:
                print('ERROR: cannot find event {}'.format(ev))
                return
        elif type(statusIn)!=dict:
            print('ERROR: status must be dict object:',type(status))
            return
        else:
            if 'meta' in statusIn:
                status=statusIn['meta']
            else:
                status=statusIn
        if not ev in self.status:
            self.status[ev]={}
        for m in status:
            self.status[ev][m]=status[m]
        # if (event):
        #     if not ev in self.status:
        #         self.status[ev]={}
        #     if 'meta' in status:
        #         for m in event['meta']:
        #             self.status[ev][m]=event['meta'][m]
        if desc!='':desctxt='[{}]'.format(desc)
        if verbose:print('update status: {} {}'.format(ev,desctxt))
        self.saveStatus()
        return

    def updateMapSrc(self,ev,verbose=False):
        """Update map source fits file and created date, and store in metadata. Load header and extract distance if possible. Called during importManual and importGraceDB.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        lmap=self.getLink(ev,'skymap-fits',verbose=verbose)
        if len(lmap)>1 and verbose:
            print('Warning: more than one skymap link for {}'.format(ev))
        if len(lmap)==0 and verbose:
            print('Warning: no skymap link for {}'.format(ev))
        if len(lmap)>0:
            l=lmap[0]
            tmpDir=os.path.join(self.dataDir,'tmp')
            try:
                hdr=fits.getheader(l['url'],ext=1)
                self.data[ev]['meta']['mapdatesrc']=hdr['DATE']
                self.data[ev]['meta']['mapurlsrc']=l['url']
                self.updateStatus(ev,verbose=verbose,desc='Map src')
                if verbose:print('Skymap loaded for {}:'.format(ev),l['url'])
                try:
                    distmu=hdr['DISTMEAN']
                    distsd=hdr['DISTSTD']
                    self.data[ev]['DL']={'best':distmu,'err':[-distsd,distsd]}
                    if verbose:print('obtained distance information for {} [{:.2f},{:.2f}]'.format(ev,distmu,distsd))
                except:
                    if verbose:print('no distance information for {}'.format(ev))
            except:
                print('WARNING: Error loading skymap for {}:'.format(ev),l)
        return

    def importManual(self,manIn,verbose=False):
        """import Manual data into catalogue. Overwrites any existing data for events if they are already present
        Inputs:
            * manIn [object]: event data
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        print('*** Importing Manual Data...')
        for g in manIn['data']:
            # get old metadata
            dmeta={}
            if g in self.data:
                if 'meta' in self.data[g]:
                    dmeta=self.data[g]['meta']
            self.data[g]=manIn['data'][g]
            # update metadata
            for m in manIn['data'][g]['meta']:
                dmeta[m]=manIn['data'][g]['meta'][m]
            self.data[g]['meta']=dmeta
            if g in manIn['links']:
                for l in manIn['links'][g]:
                    self.addLink(g,l,verbose=verbose)
            self.updateStatus(g,verbose=verbose,desc='Manual import')
        for ev in manIn['links']:
            self.updateMapSrc(ev,verbose=True)
            self.updateStatus(ev,verbose=verbose,desc='Map src')
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        if not 'manuall' in self.meta:
            self.meta['manual']={}
        for m in manIn['meta']:
            self.meta['manual'][m]=manIn['meta'][m]
        return

    def setPrecision(self,extraprec=3,verbose=False):
        """set precision for best, err, lower, upper for all events, based on precision in self.datadict
        Inputs:
            * extraprec [integer, optional]: additional precision to store. Default=3
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        from numbers import Number
        for ev in self.data:
            for p in self.data[ev]:
                newParam=self.getParameter(ev,p)
                try:
                    assert(newParam)
                    assert p in self.datadict
                    assert 'sigfig' in self.datadict[p]
                except:
                    # if verbose:print('unable to set precision for {}[{}]'.format(ev,p))
                    continue
                sigfig=self.datadict[p]['sigfig']
                if isinstance(newParam,Number) and not isinstance(newParam,bool):
                    newParam=self.setPrec(newParam,sigfig+extraprec)
                else:
                    if 'best' in newParam:
                        if isinstance(newParam['best'],Number) and not isinstance(newParam['best'],bool):
                            newbest=setPrec(newParam['best'],sigfig+extraprec)
                            if verbose:
                                if self.getParameter(ev,p)['best']!=newbest:
                                    print('Set precision {}[{}][best]: {}->{}'.format(ev,p,self.getParameter(ev,p)['best'],newbest))
                            newParam['best']=newbest
                    if 'err' in newParam:
                        if newbest!=0:
                            bprec=np.floor(np.log10(np.abs(newbest)))+1-(sigfig+extraprec)
                        else:
                            bprec=sigfig
                        mult=10**(-bprec)
                        newlow=np.round(newParam['err'][0]*mult)/mult
                        newhigh=np.round(newParam['err'][1]*mult)/mult
                        newerr=[newlow,newhigh]
                        # for e in range(len(newParam['err'])):
                        #     newerr.append(setPrec(newParam['err'][e],sigfig))
                        if verbose:
                            if self.getParameter(ev,p)['err'][0]!=newlow or self.getParameter(ev,p)['err'][1]!=newhigh:
                                print('Set precision {}[{}][err]: {}->{}'.format(ev,p,self.getParameter(ev,p)['err'],newerr))
                        newParam['err']=newerr
                    if 'lower' in newParam:
                        newParam['lower']=setPrec(newParam['lower'],sigfig+extraprec)
                    if 'upper' in newParam:
                        newParam['upper']=setPrec(newParam['upper'],sigfig+extraprec)
                self.data[ev][p]=newParam
        return

    def importGWTC(self,gwtcIn,verbose=False,devMode=False,forceOverwrite=False,catalog='GWTC'):
        """import GWTC data into catalogue. Uses gwosc.gwtc_to_cat to convert.
        Updates H5 source file.
        Convert to dataframe
        Inputs:
            * gwtcIn [object]: event data
            * verbose [boolean, optional]: set for verbose output. Default=False
            * devMode [boolean, optional]: set to fix DCC link for dev mode. Default=False.
            * forceOverwrite [boolean, optional: set to force overwrite of events with earlier version numbers
            * catalog [string, optional]: catalog to read from. Default=GWTC.
        Outputs:
            * None
        """
        print('*** Importing GWTC Catalog...')
        catData=gwosc.gwtc_to_cat(gwtcIn,self.datadict,verbose=verbose,devMode=self.devMode,catalog=catalog)
        for ev in catData['data']:
            # get old metadata
            dmeta={}
            if ev in self.data:
                verEx=self.getParameter(ev,'version')
                verNew=catData['data'][ev]['version']
                if (verNew) and (verEx):
                    if (verNew <= verEx):
                        if forceOverwrite:
                            print('Forcing replacement of  {} v{} with v{}'.format(ev,verEx,verNew))
                        else:
                            if verbose:
                                print('Not replacing {} v{} with v{}'.format(ev,verEx,verNew))
                            continue
                    else:
                        if verbose:
                            print('Updating {} v{} with v{}'.format(ev,verEx,verNew))
                if 'meta' in self.data[ev]:
                    dmeta=self.data[ev]['meta']
            self.data[ev]=catData['data'][ev]
            # update metadata
            for m in catData['data'][ev]['meta']:
                dmeta[m]=catData['data'][ev]['meta'][m]
            self.data[ev]['meta']=dmeta
            if ev in catData['links']:
                for l in catData['links'][ev]:
                    self.addLink(ev,l,verbose=verbose)
            self.updateStatus(ev,verbose=verbose,desc='{} Catalog import'.format(catalog))
        for ev in catData['links']:
            self.updateH5Src(ev,verbose=verbose)
            self.updateStatus(ev,verbose=verbose,desc='Data file src')
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        if not catalog in self.meta:
            self.meta[catalog]={}
        for m in gwtcIn['meta']:
            self.meta[catalog][m]=gwtcIn['meta'][m]

        return

    def matchGraceDB(self,verbose=False):
        """match confirmed events with GraceDB events (assuming <1s GPS timestamp difference)
        Inputs:
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        gdblist=[]
        gpslist=[]
        for ev in self.data:
            name=self.data[ev]['name']
            gps=self.getParameter(ev,'GPS')
            if name[0]=='S':
                try:
                    gps['best']
                except:
                    continue
                gdblist.append(name)
                gpslist.append(gps['best'])
        gpslist=np.array(gpslist)
        gdblist=np.array(gdblist)
        for ev in self.data:
            if self.data[ev]['name'][0]=='S':
                # is a gracedb event itself
                continue
            tgps=self.getParameter(ev,'GPS')
            try:
                tgps['best']
            except:
                continue
            idx=np.where(np.abs(gpslist-tgps['best'])<1)[0]
            if len(idx)>0:
                matchname=gdblist[idx][0]
                if verbose:
                    print('matching {} with {}'.format(ev,matchname))
                self.data[ev]['meta']['gracedb']=matchname
                # maplink=self.getLink(ev,'skymap-fits')
                # mapgdb=self.getLink(matchname,'skymap-fits')
                # if len(maplink)==0 and len(mapgdb)>0:
                #     if verbose:print('copying GraceDB map link from {} to {}: {}'.format(matchname,ev,mapgdb[0]))
                #     self.addLink(ev,mapgdb[0],verbose=verbose)

        return

    # backwards compatibility
    def importGwosc(self,gwoscIn,verbose=False):
        """OBSOLETE - inpurt passed to importGWTC
        Inputs:
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        print('***WARNING: importGwosc replaced by importGWTC***')
        self.importGWTC(gwoscIn,verbose=verbose)

    def importGraceDB(self,gracedbIn,verbose=False,forceUpdate=False,highSigOnly=False):
        """import GraceDB data into catalogue. Uses gracedb.gracedb2cat to convert.
        Remove retracted events.
        Remove low-significant events
        Update map source file.
        Export to dataframe
        Inputs:
            * gracedbIn [object]: event data
            * verbose [boolean, optional]: set for verbose output. Default=False
            * forceUpdate [boolean, optional]: set to force updates of gracedb events. Default=False (only update new events)
            * highSigOnly [boolean, optional]: set to remove low significance events. Default=False (include low significance events)
        Outputs:
            * None
        """
        print('*** Importing GraceDB...')
        evTimes=self.getTimestamps()
        gdb=gracedb.gracedb2cat(gracedbIn['data'],verbose=verbose,
            knownEvents=evTimes,forceUpdate=forceUpdate)
        for g in gdb['data']:
            if gdb['data'][g]['meta'].get('type')=='Retraction':
                if verbose:print('Skipping retracted event {}'.format(g))
                if g in self.data:
                    print('Removing data for retracted event {}'.format(g))
                    self.data.pop(g)
                if g in self.links:
                    print('Removing links for retracted event {}'.format(g))
                    self.links.pop(g)
                continue
            if highSigOnly:
                if gdb['data'][g].get('Significance')=='Low':
                    if g in self.data:
                        print('Removing data for low-significance event {}'.format(g))
                        self.data.pop(g)
                    if g in self.data:
                        print('Removing links for low-significance event {}'.format(g))
                        self.links.pop(g)
                    print('Skipping low-significance event {}'.format(g))
                    continue
            # get old metadata
            dmeta={}
            if g in self.data:
                if 'meta' in self.data[g]:
                    dmeta=self.data[g]['meta']
            self.data[g]=gdb['data'][g]
            # update metadata
            for m in gdb['data'][g]['meta']:
                dmeta[m]=gdb['data'][g]['meta'][m]
            self.data[g]['meta']=dmeta
            for l in gdb['links'][g]:
                self.addLink(g,l,verbose=verbose)
            self.updateStatus(g,verbose=verbose,desc='GraceDB import')
        for ev in gdb['links']:
            if gdb['data'][ev]['meta'].get('type')=='Retraction':
                continue
            self.updateMapSrc(ev,verbose=verbose)
            self.updateStatus(ev,verbose=verbose,desc='Map src')
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        if not 'graceDB' in self.meta:
            self.meta['graceDB']={}
        for m in gracedbIn['meta']:
            self.meta['graceDB'][m]=gracedbIn['meta'][m]
        return

    def updateH5Src(self,ev,verbose=False):
        """Update H5 source fits file, and store in metadata. Called during importGWTC.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        ldat=self.getLink(ev,'data-file',verbose=verbose)
        if len(ldat)>1 and verbose:
            print('Warning: more than one data file link for {} [{}]'.format(ev,self.data[ev]['catalog']))
        if len(ldat)==0 and verbose:
            print('Warning: no data file link for {} [{}]'.format(ev,self.data[ev]['catalog']))
        if len(ldat)>0:
            l=ldat[0]
            tmpDir=os.path.join(self.dataDir,'tmp')
            try:
                # self.data[ev]['meta']['h5datesrc']=Time.now().isot
                self.data[ev]['meta']['h5urlsrc']=l['url']
                if 'zenodo_version' in self.data[ev]['meta']:
                    stat={'h5zenodosrc':self.data[ev]['meta']['zenodo_version']}
                    self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='H5 zenodo src')
                if verbose:print('Data file loaded for {}:'.format(ev),l['url'])
            except:
                print('WARNING: Error loading data file for {}:'.format(ev),l)

    def updateH5(self,verbose=False,forceUpdate=False,forceUpdateData=False,replace=False,event=None):
        """Check whether H5 files need updating and re-download if necessary (using getH5 or getH5Local).
        Update parameters from H5 files.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
            * replace [boolean, optional]: set to replace existing parameters with H5. Default=False
            * forceUpdate [boolean, optional]: set to download all files. Default=False (only download updated files)
            * forceUpdateData [boolean, optional]: set to force update of all data. Default=False (only download updated files)
        Outputs:
            * None
        """
        # update h5 data files and parameters
        print('*** Updating data files...')

        for ev in self.data:
            if event!="" and event!=None:
                if ev!=event:
                    continue
                else:
                    if verbose: print("***PROCESSING H5 FOR SINGLE EVENT {}".format(ev))
            try:
                if not ev in self.status:
                    self.status[ev]={}
                # compare map creation dates and versions
                if 'h5datelocal' in self.status[ev] and 'h5datesrc' in self.status[ev]:
                    h5datelocal=Time(self.status[ev]['h5datelocal']).gps
                    h5datesrc=Time(self.status[ev]['h5datesrc']).gps
                else:
                    h5datelocal=-1
                    h5datesrc=-1
                if 'h5verlocal' in self.status[ev] and 'h5versrc' in self.status[ev]:
                    h5verlocal=self.status[ev]['h5verlocal']
                    h5versrc=self.status[ev]['h5versrc']
                else:
                    h5verlocal=-1
                    h5versrc=-1
                if 'h5zenodolocal' in self.status[ev] and 'h5zenodosrc' in self.status[ev]:
                    h5zenodolocal=self.status[ev]['h5zenodolocal']
                    h5zenodosrc=self.status[ev]['h5zenodosrc']
                else:
                    h5zenodolocal=-1
                    h5zenodosrc=-1
                h5file=self.status[ev].get('h5urllocal','')
                if not os.path.isfile(h5file):
                    if verbose:print('No file. Need to re-download data file for {}'.format(ev))
                    updateH5=True
                elif h5zenodosrc>0 and h5zenodolocal>0:
                    if h5zenodosrc>h5zenodolocal:
                        if verbose:print('Newer zenodo version [{}>{}]. Need to re-download data file for {}'.format(h5zenodosrc,h5zenodolocal,ev))
                        updateH5=True
                    else:
                        if verbose:print('Older zenodo version file. Do not need to re-download data file for {}'.format(ev))
                        updateH5=False
                elif h5versrc>0 and h5verlocal>0:
                    if h5versrc>h5verlocal:
                        if verbose:print('Newer version [{}>{}]. Need to re-download data filename for {}'.format(h5versrc,h5verlocal,ev))
                        updateH5=True
                    else:
                        if verbose:print('Older version file. Do not need to re-download data file for {}'.format(ev))
                        updateH5=False
                elif h5datesrc>0 and h5datelocal>0:
                    if (h5datesrc-h5datelocal)>1:
                        if verbose:
                            print('Newer file. Need to re-download data file for {}'.format(ev))
                            print('src',Time(h5datesrc,format='gps').isot)
                            print('local',Time(h5datelocal,format='gps').isot)
                        updateH5=True
                    else:
                        if verbose:print('Older file. Do not need to re-download data file for {}'.format(ev))
                        updateH5=False
                else:
                    if verbose:print('Unable to determine H5 file date/version for {}. Re-downloading'.format(ev))
                    updateH5=True
            except:
                if verbose:print('Error determining H5 file date/version for {}. Re-downloading'.format(ev))
                updateH5=True
            if updateH5 or forceUpdate:
                if verbose:print('Updating data file for {}'.format(ev))
                stat=self.getH5(ev,verbose=verbose)
                if stat!=0:
                    self.getH5Local(ev,verbose=verbose)
            else:
                if verbose:print('No data file update required for {}'.format(ev))

            # update data (if required)
            updateH5data=False
            try:
                if 'h5datelocal' in self.status[ev] and 'h5dateloaded' in self.status[ev]:
                    h5datelocal=Time(self.status[ev]['h5datelocal']).gps
                    h5dateloaded=Time(self.status[ev]['h5dateloaded']).gps
                else:
                    h5datelocal=-1
                    h5dateloaded=-1
                if 'h5verlocal' in self.status[ev] and 'h5verloaded' in self.status[ev]:
                    h5verlocal=self.status[ev]['h5verlocal']
                    h5verloaded=self.status[ev]['h5verloaded']
                else:
                    h5verlocal=-1
                    h5verloaded=-1
                if h5datelocal>0 and h5dateloaded>0:
                    if (h5datelocal-h5dateloaded)>1:
                        if verbose:
                            print('File newer than data. Need to re-load data file for {}'.format(ev))
                            print('local',Time(h5datelocal,format='gps').isot)
                            print('loaded',Time(h5dateloaded,format='gps').isot)
                        updateH5data=True
                    else:
                        if verbose:
                            print('No need to update h5 data for {}'.format(ev))
                        updateH5data=False
                elif h5verlocal>0 and h5verloaded>0:
                    if h5verlocal>h5verloaded:
                        if verbose:
                            print('File later version than data [{}>{}]. Need to re-load data file for {}'.format(ev))
                        updateH5data=True
                    else:
                        if verbose:
                            print('No need to update h5 data for {}'.format(ev))
                        updateH5data=False
                else:
                    if verbose: print('Unable to determine H5 data date/version for {}. Re-loading'.format(ev))
                    updateH5data=True
            except:
                if verbose: print('Error determining H5 data date/version for {}. Re-loading'.format(ev))
                updateH5data=True

            if updateH5data or forceUpdateData:
                try:
                    h5File=self.status[ev]['h5urllocal']
                except:
                    if verbose:print('no local data file for {}'.format(ev))
                    continue
                if 'approximant' in self.data[ev]:
                    approx=self.data[ev]['approximant']
                    if verbose:
                        print('approximant found for {}: {}'.format(ev,self.data[ev]['approximant']))
                else:
                    approx=None
                    if verbose:
                        print('no approximant found for {}'.format(ev))
                newparams=self.getH5Params(ev,approx=approx,verbose=verbose)
                if (newparams):
                    for p in newparams:
                        if (p in self.data[ev]):
                            if replace:
                                self.data[ev][p]=newparams[p]
                                if verbose:
                                    print('replacing {}[{}]:{}->{}'.format(ev,p,self.data[ev][p],newparams[p]))
                            else:
                                if verbose:
                                    print('NOT replacing {}[{}]:{}->{}'.format(ev,p,self.data[ev][p],newparams[p]))
                        else:
                            self.data[ev][p]=newparams[p]
                            print('adding {}[{}]:{}'.format(ev,p,newparams[p]))
                    stat={'h5dateloaded':Time.now().isot,
                        'h5verloaded':self.data[ev]["version"]}
                    if 'approximant' in newparams:
                        stat['approximant']=newparams['approximant']
                    self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='h5 data loaded')
                else:
                    if verbose: print('no parameters from H5 for {}'.format(ev))
            else:
                if verbose:print('No H5 data update required for {}'.format(ev))
        return

    def getH5(self,ev,verbose=False):
        """Download h5 files and update status with filenames, timestamps.
        Call extractTar if tarfile.
        Check file is valid H5 file.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * [object] dict containing extracted "h5" and "map" file locations.
        """
        # download and import remote h5 file into database
        import h5py
        h5Dir=os.path.join(self.dataDir,'h5')
        if not os.path.exists(h5Dir):
            # create directory
            os.mkdir(h5Dir)
            print('Created directory: {}'.format(h5Dir))
        ldat=self.getLink(ev,'data-file')
        if len(ldat)==0:
            if verbose: print('ERROR: no data file link for {}'.format(ev))
            return
        if len(ldat)>1:
            if verbose: print('WARNING: more than one data file link for {}'.format(ev))
        url=ldat[0]['url']
        if verbose: print('Downloading data file for {} from {}'.format(ev,url))
        if url.find("zenodo")>=0 and url.find("content")>=0:
            try:
                # New Zenodo format file
                zenurl=url.replace("/content","")
                srcfile=os.path.split(zenurl)[-1]
                if verbose: print('Zenodo src file:',zenurl,os.path.split(zenurl),srcfile)
            except:
                srcfile=os.path.split(url)[-1]
        else:
            srcfile=os.path.split(url)[-1]
        
        if self.devMode:
            h5req=self.sess.get(url)
        else:
            h5req=requests.get(url)
        if h5req.ok:
            try:
                # if url.find('.fits.gz')>=0:
                #     fileext='fits.gz'
                # else:
                #     fileext='fits'
                if srcfile.find(ev)<0:
                    h5File=os.path.join(h5Dir,'{}_{}'.format(ev,srcfile))
                    if h5File.find('.h')<0:
                        h5File=h5File+'.h5'
                else:
                    h5File=os.path.join(h5Dir,srcfile)
                fOut=open(h5File,'wb')
                fOut.write(h5req.content)
                fOut.close()
                if verbose:print('Data file downloaded')
            except:
                print('ERROR: Problem loading/saving data file:',h5req.status_code)
                return h5req
            if ldat[0]['filetype']=='tar':
                print('NEED TO EXTRACT TARBALL')
                stat={'tarurllocal':h5File,
                    'tarurlsrc':url,
                    'tardatelocal':Time.now().isot,
                    'tardatesrc':Time(os.path.getmtime(h5File),format='unix').isot,
                    'tarverlocal':self.data[ev]["version"],
                    'tarversrc':self.data[ev]["version"]}
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='tar file local')
                tarout=self.extractTar(ev,verbose=verbose)
                h5File=tarout['h5']
            try:
                testread=h5py.File(h5File)
                if verbose: print('Valid h5 file for {}: {}'.format(ev,h5File))
                # hdr=fits.getheader(fitsFile,ext=1)
                stat={'h5urllocal':h5File,
                    'h5datelocal':Time.now().isot,
                    'h5datesrc':Time(os.path.getmtime(h5File),format='unix').isot,
                    'h5verlocal':self.data[ev]["version"],
                    'h5versrc':self.data[ev]["version"]}
                if 'zenodo_version' in self.data[ev]['meta']:
                    stat['h5zenodolocal']=self.data[ev]['meta']['zenodo_version']
                    stat['h5zenodosrc']=self.data[ev]['meta']['zenodo_version']
                self.data[ev]['meta']['h5urllocal']=h5File
                self.data[ev]['meta']['h5datelocal']=Time.now().isot
                self.data[ev]['meta']['h5datesrc']=Time.now().isot
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='data file local')
            except:
                print('ERROR: Problem opening data file for {}:'.format(ev),h5File)
                return(-1)
            # self.calcAreas(ev,verbose=verbose)
        else:
            print('ERROR: Problem loading h5 file:',h5req.status_code)
            return h5req.status_code
        return(0)

    def extractTar(self,ev,statusid='tarurllocal',maponly=False,verbose=False):
        """Extract h5 and map file from tarfile
        Inputs:
            * ev [string]: event name
            * statusid [string, optional]: id in status file to identify tar file
            * maponly [boolean, optional]: set to only extract the map file. Default=False
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * [number] HTTP status code if unsuccessful. 0 if not.
        """
        import tarfile
        tarFile=self.status[ev][statusid]
        tarF=tarfile.open(tarFile)
        tarNames=tarF.getnames()
        # get first file with _comoving.h5 - for GWTC-2
        h5name=[n for n in tarNames if '_comoving.h5' in n]
        # h5name=[n for n in tarNames if '.h5' in n]
        if len(h5name)>0:
            h5name=h5name[0]
        mapname=[]
        mapfound=False
        if 'approximant' in self.status[ev]:
            tarapprox=self.status[ev]['approximant'].replace('C00','').replace('C01','')
        else:
            tarapprox='###NONE###'
        for n in tarNames:
            # GWTC-2 filename
            if n.find('PublicationSamples.fits')>=0:
                # for GWTC-2
                mapname.append(n)
            # GWTC-2.1 and GWTC-4.0 filename
            elif n.find(ev)>=0 and n.find(tarapprox)>=0:
                    mapname.append(n)
                    mapfound=True
            # GWTC-3 filename
            elif n.find(ev)>=0 and n.find('Mixed.fits')>=0:
                mapname.append(n)
                mapfound=True
            elif n.find(ev)>=0 and mapfound==False:
                # lsit all, and later find the first one it comes to!
                mapname.append(n)
            # catch for GWTC-3 zenodo filename error
            if n.find(ev.replace('GW200210_092254','GW200210_092255'))>=0 and n.find('Mixed.fits')>=0:
                mapname.append(n)

        # else:
        #     mapname=[n for n in tarNames if 'PublicationSamples.fits' in n]
        if verbose:print('map for {}:{}'.format(ev,mapname))
        out={'h5':'','map':''}

        # extract h5 to h5 directory tarfile
        if not maponly:
            h5Dir=os.path.join(self.dataDir,'h5')
            try:
                tarF.extract(h5name,path=h5Dir)
                h5File=os.path.join(h5Dir,h5name)
                if (verbose): print('Extracted h5 file to {}'.format(h5File))
                # save location of h5 file to status
                stat={'h5urllocal':h5File,
                    'h5datelocal':Time.now().isot,
                    'h5verlocal':self.data[ev]["version"],
                    'h5versrc':self.data[ev]["version"],
                    'h5urlsrc':self.status[ev]["tarurlsrc"]}
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='h5 file local')
                out['h5']=h5File
            except:
                pass

        # extract fits to fits directory
        if len(mapname)>0:
            fitsDir=os.path.join(self.dataDir,'fits')
            if not os.path.exists(fitsDir):
                # create directory
                os.mkdir(fitsDir)
                print('Created directory: {}'.format(fitsDir))

            try:
                tarF.extract(mapname[0],path=fitsDir)
                mapFile=os.path.join(fitsDir,mapname[0])
                if (verbose): print('Extracted Map to {}'.format(mapFile))

                # save location of map to status
                # save location of h5 file to status
                stat={'mapurllocal':mapFile,
                    'mapdatelocal':Time.now().isot,
                    'mapdatesrc':Time(os.path.getmtime(mapFile),format='unix').isot,
                    'mapverlocal':self.data[ev]["version"],
                    'mapversrc':self.data[ev]["version"]}
                if 'maptarurlsrc' in self.status[ev]:
                    stat['mapurlsrc']=self.status[ev]['maptarurlsrc']
                elif 'tarurlsrc' in self.status[ev]:
                    stat['mapurlsrc']=self.status[ev]['tarurlsrc']
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='map file local')
                out['map']=mapFile
                self.addLink(ev,{'url':self.status[ev]["mapurlsrc"],'text':'Sky Map',
                    'type':'skymap-fits','created':self.status[ev]["mapdatesrc"],'filetype':'tar'})
            except:
                print('WARNING: Error extracting {} map {} from tar {}'.format(ev,mapname[0],tarFile))
                tarF.extract(mapname[0],path=fitsDir)
                mapFile=os.path.join(fitsDir,mapname[0])
                if (verbose): print('Extracted Map to {}'.format(mapFile))

                # save location of map to status
                # save location of h5 file to status
                stat={'mapurllocal':mapFile,
                    'mapdatelocal':Time.now().isot,
                    'mapdatesrc':Time(os.path.getmtime(mapFile),format='unix').isot,
                    'mapverlocal':self.data[ev]["version"],
                    'mapversrc':self.data[ev]["version"]}
                if 'maptarurlsrc' in self.status[ev]:
                    stat['mapurlsrc']=self.status[ev]['maptarurlsrc']
                elif 'tarurlsrc' in self.status[ev]:
                    stat['mapurlsrc']=self.status[ev]['tarurlsrc']
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='map file local')
                out['map']=mapFile
                self.addLink(ev,{'url':self.status[ev]["mapurlsrc"],'text':'Sky Map',
                    'type':'skymap-fits','created':self.status[ev]["mapdatesrc"],'filetype':'tar'})
                pass
        else:
            print('no map found for {} (approx {}) in {}'.format(ev,tarapproxtarFile))
        # return h5file and mapfile locations?
        return out

    def getH5Local(self,ev,verbose=False):
        """Download h5 files from local file and update status with filenames, timestamps.
        Call extractTar if tarfile.
        Check file is valid H5 file.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * [object] dict containing extracted "h5" and "map" file locations.
        """
        # import local h5 file into database
        ldat=self.getLink(ev,'data-local')
        if len(ldat)==0:
            if verbose: print('ERROR: no local data file link for {} [{}]'.format(ev,self.data[ev]['catalog']))
            return
        if len(ldat)>1:
            if verbose: print('WARNING: more than one local data file link for {} [{}]'.format(ev,self.data[ev]['catalog']))
        h5File=ldat[0]['url']
        if verbose: print('Using local h5 file for {}: {}'.format(ev,h5File))
        if ldat[0]['filetype']=='tar':
            print('NEED TO EXTRACT TARBALL')
            try:
                stat={'tarurllocal':h5File,
                    'tarurlsrc':ldat[0]['url'],
                    'tardatelocal':Time.now().isot,
                    'tardatesrc':Time(os.path.getmtime(h5File),format='unix').isot,
                    'tarverlocal':self.data[ev]["version"],
                    'tarversrc':self.data[ev]["version"]}
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='tar file local')
                tarout=self.extractTar(ev,verbose=verbose)
                h5File=tarout['h5']
            except:
                print('ERROR: Problem opening local data file for {}:'.format(ev),h5File)
                return
        try:
            testread=h5py.File(h5File)
            if verbose: print('Valid h5 file for {}: {}'.format(ev,h5File))
            # hdr=fits.getheader(fitsFile,ext=1)
            stat={'h5urllocal':h5File,
                'h5datelocal':Time.now().isot,
                'h5datesrc':Time(os.path.getmtime(h5File),format='unix').isot}
            self.data[ev]['meta']['h5urllocal']=h5File
            self.data[ev]['meta']['h5datelocal']=Time.now().isot
            self.data[ev]['meta']['h5datesrc']=Time(os.path.getmtime(h5File),format='unix').isot
            self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='data file local')
        except:
            print('ERROR: Problem opening local data file for {}:'.format(ev),h5File)
            return
        return

    def getH5Params(self,ev,approx=None,verbose=False):
        """Get parameters from H5 file for GWTC-2. Uses comparison of parameter to select best approximant if default isn't available.
        Inputs:
            * ev [string]: event name
            * approx [string, optional]: approximant to try first
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * [object] dict containing new parameters. None if not GWTC2
        """
        validcats={'GWTC-2':{'approx':'PublicationSamples'},'GWTC-2.1-confident':{},'O3_Discovery_Papers':{'approx':'PublicationSamples'},'GWTC-3-confident':{'approx':'C01:Mixed'},'O4_Discovery_Papers':{},'GWTC-4.0':{}}
        if not self.getParameter(ev,'catalog') in validcats:
            if verbose:print('not valid catalogue')
            newparams=None
        else:
            if not approx:
                if 'approx' in validcats[self.getParameter(ev,'catalog')]:
                    approx=validcats[self.getParameter(ev,'catalog')]['approx']
            newparams={}
            if verbose:print('getting H5 parameters for {}'.format(ev))
            m1check=self.getParameter(ev,'M1')
            if not (m1check):
                if verbose:print('no M1 parameter for {}'.format(ev))
            else:
                m1check=m1check['best']
                if verbose:print('M1 parameter for {}:'.format(ev),m1check)
            try:
                h5File=self.status[ev]['h5urllocal']
            except:
                if verbose:print('no local data file for {}'.format(ev))
            newparams=gwosc.geth5paramsGWTC2(h5File,pcheck={'M1':m1check},datadict=self.datadict,approx=approx,verbose=verbose)
            # if verbose:print(newparams)
        return(newparams)

    def removeCandidates(self,verbose=False):
        """remove candidates that have a firm detection
        Inputs:
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        remCands=[]
        for ev in self.data:
            if 'gracedb' in self.data[ev]['meta']:
                evgdb=self.data[ev]['meta']['gracedb']
                if ev!=evgdb:
                    remCands.append(evgdb)
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        for remev in remCands:
            if remev in self.data:
                if verbose: print('Removing {}'.format(remev))
                self.data.pop(remev)
            if remev in self.links:
                self.links.pop(remev)
        return

    def updateMaps(self,verbose=False,forceUpdate=False,event=None):
        """Check whether map files need updating and re-download if necessary (using getMap).
        Update 90% area.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
            * forceUpdate [boolean, optional]: set to download all files. Default=False (only download updated files)
        Outputs:
            * None
        """
        print('*** Updating maps...')
        for ev in self.data:
            if event!="" and event!=None:
                if ev!=event:
                    continue
                else:
                    if verbose: print("***UPDATING MAPS FOR SINGLE EVENT {}".format(ev))
            if verbose: print('Checking map status for {}'.format(ev))
            # compare map creation dates
            try:
                if 'mapdatelocal' in self.status[ev] and 'mapdatesrc' in self.status[ev]:
                    mapdatelocal=Time(self.status[ev]['mapdatelocal']).gps
                    mapdatesrc=Time(self.status[ev]['mapdatesrc']).gps
                else:
                    mapdatelocal=-1
                    mapdatesrc=-1
                if 'mapverlocal' in self.status[ev] and 'mapversrc' in self.status[ev]:
                    mapverlocal=self.status[ev]['mapverlocal']
                    mapversrc=self.status[ev]['mapversrc']
                else:
                    mapverlocal=-1
                    mapversrc=-1
                if 'mapzenodolocal' in self.status[ev] and 'mapzenodosrc' in self.status[ev]:
                    mapzenodolocal=self.status[ev]['mapzenodolocal']
                    mapzenodosrc=self.status[ev]['mapzenodosrc']
                else:
                    mapzenodolocal=-1
                    mapzenodosrc=-1
                mapfile=self.status[ev].get('mapurllocal','')
                if not os.path.isfile(mapfile):
                    if verbose:print('No file. Need to re-download map for {}'.format(ev))
                    updateMap=True
                elif mapzenodosrc>0 and mapzenodolocal>0 and mapzenodosrc>mapzenodolocal:
                    if verbose:print('Newer Zenodo version. Need to re-download map for {}'.format(ev))
                    updateMap=True
                elif mapdatesrc>0 and mapdatelocal>0 and (mapdatesrc-mapdatelocal)>1:
                    if verbose:print('Newer date file. Need to re-download map for {}'.format(ev))
                    print('src',Time(mapdatesrc,format='gps').isot)
                    print('local',Time(mapdatelocal,format='gps').isot)
                    updateMap=True
                elif mapversrc>0 and mapverlocal>0 and mapversrc>mapverlocal:
                    if verbose:print('Newer version. Need to re-download map for {}'.format(ev))
                    updateMap=True
                else:
                    if verbose:print('Map up-to-date. No need to re-download map for {}'.format(ev))
                    updateMap=False
            except:
                if verbose:print('Unable to check map status. Redownloading')
                updateMap=True
            if updateMap or forceUpdate:
                # if verbose:print('{} status: {}'.format(ev,self.status[ev]))
                self.getMap(ev,verbose=verbose)
                if verbose:print('Updating map for {}'.format(ev))
                self.calcAreas(ev,verbose=verbose)
            else:
                if verbose:print('No map update required for {}'.format(ev))
                if not 'deltaOmega' in self.data[ev]:
                    if verbose:print('recalculating area')
                    self.calcAreas(ev,verbose=verbose)
            try:
                fitsFile=self.status[ev]['mapurllocal']
                self.addLink(ev,{'url':self.rel2abs(fitsFile,url=self.dataurl),'text':'Sky Map (local mirror)',
                    'type':'skymap-fits-local'})
                if verbose:print('added map link: {}'.format(fitsFile))
            except:
                if verbose:print('failed to add map link for {}: {}'.format(ev,self.status[ev].get('mapurllocal','UNKNOWN')))
                # fitsFile=self.status[ev]['mapurllocal']
                # self.addLink(ev,{'url':self.rel2abs(fitsFile,url=self.dataurl),'text':'Sky Map (local mirror)',
                #     'type':'skymap-fits-local'})

        return

    def getMap(self,ev,verbose=False):
        """Download map files and update status with filenames, timestamps.
        Calculate 90% areas.
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * [object] FITS header or HTTP status code if unsuccessful
        """
        fitsDir=os.path.join(self.dataDir,'fits')
        if not os.path.exists(fitsDir):
            # create directory
            os.mkdir(fitsDir)
            print('Created directory: {}'.format(fitsDir))
        lmap=self.getLink(ev,'skymap-fits')
        if len(lmap)==0:
            if verbose: print('ERROR: no skymap link for {}',format(ev))
            return
        if len(lmap)>1:
            if verbose: print('WARNING: more than one skymap link for {}',format(ev))
        url=lmap[0]['url']
        if verbose: print('Downloading skymap for {} from {}'.format(ev,url))
        if url.find("zenodo")>=0 and url.find("content")>=0:
            try:
                # New Zenodo format file
                zenurl=url.replace("/content","")
                srcfile=os.path.split(zenurl)[-1]
                if verbose: print('Zenodo src file:',srcfile)
            except:
                srcfile=os.path.split(url)[-1]
        else:
            srcfile=os.path.split(url)[-1]
        if srcfile.find(ev)<0 and srcfile.find('zenodo')<0 and srcfile.lower().find('skymaps.tar.gz')<0:
            fitsFile=os.path.join(fitsDir,'{}_{}'.format(ev,srcfile))
            if fitsFile.find('.fits')<0:
                fitsFile=fitsFile+'.fits'
        else:
            fitsFile=os.path.join(fitsDir,srcfile)
        if url.find('http')<0:
            # local file
            if verbose: print('copying local file from {} to {}'.format(srcfile,fitsFile))
            os.system('copy {} {}'.format(url,fitsFile))
        else:
            if self.devMode:
                mapreq=self.sess.get(url)
            else:
                mapreq=requests.get(url)
            try:
                fOut=open(fitsFile,'wb')
                fOut.write(mapreq.content)
                fOut.close()
                if verbose:print('Map downloaded from {} to {}'.format(srcfile,fitsFile))
            except:
                print('ERROR: Problem loading map:',mapreq.status_code)
                return mapreq.status_code

        if lmap[0].get('filetype','')=='tar':
            print('NEED TO EXTRACT TARBALL')
            stat={'maptarurllocal':fitsFile,
                'maptarurlsrc':url,
                'maptardatelocal':Time.now().isot,
                'maptardatesrc':Time(os.path.getmtime(fitsFile),format='unix').isot,
                'maptarverlocal':self.data[ev]["version"],
                'maptarversrc':self.data[ev]["version"]}
            self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='tar file local (map)')
            if 'maptarurllocal' in self.status[ev]:
                try:
                    tarout=self.extractTar(ev,statusid='maptarurllocal',maponly=True,verbose=verbose)
                    fitsFile=tarout['map']
                except:
                    print('problem reading maptarurllocal file: {}'.format(self.status[ev]['maptarurllocal']))
            elif 'tarurllocal' in self.status[ev]:
                try:
                    tarout=self.extractTar(ev,statusid='tarurllocal',maponly=True,verbose=verbose)
                    fitsFile=tarout['map']
                except:
                    print('problem reading tarurllocal file: {}'.format(self.status[ev]['tarurllocal']))
            else:
                print('ERROR: no status entry for tarurllocal or maptarurllocal')
        try:
            hdr=fits.getheader(fitsFile,ext=1)
            self.data[ev]['meta']['mapurllocal']=fitsFile
            self.data[ev]['meta']['mapdatelocal']=Time.now().isot
            stat={'mapurllocal':fitsFile,
                'mapdatelocal':Time.now().isot}
            if 'version' in self.data[ev]:
                stat['mapverlocal']=self.data[ev]["version"]
                stat['mapversrc']=self.data[ev]["version"]
            self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='maplocal')
        except:
            print('ERROR: Problem opening fits file for {}:'.format(ev),fitsFile)
            return
        self.calcAreas(ev,verbose=verbose)

        return hdr

    def calcAreas(self,ev,verbose=False):
        """Calculate 90% area of map (using plotloc.read_map, plotloc.getProbMap, plotloc.getArea) and save to database
        Inputs:
            * ev [string]: event name
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * [object] 90% likelihood area (in degree^2)
        """
        if 'mapurllocal' in self.status[ev]:
            fitsFile=self.status[ev]['mapurllocal']
            map=plotloc.read_map(fitsFile,verbose=verbose)

            totmap,a90=plotloc.getProbMap(map,prob=0.9,verbose=verbose)
            a50=plotloc.getArea(totmap,0.5,verbose=verbose)
            self.data[ev]['deltaOmega']={'best':round(a90)}
            self.data[ev]['skyarea(50)']={'best':round(a50)}
            if verbose:print('90% area',round(a90),'50% area',round(a50))
            return(a90)
            # except:
            #     print('WARNING: Problem calculating area for {}'.format(ev))
                # return

    def rel2abs(self,rel,data=False,url=None):
        """Convert relative to absolute URL (for links)
        Inputs:
            * rel [string]: relative URL to use
            * url [string, optional]: Absolute URL to pre-pend. Default = self.baseurl
        Outputs:
            * [string] absolute url
        """
        if url==None:
            if data:
                url=self.dataurl
            else:
                url=self.baseurl
        return(url + rel)

    def plotMapPngs(self,overwrite=False,verbose=False,logFile=None,updateLink=True,lowSigMaps=False,event=None):
        """Create maps of event localisations in various projections, zooms, styles etc.
        Save links to database.
        Inputs:
            * overwrite [boolean, optional]: set to overwrite all plots. Default=False (only output those needed)
            * verbose [boolean, optional]: set for verbose output. Default=False
            * logFile [string, optional]: set for output to logFile. Default=None (no logging)
            * updateLink [boolean, optional]: set to add/update link even if plot doesn't need making. Default=True
            * lowSigMaps [boolean, optional]: set to plot maps for events marked low significance. Default=False
            * event [string]: set to limit to one single event
        Outputs:
            * None
        """
        print('*** Updating plots...')
        if logFile:
            if os.path.exists(logFile):
                os.remove(logFile)
                print('Removing log file: {}'.format(logFile))
            else:
                print("Log file doesn't exist: {}".format(logFile))
            print('Writing Maps log to: {}'.format(logFile))
            logF=open(logFile,'a')

        pngDir=os.path.join(self.dataDir,'png/')
        gravDir=os.path.join(self.dataDir,'gravoscope/')
        dataDir=os.path.join(self.dataDir,'fits/')
        if not os.path.exists(pngDir):
            os.mkdir(pngDir)
        if not os.path.exists(gravDir):
            os.mkdir(gravDir)
        # print(self.data)
        for ev in self.data:
            if event!="" and event!=None:
                if ev!=event:
                    continue
                else:
                    if verbose: print("***PLOTTING MAPS FOR SINGLE EVENT {}".format(ev))
            print('plotting {}'.format(ev))
            # if not 'mapurlsrc' in self.status[ev]:
            #     if verbose:
            #         print('no mapurlsrc for {}'.format(ev))
            #     continue
            if not 'mapurllocal' in self.status[ev]:
                self.getMap(ev,verbose=verbose)
            if not 'mapurllocal' in self.status[ev]:
                if verbose:print('unable to get map for {}'.format(ev))
                continue
            if self.data[ev].get('Significance'):
                if verbose:print('{} significance: {}'.format(ev,self.data[ev].get('Significance')))
                if self.data[ev]['Significance']['best']=='Low':
                    if lowSigMaps:
                        if verbose:print('plotting maps for low-significance event {}'.format(ev))
                    else:
                        if verbose:print('SKIPPING maps for low-significance event {}'.format(ev))
                        continue
                else:
                    if verbose:print('plotting maps for high-significance event {}'.format(ev))
            filename=self.status[ev]['mapurllocal']
            # if verbose:print('plotting maps at {}'.format(filename))
            fitsCreated=Time(self.status[ev]['mapdatelocal'])
            if 'mapurlsrc' in self.status[ev]:
                srcfile=os.path.split(self.status[ev]['mapurlsrc'])[-1]
            else:
                srcfile=os.path.split(self.status[ev]['mapurllocal'])[-1]
            ptitle='{} [{}]'.format(ev,srcfile)
            plots={'moll':{'linktxt':'Skymap (Mollweide fullsky)'},
                'moll_pretty':{'linktxt':'Skymap (Mollweide fullsky, pretty)'},
                'moll_pretty_black':{'linktxt':'Skymap (Mollweide fullsky, pretty black)'},
                'cartzoom':{'linktxt':'Skymap (Cartesian zoomed)'},
                'cartzoom_pretty':{'linktxt':'Skymap (Cartesian zoomed, pretty)'},
                'cart':{'linktxt':'Skymap (Cartesian fullsky)'},
                'moll_rot':{'linktxt':'Skymap (Mollweide fullsky, rotated)'},
                'cart_rot':{'linktxt':'Skymap (Cartesian fullsky, rotated)'}
            }
            imgs={}
            thumbs={}
            nUpdate=0
            nUpdateLinks=0
            for p in plots:
                plots[p]['pngFileOnly']='{}_{}.png'.format(ev,p)
                plots[p]['thumbFileOnly']='{}_{}.thumb.png'.format(ev,p)
                plots[p]['pngFile']=os.path.join(pngDir,plots[p]['pngFileOnly'])
                plots[p]['thumbFile']=os.path.join(pngDir,plots[p]['thumbFileOnly'])
                plots[p]['exists']=os.path.isfile(plots[p]['pngFile'])
                if not plots[p]['exists'] or overwrite:
                    plots[p]['update']=True
                else:
                    plots[p]['update']=False
                link=self.getLink(ev,'skymap-plot',file=p)
                if len(link)>0:
                    if 'created' in link[0]:
                        if link[0]['created']<fitsCreated:
                            plots[p]['update']=True
                else:
                    nUpdateLinks+=1
                if plots[p]['update']: nUpdate+=1
                imgs[p]={'file':plots[p]['pngFileOnly'],'text':plots[p]['linktxt']}
                thumbs[p]={'file':plots[p]['thumbFileOnly'],'text':plots[p]['linktxt']}

            if nUpdate==0:
                if verbose:print('all plots exist for {}'.format(ev))
                mapread=False
                if nUpdateLinks>0 or updateLink:
                    if verbose: print('skipping plotting {} maps. Adding links'.format(ev))
                    for p in plots:
                        pp=plots[p]
                        # add links
                        # self.addLink(ev,
                        #     {'url':self.rel2abs(pp['pngFile']),'text':pp['linktxt'],
                        #     'file':pp['pngFile'],'url-loc':'skymap-base-url',
                        #     'type':'skymap-plot','created':Time.now().isot})
                        # self.addLink(ev,
                        #     {'url':self.rel2abs(pp['thumbFile']),'text':pp['linktxt'],
                        #     'file':pp['thumbFile'],'url-loc':'skymap-base-url',
                        #     'type':'skymap-thumbnail','created':Time.now().isot})
                    self.addLink(ev,
                        {'url':self.rel2abs(pngDir),'text':'Skymaps',
                        'type':'skymaps-plot','created':Time.now().isot,
                        'files':imgs})
                    self.addLink(ev,
                        {'url':self.rel2abs(pngDir),'text':'Skymaps',
                        'type':'skymaps-thumb','created':Time.now().isot,
                        'files':thumbs})
            else:
                try:
                    map=plotloc.read_map(filename,verbose=verbose)
                    mapread=True
                except:
                    print('ERROR: problem reading map at {}'.format(filename))
                    return
                if logFile:
                    logF.write(ev+'\n')
                for p in plots:
                    pp=plots[p]
                    # set defaults
                    logos=True
                    credit=True
                    title='{}'.format(ptitle)
                    plotbounds=True
                    plotlabels=True
                    plotlines=True
                    linewidth=1
                    markersize=2
                    alpha=1
                    margins=None
                    notext=False
                    bgcolor='w'
                    border=None
                    lw=None
                    fontsize=10
                    if p=='cartzoom':
                        zoomlim=0.8
                        rotmap=True
                        minzoom=20
                        proj='cart'
                    elif p=='cartzoom_pretty':
                        zoomlim=0.8
                        rotmap=True
                        minzoom=20
                        proj='cart'
                        logos=False
                        credit=False
                        notext=True
                        title=' '
                        bgcolor='k'
                        margins=[0,0,0,0]
                        lw=3
                        fontsize=20
                    elif p=='cart':
                        zoomlim=None
                        rotmap=False
                        minzoom=None
                        proj='cart'
                    elif p=='moll':
                        zoomlim=None
                        rotmap=False
                        minzoom=None
                        proj='moll'
                    elif p=='moll_pretty':
                        zoomlim=None
                        rotmap=False
                        minzoom=None
                        proj='moll'
                        logos=False
                        credit=False
                        notext=True
                        plotlines=False
                        plotlabels=False
                        plotbounds=False
                        bgcolor='w'
                        margins=[0,0,0,0]
                        title=' '
                        border='black'
                        lw=2
                    elif p=='moll_pretty_black':
                        zoomlim=None
                        rotmap=False
                        minzoom=None
                        proj='moll'
                        logos=False
                        credit=False
                        notext=True
                        plotlines=False
                        plotlabels=False
                        plotbounds=False
                        bgcolor='black'
                        margins=[0,0,0,0]
                        title=' '
                        border='white'
                        lw=2
                    elif p=='moll_rot':
                        zoomlim=None
                        rotmap=True
                        minzoom=None
                        proj='moll'
                    elif p=='cart_rot':
                        zoomlim=None
                        rotmap=True
                        minzoom=None
                        proj='cart'
                    # plot map
                    if not pp['update']:
                        if verbose: print('skipping plotting {} map. Adding links'.format(pp['linktxt']))
                        # add links
                        # self.addLink(ev,
                        #     {'url':self.rel2abs(pp['pngFile']),'text':pp['linktxt'],
                        #     'file':pp['pngFile'],'url-loc':'skymap-base-url',
                        #     'type':'skymap-plot','created':Time.now().isot})
                        # self.addLink(ev,
                        #     {'url':self.rel2abs(pp['thumbFile']),'text':pp['linktxt'],
                        #     'file':pp['thumbFile'],'url-loc':'skymap-base-url',
                        #     'type':'skymap-thumbnail','created':Time.now().isot})
                    else:
                        if verbose:print('plotting {} map to {}'.format(pp['linktxt'],pp['pngFile']))
                        plotloc.makePlot(ev=ev,mapIn=map,dirData=dataDir,
                            proj=proj,plotcont=False,smooth=0,zoomlim=zoomlim,
                            rotmap=rotmap,minzoom=minzoom,margins=margins,
                            verbose=verbose,title=title,notext=notext,cbg=bgcolor,
                            pngOut=pp['pngFile'],thumbOut=pp['thumbFile'],
                            plotbounds=plotbounds,plotlabels=plotlabels,plotlines=plotlines,
                            addCredit=credit,addLogos=logos,border=border,lw=lw,fontsize=fontsize)
                        # add links
                        # self.addLink(ev,
                        #     {'url':self.rel2abs(pp['pngFile']),'text':pp['linktxt'],
                        #     'file':pp['pngFile'],'url-loc':'skymap-base-url',
                        #     'type':'skymap-plot','created':Time.now().isot})
                        # self.addLink(ev,
                        #     {'url':self.rel2abs(pp['thumbFile']),'text':pp['linktxt'],
                        #     'file':pp['thumbFile'],'url-loc':'skymap-base-url',
                        #     'type':'skymap-thumbnail','created':Time.now().isot})
                self.addLink(ev,
                    {'url':self.rel2abs(pngDir),'text':'Skymaps',
                    'type':'skymaps-plot','created':Time.now().isot,
                    'files':imgs})
                self.addLink(ev,
                    {'url':self.rel2abs(pngDir),'text':'Skymaps',
                    'type':'skymaps-thumb','created':Time.now().isot,
                    'files':thumbs})

            gravs={'gal_8192':{'text':'Skymap (no annotations)'},
                'eq_8192':{'text':'Skymap (Equatorial, no annotations)'},
                'eq_4096':{'text':'Skymap 4096px (Equatorial, no annotations)'}
            }
            res=8
            gravNpix=int(res*1024)
            updateGrav=False
            gravFileOnly='{}_{}.png'.format(ev,gravNpix)
            gravFile=os.path.join(gravDir,gravFileOnly)
            gravs['gal_8192']['file']=gravFileOnly
            if not os.path.isfile(gravFile):
                updateGrav=True
            gravLink=self.getLink(ev,'skymap-plain',file='gal_8192')
            if len(gravLink)>0:
                if 'created' in gravLink[0]:
                    if gravLink[0]['created']<fitsCreated:
                        updateGrav=True
            if updateGrav:
                if not mapread:
                    map=plotloc.read_map(filename,verbose=verbose)
                if verbose:
                    print('plotting Gravoscope for {} ({}x{})'.format(ev,gravNpix,int(gravNpix/2)))
                plotloc.plotGravoscope(mapIn=map,pngOut=gravFile,verbose=verbose,res=res)
                self.addLink(ev,
                    {'url':self.rel2abs(gravFile),'text':gravs['gal_8192']['text'],
                    'file':gravFile,'url-loc':'skymap-base-url',
                    'type':'skymap-plain','created':Time.now().isot})

            updateGravEq=False
            gravFileOnlyEq='{}_{}_eq.png'.format(ev,gravNpix)
            gravFileEq=os.path.join(gravDir,gravFileOnlyEq)
            gravs['eq_8192']['file']=gravFileOnlyEq
            if not os.path.isfile(gravFileEq):
                updateGravEq=True
            gravLinkEq=self.getLink(ev,'skymap-plain',file='eq_8192')
            if len(gravLinkEq)>0:
                if 'created' in gravLinkEq[0]:
                    if gravLinkEq[0]['created']<fitsCreated:
                        updateGravEq=True
            if updateGravEq:
                if not mapread:
                    map=plotloc.read_map(filename,verbose=verbose)
                if verbose:
                    print('plotting Gravoscope (Equatorial) for {} ({}x{})'.format(ev,gravNpix,int(gravNpix/2)))
                plotloc.plotGravoscope(mapIn=map,pngOut=gravFileEq,verbose=verbose,res=res,coord='C')
                self.addLink(ev,
                    {'url':self.rel2abs(gravFileEq),'text':gravs['eq_8192']['text'],
                    'file':gravFileEq,'url-loc':'skymap-base-url',
                    'type':'skymaps-plain','created':Time.now().isot})

            res4096=4
            gravNpix4096=int(res4096*1024)
            updateGravEq4096=False
            gravFileOnlyEq4096='{}_{}_eq.png'.format(ev,gravNpix4096)
            gravFileEq4096=os.path.join(gravDir,gravFileOnlyEq4096)
            gravs['eq_4096']['file']=gravFileOnlyEq4096
            if not os.path.isfile(gravFileEq4096):
                updateGravEq4096=True
            gravLinkEq4096=self.getLink(ev,'skymap-plain',file='eq_4096')
            if len(gravLinkEq4096)>0:
                if 'created' in gravLinkEq4096[0]:
                    if gravLinkEq4096[0]['created']<fitsCreated:
                        updateGravEq4096=True
            if updateGravEq4096:
                if not mapread:
                    map=plotloc.read_map(filename,verbose=verbose)
                if verbose:
                    print('plotting Gravoscope (Equatorial, 4096) for {} ({}x{})'.format(ev,gravNpix4096,int(gravNpix4096/2)))
                plotloc.plotGravoscope(mapIn=map,pngOut=gravFileEq4096,verbose=verbose,res=res4096,coord='C')
                self.addLink(ev,
                    {'url':self.rel2abs(gravFileEq4096),'text':gravs['eq_4096']['text'],
                    'file':gravFileEq4096,'url-loc':'skymap-base-url',
                    'type':'skymap-plain','created':Time.now().isot})
            # self.addLink(ev,
            #     {'url':self.rel2abs(''),'text':'Skymap base url',
            #     'type':'skymap-base-url','created':Time.now().isot})
            self.addLink(ev,
                {'url':self.rel2abs(gravDir,data=True),'text':'Skymaps (plain)',
                'type':'skymaps-plain','created':Time.now().isot,
                'files':gravs})

        if logFile:
            logF.close()
        return

    def makeGravoscopeTilesPerl(self,overwrite=False,verbose=False,event=None):
        """OBSOLETE Create tile-sets of event localisation maps using perl script for use with Gravoscope.
        Save links to database.
        Inputs:
            * overwrite [boolean, optional]: set to overwrite all plots. Default=False (only output those needed)
            * verbose [boolean, optional]: set for verbose output. Default=False
            * event [string, optional]: set to plot for single event. Default=None
        Outputs:
            * None
        """

        gravDir=os.path.join(self.dataDir,'gravoscope')
        dataDir=os.path.join(self.dataDir,'fits')
        for ev in self.events:
            if event!="" and event !=None:
                if ev != event:
                    continue
                else:
                    if verbose: print("PLOTTING GRAVOSCOPE FOR SINGLE EVENT {}".format(ev))
            fitsCreated=Time(self.status[ev]['mapdatelocal'])
            res=8
            gravNpix=int(8*1024)
            updateGrav=False
            gravFile=os.path.join(gravDir,'{}_{}.png'.format(ev,gravNpix))
            tileFile=os.path.join(gravDir,'{}_{}-tiles/tt.png'.format(ev,gravNpix))

            if not os.path.isfile(gravFile):
                if verbose:print('No source PNG file for {}: {}'.format(ev,gravFile))
            else:
                if not os.path.isfile(tileFile) or overwrite:
                    if verbose:print('plotting Gravoscope files for {}: {}'.format(ev,gravFile))
                    plotloc.makeTiles(gravFile,verbose=verbose)
        return

    def makeGravoscopeTiles(self,maxres=3,overwrite=False,verbose=False,tilesurl=None,updateLink=True):
        """Create tile-sets of event localisation maps using perl script for use with Gravoscope.
        Save links to database.
        Inputs:
            * overwrite [boolean, optional]: set to overwrite all plots. Default=False (only output those needed)
            * verbose [boolean, optional]: set for verbose output. Default=False
            * tilesurl [string, optional]: base URL for tilesets. Default=None (uses self.baseurl)
            * updateLink [boolean, optional]: set to update link (regardless of whether tiles are created). Default=True
        Outputs:
            * None
        """
        gravDir=os.path.join(self.dataDir,'gravoscope')
        for ev in self.events:
            tilesDir=os.path.join(gravDir,'{}-tiles'.format(ev))
            if not os.path.exists(tilesDir):
                os.mkdir(tilesDir)
            if not 'mapdatelocal' in self.status[ev]:
                if verbose: print('no local map found for {}. Skipping.'.format(ev))
                continue
            fitsCreated=Time(self.status[ev]['mapdatelocal'])
            filename=self.status[ev]['mapurllocal']
            gravLinktxt='Gravoscope tileset'
            tileFile=os.path.join(gravDir,'{}-tiles/{}.png'.format(ev,'ttrtttttt'[0:maxres+1]))
            tilesLinktype='gravoscope-tiles'
            tilesLink=self.getLink(ev,tilesLinktype,srchtype='type')
            tilesLinkIdx=self.getLink(ev,tilesLinktype,srchtype='type',retIdx=True)
            updateTiles=False
            if not os.path.isfile(tileFile):
                if verbose: print('file {} does not exist'.format(tileFile))
                updateTiles=True
            elif overwrite:
                if verbose: print('file {} exists, but overwriting'.format(tileFile))
                updateTiles=True
            else:
                if verbose: print('file {} exists'.format(tileFile))
            if not updateTiles:
                # file exists, but need to check date
                timeTiles=Time(os.path.getmtime(tileFile),format='unix')
                if timeTiles.gps < fitsCreated.gps-1:
                    # tiles are older than map
                    if verbose: print('tiles for {} exists, but are older than map. {} < {}; {} < {}'.format(ev,timeTiles.isot,fitsCreated,timeTiles.gps,fitsCreated.gps))
                    updateTiles=True
                else:
                    # tiles are newer than map
                    if verbose: print('tiles for {} exists, and are newer than map. {} > {}; {} > {}'.format(ev,timeTiles.isot,fitsCreated,timeTiles.gps,fitsCreated.gps))
                    if updateLink:
                        # update link (but don't necessarily recreate tiles)
                        if verbose: print('replacing old link for Gravoscope tileset for {}'.format(ev))
                        self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                            'type':'gravoscope-tiles','created':timeTiles.isot})
                if len(tilesLink)>0:
                    if not 'created' in tilesLink[0]:
                        # no created date. recreate link
                        if verbose: print('No created date for tiles link for {}. Regenerating all tiles'.format(ev))
                        self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                            'type':'gravoscope-tiles','created':timeTiles.isot})
                else:
                    # no link. Need to add link
                    if verbose: print('adding tiles link for Gravoscope tileset for {}'.format(ev))
                    gravLinktxt='Gravoscope tileset'
                    self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                        'type':'gravoscope-tiles','created':fitsCreated.isot})
            res=8
            gravNpix=int(8*1024)

            if updateTiles:
                if verbose:print('plotting Gravoscope files for {}: {}'.format(ev,tilesDir))
                map=plotloc.read_map(filename,verbose=verbose)
                plotloc.makeTiles(map,dirOut=tilesDir,maxres=maxres,verbose=verbose)
                if verbose: print('adding tiles link for Gravoscope tileset for {}'.format(ev))
                self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                    'type':'gravoscope-tiles','created':fitsCreated.isot})
        return

    def getLink(self,ev,srchtxt,srchtype='type',verbose=False,retIdx=False,file=None):
        """Get link(s) for event based on string search of type or text
        Inputs:
            * ev [string]: event name
            * srchtxt [string]: text to search for (assumes regex)
            * srchtype [string, optional]: field to search. Default='type'
            * verbose [boolean, optional]: set for verbose output. Default=False
            * retIdx [boolean, optional]: set to output index of link. Default=False (output link itself)
            * file [string, optional]: file to search for within link
        Outputs:
            * If retIdx==True: [list] List of indices of matching links
            * If retIdx==False: [list] List of matching link objects
        """
        if not ev in self.links:
            return []
        lOut=[]
        if len (self.links[ev])>0:
            for ol in range(len(self.links[ev])):
                if self.links[ev][ol][srchtype]==srchtxt:
                    if retIdx:
                        lOut.append(ol)
                    else:
                        linkOut=self.links[ev][ol]
                        if file and 'files' in self.links[ev][ol]:
                            linkOut['url']=os.path.join(linkOut['url'],linkOut['files'][file]['file'])
                            if verbose:
                                print('found {} links for {} with {}=={} and file={}'.format(len(lOut),ev,srchtype,srchtxt,file))
                        else:
                            if verbose:
                                print('found {} links for {} with {}=={}'.format(len(lOut),ev,srchtype,srchtxt))
                        lOut.append(linkOut)
        return(lOut)


    def removeLink(self,ev,srchtxt,srchtype='type',verbose=False,retIdx=False,file=None):
        """Get link(s) for event based on string search of type or text
        Inputs:
            * ev [string]: event name
            * srchtxt [string]: text to search for (assumes regex)
            * srchtype [string, optional]: field to search. Default='type'
            * verbose [boolean, optional]: set for verbose output. Default=False
            * retIdx [boolean, optional]: set to output index of link. Default=False (output link itself)
            * file [string, optional]: file to search for within link
        Outputs:
            * If retIdx==True: [list] List of indices of matching links
            * If retIdx==False: [list] List of matching link objects
        """
        if not ev in self.links:
            return
        nRem=0
        idxOut=[]
        idxIn=[]
        if len (self.links[ev])>0:
            for ol in range(len(self.links[ev])):
                if self.links[ev][ol][srchtype]==srchtxt:
                    nRem=nRem+1
                    idxOut.append(ol)
                else:
                    idxIn.append(ol)
                    # self.links[ev].remove(self.links[ev][ol])
            if verbose:print('removing {} links from {}'.format(nRem,ev))
            idxOut.sort(reverse=True)
            for i in idxOut:
                self.links[ev].pop(i)
        return

    def addLink(self,ev,link,replace=True,verbose=False):
        """Add or replace link(s) for event
        Inputs:
            * ev [string]: event name
            * link [object]: json-style object to add as link
            * replace [bookean, optional]: set to replace links with same type and text. Default=True.
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        # replace: replace link with same type and text. (N.B. Always replace open-data links)
        if not ev in self.links:
            self.links[ev]=[]
        # ltype=link['type']
        exlink=None
        if verbose:print('adding {} link to {} [{}]'.format(link['type'],ev,link['url']))
        if len(self.links[ev])>0:
            for ol in range(len(self.links[ev])):
                oldlink=self.links[ev][ol]
                if replace==True:
                    # find link with same type and text
                    if oldlink['type']==link['type'] and oldlink['text']==link['text']:
                        exlink=ol
                    elif link['type']=='open-data' and oldlink['type']==link['type']:
                        # always remove open-data if it exists
                        if verbose:print('removing old {} link from {} [{}]'.format(ltype,g,self.links[ev][ol]['url']))
                        exlink=ol
        if (exlink):
            # replace link
            self.links[ev][exlink]=link
        else:
            # append new link
            self.links[ev].append(link)
        return

    def json2dataframe(self,verbose=False):
        """Convert dictionaries into pandas DataFrames
        data -> events
        datadict -> units
        links -> eventrefs
        Inputs:
            * None
        Outputs:
            * Tuple of DataFrames:
                - 0: events
                - 1: units
                - 2: eventlinks
        """

        # convert datadict into units DataFrame
        units=pd.DataFrame(self.datadict).transpose()
        self.units=units

        # convert data into events DataFrame
        data=self.data
        # if verbose:
        #     print(self.data)
        #     print(data)
        dataOut={}
        series={}
        for d in self.cols:
            # if verbose:print('col:',d)
            dataOut[d]={}
            dataOut[d+'_valtype']={}
            dataOut[d+'_errp']={}
            dataOut[d+'_errm']={}
            for e in data:
                # if verbose:
                #     print('event',e,data[e])
                if d not in data[e]:
                    continue
                if 'best' in data[e][d]:
                    dataOut[d][e]=data[e][d]['best']
                    if 'err' in data[e][d]:
                        if data[e][d]['err']=='lowerbound':
                            dataOut[d][e]=data[e][d]['best']
                            dataOut[d+'_valtype'][e]='lower'
                        elif data[e][d]['err']=='upperbound':
                            dataOut[d][e]=data[e][d]['best']
                            dataOut[d+'_valtype'][e]='upper'
                        elif type(data[e][d]['err']) == list:
                            dataOut[d+'_valtype'][e]='bestfit'
                            dataOut[d+'_errp'][e]=np.max(data[e][d]['err'])
                            dataOut[d+'_errm'][e]=np.min(data[e][d]['err'])
                    else:
                        dataOut[d+'_valtype'][e]='value'
                elif 'lim' in data[e][d]:
                    dataOut[d+'_valtype'][e]='range'
                    dataOut[d+'_errp'][e]=data[e][d]['lim'][0]
                    dataOut[d+'_errm'][e]=data[e][d]['lim'][1]
                elif 'lower' in data[e][d]:
                    dataOut[d+'_valtype'][e]='lower'
                    dataOut[d][e]=data[e][d]['lower']
                elif 'upper' in data[e][d]:
                    dataOut[d][e]=data[e][d]['upper']
                    dataOut[d+'_valtype'][e]='upper'
        # convert to series
        for dOut in dataOut:
            series[dOut]=pd.Series(dataOut[dOut],index=dataOut[dOut].keys(),dtype=object)
            # rows.append(d)
        # combine into DataFrame
        events=pd.DataFrame(series).T
        self.events=events

        # convert links into eventlinks DataFrame

        linksSeries=[]
        for ev in self.links:
            for l in range(len(self.links[ev])):
                link=self.links[ev][l]
                linkOut={'event':ev}
                for r in link:
                    linkOut[r]=link[r]
                # convert to series
                linksSeries.append(pd.Series(linkOut))
        # combine into DataFrame
        eventrefs=pd.DataFrame(linksSeries).T
        self.eventrefs=eventrefs

        return(events,units,eventrefs)

    def dataframe2json(self,dataIn,unitsIn,linksIn,mode,verbose=False):
        """Convert pandas DataFrame objects into dictionaries, merging with existing. Used to read from CSV files
        Inputs:
            * dataIn [dictionary]: data for events (merge with data)
            * unitsIn [dictionary]: units information (merge with datadict)
            * linksIn [dictionary]: links information (merge with links)
            * mode [string]: Mode to use:
                - "replace": Remove and replace entire dataset from imports
                - "update": Update existing data from imports
                - "append": Append new events/parameters, leave existing events/parameters unchanged
        Outputs:
            * Tuple of Dictionary:
                - 0: data
                - 1: datadict
                - 2: links
        """
        paramsIn=list(unitsIn.keys())
        eventsIn=list(dataIn.keys())

        # eventsInDict=dataIn.to_dict(orient='index')
        # unitsInDict=unitsIn.to_dict(orient='index')
        # linksInDict=linksIn.to_dict()

        if mode=='replace':
            # remove existing dataset
            if verbose:print('Removing existing data')
            for k in list(self.datadict.keys()):
                if verbose: print('Removing parameter {}'.format(k))
                self.datadict.pop(k,None)
            for k in list(self.data.keys()):
                if verbose: print('Removing event {}'.format(k))
                self.data.pop(k,None)
            for k in list(self.links.keys()):
                if verbose: print('Removing links {}'.format(k))
                self.links.pop(k,None)

        if verbose: print('\n*** Udating parameters ***')
        # create list of parameters
        for param in unitsIn:
            # check parameters are in current database
            if param not in self.datadict:
                if verbose: print('Adding parameter {}'.format(param))
                self.datadict[param]={}
                for k in unitsIn[param]:
                    try:
                        if np.isnan(unitsIn[param][k]):
                            continue
                    except:
                        pass
                    if not unitsIn[param][k]=='':
                        self.datadict[param][k]=unitsIn[param][k]
            elif mode=="append":
                # don't update existing parameter
                if verbose: print('Mode=append. Skipping parameter {}'.format(param))
                pass
            else:
                # replace existing parameter
                if verbose: print('Updating parameter {}'.format(param))
                self.datadict[param]={}
                for k in unitsIn[param]:
                    try:
                        if np.isnan(unitsIn[param][k]):
                            continue
                    except:
                        pass
                    if not unitsIn[param][k]=='':
                        self.datadict[param][k]=unitsIn[param][k]

        if verbose: print('\n*** Updating event data ***')
        # update events dictionary
        for ev in eventsIn:
            if ev not in self.data:
                # event is new
                if verbose: print('Adding event %s'%(ev))
                event=dataframe2jsonEvent(dataIn[ev],paramsIn,verbose=verbose)
                self.data[ev]=event
            else:
                # event exists in data
                if mode=="append":
                    # don't update existing events
                    if verbose: print('Mode=append. Skipping event {}'.format(ev))
                    pass
                else:
                    # update existing event
                    if verbose: print('Merging events {}'.format(ev))
                    event=dataframe2jsonEvent(dataIn[ev],paramsIn,verbose=verbose)
                    for el in event:
                        if el not in self.data[ev]:
                            if verbose: print ('Adding value {}'.format(el))
                            self.data[ev][el]=event[el]
                        else:
                            if verbose: print('Merging element {} for event {}'.format(el,ev))
                            if not compareElement(self.data[ev][el],event[el],verbose=verbose):
                                # if verbose: print('Updating value {}'.format(el))
                                self.data[ev][el]=event[el]
                            else:
                                pass
                                # if verbose: print('Keeping value {}'.format(el))
                    for el in list(self.data[ev].keys()):
                        if el not in event:
                            # element existed, but not in input.
                            if mode=='replace':
                                # Remove from dictionary
                                self.data[ev].pop(el,None)

        if verbose: print('\n*** Udating links ***')

        # get current links
        oldLinks=[]
        newLinks=[]
        skipLinks=[]
        for ev in self.links:
            oldLinks.append(ev)

        # update links
        for l in linksIn:
            link=linksIn[l]
            ev=link['event']
            if ev not in oldLinks and ev not in newLinks:
                # event is new in links
                if verbose: print('Adding links for event {}'.format(ev))
                self.links[ev]=[]

            if mode=="append":
                # don't update existing events
                if ev not in skipLinks:
                    if verbose: print('Skipping event {}'.format(ev))
                    skipLinks.append(ev)
                pass
            else:
                if ev not in newLinks:
                    # links haven't been replaced yet for this event
                    if verbose: print('Updating links for event {}'.format(ev))
                    self.links[ev]=[]
                    newLinks.append(ev)
                # add link to links-list for event
                newLink={}
                for key in link:
                    try:
                        if np.isnan(link[key]):
                            continue
                    except:
                        pass
                    if key!='event' and link[key]!='':
                        newLink[key]=link[key]
                self.links[ev].append(newLink)

        # put back into dataframe structures
        self.json2dataframe()

        return(self.data,self.datadict,self.links)


    def getValue(self,event,param,value):
        try:
            return self.data[event][param][value]
        except:
            print('Error finding value %s for parameter %s in event %s'%(value,param,event))
            return np.NaN

    def exportJson(self,fileout,dir='',verbose=False,contents='all'):
        """Export parameters, data and links to single JSON file
        Inputs:
            * fileout [string]: filename to write all data to
            * dir [string OPTIONAL]: directory to write files to. Default=''
            * contents [string OPTIONAL]: contents of file ('all','data','links','datadict'). Default='all'
        """
        if contents=='all':
            alldata={'meta':self.meta,'datadict':self.datadict,'data':self.data,'links':self.links}
            alldata['meta']['contents']={'data':True,'links':True,'datadict':True}
        elif contents=='data':
            alldata={'meta':self.meta,'data':self.data}
            alldata['meta']['contents']={'data':True,'links':False,'datadict':False}
        elif contents=='links':
            alldata={'meta':self.meta,'links':self.links}
            alldata['meta']['contents']={'data':False,'links':True,'datadict':False}
        elif contents=='datadict':
            alldata={'meta':self.meta,'datadict':self.datadict}
            alldata['meta']['contents']={'data':False,'links':False,'datadict':True}
        if verbose: print('Writing data to JSON: {}'.format(os.path.join(dir,fileout)))
        json.dump(alldata,open(os.path.join(dir,fileout),'w'),indent=4)


        return()

    def exportCSV(self,datafileout,dictfileout=None,linksfileout=None,dir='',verbose=False,clearcols=True):
        """Export data to CSV file(s)
        Inputs:
            * datafileout [string]: filename to write events data to
            * dictfileout [string, OPTIONAL]: filename to write data dictionary to. Default: do not export
            * linksfileout [string OPTIONAL]: filename to write references to. Default: to not export
            * clearcols [boolean OPTIONAL]: Remove columns that are empty. Default=True
            * dir [string OPTIONAL]: directory to write files to. Default=''
        """
        (dataframe,units,links) = self.json2dataframe(verbose=verbose)
        if clearcols:
            if verbose:
                print('removing empty rows')
            dataframe.dropna(axis=0,how='all',inplace=True)

        if verbose: print('Writing data to CSV: {}'.format(os.path.join(dir,datafileout)))
        dataframe.transpose().to_csv(os.path.join(dir,datafileout),encoding='utf8',index_label='Event')

        if dictfileout!=None:
            if verbose: print('Writing datadict to CSV: {}'.format(os.path.join(dir,dictfileout)))
            units.to_csv(os.path.join(dir,dictfileout),encoding='utf8',index_label='Parameter')

        if linksfileout!=None:
            if verbose: print('Writing links to CSV:{}'.format(os.path.join(dir,linksfileout)))
            links.transpose().to_csv(os.path.join(dir,linksfileout),encoding='utf8',index_label='Ref')

        return()

    def exportExcel(self,fileout,dir='',verbose=False):
        """Export datadict, events data and links to CSV file(s)
        Inputs:
            * fileout [string]: filename to write all data to
            * dir [string OPTIONAL]: directory to write files to. Default=''
        """
        (dataframe,units,links) = self.json2dataframe(verbose=verbose)

        writer=pd.ExcelWriter(os.path.join(dir,fileout),engine='xlsxwriter')

        if verbose: print('Writing data to Excel: {}'.format(os.path.join(dir,fileout)))
        dataframe.transpose().to_excel(writer,sheet_name='Events',encoding='utf8',index_label='Event')

        if verbose: print('Writing datadict to Excel: {}'.format(os.path.join(dir,fileout)))
        units.to_excel(writer,sheet_name='Parameters',encoding='utf8',index_label='Parameter')

        if verbose: print('Writing links to Excel:{}'.format(os.path.join(dir,fileout)))
        links.transpose().to_excel(writer,sheet_name='Links',encoding='utf8',index_label='Ref')

        writer.save()

        return()

    def importCSV(self,datafilein,dictfilein=None,linksfilein=None,dir='',mode=None,verbose=False):
        """Read CSV file of data and replace in database
        Inputs:
            * datafilein [string]: filename of events data file to read in
            * dictfilein [string, OPTIONAL]: filename of data dictionary to read in
            * linkfile [string OPTIONAL]: filename of references csv file to read in
            * dir [string OPTIONAL]: directory to read from. Default=''
            * mode [string, OPTIONAL]: import mode [replace,update,append]
                - replace: replace entire dataset with input data
                - update: update existing events and add new events [default]
                - append: add new events, leave existing events unchanged
        """

        # set mode to default if not set
        if mode==None:
            mode='update'

        # read CSV files
        if verbose: print('Reading data from {}'.format(os.path.join(dir,datafilein)))
        datain=pd.read_csv(os.path.join(dir,datafilein),index_col=0).to_dict(orient='index')
        if dictfilein!=None:
            if verbose: print('Reading data dict from {}'.format(os.path.join(dir,dictfilein)))
            unitsin=pd.read_csv(os.path.join(dir,dictfilein),index_col=0).to_dict(orient='index')
        else:
            unitsin=self.datadict
        if linksfilein!=None:
            if verbose: print('Reading links from {}'.format(os.path.join(dir,linksfilein)))
            linksin=pd.read_csv(os.path.join(dir,linksfilein),index_col=0).to_dict(orient='index')
        else:
            linksin=self.links

        # merge imports with existing
        self.dataframe2json(datain,unitsin,linksin,mode=mode,verbose=verbose)

        return()

    def importExcel(self,filein,dir='',sheet_events='Events',sheet_dict='Parameters',sheet_links='Links',mode='update',verbose=False):
        """Read Excel file of data and replace in database
        Inputs:
            * filein [string]: filename of data file to read in from
            * sheet_events [string, OPTIONAL]: sheet name for events data. Default="Events"
            * sheet_dict [string, OPTIONAL]: sheet name for parameters data. Default="Parameters"
            * sheet_links [string, OPTIONAL]: sheet name for links data. Default="Links"
            * dir [string OPTIONAL]: directory to read from. Default=''
            * mode [string, OPTIONAL]: import mode [replace,update,append]
                - replace: replace entire dataset with input data
                - update: update existing events and add new events [default]
                - append: add new events, leave existing events unchanged
        """

        # set mode to default if not set
        if mode==None:
            mode='update'

        # read Excel file
        if verbose: print('Reading data from {}'.format(os.path.join(dir,filein)))
        datain=pd.read_excel(os.path.join(dir,filein),sheet_name='Events',index_col=0).to_dict(orient='index')

        if verbose: print('Reading data dict from {}'.format(os.path.join(dir,filein)))
        unitsin=pd.read_excel(os.path.join(dir,filein),sheet_name='Parameters',index_col=0).to_dict(orient='index')

        if verbose: print('Reading links from {}'.format(os.path.join(dir,filein)))
        linksin=pd.read_excel(os.path.join(dir,filein),sheet_name='Links',index_col=0).to_dict(orient='index')


        # merge imports with existing
        self.dataframe2json(datain,unitsin,linksin,mode=mode,verbose=verbose)

        return()

    def addRefs(self,refsFile=None,verbose=False):
        """Add references from refs file
        Inputs:
            * refsFile [string, optional]: file to read in. Default='refs.json'
            * verbose [boolean, optional]: set for verbose output. Default=False
        Outputs:
            * None
        """
        if not (refsFile):
            refsFile='refs.json'
        fileIn=os.path.join(self.dataDir,refsFile)
        try:
            refsIn=json.load(open(fileIn))
        except:
            if verbose:
                print('error loading {}',format(fileIn))
                return
        for ev in refsIn:
            if ev in self.data:
                for r in refsIn[ev]:
                    self.addLink(ev,r,verbose=verbose)
                    if r['type']=='skymap-fits':
                        self.updateMapSrc(ev,verbose=verbose)
        return

    def makeWaveforms(self,verbose=False,overwrite=False,event=None):
        """Create waveforms for events and add to database.
        Inputs:
            * overwrite [boolean, optional]: set to overwrite all waveforms, not just new ones. Default=False
            * verbose [boolean, optional]: set for verbose output. Default=False
            * event [string, optional]: set to process single event. Default=None
        Outputs:
            * None
        """
        print('*** Updating waveforms...')
        # if os.path.exists(logFile):
        #     os.remove(logFile)
        #     print('Removing log file: {}'.format(logFile))
        # else:
        #     print("Log file doesn't exist: {}".format(logFile))
        # if logFile:
        #     print('Writing Maps log to: {}'.format(logFile))
        #     logF=open(logFile,'a')

        wfDir=os.path.join(self.dataDir,'waveforms')
        wfDirFull=os.path.join(self.dataDir,'waveforms-full')
        if not os.path.exists(wfDir):
            os.mkdir(wfDir)
        if not os.path.exists(wfDirFull):
            os.mkdir(wfDirFull)
        wfs={}
        for ev in self.events:
            if event!="" and event!=None:
                if ev != event:
                    continue
                else:
                    if verbose: print("CALCULATING WAVEFORM FOR SINGLE EVENT {}".format(ev))
            # check parameters exist
            M1Param=self.getParameter(ev,'M1')
            M2Param=self.getParameter(ev,'M2')
            MchirpParam=self.getParameter(ev,'Mchirp')
            DLParam=self.getParameter(ev,'DL')
            if (M1Param)and(M2Param)and(MchirpParam)and(DLParam):
                if verbose:print('all parameters exist')
            else:
                if verbose:print('not all parameters exist [M1,M2,Mchirp,DL]:',M1Param,M2Param,MchirpParam,DLParam)
                continue
            params={}
            if 'best' in M1Param:
                params["M1"]=M1Param['best']
            if 'best' in M2Param:
                params["M2"]=M2Param['best']
            if 'best' in MchirpParam:
                params["Mchirp"]=MchirpParam['best']
            if 'best' in DLParam:
                params["DL"]=DLParam['best']
            # check parameters exist
            if ('M1' in params)and('M2' in params)and('Mchirp' in params)and('DL' in params):
                if verbose:print('all values exist')
                wfs[ev]=params
            else:
                if verbose:print('not all parameters exist:',params)
                continue

        linktxt='Simulated waveform (compressed)'
        for ev in wfs:
            wfs[ev]['wfFile']=os.path.join(wfDir,'waveform_{}_compress.txt'.format(ev))
            wfs[ev]['wfFileFull']=os.path.join(wfDirFull,'waveform_{}_full.txt'.format(ev))
            wfs[ev]['exists']=os.path.isfile(wfs[ev]['wfFile'])
            if not wfs[ev]['exists'] or overwrite:
                wfs[ev]['update']=True
            else:
                wfs[ev]['update']=False
            # link=self.getLink(ev,linktxt,srchtype='text')
            # if len(link)>0:
            #     if 'created' in link[0]:
            #         if link[0]['created']<fitsCreated:
            #             wfs[ev]['update']=True
        for ev in wfs:
            if not wfs[ev]['update']:
                continue
            m1=wfs[ev]['M1']
            m2=wfs[ev]['M2']
            mch=wfs[ev]['Mchirp']
            K0=2.7e17 #Msun^5 s^-5

            if m1+m2 > 100:
                tres=1.0/4096
                f_lower=10
            elif m1+m2 > 67:
                tres=1.0/4096
                f_lower=20
            elif m1+m2>5:
                tres=1.0/4096
                f_lower=25
            else:
                tres=1.0/8192
                f_lower=30
            wfs[ev]['fmin']=f_lower
            f30=30
            f25=25
            tmin=K0**(1./3.) * mch**(-5./3.) * f_lower**(-8./3.)
            t30=K0**(1./3.) * mch**(-5./3.) * f30**(-8./3.)
            t25=K0**(1./3.) * mch**(-5./3.) * f25**(-8./3.)
            fitparam=[0.0029658 , 0.96112625]
            tmin=tmin*(mch*fitparam[0] + fitparam[1])
            t30=t30*(mch*fitparam[0] + fitparam[1])
            t25=t25*(mch*fitparam[0] + fitparam[1])
            wfs[ev]['tmin']=tmin
            wfs[ev]['t25']=t25
            wfs[ev]['t30']=t30
            print('processing {}: {} + {} [{}] ({} MPc) at 1/{}s resolution from {}Hz [{:.2f}s from {:.2f}Hz]'.format(ev,wfs[ev]['M1'],wfs[ev]['M2'],wfs[ev]['Mchirp'],wfs[ev]['DL'],1./tres,f_lower,t30,f30))
            # print(' t(30Hz) = {:.2}'.format(tmin))

            hp,hc = get_td_waveform(approximant="SEOBNRv3_opt_rk4",
                             mass1=wfs[ev]['M1'],
                             mass2=wfs[ev]['M2'],
                             delta_t=tres,
                             f_lower=f_lower,
                             distance=wfs[ev]['DL'])
            t= hp.sample_times
            wfs[ev]['data']=Table({'t':t,'hp':hp,'hc':hc})
            # wfs[d]['data'].write('full-data/waveform_{}.csv'.format(d),format='ascii.csv',overwrite=True)
            if verbose:print('  produced {:.2f}s from {:.2f}Hz'.format(-t[0],f_lower))

            # crop waveform
            cropt=np.where(wfs[ev]['data']['t']>-wfs[ev]['t25'])[0]
            hp=wfs[ev]['data']['hp'][cropt]
            t=wfs[ev]['data']['t'][cropt]
            if verbose:
                print('{} t30={:.2f}: {} samples'.format(ev,wfs[ev]['t30'],len(t)))
                print('compressing {}'.format(ev))
            hp2=np.where(np.abs(hp)<1e-24,0,hp*1e23)
            wfs[ev]['datacomp']=Table([t,hp2],names=['t','strain*1e23'])
            for l in range(len(wfs[ev]['datacomp'])):
                wfs[ev]['datacomp'][l]['t']=round(wfs[ev]['datacomp'][l]['t'],5)
                wfs[ev]['datacomp'][l]['strain*1e23']=round(wfs[ev]['datacomp'][l]['strain*1e23'],1)
            wfs[ev]['datacomp'].write(wfs[ev]['wfFile'],format='ascii.basic',delimiter=" ",overwrite=True)

            # add link:
            if verbose: print('adding waveform link for {}'.format(ev))
            link={'url':self.rel2abs(wfs[ev]['wfFile'],data=True),'text':linktxt,
                'type':'waveform-compressed','created':Time.now().isot,'offset':float('{:.5f}'.format(-wfs[ev]['t30'])),'tmerge':0.0}
            self.addLink(ev,link)
        return
