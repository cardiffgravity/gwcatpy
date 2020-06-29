import json
import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.io import fits
import healpy as hp
import h5py
import os
import requests
from . import gwosc
from . import gracedb
from . import plotloc
from astropy.table import Table
import astropy_healpix as ah
from pycbc.waveform import get_td_waveform
# import gracedb
# import gwosc

def json2jsonp(fileIn,fileOut=None,verbose=True):
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

def compareElement(el1,el2,verbose=False):
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
    # convert dataframe element to json
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
    if prec:
        precstr='{:.'+'{}'.format(prec)+'g}'
    else:precstr='{}'
    return float(precstr.format(v))

class GWCat(object):
    def __init__(self,fileIn='../data/events.json',statusFile='status.json',
        dataDir='data/',baseurl='https://data.cardiffgravity.org/gwcat-data/',verbose=False):
        """Initialise catalogue from input file
        Input: fileIn [string, OPTIONAL]: filename to read data from
        """
        self.dataDir=dataDir
        if baseurl[-1]!='/':baseurl=baseurl+'/'
        self.baseurl=baseurl
        self.statusFile=os.path.join(dataDir,statusFile)
        
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
        if ev in self.data:
            return(self.data[ev])
        else:
            return
    def getParameter(self,ev,param):
        pOut=None
        if ev in self.data:
            if param in self.data[ev]:
                pOut=self.data[ev][param]
        return(pOut)
    def getTimestamps(self):
        evTimes={}
        for ev in self.data:
            try:
                evTimes[ev]=Time(self.data[ev]['meta']['created_date'])
            except:
                # print('no created_date for {}'.format(ev))
                pass
        return(evTimes)

    def getStatus(self):
        try:
            self.status=json.load(open(self.statusFile))
        except:
            print("Initialising status")
            self.status={}
            json.dump(self.status,open(self.statusFile,'w'))
        return

    def saveStatus(self):
        json.dump(self.status,open(self.statusFile,'w'),indent=4)
        return

    def updateStatus(self,ev,desc='',statusIn=None,verbose=False):
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

    def importManual(self,manIn,verbose=False):
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
        
    def matchGraceDB(self,verbose=False):
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
                maplink=self.getLink(ev,'skymap-fits')
                mapgdb=self.getLink(matchname,'skymap-fits')
                if len(maplink)==0 and len(mapgdb)>0:
                    if verbose:print('copying GraceDB map link from {} to {}: {}'.format(matchname,ev,mapgdb[0]))
                    self.addLink(ev,mapgdb[0],verbose=verbose)
            
        return
        
    def importGWTC1(self,gwtc1In,verbose=False):
        print('*** Importing GWTC-1...')
        catData=gwosc.gwtc1_to_cat(gwtc1In,verbose=verbose)
        for g in catData['data']:
            # get old metadata
            dmeta={}
            if g in self.data:
                if 'meta' in self.data[g]:
                    dmeta=self.data[g]['meta']
            self.data[g]=catData['data'][g]
            # update metadata
            for m in catData['data'][g]['meta']:
                dmeta[m]=catData['data'][g]['meta'][m]
            self.data[g]['meta']=dmeta
            if g in catData['links']:
                for l in catData['links'][g]:
                    self.addLink(g,l,verbose=verbose)
            self.updateStatus(g,verbose=verbose,desc='GWTC-1 import')
        for ev in catData['links']:
            self.updateMapSrc(ev,verbose=True)
            self.updateStatus(ev,verbose=verbose,desc='Map src')
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        if not 'GWTC-1' in self.meta:
            self.meta['GWTC-1']={}
        for m in gwtc1In['meta']:
            self.meta['GWTC-1'][m]=gwtc1In['meta'][m]
        return
    
    # backwards compatibility
    def importGwosc(self,gwoscIn,verbose=False):
        print('***WARNING: importGwosc replaced by importGWTC1***')
        self.importGWTC1(gwoscIn,verbose=verbose)
    
    def importGraceDB(self,gracedbIn,verbose=False,forceUpdate=False):
        print('*** Importing GraceDB...')
        evTimes=self.getTimestamps()
        gdb=gracedb.gracedb2cat(gracedbIn['data'],verbose=verbose,
            knownEvents=evTimes,forceUpdate=forceUpdate)
        for g in gdb['data']:
            if gdb['data'][g]['meta'].get('type')=='Retraction':
                if verbose:print('Skipping retracted event {}'.format(g))
                if g in self.data:
                    if verbose:print('Removing data for retracted event {}'.format(g))
                    self.data.pop(g)
                if g in self.links:
                    if verbose:print('Removing links for retracted event {}'.format(g))
                    self.links.pop(g)
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
            self.updateMapSrc(ev)
            self.updateStatus(ev,verbose=verbose,desc='Map src')
        self.evTimes = self.getTimestamps()
        self.json2dataframe(verbose=verbose)
        if not 'graceDB' in self.meta:
            self.meta['graceDB']={}
        for m in gracedbIn['meta']:
            self.meta['graceDB'][m]=gracedbIn['meta'][m]
        return

    def updateH5Src(self,ev,verbose=False):
        ldat=self.getLink(ev,'data-file',verbose=verbose)
        if len(ldat)>1 and verbose:
            print('Warning: more than one data file link for {}'.format(ev))
        if len(ldat)==0 and verbose:
            print('Warning: no data file link for {}'.format(ev))
        if len(ldat)>0:
            l=ldat[0]
            tmpDir=os.path.join(self.dataDir,'tmp')
            try:
                # self.data[ev]['meta']['h5datesrc']=Time.now().isot
                self.data[ev]['meta']['h5urlsrc']=l['url']
                if verbose:print('Data file loaded for {}:'.format(ev),l['url'])
            except:
                print('WARNING: Error loading data file for {}:'.format(ev),l)

    def updateH5(self,verbose=False,forceUpdate=False):
        print('*** Updating data files...')
        for ev in self.data:
            # compare map creation dates
            try:
                h5datelocal=Time(self.status[ev]['h5datelocal'])
                h5datesrc=Time(self.status[ev]['h5datesrc'])
                h5file=self.status[ev]['h5urllocal']
                if not os.path.isfile(h5file):
                    if verbose:print('No file. Need to re-download data file for {}'.format(ev))
                    updateH5=True
                elif h5datesrc>h5datelocal:
                    if verbose:
                        print('Newer file. Need to re-download data file for {}'.format(ev))
                        print('src',h5datesrc)
                        print('local',h5datelocal)
                    updateH5=True
                else:
                    if verbose:print('Older file. Do not need to re-download data file for {}'.format(ev))
                    updateH5=False
            except:
                updateH5=True
            if updateH5 or forceUpdate:
                if verbose:print('Updating data file for {}'.format(ev))
                stat=self.getH5(ev,verbose=verbose)
                if stat!=0:
                    self.getH5Local(ev,verbose=verbose)
            else:
                if verbose:print('No data file update required for {}'.format(ev))
            
            try:
                h5File=self.status[ev]['h5urllocal']
            except:
                if verbose:print('no local data file for {}'.format(ev))
                continue
            newparams=self.getH5Params(ev,verbose=verbose)
            for p in newparams:
                if verbose:
                    if (p in self.data[ev]):
                        print('replacing {}[{}]:{}->{}'.format(ev,p,self.data[ev][p],newparams[p]))
                else:
                    print('adding {}[{}]:{}'.format(ev,p,newparams[p]))
                self.data[ev][p]=newparams[p]
        return

    def getH5(self,ev,verbose=False):
        import h5py
        h5Dir=os.path.join(self.dataDir,'h5')
        if not os.path.exists(h5Dir):
            # create directory
            os.mkdir(h5Dir)
            print('Created directory: {}'.format(h5Dir))
        ldat=self.getLink(ev,'data-file')
        if len(ldat)==0:
            if verbose: print('ERROR: no data file link for {}',format(ev))
            return
        if len(ldat)>1:
            if verbose: print('WARNING: more than one data file link for {}',format(ev))
        url=ldat[0]['url']
        if verbose: print('Downloading data file for {} from {}'.format(ev,url))
        srcfile=os.path.split(url)[-1]
        h5req=requests.get(url)
        if h5req.ok:
            try:
                # if url.find('.fits.gz')>=0:
                #     fileext='fits.gz'
                # else:
                #     fileext='fits'
                if srcfile.find(ev)<0:
                    h5File=os.path.join(h5Dir,'{}_{}'.format(ev,srcfile))
                else:
                    h5File=os.path.join(h5Dir,srcfile)
                fOut=open(h5File,'wb')
                fOut.write(h5req.content)
                fOut.close()
                if verbose:print('Data file downloaded')
            except:
                print('ERROR: Problem loading/saving data file:',h5req.status_code)
                return h5req
            try:
                testread=h5py.File(h5File)
                if verbose: print('Valid h5 file for {}: {}'.format(ev,h5File))
                # hdr=fits.getheader(fitsFile,ext=1)
                stat={'h5urllocal':h5File,
                    'h5datelocal':Time.now().isot,
                    'h5datesrc':Time.now().isot}
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
        
    def getH5Local(self,ev,verbose=False):
        ldat=self.getLink(ev,'data-local')
        if len(ldat)==0:
            if verbose: print('ERROR: no local data file link for {}',format(ev))
            return
        if len(ldat)>1:
            if verbose: print('WARNING: more than one local data file link for {}',format(ev))
        h5File=ldat[0]['url']
        if verbose: print('Using local h5 file for {}: {}'.format(ev,h5File))
        try:
            testread=h5py.File(h5File)
            if verbose: print('Valid h5 file for {}: {}'.format(ev,h5File))
            # hdr=fits.getheader(fitsFile,ext=1)
            stat={'h5urllocal':h5File,
                'h5datelocal':Time.now().isot,
                'h5datesrc':Time.now().isot}
            self.data[ev]['meta']['h5urllocal']=h5File
            self.data[ev]['meta']['h5datelocal']=Time.now().isot
            self.data[ev]['meta']['h5datesrc']=Time.now().isot
            self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='data file local')
        except:
            print('ERROR: Problem opening local data file for {}:'.format(ev),h5File)
            return
        return
        
    def getH5Params(self,ev,verbose=False):
        if verbose:print('getting H5 parameters for {}'.format(ev))
        m1check=self.getParameter(ev,'M1')
        approx=self.getParameter(ev,'approximant')
        if not (m1check):
            if verbose:print('no M1 parameter for {}'.format(ev))
            return({})
        else:
            m1check=m1check['best']
            if verbose:print('M1 parameter for {}:'.format(ev),m1check)
        try:
            h5File=self.status[ev]['h5urllocal']
        except:
            if verbose:print('no local data file for {}'.format(ev))
            return({})
        newparams=gwosc.geth5params(h5File,pcheck={'M1':m1check},approx=approx,datadict=self.datadict,verbose=verbose)
        if verbose:print(newparams)
        return(newparams)
        
    def removeCandidates(self,verbose=False):
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
                if verbose: print('Removing {}'.format(evgdb))
                self.data.pop(remev)
            if remev in self.links:
                self.links.pop(remev)
        return

    def updateMaps(self,verbose=False,forceUpdate=False):
        print('*** Updating maps...')
        for ev in self.data:
            # compare map creation dates
            try:
                mapdatelocal=Time(self.status[ev]['mapdatelocal'])
                mapdatesrc=Time(self.status[ev]['mapdatesrc'])
                mapfile=self.status[ev]['mapurllocal']
                if not os.path.isfile(mapfile):
                    if verbose:print('Need to re-download map for {}'.format(ev))
                    updateMap=True
                elif mapdatesrc>mapdatelocal:
                    updateMap=True
                else:
                    updateMap=False
            except:
                updateMap=True
            if updateMap or forceUpdate:
                self.getMap(ev,verbose=verbose)
                if verbose:print('Updating map for {}'.format(ev))
                self.calcAreas(ev,verbose=verbose)
            else:
                if verbose:print('No map update required for {}'.format(ev))
                if not 'deltaOmega' in self.data[ev]:
                    self.calcAreas(ev,verbose=verbose)
            # fitsFile=self.status[ev]['mapdatelocal']
        return

    def getMap(self,ev,verbose=False):
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
        srcfile=os.path.split(url)[-1]
        mapreq=requests.get(url)
        if mapreq.ok:
            try:
                # if url.find('.fits.gz')>=0:
                #     fileext='fits.gz'
                # else:
                #     fileext='fits'
                if srcfile.find(ev)<0:
                    fitsFile=os.path.join(fitsDir,'{}_{}'.format(ev,srcfile))
                else:
                    fitsFile=os.path.join(fitsDir,srcfile)
                fOut=open(fitsFile,'wb')
                fOut.write(mapreq.content)
                fOut.close()
                if verbose:print('Map downloaded')
            except:
                print('ERROR: Problem loading/saving map:',mapreq.status_code)
                return mapreq
            try:
                hdr=fits.getheader(fitsFile,ext=1)
                stat={'mapurllocal':fitsFile,
                    'mapdatelocal':hdr['DATE']}
                self.data[ev]['meta']['mapurllocal']=fitsFile
                self.data[ev]['meta']['mapdatelocal']=hdr['DATE']
                self.updateStatus(ev,statusIn=stat,verbose=verbose,desc='maplocal')
            except:
                print('ERROR: Problem opening fits file for {}:'.format(ev),fitsFile)
                return
            self.calcAreas(ev,verbose=verbose)
        else:
            print('ERROR: Problem loading map:',mapreq.status_code)
            return mapreq.status_code
        return hdr

    def calcAreas(self,ev,verbose=False):
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

    def rel2abs(self,rel,url=None):
        if url==None:
            url=self.baseurl
        return(url + rel)

    def plotMapPngs(self,overwrite=False,verbose=False,logFile=None):
        print('*** Updating plots...')
        if os.path.exists(logFile):
            os.remove(logFile)
            print('Removing log file: {}'.format(logFile))
        else:
            print("Log file doesn't exist: {}".format(logFile))
        if logFile:
            print('Writing Maps log to: {}'.format(logFile))
            logF=open(logFile,'a')

        pngDir=os.path.join(self.dataDir,'png')
        gravDir=os.path.join(self.dataDir,'gravoscope')
        dataDir=os.path.join(self.dataDir,'fits')
        if not os.path.exists(pngDir):
            os.mkdir(pngDir)
        if not os.path.exists(gravDir):
            os.mkdir(gravDir)
        for ev in self.events:
            if not 'mapurlsrc' in self.status[ev]:
                if verbose:
                    print('no mapurlsrc for {}'.format(ev))
                continue
            if not 'mapurllocal' in self.status[ev]:
                self.getMap(ev,verbose=verbose)
            filename=self.status[ev]['mapurllocal']
            # if verbose:print('plotting maps at {}'.format(filename))
            fitsCreated=Time(self.status[ev]['mapdatelocal'])
            srcfile=os.path.split(self.status[ev]['mapurlsrc'])[-1]
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
            nUpdate=0
            for p in plots:
                plots[p]['pngFile']=os.path.join(pngDir,'{}_{}.png'.format(ev,p))
                plots[p]['thumbFile']=os.path.join(pngDir,'{}_{}.thumb.png'.format(ev,p))
                plots[p]['exists']=os.path.isfile(plots[p]['pngFile'])
                if not plots[p]['exists'] or overwrite:
                    plots[p]['update']=True
                else:
                    plots[p]['update']=False
                link=self.getLink(ev,plots[p]['linktxt'],srchtype='text')
                if len(link)>0:
                    if 'created' in link[0]:
                        if link[0]['created']<fitsCreated:
                            plots[p]['update']=True
                if plots[p]['update']: nUpdate+=1

            if nUpdate==0:
                if verbose:print('all plots exist for {}'.format(ev))
                mapread=False
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
                        if verbose: print('skipping plotting {} map'.format(pp['linktxt']))
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
                        self.addLink(ev,{'url':self.rel2abs(pp['pngFile']),'text':pp['linktxt'],
                            'type':'skymap-plot','created':Time.now().isot})
                        self.addLink(ev,{'url':self.rel2abs(pp['thumbFile']),'text':pp['linktxt'],
                            'type':'skymap-thumbnail','created':Time.now().isot})

            res=8
            gravNpix=int(8*1024)
            updateGrav=False
            gravFile=os.path.join(gravDir,'{}_{}.png'.format(ev,gravNpix))
            if not os.path.isfile(gravFile):
                updateGrav=True
            gravLinktxt='Skymap (no annotations)'
            gravLink=self.getLink(ev,gravLinktxt,srchtype='text')
            if len(gravLink)>0:
                if 'created' in gravLink[0]:
                    if link[0]['created']<fitsCreated:
                        updateGrav=True
            if updateGrav:
                if not mapread:
                    map=plotloc.read_map(filename,verbose=verbose)
                if verbose:
                    print('plotting Gravoscope for {} ({}x{})'.format(ev,gravNpix,int(gravNpix/2)))
                plotloc.plotGravoscope(mapIn=map,pngOut=gravFile,verbose=verbose,res=res)
                self.addLink(ev,{'url':self.rel2abs(gravFile),'text':gravLinktxt,
                    'type':'skymap-plain','created':Time.now().isot})

        if logFile:
            logF.close()
        return

    def makeGravoscopeTilesPerl(self,overwrite=False,verbose=False):

        gravDir=os.path.join(self.dataDir,'gravoscope')
        dataDir=os.path.join(self.dataDir,'fits')
        for ev in self.events:
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

        gravDir=os.path.join(self.dataDir,'gravoscope')
        for ev in self.events:
            tilesDir=os.path.join(gravDir,'{}-tiles'.format(ev))
            if not os.path.exists(tilesDir):
                os.mkdir(tilesDir)
            fitsCreated=Time(self.status[ev]['mapdatelocal'])
            filename=self.status[ev]['mapurllocal']
            gravLinktxt='Gravoscope tileset'
            tilesLinktype='gravoscope-tiles'
            tilesLink=self.getLink(ev,tilesLinktype,srchtype='type')
            tilesLinkIdx=self.getLink(ev,tilesLinktype,srchtype='type',retIdx=True)
            updateTiles=False
            if len(tilesLink)>0:
                # link is present
                if 'created' in tilesLink[0]:
                    timeTiles=Time(tilesLink[0]['created'])
                    if timeTiles.gps < fitsCreated.gps-1:
                        # tiles are older than map
                        if verbose: print('tiles for {} exists, but are older than map. {} < {}; {} < {}'.format(ev,timeTiles,fitsCreated,timeTiles.gps,fitsCreated.gps))
                        updateTiles=True
                    else:
                        # tiles are newer than map
                        if verbose: print('tiles for {} exists, and are newer than map. {} > {}; {} > {}'.format(ev,timeTiles,fitsCreated,timeTiles.gps,fitsCreated.gps))
                else:
                    # no created date. remove link and regenerate tiles (just in case)
                    if verbose: print('No created date for tiles link for {}. Regenerating all tiles'.format(ev))
                    self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                        'type':'gravoscope-tiles','created':fitsCreated.isot})
                    updateTiles=True
                if updateLink:
                    # update link (but don't necessarily recreate tiles)
                    if verbose: print('replacing old link for Gravoscope tileset for {}'.format(ev))
                    self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                        'type':'gravoscope-tiles','created':fitsCreated.isot})
            else:
                # no link. Need to add link and regenerate tiles
                if verbose: print('adding tiles link for Gravoscope tileset for {}'.format(ev))
                self.addLink(ev,{'url':self.rel2abs(tilesDir,url=tilesurl),'text':gravLinktxt,
                    'type':'gravoscope-tiles','created':fitsCreated.isot})
                updateTiles=True
            tileFile=os.path.join(gravDir,'{}-tiles/{}.png'.format(ev,'ttrtttttt'[0:maxres+1]))
            if not os.path.isfile(tileFile):
                if verbose: print('file {} does not exist'.format(tileFile))
                updateTiles=True
            elif overwrite:
                if verbose: print('file {} exists, but overwriting'.format(tileFile))
                updateTiles=True
            else:
                if verbose: print('file {} exists'.format(tileFile))
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

    def getLink(self,ev,srchtxt,srchtype='type',verbose=False,retIdx=False):
        if not ev in self.links:
            return []
        lOut=[]
        if len (self.links[ev])>0:
            for ol in range(len(self.links[ev])):
                if self.links[ev][ol][srchtype]==srchtxt:
                    if retIdx:
                        lOut.append(ol)
                    else:
                        lOut.append(self.links[ev][ol])
        if verbose:
            print('found {} links for {} with {}=={}'.format(len(lOut),ev,srchtype,srchtxt))
        return(lOut)

    def addLink(self,ev,link,replace=True,verbose=False):
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
            series[dOut]=pd.Series(dataOut[dOut],index=dataOut[dOut].keys())
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

    def exportJson(self,fileout,dir='',verbose=False):
        """Export parameters, data and links to single JSON file
        Inputs:
            * fileout [string]: filename to write all data to
            * dir [string OPTIONAL]: directory to write files to. Default=''
        """
        alldata={'meta':self.meta,'datadict':self.datadict,'data':self.data,'links':self.links}
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

    def addrefs(self,verbose=False):
        fileIn=os.path.join(self.dataDir,'refs.json')
        try:
            refsIn=json.load(open(fileIn))
        except:
            if verbose:
                print('error loading {}',format(fileIn))
                return
        for ev in refsIn:
            if ev in self.events:
                for r in refsIn[ev]:
                    self.addLink(ev,r,verbose=verbose)
        return
    
    def makeWaveforms(self,verbose=False,overwrite=False):
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

            if m1+m2 > 67:
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
            link={'url':self.rel2abs(wfs[ev]['wfFile']),'text':linktxt,
                'type':'waveform-compressed','created':Time.now().isot,'offset':float('{:.5f}'.format(-wfs[ev]['t30'])),'tmerge':0.0}
            self.addLink(ev,link)
        return