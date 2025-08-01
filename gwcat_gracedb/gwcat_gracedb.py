import os
from ligo.gracedb.rest import GraceDb
from numpy import argmax
import json
import requests
from bs4 import BeautifulSoup
from astropy.time import Time
from astropy import units as un

def gps2obsrun(GPS):
    """Convert GPS to observing run
    Inputs:
        * GPS [float]: sGPS timestamp
    Outputs:
        * [string] observing run
    """
    obsruns={'O1':[1126051217,1137254417],
        'O2':[1164556817,1187733618],
        'O3a':[1238112018,1253977218],
        'O3b':[1256655618,1269363618],
        'O4a':[1366556418,1389452418],
        'O4b':[1396796418,1427414418],
        'O4c':[1433548818,1447545617]}
    obs=''
    for o in obsruns:
        if GPS > obsruns[o][0] and GPS < obsruns[o][1]:
            obs=o
    return obs

def gracedb2cat(gdb,verbose=False,knownEvents={},forceUpdate=False):
    """Convert GraceDB data to gwcat format.
    Inputs:
        * gwtcdata [object]: GWTC data (from getGWTC)
        * verbose [boolean, optional]: set for verbose output. Default=False
        * knownEvents [object]: timestamps of existing events. Default=empty (ignore known events)
        * forceUpdate [boolean, optional]: set to force output of all events, not just updated ones. Default=False
    Outputs:
        * [object] object containing 'data' and 'links'
    """
    catOut={}
    linksOut={}
    if 'data' in gdb:
        gdbIn=gdb['data']
    else:
        gdbIn=gdb
    for g in gdbIn:
        if 'validXML' in gdbIn[g]['meta']:
            if not gdbIn[g]['meta']['validXML']:
                if verbose:print('skipping {} due to invalid XML'.format(g))
                continue
        if g in knownEvents:
            tOld=Time(knownEvents[g])
            try:
                tNew=Time(gdbIn[g]['meta']['created_date'])
            except:
                tNew=Time.now()
            if tNew <= tOld:
                if forceUpdate:
                    if verbose:print('forcing update for {}: [{}<={}]'.format(g,tNew.isot,tOld.isot))
                else:
                    if verbose:print('no update needed for {}: [{}<={}]'.format(g,tNew.isot,tOld.isot))
                    continue
            else:
                if verbose:print('updating {}: [{}>{}]'.format(g,tNew.isot,tOld.isot))
        if verbose: print('importing GraceDB event {}'.format(g))
        # check for Retraction
        if 'xml' in gdbIn[g]:
            xml=gdbIn[g]['xml']
            if 'AlertType' in xml:
                if xml['AlertType']=='Retraction':
                    if verbose: print('skipping Retraction {}'.format(g))
                    catOut[g]={'Retracted':True}
                    continue
            if 'Significant' in xml:
                if xml['Significant']=='0':
                    if verbose: print('skipping Low-significance {}'.format(g))
                    # catOut[g]={'Significance':'Low'}
                    continue

        catOut[g]={}
        linksOut[g]=[]
        if 'superevent_id' in gdbIn[g]: catOut[g]['name']=gdbIn[g]['superevent_id']
        catOut[g]['detType']={'best':'Candidate'}
        catOut[g]['conf']={'best':'Candidate'}
        catOut[g]['catalog']='gracedb'
        if 't_0' in gdbIn[g]:
            catOut[g]['GPS']={'best':gdbIn[g]['t_0']}
            dtIn=Time(gdbIn[g]['t_0'],format='gps')
            dtOut=Time(dtIn,format='iso').isot
            catOut[g]['UTC']={'best':dtOut}
            catOut[g]['obsrun']={'best':gps2obsrun(gdbIn[g]['t_0'])}
        if 'far' in gdbIn[g]:
            catOut[g]['FAR']={'best':gdbIn[g]['far']*un.year.to('s')}

        # if 'hdr' in gdbIn[g]:
        #     hdr=gdbIn[g]['hdr']
        #     if 'DISTMEAN' in hdr and 'DISTSTD' in hdr:
        #         catOut[g]['DL']={
        #             'best':float(hdr['DISTMEAN']['value']),
        #             'err':[-float(hdr['DISTSTD']['value']),float(hdr['DISTSTD']['value'])]
        #         }
            # if 'DATE' in hdr:

        if 'xml' in gdbIn[g]:
            xml=gdbIn[g]['xml']
            if 'Pipeline' in xml:
                catOut[g]['FAR']['fartype']=xml['Pipeline']
            if 'Instruments' in xml:
                instOut=''
                inst=xml['Instruments']
                if inst.find('H1')>=0: instOut=instOut+'H'
                if inst.find('L1')>=0: instOut=instOut+'L'
                if inst.find('V1')>=0: instOut=instOut+'V'
                catOut[g]['net']={'best':instOut}
            if 'Classification' in xml:
                catOut[g]['objType']={'prob':{}}
                bestObjType=''
                bestProb=0
                for t in xml['Classification']:
                    prob=float(xml['Classification'][t])
                    catOut[g]['objType']['prob'][t]=prob
                    if prob > bestProb:
                        bestObjType=t
                        bestProb=float(xml['Classification'][t])
                catOut[g]['objType']['best']=bestObjType
            if 'Significant' in xml:
                if float(xml['Significant'])==1:
                    sig='High'
                else:
                    sig='Low'
                catOut[g]['Significance']={'best':sig}
        if 'meta' in gdbIn[g]:
            catOut[g]['meta']=gdbIn[g]['meta']
        # update links
        if 'links' in gdbIn[g]:
            if 'self' in gdbIn[g]['links']:
                opendata={'url':gdbIn[g]['links']['self'],
                    'text':'GraceDB data',
                    'type':'open-data'}
                linksOut[g].append(opendata)
        if 'mapfile' in gdbIn[g]:
            skymap={'url':gdbIn[g]['mapfile'][1],
                'text':'Sky Map',
                'type':'skymap-fits'}
            linksOut[g].append(skymap)

    return {'data':catOut,'links':linksOut}

def getSuperevents(export=False,dirOut=None,fileOut=None,indent=2,verbose=False,knownEvents={},forceUpdate=False,datelim=999,logFile=None,highSigOnly=False):
    """Get GWTC from GWOSC json files and add parameters.
    Inputs:
        * export [boolean, optional]: set for export to JSON file output. Default=False
        * dirOut [string, optional]: directory to export to (if export=True). Default='data/'
        * fileOut [string, optional]: file to export to (if export=True). Default='GWTC.json'
        * indent [integer, optional]: json indent for exported file. Default=2
        * verbose [boolean, optional]: set for verbose output. Default=False
        * knownEvents [object]: timestamps of existing events. Default=empty (ignore known events)
        * forceUpdate [boolean, optional]: set to force output of all events, not just updated ones. Default=False
        * datelim [integer, optional]: number of days old to skip events for. Default=999
        * lodGile [string, optional]: logFile to output logging to. Default=None (no logging)
        * highSigOnly [boolean, optional]: set to remove low significance events. Default=False (include low significance events)
    Outputs:
        * [object] object containing data (can be read by gwcatpy.importGWTC)
    """
    service_url = 'https://gracedb.ligo.org/api/'
    if verbose: print('Retrieving GraceDB data from {}'.format(service_url))
    client = GraceDb(service_url,force_noauth=True)
    if verbose: print('Limiting to {} days'.format(datelim))
    # Retrieve an iterator for events matching a query.
    events = client.superevents('far < 1.0e-6')
    # if verbose: print('retrieved {} events'.format(len(events)))
    # For each event in the search results, add the graceid
    # and chirp mass to a dictionary.
    results = {}
    links = {}

    if logFile:
        if os.path.exists(logFile):
            os.remove(logFile)
            print('Removing log file: {}'.format(logFile))
        else:
            print("Log file doesn't exist: {}".format(logFile))
        print('Writing GraceDB log to: {}'.format(logFile))
        logF=open(logFile,'a')

    for event in events:
        sid = event['superevent_id']
        tEvent=Time(event['t_0'],format='gps')
        tNow=Time.now()
        dtEvent=(tNow-tEvent).jd
        if dtEvent>datelim:
            print('Too old ({} days). Skipping {}'.format(dtEvent,sid))
            continue
        # if sid in knownEvents:
        #     tOld=Time(knownEvents[sid])
        #     try:
        #         tNew=Time(cdate)
        #     except:
        #         tNew=Time.now()
        #     if tNew <= tOld:
        #         if verbose:print('no import needed for {}: [{}<={}]'.format(sid,tNew.isot,tOld.isot))
        #         continue
        #     else:
        #         if verbose:print('importing for {}: [{}>{}]'.format(sid,tNew.isot,tOld.isot))
        evOut=event
        evOut['meta']={'retrieved':Time.now().isot,'src':service_url}

        voreq=client.voevents(sid)
        voevents=json.loads(voreq.content)
        evOut['voevents']=voevents

        volist=voevents['voevents']
        good_voevents = [voe for voe in volist if voe['voevent_type'] != 'RE']
        retraction_list = [voe for voe in volist if voe['voevent_type'] == 'RE']
        if len(retraction_list)>0:
            print('Event {} retracted. Skipping'.format(sid))
            cdate=Time(' '.join(retraction_list[-1]['created'].split(' ')[0:2]))
            evOut['meta']['type']='Retraction'
        else:
            evOut['meta']['type']='Candidate'

            Ngood=len(good_voevents)
            validXML=False
            nvo=Ngood
            while validXML==False and nvo>0:
                thisvo=good_voevents[nvo-1]
                cdate=Time(' '.join(thisvo['created'].split(' ')[0:2]))
                if sid in knownEvents:
                    tOld=Time(knownEvents[sid])
                    tNew=cdate
                    if tNew <= tOld:
                        if forceUpdate:
                            if verbose:print('forcing update for {}: [{}<={}]'.format(sid,tNew.isot,tOld.isot))
                            update=True
                        else:
                            if verbose:print('no update needed for {}: [{}<={}]'.format(sid, tNew.isot,tOld.isot))
                            update=False
                            validXML=True
                    else:
                        if verbose:print('getting files for {}: [{}>{}]'.format(sid,tNew.isot,tOld.isot))
                        update=True
                else:
                    update=True

                thisvoe_url = thisvo['links']['file']
                vonum=thisvo['N']

                evOut['xmlfile']=[os.path.split(thisvoe_url)[-1],thisvoe_url]
                if update:
                    # parse XML
                    xml={}
                    if verbose: print('  parsing {}'.format(evOut['xmlfile'][0]))
                    xmlurl=evOut['xmlfile'][1]
                    xmlreq=requests.get(xmlurl)
                    soup=BeautifulSoup(xmlreq.text,'lxml')
                    try:
                        params=soup.what.find_all('param',recursive=False)
                        validXML=True
                    except:
                        print('problem with {}: {}'.format(sid,evOut['xmlfile'][0]))
                    evOut['validXML']=validXML
                    if validXML:
                        for p in params:
                            xml[p.attrs['name']]=p.attrs['value']
                        groups=soup.what.find_all('group',recursive=False)
                        for g in groups:
                            gt=g.attrs['type']
                            xml[gt]={}
                            gparams=g.find_all('param',recursice=False)
                            for gp in gparams:
                                xml[gt][gp.attrs['name']]=gp.attrs['value']
                        if 'GW_SKYMAP' in xml:
                            if 'skymap_fits' in xml['GW_SKYMAP']:
                                mapfile=xml['GW_SKYMAP']['skymap_fits']
                                evOut['mapfile']=[os.path.split(mapfile)[-1],mapfile]
                        if 'Significant' in xml:
                            if xml.get('Significant')=="0":
                                if highSigOnly:
                                    evOut['meta']['include']=False
                                    print('Not including low-significance event {} [{}].'.format(sid,xml.get('Significant')))
                                else:
                                    print('Including low-significance event {} [{}].'.format(sid,xml.get('Significant')))
                                    evOut['meta']['include']=True
                            else:
                                print('Including high-significance event {} [{}].'.format(sid,xml.get('Significant')))
                                evOut['meta']['include']=True
                        else:
                            print('Including unknown significance event {} [{}].'.format(sid,xml.get('Significant')))
                            evOut['meta']['include']=True
                        if evOut['meta']['include']:
                            if logFile:
                                logF.write(sid+'\n')
                        evOut['xml']=xml
                nvo-=1
            evOut['meta']['validXML']=validXML
                        # create meta data

        evOut['meta']['created_date']=cdate.isot

        results[sid]=evOut

    if logFile:
        logF.close()
    if verbose: print('Retrieved data for {} events'.format(len(results)))

    cat={'meta':{'retrieved':Time.now().isot,'src':service_url},'data':results}

    if export:
        if dirOut==None:
            dirOut='../../data/'
        if fileOut==None:
            fileOut='gracedb.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(cat,fOut,indent=indent)
        fOut.close()

    return(cat)

# def getMap(sid,fileout=None,dirOut=None,getLAL=False):
#     if getLAL:
#         filename = 'bayestar.fits'
#     else:
#         filename = 'LALInference.fits'
#     if fileout==None:
#         outFilename = '{}_{}'.format(sid,filename)
#     if dirOut==None:
#         dirOut='../../data/'
#     print('downloading {} for superevent {}'.format(filename,sid))
#     clFits=GraceDbBasic(service_url)
#     fout=open(os.path.join(dirOut,outFilename),'wb')
#     r = clFits.files(sid,filename)
#     fout.write(r.read())
#     fout.close()