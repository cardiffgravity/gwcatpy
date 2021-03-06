from astropy.time import Time
import os
import numpy as np

def gps2obsrun(GPS):
    """Convert GPS to observing run
    Inputs:
        * GPS [float]: sGPS timestamp
    Outputs:
        * [string] observing run
    """
    obsruns={'O1':[1126051217,1137254417],
        'O2':[1164556817,1187733618],
        'O3':[1238112018,1269363618]}
    obs=''
    for o in obsruns:
        if GPS > obsruns[o][0] and GPS < obsruns[o][1]:
            obs=o
    return obs

def gwosc2cat(gwosc,verbose=False):
    """OBSOLETE - passed to gwtc_to_cat
    """
    print('***WARNING gwosc2cat replaced by gwtc_to_cat***')
    out=gwct_to_cat(gwosc,verbose=False)
    return(out)

def getGwosc(url='',verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False):
    """OBSOLETE - passed to getGWTC
    """
    print('***WARNING getGwosc replaced by getGWTC1***')
    out=getGWTC(url=url,verbose=verbose,export=export,dirOut=dirOut,fileOut=fileOut,indent=indent,triggers=triggers)
    return(out)

def getGWTC(url='',useLocal=False,verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False,devMode=False,sess=None):
    """Get GWTC from GWOSC json files and add parameters.
    Inputs:
        * url [string, optional]: URL to download from. Default=https://www.gw-openscience.org/eventapi/json/GWTC/
        * useLocal [boolean, optional]: Set to use local files (default URL=data/local-mirror/GWTC-dev.json)
        * verbose [boolean, optional]: set for verbose output. Default=False
        * export [boolean, optional]: set for export to JSON file output. Default=False
        * dirOut [string, optional]: directory to export to (if export=True). Default='data/'
        * fileOut [string, optional]: file to export to (if export=True). Default='GWTC.json'
        * indent [integer, optional]: json indent for exported file. Default=2
        * triggers [boolean, optional]: OBSOLETE
        * devMode [boolean, optional]: set to ciecplib to access data, and use dev link as default (https://openscience-dev.ligo.caltech.edu/eventapi/json/GWTC/)
        * sess [ciecplib.Session, optional]: ciecplib session to use in dev mode. Default=None (creates new one)
    Outputs:
        * [object] object containing data (can be read by gwcatpy.importGWTC)
    """
    import requests
    import json
    import h5py
    import ciecplib
    if url=='':
        if devMode:
            url='https://openscience-dev.ligo.caltech.edu/eventapi/json/GWTC/'
        elif useLocal:
            url='data/local-mirror/GWTC-dev.json'
        else:
            url='https://www.gw-openscience.org/eventapi/json/GWTC/'
    if verbose: print('Retrieving GWTC data from {}'.format(url))
    if devMode:
        if not sess:
            sess=ciecplib.Session("LIGO")
        gwtcresp = sess.get(url)
        gwtcread=json.loads(gwtcresp.content)
    elif url[0:4]=='http':
        gwtcresp = requests.get(url)
        gwtcread=json.loads(gwtcresp.text)
    else:
        gwtcread=json.load(open(url,'r'))
        
    gwtcin={'meta':{'retrieved':Time.now().isot,'src':url}}
    for s in gwtcread:
        gwtcin[s]=gwtcread[s]
    
    gwtcdata={'meta':{'retrieved':Time.now().isot,'src':url}}
    evvers={}
    for entry in gwtcin['events']:
        evname=gwtcin['events'][entry]['commonName']
        ver=gwtcin['events'][entry]['version']
        if not evname in evvers:
            evvers[evname]={'vers':{},'latest':0}
        evvers[evname]['vers'][ver]=entry
        if ver>evvers[evname]['latest']:
            evvers[evname]['latest']=ver
    # print(evvers)
    gwtcdata['data']={}
    for ev in evvers:
        print(ev,evvers[ev],evvers[ev]['vers'][evvers[ev]['latest']])
        gwtcdata['data'][ev]=gwtcin['events'][evvers[ev]['vers'][evvers[ev]['latest']]]
        # print(gwtcdata['data'][ev])
        if 'jsonurl' in gwtcdata['data'][ev]:
            jsonurl=gwtcdata['data'][ev]['jsonurl']
            # ***replace with local version***
            jsonurllocal='data/local-mirror/{}-v{}.json'.format(gwtcdata['data'][ev]['commonName'],gwtcdata['data'][ev]['version'])
            gwtcdata['data'][ev]['jsonurl_local']=jsonurllocal
            if useLocal: 
                if gwtcdata['data'][ev]['catalog.shortName']=='GWTC-1-confident':
                    jsonurl=jsonurl.replace('https://openscience-dev.ligo.caltech.edu/','https://www.gw-openscience.org/')
                if gwtcdata['data'][ev]['catalog.shortName']=='GWTC-2':
                    # replace GWTC-2 with local versions
                    jsonurl=jsonurllocal
            if verbose: print('Retrieving JSON {} data from {}'.format(ev,jsonurl))
            
            try:
                if devMode:
                    jsonresp = sess.get(jsonurl)
                    jsondata=json.loads(jsonresp.content)
                elif jsonurl[0:4]=='http':
                    if verbose:print('reading remote',jsonurl)
                    jsonresp = requests.get(jsonurl)
                    jsondata=json.loads(jsonresp.text)
                else:
                    if verbose:print('reading local',jsonurl)
                    jsondata=json.load(open(jsonurl,'r'))
            except:
                jsondata=None
                if verbose:print('unable to read data for {}'.format(ev))
                gwtcdata['data'][ev]['parameters']={}
                continue
            eventdata=jsondata['events']
            evdata=eventdata[next(iter(eventdata.keys()))]
            evname=evdata['commonName']
            if 'parameters' in evdata:
                gwtcdata['data'][ev]['parameters']=evdata['parameters']
                if evdata['catalog.shortName']=='GWTC-1-confident':
                    petag='gwtc1_pe_{}'.format(evdata['commonName'])
                elif evdata['catalog.shortName']=='GWTC-2':
                    petag='gwtc2_pe_{}'.format(evdata['commonName'])
                else:
                    petag='UNKNOWN'
                if petag in evdata['parameters']:
                    params=evdata['parameters'][petag]
                    if 'data_url' in params:
                        data_url=params['data_url']
                        gwtcdata['data'][ev]['data_link']=data_url
                        gwtcdata['data'][ev]['data_link_local']='data/local-mirror/{}-v{}.{}'.format(gwtcdata['data'][ev]['commonName'],gwtcdata['data'][ev]['version'],data_url.split('.')[-1])
            if 'strain' in evdata:
                gwtcdata['data'][ev]['strain']=evdata['strain']
    
    if export:
        if dirOut==None:
            dirOut='data/'
        if fileOut==None:
            fileOut='gwtc.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(gwtcdata,fOut,indent=indent)
        fOut.close()
    
    if verbose: print('Retrieved data for {} events'.format(len(gwtcdata['data'])))
    return gwtcdata

def gwtc_to_cat(gwtcdata,datadict,verbose=False,devMode=False):
    """Convert GWTC data to gwcat format.
    Inputs:
        * gwtcdata [object]: GWTC data (from getGWTC)
        * datadict [object]: Data dictionary (i.e. parameters)
        * verbose [boolean, optional]: set for verbose output. Default=False
        * devMode [boolean, optional]: set to use dev mode; replace "public/" in DCC links with "DocDB"
    Outputs:
        * [object] object containing 'data' and 'links'
    """
    import requests
    import json
    if 'data' in gwtcdata:
        gwtcin=gwtcdata['data']
    else:
        gwtcin=gwtcdata
    try:
        url=['meta']['url']
    except:
        url='https://openscience-dev.ligo.caltech.edu/eventapi/json/GWTC/'
    catOut={}
    linksOut={}
    conv={
        'mass_1_source':'M1',
        'mass_2_source':'M2',
        'luminosity_distance':'DL',
        'redshift':'z',
        'total_mass_source':'Mtotal',
        'chirp_mass_source':'Mchirp',
        'final_mass_source':'Mfinal',
        'chi_eff':'chi',
        'mass_ratio':'Mratio',
        'final_spin':'af',
        'a_final':'af',
        'L_peak':'lpeak',
        'E_rad':'Erad',
        'sky_size':'deltaOmega'
    }
    convsnr={
        'GPS':'GPS',
        'network_matched_filter_snr':'rho',
        'far':'FAR'
    }
    for e in gwtcin:
        catOut[e]={}
        linksOut[e]=[]
        if verbose: print('converting event',e)
        catOut[e]['name']=gwtcin[e]['commonName']
        catOut[e]['catalog']=gwtcin[e]['catalog.shortName']
        catOut[e]['version']=gwtcin[e]['version']
        if 'GPS' in gwtcin[e]:
            catOut[e]['GPS']={'best':gwtcin[e]['GPS']}
        # basic parameter conversions
        # conv = gwosc-name : gwcat-name
        print(catOut[e])
        catOut[e]['jsonURL']=gwtcin[e]['jsonurl']
        
        for c in convsnr:
            pdict=None
            if convsnr[c] in datadict:
                pdict=datadict[convsnr[c]]
            param=paramConv(gwtcin[e],c,pdict,verbose=verbose)
            if (param):
                catOut[e][convsnr[c]]=param
                # catOut[e][convsnr[c]]['src']=psnrname
        
        if gwtcin[e]['catalog.shortName']=='GWTC-1-confident':
            pname='gwtc1_pe_{}'.format(catOut[e]['name'])
            pycbcname='gwtc1_pycbc_{}'.format(catOut[e]['name'])
            gstlalname='gwtc1_gstlal_{}'.format(catOut[e]['name'])
            cwbname='gwtc1_cwb_{}'.format(catOut[e]['name'])
        elif gwtcin[e]['catalog.shortName']=='GWTC-2':
            pname='gwtc2_pe_{}'.format(catOut[e]['name'])
            pycbcname='gwtc2_pycbc_allsky_{}'.format(catOut[e]['name'])
            gstlalname='gwtc2_gstlal_allsky_{}'.format(catOut[e]['name'])
            cwbname='gwtc2_cwb_allsky_{}'.format(catOut[e]['name'])
        else:
            pname='UNKNOWN'
            pycbcname='UNKNOWN'
            gstlalname='UNKNOWN'
            cwbname='UNKNOWN'
        if pycbcname in gwtcin[e]['parameters']:
            psnrname=pycbcname
        elif gstlalname in gwtcin[e]['parameters']:
            psnrname=gstlalname
        elif cwbname in gwtcin[e]['parameters']:
            psnrname=cwbname
        else:
            psnrname='UNKNOWN'
        if pname in gwtcin[e]['parameters']:
            for c in conv:
                pdict=None
                if conv[c] in datadict:
                    pdict=datadict[conv[c]]
                param=paramConv(gwtcin[e]['parameters'][pname],c,pdict,verbose=verbose)
                if (param):
                    catOut[e][conv[c]]=param
        # if psnrname in gwtcin[e]['parameters']:
        #     for c in convsnr:
        #         pdict=None
        #         if conv[c] in datadict:
        #             pdict=datadict[conv[c]]
        #         param=paramConv(gwtcin[e]['parameters'][psnrname],c,pdict,verbose=verbose)
        #         if (param):
        #             catOut[e][conv[c]]=param
        #             catOut[e][conv[c]]['src']=psnrname
            # if 'far' in gwtcin[e]['parameters'][pycbcname]:
            #     catOut[e]['FAR']={'best':gwtcin[e]['parameters'][pycbcname]['far']}
            # if 'far' in gwtcin[e]['parameters'][pycbcname]:
            #     catOut[e]['FAR']={'best':gwtcin[e]['parameters'][pycbcname]['far']}
        catOut[e]['obsrun']={'best':gps2obsrun(gwtcin[e]['GPS'])}
        if verbose:print('GPS={} ; obsrun={}'.format(gwtcin[e]['GPS'],catOut[e]['obsrun']))
        catOut[e]['jsonURL']=gwtcin[e]['jsonurl']
        # get detectors:
        dets={'H1':False,'L1':False,'V1':False}
        if 'strain' in gwtcin[e]:
            if verbose:print('getting detector network')
            strains=gwtcin[e]['strain']
            for s in strains:
                if s['detector'] in dets:
                    dets[s['detector']]=True
            net=''
            for d in dets:
                if (dets[d]):
                    net+=d[0]
            if (net):
                catOut[e]['net']={'best':net}
            if verbose:print('{} Network:'.format(e),dets)
        else:
            if verbose:
                print('no "strain" found for {}:'.format(e),gwtcin[e])
        if 'jsonurl' in gwtcin[e]:
            jsonurllink={'url':gwtcin[e]['jsonurl'],
                'text':'JSON file',
                'type':'json-url'}
            linksOut[e].append(jsonurllink)
            htmlurllink={'url':gwtcin[e]['jsonurl'].replace('/json/','/html/'),
                'text':'GWTC link',
                'type':'open-data'}
            linksOut[e].append(htmlurllink)
        if 'jsonurl_local' in gwtcin[e]:
            jsonurllinklocal={'url':gwtcin[e]['jsonurl_local'],
                'text':'JSON file (local)',
                'type':'json-url-local'}
            linksOut[e].append(jsonurllinklocal)
        if 'data_link' in gwtcin[e]:
            datalink={'url':gwtcin[e]['data_link'],
                'text':'Data file',
                'type':'data-file'}
            if devMode:
                datalink['url']=datalink['url'].replace('public/','DocDB/')
            if datalink['url'].find('.tar')>0:
                datalink['filetype']='tar'
            elif datalink['url'].find('.h5')>0 or datalink['url'].find('.hdf5'):
                datalink['filetype']='h5'
            else:
                datalink['filetype']='unknown'
            linksOut[e].append(datalink)
        if 'data_link_local' in gwtcin[e]:
            datalinklocal={'url':gwtcin[e]['data_link_local'],
                'text':'Data file (local)',
                'type':'data-local'}
            if datalinklocal['url'].find('.tar')>0:
                datalinklocal['filetype']='tar'
            elif datalinklocal['url'].find('.h5')>0 or datalinklocal['url'].find('.hdf5'):
                datalinklocal['filetype']='h5'
            else:
                datalinklocal['filetype']='unknown'
            linksOut[e].append(datalinklocal)
            
        if catOut[e]['name'][0]=='G':
            catOut[e]['detType']={'best':'GW'}
            catOut[e]['conf']={'best':'GW'}
        dtIn=Time(gwtcin[e]['GPS'],format='gps')
        dtOut=Time(dtIn,format='iso').isot
        catOut[e]['UTC']={'best':dtOut}
        catOut[e]['meta']={'retrieved':Time.now().isot,'src':url}
        
    
    return ({'data':catOut,'links':linksOut})
    
def geth5paramsGWTC2(h5File,pcheck={},datadict={},verbose=False):
    """Extract parameters from GWTC2 HDF (.h5) files using pesummary
    Inputs:
        * hfFile [string]: filename to read
        * pcheck [object, optional]: Object containing parameter:value (used to identify closest approximant)
        * datadict [object]: Object containing data dictionary
        * verbose [boolean, optional]: set for verbose output. Default=False
    Outputs:
        * [object] object containing updated parameters
    """
    from pesummary.gw.file.read import read
    paramsOut={}
    conv={'mass_1_source':'M1',
        'mass_2_source':'M2',
        'luminosity_distance':'DL',
        'redshift':'z',
        'total_mass_source':'Mtotal',
        'chirp_mass_source':'Mchirp',
        'final_mass_source':'Mfinal',
        'geocent_time':'GPS',
        'chi_eff':'chi',
        'mass_ratio':'Mratio',
        'final_spin':'af',
        'peak_luminosity':'lpeak',
        'radiated_energy':'Erad'
    }
    def getParam(samps,param):
        if param in samps:
            best=samps.median[param][0]
            lower=np.percentile(samps[param],5)
            upper=np.percentile(samps[param],95)
            return({'best':best,'err':[lower-best,upper-best]})
        else:
            return(None)
    
    try:
        h5dat=read(h5File)
    except:
        if verbose:print('error reading data file: {}'.format(h5File))
        return({})
    if verbose:print('getting parameters from data file: {}'.format(h5File))
    approximants=h5dat.labels
    thisapprox='PublicationSamples'
    if not thisapprox in approximants:
        if len(pcheck)>0:
            pcheckl=next(iter(pcheck.keys()))
            pcheckv=pcheck[pcheckl]
            if verbose:print('finding best match for {}={}'.format(pcheckl,pcheckv))
            for c in conv:
                if conv[c]==pcheckl:
                    pcheckl=c
            diff=[]
            for a in approximants:
                checka=h5dat.samples_dict[a].median[pcheckl][0]
                diffa=np.abs(checka-pcheckv)
                diff.append(diffa)
                if verbose:print(a,diffa)
            thisapprox=approximants[np.argmin(diff)]
        else:
            # just use first
            thisapprox=approximants[0]
    if verbose:print('using approximant {}'.format(thisapprox))
    
    # get more params
    h5samp=h5dat.samples_dict[thisapprox]
    params={}
    for c in conv:
        if c in h5samp.parameters:
            params[c]=h5samp.median[c][0]
            params[c+'_lower']=np.percentile(h5samp[c],5)-params[c]
            params[c+'_upper']=np.percentile(h5samp[c],95)-params[c]
    
    for c in conv:
        pdict=None
        if conv[c] in datadict:
            pdict=datadict[conv[c]]
        param=paramConv(params,c,pdict,verbose=verbose)
        if (param):
            paramsOut[conv[c]]=param
    
    if (h5dat.detectors):
        dets=''
        for d in h5dat.detectors[0]:
            dets=dets+d[0]
        paramsOut['net']={'best':dets}
    return(paramsOut)

def paramConv(evdat,param,paramdict,verbose=False):
    """Convert parameters from GWTC to gwcat formats, and truncate to sigfigs
    Inputs:
        * evdat [object]: event data in GWTC format
        * param [object, optional]: GWTC parameter name
        * paramdict [object]: data dictionary object for parameter (used to check significant figures)
        * verbose [boolean, optional]: set for verbose output. Default=False
    Outputs:
        * [object] object containing updated parameters
    """
    # convert from gwosc O3 parameter format to gwcat parameter format
    if not param in evdat:
        if verbose: print('no parameter {} found'.format(param))
        return
    if evdat[param]==None:
        if verbose: print('null parameter {} found'.format(param))
        return
    sigfig=None
    if (paramdict):
        if 'sigfig' in paramdict:
            sigfig=paramdict['sigfig']
    
    plo=param+'_lower'
    phi=param+'_upper'
    if plo in evdat and phi in evdat:
        haserr=True
    else:
        haserr=False
    
    pOut={'best':evdat[param]}
    if haserr:
        if (evdat[plo]) and (evdat[phi]):
            pOut['err']=[evdat[plo],evdat[phi]]

    return(pOut)
