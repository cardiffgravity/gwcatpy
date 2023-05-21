from astropy.time import Time
import os
import numpy as np
import zenodo_get as zen
import re

catlist={'GWTC':{'type':'confident'},
    'GWTC-1-confident':{'type':'confident'},
    'GWTC-2':{'type':'confident'},
    'GWTC-2.1-confident':{'type':'confident','zenodo':'5117703'},
    # 'GWTC-2.1-confident':{'type':'confident'},
    'GWTC-3-confident':{'type':'confident','zenodo':'5546663'},
    'O3_Discovery_Papers':{'type':'confident'},
    'GWTC-1-marginal':{'type':'marginal'},
    'GWTC-2-marginal':{'type':'marginal'},
    'GWTC-2.1-marginal':{'type':'marginal'},
    'GWTC-2.1-auxiliary':{'type':'marginal'},
    'GWTC-3-marginal':{'type':'marginal'},
}

def cat2url(cat,devMode=False):
    roots={'dev':'https://openscience-dev.ligo.caltech.edu/eventapi/json/',
        'public':'https://www.gw-openscience.org/eventapi/json/'
    }
    if not cat in catlist:
        print('WARNING: catalog {} not valid.'.format(cat))
        url=''
    else:
        if devMode:
            url=roots['dev']+cat
        else:
            url=roots['public']+cat
    return url

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
        'O4':[1366556418,1419724817]}
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

def getGWTC(url='',useLocal=False,verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False,devMode=False,catalog='GWTC',sess=None):
    """Get GWTC from GWOSC json files and add parameters.
    Inputs:
        * url [string, optional]: URL to download from. Default=https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/
        * useLocal [boolean, optional]: Set to use local files (default URL=data/local-mirror/GWTC-dev.json)
        * verbose [boolean, optional]: set for verbose output. Default=False
        * export [boolean, optional]: set for export to JSON file output. Default=False
        * dirOut [string, optional]: directory to export to (if export=True). Default='data/'
        * fileOut [string, optional]: file to export to (if export=True). Default='GWTC.json'
        * indent [integer, optional]: json indent for exported file. Default=2
        * triggers [boolean, optional]: OBSOLETE
        * devMode [boolean, optional]: set to ciecplib to access data, and use dev link as default (https://openscience-dev.ligo.caltech.edu/eventapi/json/GWTC/)
        * catalog [string, optional]: catalog to read from. Default=GWTC.
        * sess [ciecplib.Session, optional]: ciecplib session to use in dev mode. Default=None (creates new one)
    Outputs:
        * [object] object containing data (can be read by gwcatpy.importGWTC)
    """
    import requests
    import json
    import h5py
    import ciecplib
    if url=='':
        url=cat2url(catalog,devMode=devMode)
        # if devMode:
        #     url='https://openscience-dev.ligo.caltech.edu/eventapi/json/GWTC/'
        # elif useLocal:
        #     url='data/local-mirror/GWTC-dev.json'
        # else:
        #     url='https://www.gw-openscience.org/eventapi/json/GWTC/'
    if verbose: print('Retrieving {} data from {}'.format(catalog,url))
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

        # get zenodo
        evcat=gwtcdata['data'][ev]['catalog.shortName']
        zenFiles=[]
        if 'zenodo' in catlist[evcat]:
            if not 'zenFiles' in catlist[evcat]:
                # (re)download zenodo file list
                zenFileList=os.path.join(dirOut,'{}_zenodo-filelist.txt'.format(evcat))
                try:
                    # download new file
                    zen.zenodo_get(['--wget={}'.format(zenFileList),catlist[evcat]['zenodo']])
                    gwtcdata['meta']['zenodoFileList']=zenFileList
                    fzen=open(zenFileList,'r')
                    catlist[evcat]['zenFiles']=fzen.readlines()
                    fzen.close()
                    catlist[evcat]['zenodoLoaded']=True
                    if verbose:print('Downloaded {} Zenodo file list from {} to {}'.format(evcat,catlist[evcat]['zenodo'],zenFileList))
                    zenFiles=catlist[evcat]['zenFiles']
                except:
                    print('WARNING: unable to save {} Zenodo file list from {} to {}'.format(evcat,catlist[evcat]['zenodo'],zenFileList))
                    try:
                        fzen=open(zenFileList,'r')
                        catlist[evcat]['zenFiles']=fzen.readlines()
                        fzen.close()
                        catlist[evcat]['zenodoLoaded']=True
                        if verbose:print('Using existing Zenodo file for {}: {}'.format(evcat,zenFileList))
                        zenFiles=catlist[evcat]['zenFiles']
                    except:
                        print('WARNING: no existing Zenodo file to use')
            else:
                # zenFiles list already exists
                zenFiles=catlist[evcat]['zenFiles']

        if 'jsonurl' in gwtcdata['data'][ev]:
            jsonurl=gwtcdata['data'][ev]['jsonurl']
            if devMode:
                jsonurl=jsonurl.replace('https://www.gw-openscience.org','https://openscience-dev.ligo.caltech.edu')
            # ***replace with local version***
            jsonurllocal='data/local-mirror/{}-v{}.json'.format(gwtcdata['data'][ev]['commonName'],gwtcdata['data'][ev]['version'])
            if useLocal:
                # fudge for local mirrors
                gwtcdata['data'][ev]['jsonurl_local']=jsonurllocal
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
                petag='UNKNOWN'
                srchtag='UNKNOWN'
                #find default PE and search tags
                for par in evdata['parameters']:
                    if evdata['parameters'][par]['is_preferred']:
                        if evdata['parameters'][par]['pipeline_type']=='pe':
                            petag=par
                        elif evdata['parameters'][par]['pipeline_type']=='search':
                            srchtag=par
                # find PE tag that matches if not found
                if not petag in evdata['parameters']:
                    petagbest=getBestParam(evdata,{'M1':None},{'mass_1_source':'M1'},verbose=verbose)
                    if petagbest:
                        petag=petagbest
                if petag in evdata['parameters']:
                    params=evdata['parameters'][petag]
                    gwtcdata['data'][ev]['parameter_tag']=petag
                    if 'data_url' in params:
                        data_url=params['data_url']
                        gwtcdata['data'][ev]['data_link']=data_url
                        if useLocal:
                            gwtcdata['data'][ev]['data_link_local']='data/local-mirror/{}-v{}.{}'.format(gwtcdata['data'][ev]['commonName'],gwtcdata['data'][ev]['version'],data_url.split('.')[-1])
                if not srchtag in evdata['parameters']:
                    srchtagbest=getBestParam(evdata,{'rho':None},{'network_matched_filter_snr':'rho'},verbose=verbose,search=True)
                    if srchtagbest:
                        srchtag=srchtagbest
                if srchtag in evdata['parameters']:
                    gwtcdata['data'][ev]['search_parameter_tag']=srchtag
            if 'zenodo' in catlist[evcat]:
                # find h5 file
                if ev=='GW200210_092254':
                    # to catch zenodo filename typo!
                    zenev=ev.replace('GW200210_092254','GW200210_092255')
                else:
                    zenev=ev
                for zenf in zenFiles:
                    if zenf.find('nocosmo.h5')<0 and zenf.find(zenev)>=0:
                        # make sure not to include nocosmo files from GWTC-3-confident
                        if verbose:print('zenodo data link for {}:{}'.format(ev,zenf))
                        if not 'data_link' in gwtcdata['data'][ev]:
                            gwtcdata['data'][ev]['data_link']=zenf.replace('\n','')
                            if verbose:print('using zenodo data_link')
                        else:
                            if verbose:print('using old data_link: {}'.format(gwtcdata['data'][ev]['data_link']))
                        gwtcdata['data'][ev]['zenodo_version']=int(re.match('.*\/record\/([0-9]*)\/',zenf).groups()[0])
                        if verbose:print('zenodo data link for {}:{}'.format(ev,zenf))
                    if zenf.find('skymaps.tar.gz')>=0 or zenf.find('SkyMaps.tar.gz')>=0:
                        gwtcdata['data'][ev]['map_link']=zenf.replace('\n','')
                        if verbose:print('map link for {}:{}'.format(ev,zenf))
            if 'strain' in evdata:
                gwtcdata['data'][ev]['strain']=evdata['strain']

    if export:
        if dirOut==None:
            dirOut='data/'
        if fileOut==None:
            fileOut=catalog+'.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(gwtcdata,fOut,indent=indent)
        fOut.close()

    if verbose: print('Retrieved data for {} events'.format(len(gwtcdata['data'])))
    return gwtcdata

def gwtc_to_cat(gwtcdata,datadict,verbose=False,devMode=False,catalog='GWTC'):
    """Convert GWTC data to gwcat format.
    Inputs:
        * gwtcdata [object]: GWTC data (from getGWTC)
        * datadict [object]: Data dictionary (i.e. parameters)
        * verbose [boolean, optional]: set for verbose output. Default=False
        * devMode [boolean, optional]: set to use dev mode; replace "public/" in DCC links with "DocDB"
        * catalog [string, optional]: catalog to read from. Default=GWTC.
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
        url=gwtcdata['meta']['src']
    except:
        url=cat2url(catalog,devMode=devMode)
        # if devMode:
        #     url='https://openscience-dev.ligo.caltech.edu/eventapi/json/GWTC/'
        # else:
        #     url='https://gw-openscience.org/eventapi/json/GWTC/'
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
        'far':'FAR',
        'p_astro':'p_astro'
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

        # srchname=gwtcin[e].get('search_parameter_tag','UNKNOWN')
        # if srchname in gwtcin[e]['parameters']:
        #     if verbose:
        #         print('reading search parameters from {}'.format(srchname))
        #         print(gwtcin[e]['parameters'][srchname])
        #     srchparamIn=gwtcin[e]['parameters'][srchname]
        # else:
        #     if verbose:
        #         print('cannot find {} parameters. Using root parameters'.format(srchname))
        srchparamIn=gwtcin[e]
        for c in convsnr:
            srchpdict=None
            if convsnr[c] in datadict:
                srchpdict=datadict[convsnr[c]]
            else:
                if verbose:print('{} not in datadict'.format(convsnr[c]))
            srchparam=paramConv(srchparamIn,c,srchpdict,verbose=verbose)
            if (srchparam):
                catOut[e][convsnr[c]]=srchparam
                # catOut[e][convsnr[c]]['src']=psnrname

        pname=gwtcin[e].get('parameter_tag','UNKNOWN')
        if pname in gwtcin[e]['parameters']:
            if verbose:
                print('reading parameters from {}'.format(pname))
            paramIn=gwtcin[e]['parameters'][pname]
        else:
            if verbose:
                print('cannot find {} parameters. Using root parameters'.format(pname))
            paramIn=gwtcin[e]
        for c in conv:
            pdict=None
            if conv[c] in datadict:
                pdict=datadict[conv[c]]
            param=paramConv(paramIn,c,pdict,verbose=verbose)
            if (param):
                catOut[e][conv[c]]=param
            if 'waveform_family' in paramIn:
                catOut[e]['approximant']=paramIn['waveform_family']

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
            if gwtcin[e]['data_link']!='':
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
                if datalink['url'].find('gw-openscience')>0:
                    datalink['src']='gwosc'
                elif datalink['url'].find('zenodo')>0:
                    datalink['src']='zenodo'
                linksOut[e].append(datalink)
        if 'map_link' in gwtcin[e]:
            if gwtcin[e]['map_link']!='':
                maplink={'url':gwtcin[e]['map_link'],'text':'Sky Map',
                        'type':'skymap-fits'}
                if maplink['url'].find('tar')>0:
                    maplink['filetype']='tar'
                linksOut[e].append(maplink)
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

        if catlist[catOut[e]['catalog']]['type']=='marginal':
            if catOut[e]['name'][0]=='G':
                catOut[e]['detType']={'best':'Marginal'}
                catOut[e]['conf']={'best':'Marginal'}
            else:
                catOut[e]['detType']={'best':'Marginal'}
                catOut[e]['conf']={'best':'Candidate'}
        elif catlist[catOut[e]['catalog']]['type']=='confident':
            catOut[e]['detType']={'best':'GW'}
            catOut[e]['conf']={'best':'GW'}
        elif catOut[e]['name'][0]=='G':
            catOut[e]['detType']={'best':'GW'}
            catOut[e]['conf']={'best':'GW'}
        else:
            catOut[e]['detType']={'best':'Candidate'}
            catOut[e]['conf']={'best':'Candidate'}
        dtIn=Time(gwtcin[e]['GPS'],format='gps')
        dtOut=Time(dtIn,format='iso').isot
        catOut[e]['UTC']={'best':dtOut}
        catOut[e]['meta']={'retrieved':Time.now().isot,'src':url}
        if 'date_added' in paramIn:
            catOut[e]['meta']['data_date_added']=paramIn['date_added']
        if 'zenodo_version' in gwtcin[e]:
            catOut[e]['meta']['zenodo_version']=gwtcin[e]['zenodo_version']
        if 'parameter_tag' in gwtcin[e]:
            catOut[e]['meta']['parameter_tag']=gwtcin[e]['parameter_tag']


    return ({'data':catOut,'links':linksOut})

def geth5paramsGWTC2(h5File,pcheck={},datadict={},approx=None,verbose=False):
    """Extract parameters from GWTC2 HDF (.h5) files using pesummary
    Inputs:
        * hfFile [string]: filename to read
        * pcheck [object, optional]: Object containing parameter:value (used to identify closest approximant)
        * datadict [object]: Object containing data dictionary
        * approx [string, optional]: Approximant to try first (Detault="PublicationSamples")
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
    if approx:
        thisapprox=approx
    else:
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
    paramsOut['approximant']=thisapprox
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

def getBestParam(gwtcin_e,datadict={'M1':None},conv=None,verbose=False,search=False):
    # find best-fit parameters
    # if verbose:
    #     if search:
    #         print('searching for best search parameter')
    #     else:
    #         print('searching for best PE parameter')
    if not conv:
        if search:
            conv={'network_matched_filter_snr':'rho'}
        else:
            conv={'mass_1_source':'M1'}
    pnbest=''
    for col in conv:
        c_chk=col
    # if verbose:
        # print('checking {} against {}'.format(c_chk,conv[c_chk]))
    # check whether dataset has is_preferred set and mass_1_source
    for pn in gwtcin_e['parameters']:
        # if verbose:print('checking {}'.format(pn))
        if c_chk in gwtcin_e['parameters'][pn]:
            # if verbose:print('{} in {}'.format(c_chk,pn),gwtcin_e['parameters'][pn])
            if gwtcin_e['parameters'][pn].get('is_preferred','')==True:
                # if verbose:print('{} preferred'.format(pn))
                pnbest=pn
                if verbose:
                    if search:
                        print('found preferred search parameters in {}'.format(pn))
                    else:
                        print('found preferred PE parameters in {}'.format(pn))
    if (pnbest):
        # preferred parameter set is found
        return pnbest

    # compare parameter values and compare with "root" values
    chkroot=paramConv(gwtcin_e,c_chk,datadict[conv[c_chk]],verbose=False)
    if chkroot:
        chkroot=chkroot['best']
    else:
        chkroot=np.nan
    pnames=[]
    pchks=[]
    peparam=[]
    srchparam=[]
    pdates=[]

    for pn in gwtcin_e['parameters']:
        pnames.append(pn)
        # figure out if parameters are srch/PE parameters
        if search:
            srchparam.append(c_chk in gwtcin_e['parameters'][pn])
            peparam.append(np.nan)
        else:
            ispe=pn.find('pe')>=0
            if not ispe:
                ispe=(c_chk in gwtcin_e['parameters'][pn])
            peparam.append(ispe)
            srchparam.append(np.nan)
        # get date_added from parameters
        if 'date_added' in gwtcin_e['parameters'][pn]:
            pdates.append(Time(gwtcin_e['parameters'][pn]['date_added']).gps)
        else:
            pdates.append(np.nan)
        # get check parameter
        paramch=paramConv(gwtcin_e['parameters'][pn],c_chk,datadict[conv[c_chk]],verbose=False)
        if paramch:
            paramch=paramch['best']
        else:
            paramch=np.nan
        pchks.append(paramch)
    pnames=np.array(pnames)
    pdates=np.array(pdates)
    peparam=np.array(peparam)
    srchparam=np.array(srchparam)
    pchks=np.array(pchks)
    finchk=np.isfinite(pchks)
    findates=np.isfinite(pdates)

    if search:
        if len(np.argwhere(srchparam))==1:
            pnbest=pnames[np.argwhere(srchparam)][0][0]
            if verbose:
                print('found only search parameters in {}'.format(pnbest))
        elif len(pdates[findates])>0 and len(np.argwhere(srchparam))>0:
            pnbest=pnames[srchparam][np.argmax(pdates[srchparam])]
            if verbose:
                print('found latest search parameters in {}'.format(pnbest))
    else:
        if len(np.argwhere(peparam))==1:
            pnbest=pnames[np.argwhere(peparam)][0][0]
            if verbose:
                print('found only PE parameters in {}'.format(pnbest))
        elif len(pdates[findates])>0 and len(np.argwhere(peparam))>0:
            pnbest=pnames[peparam][np.argmin(pdates[peparam])]
            if verbose:
                print('found latest PE parameters in {}'.format(pnbest))
    if not pnbest:
        if len(pchks[finchk])>0:
            bestidx=np.argmin(np.abs(pchks[finchk]-chkroot))
            pnbest=pnames[finchk][bestidx]
            if verbose:
                print('found closest-matching parameters in {} [{} vs {}]'.format(pnbest,pchks[finchk],chkroot))
    return pnbest
