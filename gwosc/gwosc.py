from astropy.time import Time
import os
import numpy as np
def paramConv(param):
    # convert from gwosc parameter format to gwcat parameter format
    if not 'best' in param:
        pOut=param
    elif not 'err' in param:
        pOut={'best':param['best']}
    elif type(param['err'])==str:
        if param['err']=='lowerbound':
            pOut={'lower':param['best']}
        elif param['err']=='upperbound':
            pOut={'upper':param['best']}
    else:
        pOut={'best':param['best'],'err':param['err']}
    return(pOut)


def gwtc1_to_cat(gwtc1,verbose=False):
    # convert dataset from gwosc to gwcat format
    if 'data' in gwtc1:
        gwtc1In=gwtc1['data']
    else:
        gwtc1In=gwtc1
    try:
        url=gwtc1['meta']['url']
    except:
        url='https://www.gw-openscience.org/catalog/GWTC-1-confident/filelist/'
    catOut={}
    linksOut={}
    for e in gwtc1In:
        if verbose: print('converting event',e)
        catOut[e]={}
        # basic parameter conversions
        # conv = gwosc-name : gwcat-name
        conv={
            'mass1':'M1',
            'mass2':'M2',
            'E_rad':'Erad',
            'L_peak':'lpeak',
            'distance':'DL',
            'redshift':'z',
            'sky_size':'deltaOmega',
            'mfinal':'Mfinal',
            'mchirp':'Mchirp',
            'chi_eff':'chi',
            'a_final':'af',
            'tc':'GPS'
        }
        for c in conv:
            if c in gwtc1In[e]:
                param=gwtc1In[e][c]
                catOut[e][conv[c]]=paramConv(param)

        # convert parameters in 'files'
        if 'files' in gwtc1In[e]:
            if 'ObsRun' in gwtc1In[e]['files']:
                catOut[e]['obsrun']={'best':gwtc1In[e]['files']['ObsRun']}
            if 'eventName' in gwtc1In[e]['files']:
                catOut[e]['name']=paramConv(gwtc1In[e]['files']['eventName'])
                if catOut[e]['name'][0]=='G':
                    catOut[e]['detType']={'best':'GW'}
                    catOut[e]['conf']={'best':'GW'}
                else:
                    catOut[e]['detType']={'best':'LVT'}
                    catOut[e]['conf']={'best':'LVT'}
            # convert UTC
            if 'PeakAmpGPS' in gwtc1In[e]['files']:
                dtIn=Time(gwtc1In[e]['files']['PeakAmpGPS'],format='gps')
                dtOut=Time(dtIn,format='iso').isot
                catOut[e]['UTC']={'best':dtOut}
            if 'OperatingIFOs' in gwtc1In[e]['files']:
                netOut=''
                ifos=gwtc1In[e]['files']['OperatingIFOs']
                if ifos.find('H1')>=0: netOut=netOut+'H'
                if ifos.find('L1')>=0: netOut=netOut+'L'
                if ifos.find('V1')>=0: netOut=netOut+'V'
                catOut[e]['net']={'best':netOut}
        # convert object type based on masses
        if 'M1' in catOut[e] and 'M2' in catOut[e]:
            if catOut[e]['M1']['best']<3:
                catOut[e]['objType']={'best':'BNS'}
            elif catOut[e]['M1']['best']<5:
                catOut[e]['objType']={'best':'MassGap'}
            else:
                if catOut[e]['M2']['best']<3:
                    catOut[e]['objType']={'best':'NSBH'}
                elif catOut[e]['M2']['best']<5:
                    catOut[e]['objType']={'best':'MassGap'}
                else: catOut[e]['objType']={'best':'BBH'}

        # convert FAR and SNR (priority pycbc - gstlal - cwb)
        if 'far_pycbc' in gwtc1In[e]:
            if gwtc1In[e]['far_pycbc']['best']!='NA':
                catOut[e]['FAR']={'best':gwtc1In[e]['far_pycbc']['best'],'fartype':'pycbc'}
            elif gwtc1In[e]['far_gstlal']['best']!='NA':
                catOut[e]['FAR']={'best':gwtc1In[e]['far_gstlal']['best'],'fartype':'gstlal'}
            elif gwtc1In[e]['far_cwb']['best']!='NA':
                catOut[e]['FAR']={'best':gwtc1In[e]['far_cwb']['best'],'fartype':'cwb'}

        if 'snr_pycbc' in gwtc1In[e]:
            if gwtc1In[e]['snr_pycbc']['best']!='NA':
                catOut[e]['rho']={'best':gwtc1In[e]['snr_pycbc']['best'],'snrtype':'pycbc'}
            elif gwtc1In[e]['snr_gstlal']['best']!='NA':
                catOut[e]['rho']={'best':gwtc1In[e]['snr_gstlal']['best'],'snrtype':'gstlal'}
            elif gwtc1In[e]['snr_cwb']['best']!='NA':
                catOut[e]['rho']={'best':gwtc1In[e]['snr_cwb']['best'],'snrtype':'cwb'}
        elif 'snr' in gwtc1In[e]:
            catOut[e]['rho']={'best':gwtc1In[e]['snr']['best']}
            if 'pipeline' in gwtc1In[e]:
                catOut[e]['rho']['snrtype']=gwtc1In[e]['pipeline']['best']

        catOut[e]['meta']={'retrieved':Time.now().isot,'src':url}
    return({'data':catOut,'links':linksOut})

def gwosc2cat(gwosc,verbose=False):
    print('***WARNING gwosc2cat replaced by gwtc1_to_cat***')
    out=gwct1_to_cat(gwosc,verbose=False)
    return(out)
    
def getGWTC1(url='',verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False):
    import requests
    import json
    if url=='':
        if triggers:
            url='https://www.gw-openscience.org/catalog/GWTC-1-marginal/filelist/'
        else:
            url='https://www.gw-openscience.org/catalog/GWTC-1-confident/filelist/'
    if verbose: print('Retrieving GWTC-1 data from {}'.format(url))
    gwtc1resp = requests.get(url)
    gwtc1read=json.loads(gwtc1resp.text)
    gwtc1data={'meta':{'retrieved':Time.now().isot,'src':url}}
    for s in gwtc1read:
        gwtc1data[s]=gwtc1read[s]

    if verbose: print('Retrieved data for {} events'.format(len(gwtc1data['data'])))

    if export:
        if dirOut==None:
            dirOut='../../data/'
        if fileOut==None:
            if triggers:
                fileOut='gwtc1-marginal.json'
            else:
                fileOut='gwtc1.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(gwtc1data,fOut,indent=indent)
        fOut.close()
    return gwtc1data
    
def getGwosc(url='',verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False):
    print('***WARNING getGwosc replaced by getGWTC1***')
    out=getGWTC1(url=url,verbose=verbose,export=export,dirOut=dirOut,fileOut=fileOut,indent=indent,triggers=triggers)
    return(out)
    
        
def getGwoscO3(url='',verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False):
    import requests
    import json
    import h5py
    if url=='':
        url='../../data/O3cat.json'
    if verbose: print('Retrieving O3 data from {}'.format(url))
    if url[0:4]=='http':
        o3resp = requests.get(url)
        o3read=json.loads(o3resp.text)
    else:
        o3read=json.load(open(url,'r'))

    o3in={'meta':{'retrieved':Time.now().isot,'src':url}}
    for s in o3read:
        o3in[s]=o3read[s]

    o3data={'meta':{'retrieved':Time.now().isot,'src':url}}
    evvers={}
    for entry in o3in['events']:
        evname=o3in['events'][entry]['commonName']
        ver=o3in['events'][entry]['version']
        if not evname in evvers:
            evvers[evname]={'vers':{},'latest':0}
        evvers[evname]['vers'][ver]=entry
        if ver>evvers[evname]['latest']:
            evvers[evname]['latest']=ver
    o3data['data']={}
    for ev in evvers:
        o3data['data'][ev]=o3in['events'][evvers[ev]['vers'][evvers[ev]['latest']]]
        if 'jsonurl' in o3data['data'][ev]:
            jsonurl=o3data['data'][ev]['jsonurl']
            # ***replace with local version***
            jsonurllocal='data/json/{}_v{}.json'.format(o3data['data'][ev]['commonName'],o3data['data'][ev]['version'])
            o3data['data'][ev]['jsonurl_local']=jsonurllocal
            jsonurl=jsonurllocal
            if verbose: print('Retrieving O3 {} data from {}'.format(ev,jsonurl))
            if jsonurl[0:4]=='http':
                if verbose:print('reading remote',jsonurl)
                jsonresp = requests.get(jsonurl)
                jsondata=json.loads(jsonresp.text)
            else:
                if verbose:print('reading local',jsonurl)
                jsondata=json.load(open(jsonurl,'r'))
            eventsdata=jsondata[next(iter(jsondata.keys()))]
            evdata=eventsdata[next(iter(eventsdata.keys()))]
            o3data['data'][ev]=evdata
            if 'parameters' in evdata:
                params=evdata['parameters']
                if len(params)>0:
                    param0=params[next(iter(params.keys()))]
                    if 'data_url' in param0:
                        data_url=param0['data_url']
                        o3data['data'][ev]['data_link']=data_url
                        o3data['data'][ev]['data_link_local']='data/h5/{}_v{}.h5'.format(o3data['data'][ev]['commonName'],o3data['data'][ev]['version'])
        # if 'data_link_local' in o3data['data'][ev]:
        #     h5=h5py.File(o3data['data'][ev]['data_link_local'],'r')
        # print(ev,evvers[ev],evvers[ev]['vers'][evvers[ev]['latest']])

    if export:
        if dirOut==None:
            dirOut='../../data/'
        if fileOut==None:
            fileOut='o3data.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(o3data,fOut,indent=indent)
        fOut.close()

    if verbose: print('Retrieved data for {} events'.format(len(o3data['data'])))
    return o3data

def geth5params(h5File,pcheck={},datadict={},verbose=False):
    from pesummary.gw.file.read import read
    paramsOut={}
    conv={'mass_1_source':'M1',
        'mass_2_source':'M2',
        'luminosity_distance':'DL',
        'redshift':'z',
        'total_mass_source':'Mtotal',
        'chirp_mass_source':'Mchirp',
        'final_mass_source':'Mfinal',
        't0':'GPS',
        'network_matched_filter_snr':'rho',
        'chi_eff':'chi',
        'mass_ratio':'Mratio',
        'final_spin':'af',
        'peak_luminosity':'lpeak',
        'radiated energy':'Erad'
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
        param=paramConvO3(params,c,pdict,verbose=verbose)
        if (param):
            paramsOut[conv[c]]=param
    
    if (h5dat.detectors):
        dets=''
        for d in h5dat.detectors[0]:
            dets=dets+d[0]
        paramsOut['net']={'best':dets}
    return(paramsOut)

def paramConvO3(evdat,param,paramdict,verbose=False):
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
        pOut['err']=[evdat[plo],evdat[phi]]

    return(pOut)
