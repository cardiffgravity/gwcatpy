from astropy.time import Time
import os
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


def gwosc2cat(gwosc,verbose=False):
    # convert dataset from gwosc to gwcat format
    if 'data' in gwosc:
        gwoscIn=gwosc['data']
    else:
        gwoscIn=gwosc
    try:
        url=gwosc['meta']['url']
    except:
        url='https://www.gw-openscience.org/catalog/GWTC-1-confident/filelist/'
    catOut={}
    linksOut={}
    for e in gwoscIn:
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
            if c in gwoscIn[e]:
                param=gwoscIn[e][c]
                catOut[e][conv[c]]=paramConv(param)

        # convert parameters in 'files'
        if 'files' in gwoscIn[e]:
            if 'ObsRun' in gwoscIn[e]['files']:
                catOut[e]['obsrun']={'best':gwoscIn[e]['files']['ObsRun']}
            if 'eventName' in gwoscIn[e]['files']:
                catOut[e]['name']=paramConv(gwoscIn[e]['files']['eventName'])
                if catOut[e]['name'][0]=='G':
                    catOut[e]['detType']={'best':'GW'}
                    catOut[e]['conf']={'best':'GW'}
                else:
                    catOut[e]['detType']={'best':'LVT'}
                    catOut[e]['conf']={'best':'LVT'}
            # convert UTC
            if 'PeakAmpGPS' in gwoscIn[e]['files']:
                dtIn=Time(gwoscIn[e]['files']['PeakAmpGPS'],format='gps')
                dtOut=Time(dtIn,format='iso').isot
                catOut[e]['UTC']={'best':dtOut}
            if 'OperatingIFOs' in gwoscIn[e]['files']:
                netOut=''
                ifos=gwoscIn[e]['files']['OperatingIFOs']
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
        if 'far_pycbc' in gwoscIn[e]:
            if gwoscIn[e]['far_pycbc']['best']!='NA':
                catOut[e]['FAR']={'best':gwoscIn[e]['far_pycbc']['best'],'fartype':'pycbc'}
            elif gwoscIn[e]['far_gstlal']['best']!='NA':
                catOut[e]['FAR']={'best':gwoscIn[e]['far_gstlal']['best'],'fartype':'gstlal'}
            elif gwoscIn[e]['far_cwb']['best']!='NA':
                catOut[e]['FAR']={'best':gwoscIn[e]['far_cwb']['best'],'fartype':'cwb'}

        if 'snr_pycbc' in gwoscIn[e]:
            if gwoscIn[e]['snr_pycbc']['best']!='NA':
                catOut[e]['rho']={'best':gwoscIn[e]['snr_pycbc']['best'],'snrtype':'pycbc'}
            elif gwoscIn[e]['snr_gstlal']['best']!='NA':
                catOut[e]['rho']={'best':gwoscIn[e]['snr_gstlal']['best'],'snrtype':'gstlal'}
            elif gwoscIn[e]['snr_cwb']['best']!='NA':
                catOut[e]['rho']={'best':gwoscIn[e]['snr_cwb']['best'],'snrtype':'cwb'}
        elif 'snr' in gwoscIn[e]:
            catOut[e]['rho']={'best':gwoscIn[e]['snr']['best']}
            if 'pipeline' in gwoscIn[e]:
                catOut[e]['rho']['snrtype']=gwoscIn[e]['pipeline']['best']

        catOut[e]['meta']={'retrieved':Time.now().isot,'src':url}
    return({'data':catOut,'links':linksOut})

def getGwosc(url='',verbose=True,export=False,dirOut=None,fileOut=None,indent=2,triggers=False):
    import requests
    import json
    if url=='':
        if triggers:
            url='https://www.gw-openscience.org/catalog/GWTC-1-marginal/filelist/'
        else:
            url='https://www.gw-openscience.org/catalog/GWTC-1-confident/filelist/'
    if verbose: print('Retrieving GWOSC data from {}'.format(url))
    gwoscresp = requests.get(url)
    gwoscread=json.loads(gwoscresp.text)
    gwoscdata={'meta':{'retrieved':Time.now().isot,'src':url}}
    for s in gwoscread:
        gwoscdata[s]=gwoscread[s]

    if verbose: print('Retrieved data for {} events'.format(len(gwoscdata['data'])))

    if export:
        if dirOut==None:
            dirOut='../../data/'
        if fileOut==None:
            if triggers:
                fileOut='gwosc-marginal.json'
            else:
                fileOut='gwosc.json'
        if verbose: print('Exporting to {}'.format(os.path.join(dirOut,fileOut)))
        fOut=open(os.path.join(dirOut,fileOut),'w')
        json.dump(gwoscdata,fOut,indent=indent)
        fOut.close()
    return gwoscdata