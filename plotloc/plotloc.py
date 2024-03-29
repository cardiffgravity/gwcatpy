from ligo.gracedb.rest import GraceDb
import json
import healpy as hp
import numpy as np
import requests
import os
from matplotlib import pyplot as plot
from matplotlib import cm
from astropy.io import fits
from astropy.table import Table, Column
import astropy_healpix as ah
# print ('plotloc file',__file__,os.path.dirname(__file__))
# plot.ion()

def read_map(fileIn,verbose=False,force=False,fullout=False):
    """Read FITS map and convert to Healpix if required
    Inputs:
        * fileIn [string]: fits filename
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * force [string, optional]: set to force re-conversion of map (rather than only if converted file can't be read)
        * fullout [boolean, optional]: set to output more variables. Default=False
    Outputs:
        * if fullout==True: [list] [output map, multi-order maps, reordered multi-ordered maps]
        * if fullout==False: [array] output map
    """
    hdr=fits.getheader(fileIn,ext=1)
    if hdr.get('ORDERING')=='RING' or hdr.get('ORDERING')=='NESTED':
        # read as normal map
        hpmap=hp.read_map(fileIn)
        return(hpmap)
    maps={}
    maps_re={}
    maptot=[]
    hpmap=[]
    if hdr.get('ORDERING')=='NUNIQ':
        if fileIn.find('multiorder')>=0:
            fileOut=fileIn.replace('multiorder','healpix')
        else:
            fileOut=fileIn.replace('.fits','.healpix.fits')
        fileRead=False
        if os.path.isfile(fileOut):
            try:
                hpmap=hp.read_map(fileOut)
                if verbose:
                    print('Reading converted map: {}'.format(fileOut))
                fileRead=True
                return hpmap
            except:
                fileRead=False
        if not fileRead or force:
            if verbose:
                print('Converting multiorder map: {}'.format(fileIn))
            skymap=Table.read(fileIn)
            # get pixel number for each map
            level,ipix=ah.uniq_to_level_ipix(skymap['UNIQ'])
            nside=ah.level_to_nside(level)
            skymap.add_columns([Column(nside),Column(ipix)],names=['NSIDE','IPIX'])
            #
            maxns=np.max(nside)
            pixarea=hp.nside2pixarea(maxns)
            npixtot=hp.nside2npix(maxns)
            # create maps at each level
            for row in skymap:
                ns=row['NSIDE']
                nsstr='{}'.format(ns)
                pa=hp.nside2pixarea(ns)
                if not nsstr in maps:
                    maps[nsstr]=np.zeros(hp.nside2npix(ns))
                maps[nsstr][row['IPIX']]=row['PROBDENSITY']*pa
            # reorder maps
            hpmap=np.zeros(npixtot)
            for n in maps:
                maps_re[n]=hp.ud_grade(maps[n],maxns,order_in='NESTED',order_out='RING',power=-2)
                hpmap += maps_re[n]
            hp.write_map(fileOut,hpmap,overwrite=True)
            if verbose:
                print('Writing converted map: {}'.format(fileOut))
    if fullout:
        return(hpmap,maps,maps_re,maptot)
    else:
        return hpmap

def getSuperevents(verbose=False):
    """get list of superevents from GraceDB (used for plotall)
    Inputs:
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [list] list of events
    """
    # get events list from GraceDB
    superevents=[]
    service_url = 'https://gracedb.ligo.org/api/'
    if verbose: print('Retrieving GraceDB data from {}'.format(service_url))
    client = GraceDb(service_url,force_noauth=True)
    events = client.superevents('far < 1.0e-4')
    results=[]
    for event in events:
        ev = event['superevent_id']
        results.append(ev)
    return results

def getSuperevent(ev,verbose=False):
    """get specific superevent from GraceDB
    Inputs:
        * ev [string]: string name
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [list] list of events
    """
    # get event properties from GraceDB
    service_url = 'https://gracedb.ligo.org/api/'
    client = GraceDb(service_url,force_noauth=True)
    evreq=client.superevent(ev)
    event=json.loads(evreq.read())
    if verbose: print('getting files for {}'.format(ev))
    freq=client.files(ev)
    files=json.loads(freq.read())
    event['files'] = files

    mapsrch=['LALInference1','LALInference','bayestar1','bayestar']
    mapFound=False
    m=0
    while not mapFound and m<len(mapsrch):
        mapfile='{}.fits.gz'.format(mapsrch[m])
        if mapfile in files:
            event['mapfile']=[mapfile,event['files'][mapfile]]
            mapFound=True
            # if verbose: print('  found {}'.format(mapfile))
        m+=1
    return event

def getMapFile(urlIn,dirOut='data/',fileOut='',verbose=True,overwrite=False):
    """Get mapFile from location
    Inputs:
        * urlIn [string]: url to read map from (must be HTTP)
        * dirOut [string, optional]: directory to output map to. Default='data/'
        * fileOut [string, optional]: file to output map to. Default='dkyloc.fits'
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * overwrite [boolean, optional]: set to for overwrite file, even if it exists. Default=False
    Outputs:
        * HTTP status code
    """
    # download map file from GraceDB (if not already present)
    if fileOut=='':
        fileOut='skyloc.fits'.format()
    fOut=os.path.join(dirOut,fileOut)
    # check if file exists
    fExists=os.path.isfile(fOut)
    if fExists:
        if overwrite:
            if verbose:print('Overwriting file: {}'.format(fOut))
        else:
            if verbose:print('File already exists: {}'.format(fOut))
            return()
    if verbose:print('Downloading file from {}'.format(urlIn))
    mapreq=requests.get(urlIn)
    if verbose: print('status code:',mapreq.status_code)
    if verbose: print('Saving file to {}'.format(fOut))
    f=open(fOut, 'wb')
    f.write(mapreq.content)
    f.close()
    return(mapreq.status_code)


def getConstBounds(fIn=None,verbose=False):
    """Read constellation boundaries from file
    Inputs:
        * fIn [string, optional]: file to read map from. Default='constdata/constellations.bounds.json'
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [object] object containing constellation boundaries
    """
    # get constellation boundaries from file
    if not fIn:
        fIn=os.path.join(os.path.dirname(__file__),'constdata/constellations.bounds.json')
    cBounds=json.load(open(fIn))['features']
    const={}
    for c in cBounds:
        coords=c['geometry']['coordinates'][0]
        ra=[]
        dec=[]
        for xy in coords:
            ra.append(xy[0])
            dec.append(xy[1])
        ra.append(coords[0][0])
        dec.append(coords[0][1])
        const[c['id']]={'coord':coords,'ra':ra,'dec':dec}
    return(const)

def plotConstBounds(color='w',alpha=1,verbose=False,lw=1):
    """Plot constellation boundaries using matplotlib
    Inputs:
        * color [string, optional]: colour to use when plotting lines. Default='w'
        * alpha [float, optional]: alpha opacity to use when plotting lines. Default=1
        * lw [float, optional]: Line width to use when plot lines. Default=1
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * None
    """
    # plot constellation boundaries
    cb=getConstBounds(verbose=verbose)
    for c in cb:
        const=cb[c]
        nra=len(const['ra'])
        for i in range(nra-1):
            line=hp.projplot(const['ra'][i:i+2],const['dec'][i:i+2],lonlat=True,color=color,linewidth=lw,alpha=alpha)
            if len(line)>1:
                # catch lines that don't plot properly
                # print(c,i,nra,len(line))
                # for l in range(len(line)):
                #     print(l,line[l].get_xydata())
                xy0=line[1].get_xydata()[0]
                xy1=line[2].get_xydata()[0]
                if np.sign(xy0[0])==np.sign(xy1[0]):
                    plot.plot([xy0[0],xy1[0]],[xy0[1],xy1[1]],color=color,linewidth=lw,alpha=alpha)
                # else:
                #     print('***')


def getConstLines(fIn=None,verbose=False):
    """Read constellation lines from file
    Inputs:
        * fIn [string, optional]: file to read map from. Default='constdata/constellations.lines.json'
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [object] object containing constellation lines
    """
    # get constellation lines from files
    if not fIn:
        fIn=os.path.join(os.path.dirname(__file__),'constdata/constellations.lines.json')
    cLines=json.load(open(fIn))['features']
    const={}
    for c in cLines:
        coords=c['geometry']['coordinates']
        ra=[]
        dec=[]
        for i in range(len(coords)):
            ra.append([])
            dec.append([])
            for xy in coords[i]:
                ra[i].append(xy[0])
                dec[i].append(xy[1])
        const[c['id']]={'coord':coords,'ra':ra,'dec':dec}
    return(const)

def plotConstLines(color='y',colorstars='y',alpha=1,plotstars=True,verbose=False,lw=1,ms=2):
    """Plot constellation lines using matplotlib
    Inputs:
        * color [string, optional]: colour to use when plotting lines. Default='y'
        * colorstars [string, optional]: colour to use when plotting stars. Default='y'
        * alpha [float, optional]: alpha opacity to use when plotting lines. Default=1
        * plotstars [boolean, optional]: set to plot star locations. Default=False
        * lw [float, optional]: Line width to use when plotting lines. Default=1
        * ms [float, optional]: Marker size to use when plotting stars. Default=1
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * None
    """
    # plot constellation lines
    cb=getConstLines(verbose=verbose)
    for c in cb:
        const=cb[c]
        # if verbose:print('plotting {}:'.format(c),const['ra'],const['dec'])
        for i in range(len(const['ra'])):
            hp.projplot(const['ra'][i],const['dec'][i],lonlat=True,color=color,linewidth=lw,alpha=alpha)
            if plotstars:
                hp.projplot(const['ra'][i],const['dec'][i],'o',lonlat=True,color=colorstars,markersize=ms,alpha=alpha)
    return

def plotGrid(dRA=45,dDec=30,color='w',alpha=0.5,ls=':'):
    """Plot grid on map using matplotlib
    Inputs:
        * dRA [string, optional]: RA spacing of grid lines [in degrees]. Default=45
        * dDec [string, optional]: Dec spacing of grid lines [in degrees]. Default=30
        * color [string, optional]: colour to use when plotting lines. Default='w'
        * alpha [float, optional]: alpha opacity to use when plotting lines. Default=0.5
        * ls [string, optional]: Line style to use when plotting lines. Default=':'
    Outputs:
        * None
    """
    print('plotting grid')
    for ra in range(-180,180,dRA):
        print(ra)
        decl=np.arange(-90,90,1)
        ral=[ra]*len(decl)
        line=hp.projplot(ral,decl,lonlat=True,color=color,alpha=alpha,ls=ls)
        if len(line)>1:
            xy0=line[1].get_xydata()[0]
            xy1=line[2].get_xydata()[0]
            if np.sign(xy0[0])==np.sign(xy1[0]):
                plot.plot([xy0[0],xy1[0]],[xy0[1],xy1[1]],color=color,linewidth=1,alpha=alpha,ls=ls)
    for dec in range(-90,90,dDec):
        print(dec)
        ral=np.arange(-180,180,1)
        decl=[dec]*len(ral)
        line=hp.projplot(ral,decl,lonlat=True,color=color,alpha=alpha,ls=ls)
        if len(line)>1:
            xy0=line[1].get_xydata()[0]
            xy1=line[2].get_xydata()[0]
            if np.sign(xy0[0])==np.sign(xy1[0]):
                plot.plot([xy0[0],xy1[0]],[xy0[1],xy1[1]],color=color,linewidth=1,alpha=alpha,ls=ls)
    return

def getConstLabs(fIn=None,verbose=False):
    """Read constellation labels from file
    Inputs:
        * fIn [string, optional]: file to read map from. Default='constdata/constellations.json'
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [object] object containing constellation lines
    """
    # get constellation labels from input file
    if not fIn:
        fIn=os.path.join(os.path.dirname(__file__),'constdata/constellations.json')
    constellations=json.load(open(fIn))['features']
    const={}
    for c in constellations:
        const[c['id']]={}
        for p in c['properties']:
            const[c['id']][p]=c['properties'][p]
            const[c['id']]['coord']=c['geometry']['coordinates']
    return(const)

def plotConstLabs(color='w',alpha=1,verbose=False,radeclim=None,maxdist=None,plotcentre=[0,0],fontsize=10):
    """Plot constellation labels using matplotlib
    Inputs:
        * color [string, optional]: colour to use when plotting text. Default='y'
        * alpha [float, optional]: alpha opacity to use when plotting lines. Default=1
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * radeclim [list, optional]: list of RA/Dec limits of plot [RAmin,RAmax,Decmin,Decmax] (in degrees). default=[-180,180,-90,90]
        * maxdist [float, optional]: max distance at which to plot labels (in degrees)]. Default=180
        * plotcentre [list, optional]: coordinates of plot centre [RA,Dec]. Default=[0,0]
        * fontsize [string, optional]: fontsize to use when plotting text. Default=10
    Outputs:
        * None
    """
    # plot constellation labels, restricted to RA/Dec range if set
    if not radeclim:
        radeclim=[-180,180,-90,90]
    if maxdist==None:
        print('setting maxdist to 180')
        maxdist=180
        limplot=False
    else:
        limplot=True
    if plotcentre==None:
        plotcentre=[0,0]
    cl=getConstLabs(verbose=verbose)
    for c in cl:
        const=cl[c]
        ra=const['coord'][0]
        dec=const['coord'][1]
        angdist=hp.rotator.angdist([ra,dec],plotcentre,lonlat=True)*180/np.pi
        angdist=np.min([angdist[0],90])
        # if ra>radeclim[0] and ra<radeclim[1] and dec>radeclim[2] and dec<radeclim[3]:
        if not limplot or angdist<maxdist:
            # point is in range
            # print(ra,dec,angdist,maxdist)
            # if verbose:print('printing {:s} [{:.2f},{:2f} : {:2f}<{:.2f}]'.format(const['name'],ra,dec,angdist,maxdist))
            hp.projtext(ra,dec,const['name'],lonlat=True,color=color,alpha=alpha,ha='center',va='center',fontsize=fontsize)
        # elif verbose:print('not printing {:s} [{:.2f},{:2f} : {:2f}>{:.2f}]'.format(const['name'],ra,dec,angdist,maxdist))
    return

def smoothMap(map,fwhm=None,verbose=False):
    """Smooth map (NOT USED)
    Inputs:
        * map [array]: array of map to smooth
        * fwhm [float, optional]: FWHM of smoothing (in degrees). Default=1.
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [array] smoothed map
    """
    # smooth map do specified FWHM
    if fwhm==None:
        fwhm=1
    fwhmrad=fwhm*np.pi/180
    if verbose: print('smoothing map by {} deg ({:.2f} radians)'.format(fwhm,fwhmrad))
    mapsm=hp.smoothing(map,fwhm=fwhmrad)
    return(mapsm)

def readMap(fileIn,dirIn='data/',smooth=0,overwrite=False,verbose=False):
    """read msp from file. Smooth if necessary and export.
    Inputs:
        * fileIn [string]: filename of map to read in
        * dirIn [string, optional]: directory containing map
        * smooth [float, optional]: FWHM of smoothing (in degrees). Default=0 (no smoothing)
        * overwrite [boolean, optional]: set to overwrite smoothed map if it exists
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [array] map (smoothed if required)
    """
    # read in map file, and smooth if requested
    # writes out smoothed file to disc
    if smooth==0:
        map=hp.read_map(os.path.join(dirIn,fileIn))
        return(map)
    else:
        smoothFile=os.path.join(dirIn,fileIn.split('.fits')[0]+'_sm{}.fits'.format(smooth))
        fwhmrad=smooth*np.pi/180
        fExists=os.path.isfile(smoothFile)
        if fExists and not overwrite:
            if verbose:print('Loading file: {}'.format(smoothFile))
            mapsm=hp.read_map(smoothFile)
            return(mapsm)
        else:
            map=hp.read_map(os.path.join(dirIn,fileIn))
            if verbose: print('smoothing map by {} deg ({:.2f} radians)'.format(smooth,fwhmrad))
            mapsm=hp.smoothing(map,fwhm=fwhmrad)
            if verbose:print('Writing smoothed file: {}'.format(smoothFile))
            hp.write_map(smoothFile,mapsm)
            return(mapsm)

def getProbMap(map,prob=0.9,verbose=False):
    """Get cumulative probability map and calculate area containing given localisation probability. Map pixels contains integrated probability from most probable pixel to least
    Inputs:
        * map [array]: map to read in
        * prob [float, optional]: Probabilty to use to calculate area. Default=0.9
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [list] (map, area within which prob=90% [deg^2]
    """
    # get map of area > prob
    # returns area and area within limit
    mapP=np.ones_like(map)
    # mapP[::]=hp.UNSEEN
    nside=hp.get_nside(map)
    pixarea=hp.nside2pixarea(nside,degrees=True)
    isort=np.argsort(map)[::-1]
    probtot=0
    foundP=False
    i=0
    while not foundP:
        # print(i,isort[i],map[isort[i]])
        probtot += map[isort[i]]
        if probtot>prob:
            if not foundP:
                iP=i
                if verbose: print('{} pixels in {:d}% area'.format(i,int(prob*100)))
                foundP=True
        i += 1
        mapP[isort[i]]=probtot
    areaP=(iP+1)*pixarea

    if verbose: print('{:d}% area={:d} deg^2'.format(int(prob*100),int(areaP)))

    return(mapP,areaP)

def getArea(maptot,prob=0.9,verbose=False):
    """Get area within which likelihood > given probability
    Inputs:
        * maptot [array]: total probability map
        * prob [float, optional]: Probabilty to use to calculate area. Default=0.9
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [float] area within which prob=90% [deg^2]
    """
    npix=len(np.where(maptot<=prob)[0])
    pixarea=hp.nside2pixarea(hp.get_nside(maptot),degrees=True)
    areaP=npix*pixarea
    if verbose: print('{:d}% area={:d} deg^2'.format(int(prob*100),int(areaP)))
    return(areaP)

def getRaDecRange(map,lim=0.5,ltype='max',border=0,verbose=False):
    """Get RA/Dec range with cumulative probability > limit
    Inputs:
        * map [array]: total probability map
        * lim [float, optional]: Probabilty to use to calculate area. Default=0.5
        * ltype [string, optional]: 'min' or 'max' depending on looking for min or max. Default='max'
        * border [float, optional]: UNUSED
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [float] area within which prob=90% [deg^2]
    """
    # get RA/Dec limits of all pixels which satisfy criteria
    # returns min(RA), max(RA), min(Dec), max(Dec)
    if ltype=='max':
        # get all pixels less than limit
        pixlim=np.where(map<lim)[0]
    elif ltype=='min':
        # get all pixels greater than limit
        pixlim=np.where(map>lim)[0]
    # convert HP pixels to RA/Dec
    nside=hp.get_nside(map)
    ra,dec=hp.pix2ang(nside, pixlim, lonlat=True)
    # keep in range [-180,180]
    ra=np.where(ra>180,ra-360,ra)
    # find minima
    decMin=np.min(dec)
    decMax=np.max(dec)
    raMin=np.min(ra)
    raMax=np.max(ra)
    # if raMax-raMin > 359:
    #     # crosses meridian
    #     # ra=np.where(ra>180,ra-360,ra)
    #     raMin=np.min(ra)
    #     raMax=np.max(ra)

    # get ra/dec of peak of map
    raPeak,decPeak=getPeak(map,getmin=True,verbose=verbose)
    # calculate angular distance to corners/edges of map and find max distance
    maxdist=0
    for ra in [raMin,raMax,np.mean([raMin,raMax])]:
        for dec in [decMin,decMax,np.mean([decMin,decMax])]:
            dist=hp.rotator.angdist([ra,dec],[raPeak,decPeak],lonlat=True)
            try:
                maxdist=np.max([dist,maxdist])
            except:
                print('WARNING: error calculating max of [distance {}, maxdistance {}] for ra {} dec {}'.format(dist,maxdist,ra,dec))
    maxdist=np.min([maxdist*180./np.pi,90])

    if verbose: print('RA range=[{:.1f},{:.1f}] ; Dec range=[{:.1f},{:.1f}]'.format(raMin,raMax,decMin,decMax))
    if verbose: print('maximum angular distance: {:.1f}'.format(maxdist))
    return([raMin,raMax,decMin,decMax,maxdist])

def makeGrid(radecLim,nX=4096,nY=2048,verbose=False):
    """Create RA/Dec grid over area (for contour plotting)
    Inputs:
        * radecLim [list]: Values to use for grid [raMin,raMax,decMin,decMax,maxdist]
        * nX [float, optional]: Number of grid points in X dimension across whole RA range. Default=4096
        * nY [float, optional]: Number of grid points in Y dimension across whole Dec rance. Default=4096
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [object] [x grid locations, y grid locations
    """
    # make RA/Dec grid over spcified RA/Dec area
    # returns x,y coords of grid
    raMin,raMax,decMin,decMax,maxdist=radecLim

    # get grid of RA, Dec within range
    xxAll=np.arange(-180,180,360/nX)+0.5
    yyAll=np.arange(-90,90,180/nY)+0.5

    # restrict to RA/Dec range
    xx=xxAll[np.where((xxAll<=raMax) & (xxAll>=raMin))]
    yy=yyAll[np.where((yyAll<=decMax) & (yyAll>=decMin))]
    # if len(np.where(xx>180)[0]):
    #     xx2=np.zeros_like(xx)
    #     no180=len(np.where(xx>180)[0])
    #     xx2[0:no180]=xx[np.where(xx>180)]-360
    #     xx2[no180:]=xx[np.where(xx<=180)]
    if verbose: print('Created {} x {} grid'.format(len(xx),len(yy)))

    return(xx,yy)

def map2grid(map,xx,yy,verbose=False):
    """Convert map to x-y grid (for contour plotting). Warning - this is SLOW!
    Inputs:
        * map [array]: map to grid
        * xx [array]: x grid-points
        * yy [array]: y grid-points
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [array] gridded map of size (len(xx),len(yy))
    """
    # transform map to x-y grid
    # returns map regridded onto grid
    if verbose: print('Populating grid...')
    nside=hp.get_nside(map)
    nX=len(xx)
    nY=len(yy)
    zz=np.zeros((nY,nX))
    zz[::]=np.nan
    for i in range(nX):
        x=xx[i]
        for j in range(nY):
            y=yy[j]
            ipix=hp.ang2pix(nside,x,y,lonlat=True)
            zz[j,i]=map[ipix]

    return(zz)

def getContLines(xx,yy,zz,level=0.9,ax=None,verbose=False):
    """get contour lines for a gridded map
    Inputs:
        * xx [array]: x grid-points
        * yy [array]: y grid-points
        * zz [array]: gridded map values
        * level [float, optional]: Level of contour to calculate. Default=0.9
        * ax [Axes, optional]: Mapplotlib Axes object. Default = current aces [plot.gca()]
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [object] [RA contour points, Dec contour points, contour object]
    """
    # get contour lines and return as RA/Dec points
    if ax==None:
        ax=plot.gca()
    if verbose:print('Getting contour lines for level:',level)
    cont=ax.contour(xx,yy,zz,levels=level,alpha=0)
    coll=cont.collections[0]
    segs=coll.get_segments()
    # nseg=len(segs)
    raCont=[]
    decCont=[]
    ns=0
    ncp=0
    for seg in segs:
        raCont.append([])
        decCont.append([])
        for s in range(len(seg)):
            raCont[ns].append(seg[s][0])
            decCont[ns].append(seg[s][1])
            ncp+=1
        ns+=1
    if verbose:print('Calculate {} points in contour'.format(ncp))

    return(raCont,decCont,cont)

def plotMap(map,proj='moll',fig=None,pmax=None,pmin=None,rot=None,cmap=None,cbg=None,verbose=False,half_sky=False,zoomrng=None,title=None,margins=None,notext=None,border=None):
    """plot probabilty map with set options
    Inputs:
        * map [array]: map to plot
        * proj [string, optional]: projection for healpix [moll|cart|orth|gnom]. Default='moll'
        * fig [Figure, optional]: matplotlib Figure object. Default=None (make new figure)
        * pmax [float, optional]: Max value to plot. Default=None (plot all values)
        * pmin [float, optional]: Min value to plot. Default=None (plot all values)
        * rot [array, optional]: Coordinates to rotate to centre. Default=None (no rotation)
        * cmap [cmap object, optional]: Colourmap to use. Default=cmap.gray
        * cbg [string, optional]: Background colour to use. Default='w'
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * half_sky [boolean, optional]: set to only plot half the sky in orthview. Default=False
        * zoomrng [array, optional]: coordinates of zoom range [RAmin,RAmax,Decmin,Decmax,maxdist]. Default=None (no zooming)
        * title [string, optional]: title to add to plot. Default=None (no title)
        * margins [array, optional]: margins to add aroung plot [left,right,top,bottom], as passed to healpy ????view. Default=None (default mardins)
        * notext [boolean, optional] set to not plot any text. Default=None (healpy default)
        * border [string, optional] colour of border around plot. Default=None (no border)
    Outputs:
        * [Figure] matplotlib Figure object containing figure
    """
    # plot map based on options specified
    if not cmap:
        cmap=cm.gray
    if not cbg:
        cmap.set_under('w')
        cmap.set_over('w')
        cmap.set_bad('w')
    else:
        cmap.set_under(cbg)
        cmap.set_over(cbg)
        cmap.set_bad(cbg)
    if fig==None:
        if proj=='cart' and zoomrng!=None:
            figsize=(10,10)
        else:
            figsize=(20,10)
        fig=plot.figure(figsize=figsize)

    if proj=='cart':
        if zoomrng==None:
            # plot full-sky Cartesian view (rotated to centre if specified)
            fig=hp.cartview(map,cmap=cmap,max=pmax,min=pmin,cbar=False,rot=rot,fig=fig.number,title=title,notext=notext,bgcolor=cbg,margins=margins)
        else:
            # get the centre of the plot averate of the ra and dec limits
            ramean,decmean=getPeak(map,verbose=verbose)
            maxdist=np.min([zoomrng[4],90])
            rarng=[-maxdist,maxdist]
            decrng=[-maxdist,maxdist]
            #
            # ramean=np.mean(zoomrng[0:1])
            # decmean=np.mean(zoomrng[2:3])
            # # calculate RA/Dec ranges from centre and add 10deg border around edge
            # rawid=zoomrng[1]-zoomrng[0]
            # decht=zoomrng[3]-zoomrng[2]
            # rarng=[np.max([zoomrng[0]-ramean-10,-180]),np.min([zoomrng[1]-ramean+10,180])]
            # decrng=[np.max([zoomrng[2]-decmean-10,-90]),np.min([zoomrng[3]-decmean+10,90])]
            # ralim=[rarng[0]+ramean,rarng[1]+ramean]
            # declim=[decrng[0]+decmean,decrng[1]+decmean]


            # plot Cartesian view centred and zoomed on RA/Dec range
            fig=hp.cartview(map,cmap=cmap,max=pmax,min=pmin,cbar=False,rot=[ramean,decmean],fig=fig.number,
                lonra=rarng,latra=decrng,title=title,notext=notext,bgcolor=cbg,margins=margins)
    elif proj=='moll':
        # plot full-sky Mollweide view (rotated to centre if specified)
        fig=hp.mollview(map,cmap=cmap,max=pmax,min=pmin,cbar=False,rot=rot,fig=fig.number,title=title,
            notext=notext,bgcolor=cbg,margins=margins)
        if border!=None:
            theta= np.arange(0, 181) * np.pi/180 
            hp.projplot(theta, theta * 0 + 0.9999 * np.pi, ls="-",color=border, lw=3, direct=True)
            hp.projplot(theta, theta * 0 - np.pi, ls="-",color=border, lw=3, direct=True)
    elif proj=='orth':
        # plot full-sky Mollweide view (rotated to centre if specified)
        fig=hp.orthview(map,cmap=cmap,max=pmax,min=pmin,cbar=False,rot=rot,half_sky=half_sky,fig=fig.number,title=title,
        notext=notext,bgcolor=cbg,margins=margins)
    elif proj=='gnom':
        fig=hp.gnomview(map,cmap=cmap,max=pmax,min=pmin,cbar=False,rot=rot,half_sky=half_sky,fig=fig.number,title=title,
        notext=notext,bgcolor=cbg,margins=margins)
    else:
        print('unknown projection:',proj)
    return fig

def getPeak(map,verbose=False,getmin=False):
    """get peak of map
    Inputs:
        * map [array]: map to plot
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * getmin [boolean, optional]: set to return minimum location. Default=False (find maximum)
    Outputs:
        * [list] [RA, Dec] (in degrees)
    """
    # get RA/Dec of peak
    # set getmin to find minimum of map, otherwise fint max [default=max]
    # returns RA and Dec
    if getmin:
        iPeak=np.argmin(map)
        pval=np.min(map)
    else:
        iPeak=np.argmax(map)
        pval=np.max(map)
    ns=hp.get_nside(map)
    raPeak,decPeak=hp.pix2ang(ns,iPeak,lonlat=True)
    if verbose: print('peak [{:.2f}] found at ({:.2f},{:.2f})'.format(pval,raPeak,decPeak))
    return(raPeak,decPeak)

def plotLab(map,txt,color='r',verbose=False):
    """plot label at peak of map (on current axes)
    Inputs:
        * map [array]: map to plot
        * txt [string]: text to plot
        * color [string, optional]: colour of text to plot. Default='r'
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * None
    """
    # plot label at Peak of map
    # Depracated
    raPeak,decPeak=getPeak(map,verbose=verbose)
    hp.projtext(raPeak,decPeak,txt,lonlat=True,color=color)
    return

def plotContours(map,level=0.9,color='w',alpha=0.5,linestyle='-',linewidth=2,verbose=False):
    """plot contours at level on map (on current axes). Warning - this is SLOW
    Inputs:
        * map [array]: map to plot contours on
        * level [float, optional]: contour level to plot
        * color [string, optional]: colour of contour to plot. Default='w'
        * alpha [float, optional]: alpha opacity of contour line. Default=0.5
        * linestyle [string, optional]: linestyle of contour line. Default='-'
        * linewidth [string, optional]: linewidth of contour line. Default=2
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * [object] contour object
    """
    # plot contours of map
    if verbose:print('plotting contours for level:',level)
    radec95=getRaDecRange(map,lim=0.92,ltype='max',verbose=verbose)
    xx95,yy95=makeGrid(radec95,verbose=verbose)
    zz95=map2grid(map,xx95,yy95,verbose=verbose)
    raCont,decCont,cont=getContLines(xx95,yy95,zz95,level=level,verbose=verbose)
    nseg=len(raCont)
    for s in range(nseg):
        hp.projplot(raCont[s],decCont[s],lonlat=True,color=color,alpha=alpha,linestyle=linestyle,linewidth=linewidth)
    return(cont)

def plotGravoscope(mapIn,fileIn='',cmap=cm.gray,pngOut='',jpgOut='',res=4,verbose=False,coord='G'):
    """plot high-res Cartesian map
    Inputs:
        * mapIn [array, optional]: map create tiles of
        * fileIn [string, optional]: file to read map from (if map not valid healpy map).
        * cmap [cmap, optional]: Matplotlib Colour map to use. Default=cmap.gray
        * pngOut [string, optional]: file to output PNG map to. Default='' (no PNG file)
        * jpgOut [string, optional]: file to output JPG map to. Default='' (no PNG file)
        * res [integer, optional]: resolution to plot, with image 1024*res wide. Default=4
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * coord [string, optional]: coordinate system to plot, converted from 'C' (passed to cartview). Default='G' (Galactic)
    Outputs:
        * None
    """
    try:
        nside=hp.get_nside(mapIn)
        T=mapIn
    except:
        if verbose:print('reading map from {}'.format(fileIn))
        T = hp.read_map(fileIn)
        nside=hp.get_nside(T)

    sky = np.zeros((int(1024*res),int(1024*res/2)),dtype=np.float32)
    tmp = np.zeros((1024,1024),dtype=np.float32)

    dlon=(360/res)
    dlat=dlon
    for i in np.arange(res):
        lon_off = (i-(res/2)+0.5)*dlon
        for j in np.arange(int(res/2)):
            lat0=-90 + j*dlon
            lat1=lat0+dlon
            if verbose:print('plotting {}x{} of {}x{} : [{},{}] - [{}:{}]'.format(
                i+1,j+1,res,int(res/2),-lon_off-dlon/2,lat0,-lon_off+dlon/2,lat1))
            tmp[...] = np.transpose(hp.cartview(T,coord=['C',coord],return_projected_map=1,
                xsize=1024,ysize=1024,lonra=[-dlon/2,dlon/2],latra=[lat0,lat1],rot=[-lon_off,0]))
            # print(i,j,i*1024,(i+1)*1024,j*1024,(j+1)*1024)
            sky[i*1024:(i+1)*1024,j*1024:(j+1)*1024] = tmp

        plot.close('all')
    plot.figure(figsize=(10.24*res,10.24*res/2))
    plot.figimage(np.flipud(np.transpose(sky)),cmap=cmap)

    if pngOut!='':
        if verbose:print('saving to',pngOut)
        plot.savefig(pngOut)
    if jpgOut!='':
        # jpgOut=pngOut.replace('.png','.jpg')
        if verbose:print('saving to',jpgOut)
        plot.savefig(jpgOut)

    return

def makeTilesPerl(fileIn,verbose=False):
    """make tileset using perl script (slow)
    Inputs:
        * fileIn [string, optional]: filename to load map from
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * None
    """
    cutterfile=os.path.join(os.path.dirname(__file__),'cutter.pl')
    if verbose:
        v=1
    else:
        v=0
    print('verbose',verbose,v)
    os.system('perl {} cutter.pl file="{}" minzoom=3 maxzoom=6 ext="png" verbose={}'.format(cutterfile,fileIn,v))

    return

def makeTiles(map,dirOut='data/gravoscope/test',cmap=None,vmin=None,vmax=None,
        maxres=3,minres=2,verbose=False):
    """make tileset using python script
    Inputs:
        * map [array]: map to plot
        * dirOut [string, optional]: directory to output tiles to. Default='data/gravoscope/test'
        * cmap [cmap object, optional]: matplotlibt colormap to use. Default=cmap.gray
        * vmin [float, optional]: minimum value to plot. Default=None (no minimum value)
        * vmax [float, optional]: maximum value to plot. Default=None (no maximumvalue)
        * maxres [integer, optional]: maximum resolution to output tiles to. Default=3
        * minres [integer, optional]: minimum resolution to output tiles to. Default=2
        * verbose [boolean, optional]: set to for verbose output. Default=False
    Outputs:
        * None
    """
    if not os.path.exists(dirOut):
        os.mkdir(dirOut)
    if vmax==None:
        vmax=np.max(map)
    if vmin==None:
        vmin=np.min(map)
    if verbose:print('min={}; max={}'.format(vmin,vmax))
    if cmap==None:
        cmap=cm.gray
    loc=[0]*maxres
    name=['']*maxres
    qrst=['q','r','s','t']
    loc[0]=0
    name[0]='t'

    for res in range(1,maxres+1):
        nimg=4**(res)
        loc=[0]*(res+1)
        name=['']*(res+1)
        loc[0]=0
        name[0]='t'
        dra=360/(2**res)
        ddec=360/(2**res)
        if verbose:print('level {} ({} tiles) [{}x{}]'.format(res,nimg,dra,ddec))
        for x in range(nimg):
            xb4=np.base_repr(x,4)
            if len(xb4) < res:
                xb4=''.join(['0']*(res-len(xb4)))+xb4
            ra0=-180
            dec0=-180
            for y in range(len(xb4)):
                dy=360/(2**(y+1))
                loc[y+1]=int(xb4[y])
                name[y+1]=qrst[loc[y+1]]
                ra0=ra0+dy*np.mod(int((loc[y+1]+1)/2),2)
                dec0=dec0+dy*(1-int(loc[y+1]/2))
            fileOut=''.join(name)+'.png'
            if dec0>=-90 and dec0+ddec<=90:
                # if verbose:print('{} ({}) [{},{}]'.format(x,fileOut,ra0,dec0))
                img=np.transpose(hp.cartview(map,coord=['C','G'],return_projected_map=1,
                    xsize=256,ysize=256,lonra=[-dra,0],latra=[dec0,dec0+ddec],rot=[-ra0,0]))
                plot.figure(figsize=(2.56,2.56))
                plot.figimage(np.flipud(np.transpose(img)),cmap=cmap,vmin=vmin,vmax=vmax)
                plot.savefig(os.path.join(dirOut,fileOut))
                plot.close('all')
            # else:
            #     if verbose:print('[{} ({}) [{},{}]]'.format(x,fileOut,ra0,dec0))
    return

def makePlot(ev='S190412m',mapIn=None,proj='moll',plotcont=False,smooth=0.5,zoomlim=0.92,rotmap=True,
    half_sky=False,pngOut=None,verbose=False,cbg=None,
    dirData='data/',minzoom=10,pngSize=3000,thumbOut=None,margins=None,
    thumbSize=300,title=None,RAcen=0,grid=False,addCredit=True,addLogos=False,notext=True,
    plotbounds=True,plotlines=True,plotlabels=True,border=None,lw=1,fontsize=10):
    """make tileset using python script
    Inputs:
        * ev [string, optional]: superevent ID to read if no map provided. Default='S190412m'
        * mapIn [string or array, optional]: map to read in (filename [string] or HEALPix map). If not provided, tries to get event with superevent ID provided in ev from GraceDB. Default=None (load superevent)
        * proj [string, optional]: projection (moll=Mollweide [Default], cart=Cartesian)
        * plotcont [boolean, optional]: set to plot contours (default=False)
        * smooth [float, optional]: degrees to smooth probability densith map to get contours. 0=no smoothing. [Default=0.5deg]
        * zoomlim [float, optional]: probability to zoom map in to (cartview only). Default=0.92
        * rotmap [boolean, optional]: rotate map to centre on peak value [Default=False]
        * half_sky [boolean, optional]: only show half the sky (orthographic only). Default=False
        * pngOut [string, optional]: filename to export to PNG. Default=None (no export)
        * pngSize [integer, optional]: width of output image in pixels [default=3000]
        * thumbOut [string, optional]: filename to export to thumbnail PNG. Default=None (-> no export)
        * thumbSize [integer, optional]: width of output thumbnail image in pixels. Default=100
        * verbose [boolean, optional]: set to for verbose output. Default=False
        * cbg [string, optional]: background colour [Default=black]
        * dirData [string, optional]: directory to load files from [Default= 'data/']
        * minzoom [float, optional]: minimum map radius, in degrees (proj=cart & zoomlim!=None only) [Default=10]
        * title [string, optional]: Title to plot on image (optional string) [Default = <ev>]
        * RAcen [string, optional]: Centre RA (Default=0; applies only if rotmap not set)
        * grid [boolean, optional]: set to plot grid (Default=False)
        * addCredit [boolean, optional]: set to plot credit line (Default=True)
        * plotbounds [boolean, optional]: set to plot constellation boundaries (default=True)
        * plotlabels [boolean, optional]: set to plot constellation labels (default=True)
        * plotlines [boolean, optional]: set to plot constellation liness (default=True)
    Outputs:
        * [array] map that was plotted
    """
    if type(mapIn)==type(None):
        event=getSuperevent(ev,verbose=verbose)
        event['mapfile_local']=fileOut='{}_{}'.format(ev,event['mapfile'][0])
        getMapFile(event['mapfile'][1],fileOut=event['mapfile_local'],verbose=verbose,dirOut=dirData)
        map=readMap(event['mapfile_local'],smooth=smooth,verbose=verbose,dirIn=dirData)
    elif type(mapIn)==str:
        map=readMap(mapIn,smooth=smooth,verbose=verbose,dirIn=dirData)
    else:
        map=mapIn
    maptot,area95=getProbMap(map,0.95,verbose=verbose)

    if title==None:
        title=ev

    blankmap=np.zeros_like(maptot)
    if rotmap:
        raPeak,decPeak=getPeak(map,verbose=verbose)
        rot=[raPeak,decPeak]
    elif RAcen!=0:
        print('centering on {}'.format(RAcen))
        rot=[RAcen,0]
    else:
        rot=None
    if minzoom==None:
        minzoom=180
    if proj=='cart' and zoomlim!=None:
        radeczoom=getRaDecRange(maptot,lim=zoomlim,ltype='max',verbose=verbose)
        # apply min zoom
        if radeczoom[4] < minzoom:
            print('setting min zoom [{}->{}]'.format(radeczoom[4],minzoom))
            radeczoom[4]=minzoom
        maxdist=radeczoom[4]
        # print('zoom range',radeczoom)
        # radeclab=getRaDecRange(maptot,lim=zoomlim,ltype='min',verbose=verbose,border=0)
        # ralim=[radeclab[0],radeclab[1]]
        # declim=[radeclab[2],radeclab[3]]
    else:
        radeczoom=None
        maxdist=None
        ralim=[-180,180]
        declim=[-90,90]
    # print (radec95)
    fig=plotMap(map,cmap=cm.hot,proj=proj,rot=rot,verbose=verbose,zoomrng=radeczoom,title=title,half_sky=half_sky,cbg=cbg,notext=notext,margins=margins,border=border)

    if plotcont:
        cont90=plotContours(maptot,level=0.9,color='w',alpha=0.5,linestyle='-',linewidth=2,verbose=verbose)

    if proj=='cart' and zoomlim!=None:
        alphaLab=1
    else:
        alphaLab=0.5

    if grid:
        plotGrid(dRA=45, dDec=30)
    if plotbounds:
        plotConstBounds(color=(0.5,0.5,0.5),verbose=verbose,alpha=0.5,lw=lw)
    if plotlabels:
        plotConstLabs(color=(0,0.7,0.7),verbose=verbose,alpha=alphaLab,maxdist=maxdist,plotcentre=rot,fontsize=fontsize)
    if plotlines:
        plotConstLines(color=(0,0.7,0.7),verbose=verbose,alpha=0.5,lw=lw)

    if addCredit:
        credit='Credit: LIGO-Virgo/Cardiff Uni./C. North'
        text=plot.figtext(0.99,0.01,credit,color='black',ha='right')
    if addLogos:
        # get height/width
        figsize=plot.gcf().get_size_inches()
        aspect=figsize[1]/figsize[0]
        logosize=[0.05*aspect,0.05]
        logo1=os.path.join(os.path.dirname(__file__),'logos/cardiffuniversitylogo-spot_500px.jpg')
        im1=plot.imread(logo1)
        axlogo1=plot.gcf().add_axes([0.01,0.01,logosize[0],logosize[1]],anchor='SW')
        axlogo1.imshow(im1)
        axlogo1.axis('off')
        logo2=os.path.join(os.path.dirname(__file__),'logos/LigoVirgo14.jpg')
        im2=plot.imread(logo2)
        axlogo2=plot.gcf().add_axes([logosize[0]+0.02,0.01,logosize[0],logosize[1]],anchor='SW')
        axlogo2.imshow(im2)
        axlogo2.axis('off')
    if pngOut:
        plot.savefig(pngOut,dpi=pngSize/10)
    if thumbOut:
        plot.savefig(thumbOut,dpi=thumbSize/10)

    return map

###########################
def plotall():
    """test script to plot lots of maps from superevents
    Inputs:
        * None
    Outputs:
        * None
    """
    evs=['S190421ar']
    evs=[]
    try:
        evs[0]
    except:
        evs=getSuperevents()
    print ('plotting {} event(s):'.format(len(evs)),evs)
    # dirData='../data/'
    for ev in evs:
        map=makePlot(ev=ev,proj='cart',plotcont=False,smooth=0.5,zoomlim=0.8,rotmap=True,verbose=True,
            pngOut='png/{}_cartzoom.png'.format(ev))
        makePlot(ev=ev,mapIn=map,proj='cart',plotcont=False,smooth=0.5,zoomlim=None,rotmap=False,verbose=True,pngOut='png/{}_cart.png'.format(ev))
        makePlot(ev=ev,mapIn=map,proj='moll',plotcont=False,smooth=0.5,zoomlim=None,rotmap=False,verbose=True,pngOut='png/{}_moll.png'.format(ev),cbg='w')
        # map = makePlot(ev=ev,proj='orth',plotcont=False,smooth=0.5,zoomlim=None,rotmap=True,half_sky=True,verbose=True,pngOut='png/{}_orth.png'.format(ev))
    plot.show()
