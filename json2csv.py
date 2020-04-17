# -*- coding: utf-8 -*-
import json
import argparse
import pandas as pd

parser=argparse.ArgumentParser(prog="plotApps", \
    description="Converting JSON to CSV file")
parser.add_argument('-i','--in', dest='filein', type=str, default='events.json', help='Input file (default=../json/events.json)')
parser.add_argument('-o','--out', dest='fileout', type=str, default='events.csv', help='Input file (default=../csv/events.csv)')
parser.add_argument('-d','--datadict', dest='datadict', type=str, default=None, help='Datadict file (default: "" (included in input file))')


args=parser.parse_args()
filein=args.filein
fileout=args.fileout
datadict=args.datadict

eventsIn=json.load(open(filein))
data=eventsIn['data']
if datadict:
    ddIn=json.load(open(datadict))
else:
    ddIn=eventsIn['datadict']

dOut={'unit':{},'description':{}}
series={}
cols=['Mchirp','M1','M2','Mtotal','Mfinal','Mratio','DL','z','deltaOmega','chi','af','Erad','lpeak','rho','FAR','UTC','GPS']
for e in data:
    dOut[e]={}
    for d in cols:
        if ddIn[d].has_key('unit_en'):
            dOut['unit'][d]=ddIn[d]['unit_en']
        if ddIn[d].has_key('name_en'):
            dOut['description'][d]=ddIn[d]['name_en']
        if not data[e].has_key(d):
            continue
        if data[e][d].has_key('best'):
            dOut[e][d]=data[e][d]['best']
            if data[e][d].has_key('err'):
                dOut[e][d+'_errtype']='+/- (90%)'
                dOut[e][d+'_errp']=data[e][d]['err'][0]
                dOut[e][d+'_errm']=data[e][d]['err'][1]
            elif data[e][d].has_key('lim'):
                dOut[e][d+'_errtype']='+/- (range)'
                dOut[e][d+'_errp']=data[e][d]['lim'][0]
                dOut[e][d+'_errm']=data[e][d]['lim'][1]
        elif data[e][d].has_key('lower'):
            dOut[e][d]=data[e][d]['lower']
            dOut[e][d+'_errtype']='lower_limit'
        elif data[e][d].has_key('upper'):
            dOut[e][d]=data[e][d]['upper']
            dOut[e][d+'_errtype']='upper'
series['unit']=pd.Series(dOut['unit'],index=dOut['unit'].keys())
series['description']=pd.Series(dOut['description'],index=dOut['description'].keys())

rows=['description','unit']
for e in data:
    series[e]=pd.Series(dOut[e],index=dOut[e].keys())
    rows.append(e)

df=pd.DataFrame(series)
dfT=df.T

colsAll=[]
for c in cols:
    # print c
    try:
        dfT[c]
        colsAll.append(c)
    except:
        pass
    try:
        dfT[c+'_errtype']
        colsAll.append(c+'_errtype')
    except:
        pass
    try:
        dfT[c+'_errp']
        colsAll.append(c+'_errp')
    except:
        pass
    try:
        dfT[c+'_errm']
        colsAll.append(c+'_errm')
    except:
        pass
# print colsAll
pd.DataFrame(dfT,index=rows).to_csv(fileout,columns=colsAll,encoding='utf8')