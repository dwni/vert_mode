#equation

import numpy as np
import datetime
from datetime import datetime, timedelta
from scipy import stats


#function round time
def roundTime(dt=None, roundTo=60):
    """Round a datetime object to any time lapse in seconds
    dt : datetime.datetime object, default now.
    roundTo : Closest number of seconds to round to, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me."""

    if dt == None : dt = datetime.datetime.now()
    seconds = (dt.replace(tzinfo=None) - dt.min).seconds
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + timedelta(0,rounding-seconds,-dt.microsecond)


#time to datetime from xarray
def t2dt(time):
#time = xrs.time
    nt = len(time.values)
    date    = np.zeros((nt), dtype='object')
    for ix in range(0,nt):
        date[ix] = datetime(time[ix].dt.year, time[ix].dt.month, time[ix].dt.day, time[ix].dt.hour, time[ix].dt.minute, time[ix].dt.second)
    return date[:]

def time2slice(date, avg):
    nt=np.arange(0,60,avg)
    d=[]
    date_st = int(date.dt.day[0].values)
    date_et = int(date.dt.day[-1].values)+1
    
    st=int(date.dt.hour[0])
    if date.dt.day[-1]>date.dt.day[0]:
        et=int(date.dt.hour[-1])+1
    elif (date.dt.minute[-1].values)==0:
        et = date.dt.hour[-1].values-1
    else:    
#         et = date.dt.hour[-1].values
        et = 24
    for k in range(date_st,date_et):
        for i in range(st,et):
        # for i in range(24):
            #a=np.where(date.dt.hour==i)
            for j in (nt):
                b=np.where((date.dt.day == k) & (date.dt.hour==i) & (date.dt.minute==j))
                if len(b[0])!=0:
                    #c=a[0][b[0][0]]
                    c=b[0][0]
                    d.append((c))
                else:
                    b=np.where((date.dt.day == k) & (date.dt.hour==i) & (date.dt.minute//j==1))
                    print(b)

                    if len(b[0])!=0:
                        c=b[0][0]
                        d.append((c))
                    else:
                        c=0
                        d.append((c))
    #print(d)
    if (d[0]!=0):
        d=np.append(0,d)
    if (d[-1]!=(len(date)-1)):
        d=np.append(d,len(date)-1)

    return d

