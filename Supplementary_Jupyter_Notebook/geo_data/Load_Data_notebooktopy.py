#!/usr/bin/env python
# coding: utf-8

# # Load data

# In[1]:


from os import path
import requests


# In[14]:


# def loadINGV(minMagnitude=None, onlyFocal = False, 
# #              fname = 'LaqAll.final', # on 04/02/2019, adding 49 missing events 'Brennan_Aquila2009_51k_FocMec_CMT.dat.final',
#              fname = 'stored_variables/LaqAlldf.final', # on 04/03/2019. Added slip and radius calculated from "creating_final_catalog"
#              normalRet = True):
#     if minMagnitude is not None:
#         print("Not ready for minMagnitude in loadINGV (after getting new dataset)")
#         raise NotImplementedError()
#     if onlyFocal != False:
#         print("Not ready for onlyFocal in loadINGV (after getting new dataset)")
#         raise NotImplementedError()
    
    
#     # load data from file into panda dataframe. Extract data from there. 
#     panDat = pd.read_csv(fname, delimiter= '\t', skipinitialspace=True)
#     values = panDat.values.T
    
#     tt        = values[0 ].astype('datetime64')
#     latitude  = values[1 ].astype('float64')
#     longitude = values[2 ].astype('float64')
#     depth     = values[3 ].astype('float64') * 1000 # this becomes negative only in selectedINGV function
#     ML        = values[4 ].astype('float64')                        #calling ml for now, not sure if mw
#     ID        = values[5 ].astype('int'    )
#     ST1       = values[6 ].astype('float64') * np.pi / 180
#     DIP1      = values[7 ].astype('float64') * np.pi / 180
#     RK1       = values[8 ].astype('float64') * np.pi / 180
#     ST2       = values[9 ].astype('float64') * np.pi / 180
#     DIP2      = values[10].astype('float64') * np.pi / 180
#     RK2       = values[11].astype('float64') * np.pi / 180
#     fty       = values[12].astype('str'    )
#     td        = values[13].astype('float64')
#     tla       = values[14].astype('float64')
#     tlo       = values[15].astype('float64')
#     radius    = values[16].astype('float64')
#     SLIP      = values[17].astype('float64')
    
#     # Get CMT from Zenodo file, work it into the Valoroso dataset. Also take the magnitudes from CMT set. 
#     st1, dip1, rk1, st2, dip2, rk2, id, ml = loadCMT()
#     for i, (_st1, _dip1, _rk1, _st2, _dip2, _rk2, _id, _ml) in enumerate(
#         zip(st1, dip1, rk1, st2, dip2, rk2, id, ml)):
#         this = ID == _id 
#         ST1 [this] = _st1 
#         DIP1[this] = _dip1 
#         RK1 [this] = _rk1 
#         ST2 [this] = _st2 
#         DIP2[this] = _dip2
#         RK2 [this] = _rk2
#         ML  [this] = _ml
    
#     # strike of plane 1 is actually dip direction when fty=='MP'. This is FPFIT convention?
#     changeStrike = fty == 'MP'
#     ST1[changeStrike] += - np.pi / 2
#     ST2[changeStrike] += - np.pi / 2
    
# #     #can run the following test to show that if and only if
# #     #    all focal data == 0, there is no focal mechanism
# #     Foc1 = fty!='0'
# #     Foc2 = ~( (ST1==0) * (ST2==0) * (RK1==0) * (RK2==0) * (DIP1==0) * (DIP2==0) )
# #     print(np.sum(Foc1!=Foc2)) 
    
#     if normalRet:
#         return [latitude, longitude, depth, ML, radius, ST1, 
#             DIP1, RK1, ST2, DIP2, RK2, tt, SLIP]
#     else:
#         return [latitude, longitude, depth, ML, radius, ST1, 
#             DIP1, RK1, ST2, DIP2, RK2, tt, SLIP, fty]
def loadINGV(minMagnitude=None, onlyFocal = False, 
#              fname = 'LaqAll.final', # on 04/02/2019, adding 49 missing events 'Brennan_Aquila2009_51k_FocMec_CMT.dat.final',
#              fname = 'stored_variables/LaqAlldf.final', # on 04/03/2019. Added slip and radius calculated from "creating_final_catalog"
                          fname = '../stored_variables/Valoroso_et_al_2013.csv', # on 04/03/2019. Added slip and radius calculated from "creating_final_catalog"
             normalRet = True, use_MW = False):
    
    # Download data if it doesn't exist. 
    if not path.exists(fname): 
        url = 'https://zenodo.org/record/4036248/files/Valoroso_et_al_2013.csv?download=1' 
        r = requests.get(url, allow_redirects=True)
        open(fname, 'wb').write(r.content)
    
    
    if minMagnitude is not None:
        print("Not ready for minMagnitude in loadINGV (after getting new dataset)")
        raise NotImplementedError()
    if onlyFocal != False:
        print("Not ready for onlyFocal in loadINGV (after getting new dataset)")
        raise NotImplementedError()
    
    
    # load data from file into panda dataframe. Extract data from there. 
    panDat = pd.read_csv(fname, delimiter= '\;', skipinitialspace=True)
    values = panDat.values.T

    tt        = values[1 ].astype('datetime64')
    latitude  = values[2 ].astype('float64')
    longitude = values[3 ].astype('float64')
    depth     = values[4 ].astype('float64') * 1000 # this becomes negative only in selectedINGV function
    ML        = values[5 ].astype('float64')                        #calling ml for now, not sure if mw
    ID        = values[0 ].astype('int'    )
    flag      = values[6 ].astype('int'    )
    
    reject_flag0=True
    if reject_flag0:
        keep = flag == 1 
        
        tt        = tt       [keep]
        latitude  = latitude [keep]
        longitude = longitude[keep] 
        depth     = depth    [keep]
        ML        = ML       [keep] 
        ID        = ID       [keep] 
        flag      = flag     [keep]
                
    ST1 = np.zeros(tt.shape)
    DIP1 = np.zeros(tt.shape)
    RK1 = np.zeros(tt.shape)
    ST2 = np.zeros(tt.shape)
    DIP2 = np.zeros(tt.shape)
    RK2 = np.zeros(tt.shape)
    
    fty = np.zeros(tt.shape, dtype = object) 
    radius = np.zeros(tt.shape) # Can calculate this here if desired; not doing now. 
    SLIP = np.zeros(tt.shape) # Can calculate this here if desired; not doing now. 
    
    # Get CMT from Zenodo file, work it into the Valoroso dataset. Also take the magnitudes from CMT set. 
    st1, dip1, rk1, st2, dip2, rk2, id, ml, ftycmt = loadCMT(use_MW = use_MW)
    for i, (_st1, _dip1, _rk1, _st2, _dip2, _rk2, _id, _ml, _fty) in enumerate(
        zip(st1, dip1, rk1, st2, dip2, rk2, id, ml, ftycmt)):
        this = ID == _id 
        ST1 [this] = _st1 
        DIP1[this] = _dip1 
        RK1 [this] = _rk1 
        ST2 [this] = _st2 
        DIP2[this] = _dip2
        RK2 [this] = _rk2
        ML  [this] = _ml
        fty [this] = _fty
    
    # strike of plane 1 is actually dip direction when fty=='MP'. This is FPFIT convention?
    changeStrike = fty == 'MP'
    ST1[changeStrike] += - np.pi / 2
    ST2[changeStrike] += - np.pi / 2
        
#     #can run the following test to show that if and only if
#     #    all focal data == 0, there is no focal mechanism
#     Foc1 = fty!='0'
#     Foc2 = ~( (ST1==0) * (ST2==0) * (RK1==0) * (RK2==0) * (DIP1==0) * (DIP2==0) )
#     print(np.sum(Foc1!=Foc2)) 
    if normalRet:
        return [latitude, longitude, depth, ML, radius, ST1, 
            DIP1, RK1, ST2, DIP2, RK2, tt, SLIP]
    else:
        return [latitude, longitude, depth, ML, radius, ST1, 
            DIP1, RK1, ST2, DIP2, RK2, tt, SLIP, fty]
# loadINGV()


# In[13]:


import pandas as pd
import numpy as np
def loadCMT(
            minMagnitude=None,
            onlyFocal = False,
            fname = '../stored_variables/Brennan_Aquila2009_FocMec_CMT.dat.final.mw.zenodo',
            normalRet = True, 
            use_MW = False):
    
    # Download data if it doesn't exist. 
    if not path.exists(fname): 
        try: 
            url = '' # Add url here. 
            r = requests.get(url, allow_redirects=True)
            open(fname, 'wb').write(r.content)
        except: 
            raise(Exception('Note from Brennan: the Focal mechanism dataset is not downloaded! Please go to the function loadCMT and add the appropriate Zenodo url for that dataset'))

    if minMagnitude is not None:
        print("Not ready for minMagnitude in loadINGV (after getting new dataset)")
        raise NotImplementedError()
    if onlyFocal != False:
        print("Not ready for onlyFocal in loadINGV (after getting new dataset)")
        raise NotImplementedError()

    # load data from file into panda dataframe. Extract data from there. 
    panDat = pd.read_csv(fname, delimiter= ' ', skipinitialspace=True)
    values = panDat.values.T

    spaces = values[0].astype(object).copy()
    spaces[:] = ' '
    tt        = (values[0] + spaces + values[1]).astype('datetime64')
    latitude  = values[1 +1].astype('float64')
    longitude = values[2 +1].astype('float64')
    depth     = values[3 +1].astype('float64') * 1000 # this becomes negative only in selectedINGV function
    ML        = values[4 +1].astype('float64')                        #calling ml for now, not sure if mw
    ID        = values[5 +1].astype('int'    )
    ST1       = values[6 +1].astype('float64') * np.pi / 180
    DIP1      = values[7 +1].astype('float64') * np.pi / 180
    RK1       = values[8 +1].astype('float64') * np.pi / 180
    ST2       = values[9 +1].astype('float64') * np.pi / 180
    DIP2      = values[10+1].astype('float64') * np.pi / 180
    RK2       = values[11+1].astype('float64') * np.pi / 180
    fty       = values[12+1].astype('str'    )
    MW        = values[13+1]; MW[MW=='Null'] = 0; MW = MW.astype('float64')
    
    
    if use_MW:
        ML = MW # Not good coding but works. 
        
    return ST1, DIP1, RK1, ST2, DIP2, RK2, ID, ML, fty

# loadCMT()


# In[11]:


# # This both loads the data from the INGV file,
# # and it calculates slip from the provided radius and mL

# def loadINGV(minMagnitude=0, onlyFocal = False): 
#     """Returns: (latitude, longitude, depth,
#     mL, radius, ST1, DIP1, RK1, ST2, DIP2, RK2, time, SLIP)
    
#     Provide maxMagnitude to load only earthquake with higher than that magnitude
#     Only focal is boolean to indicate if earthquake with no focal mechanisms are desired
    
#     Angles are returned in radians
#     slip and depth are in meters
#     """

#     # We've had some problems with the address of the file. This should fix it.
#     try:
#         file=open('LAquila_2009_ALLinONE_unpub.out')
#     except:
#         url = '/work/Course'#= os.getcwd()
#         file=open(url+'/LAquila_2009_ALLinONE_unpub.out')
        
#     maxLines = sum(1 for lines in file)
#     file.seek(0)
#     splitlines = []
#     for i in range(maxLines):
#         if i == 0:
#             next(file)
#         else:
#             splitlines.append(file.readline().split())
            
#     dataSize=len(splitlines)
#     year=np.zeros(dataSize,np.int)
#     month=np.zeros(dataSize,np.int)    
#     day=np.zeros(dataSize,np.int)  
#     hours=np.zeros(dataSize,np.int)      
#     minutes=np.zeros(dataSize,np.int)    
#     seconds=np.zeros(dataSize,np.int)    
#     latitude=np.zeros(dataSize,np.float)    
#     longitude=np.zeros(dataSize,np.float)  
#     depth=np.zeros(dataSize,np.float)      
#     ML=np.zeros(dataSize,np.float)    
#     err=np.zeros(dataSize,np.float)    
#     radius=np.zeros(dataSize,np.float)    
#     ID=np.zeros(dataSize,np.int)    
#     ST1=np.zeros(dataSize,np.float)   
#     DIP1=np.zeros(dataSize,np.float)   
#     RK1=np.zeros(dataSize,np.float)   
#     ST2=np.zeros(dataSize,np.float)     
#     DIP2=np.zeros(dataSize,np.float) 
#     RK2=np.zeros(dataSize,np.float)  
#     DS=np.zeros(dataSize,np.float)
#     tt=np.zeros(dataSize,dtype='datetime64[s]')
    
#     for data in np.arange(dataSize):
#         latitude[data]=float(splitlines[data][6])
#         longitude[data]=float(splitlines[data][7])    
#         depth[data]=float(splitlines[data][8]) * 1000 # meters
#         ML[data]=float(splitlines[data][9])
#         radius[data]=float(splitlines[data][11]) # meters
#         ST1[data]=float(splitlines[data][13]) * np.pi / 180 
#         DIP1[data]=float(splitlines[data][14]) * np.pi / 180 
#         RK1[data]=float(splitlines[data][15]) * np.pi / 180 
#         ST2[data]=float(splitlines[data][16]) * np.pi / 180 
#         DIP2[data]=float(splitlines[data][17]) * np.pi / 180 
#         RK2[data]=float(splitlines[data][18])  * np.pi / 180 
#         year[data]=float(splitlines[data][0])
#         month[data]=float(splitlines[data][1])
#         day[data]=float(splitlines[data][2])
#         hours[data]=float(splitlines[data][3])
#         minutes[data]=float(splitlines[data][4])
#         seconds[data]=float(splitlines[data][5])
#         d = date(year[data], month[data], day[data])
#         t = time(hours[data], minutes[data])
#         dt = datetime.combine(d, t)
#         tt[data] = dt

#     #Here we calculate displacement based on magnitude and radius. 
#     #I think we calculated to go from a circular fault to a rectangular fault?
#     Mo=10**(1.5*ML+16.1)/1e7
#     SLIP = Mo/(2/3*3e10*radius**2*3.14) #slip is displacement in m
#     #
    
#     #boolean array that selects only the quakes with a rupture mechanism
#     if onlyFocal:
#         selectQuakesST1 = ~( (ST1 == 0) * (ST2==0) * (DIP1==0) * (DIP2==0) * (RK1==0) * (RK2==0) )
#     else:
#         selectQuakesST1 = True
        
#     selectQuakesML = ML>minMagnitude#max magnitude is the highest magnitude we want to load
#     selectQuakes = selectQuakesST1 * selectQuakesML # selectQuakes is boolean for events we want to load
        
#     return[ latitude[selectQuakes].copy(), longitude[selectQuakes].copy(), depth[selectQuakes].copy(),
#             ML[selectQuakes].copy(), radius[selectQuakes].copy(), ST1[selectQuakes].copy(),
#             DIP1[selectQuakes].copy(), RK1[selectQuakes].copy(), ST2[selectQuakes].copy(), DIP2[selectQuakes].copy(),
#             RK2[selectQuakes].copy(), tt[selectQuakes].copy(), SLIP[selectQuakes].copy() ]


# In[52]:


# The point here is to load or save INGV data as .npy for quick data loading and manimuplating

def selectedINGV(minMag = None, onlyFocal = False, reRun = False, path=None):
    """returns (latFault, lonFault, depthFault, mL, rad, ST1, DIP1, RK1, ST2, DIP2, RK2, tt, slip)
    
    minMag: minimum magnitude that is loaded
    onlyFocal: Return only earthquakes with a known focal mechanism
    
    Will access the loadINGV function if .npy file is not already saved, and it will create the .npy file.
    Otherwise, it simply loads the data from .npy file. This should save time.  """

    if not reRun:
        try: 
            mostINGV = loadArrays('INGV_data', path)
            time = loadSerialized('INGV_event_times', path=path, dtype='datetime64[s]')
        except: 
           print("""INGV_data and INGV_event_times were not found in specified path.
               Loading data directly from LAquila_2009_ALLinONE_unpub. """)
           reRun = True
        
    if reRun: 
        mostINGV = loadINGV()
        time = mostINGV.pop(-2)
        saveArrays(mostINGV, 'INGV_data', path)
        saveSerialized(time, 'INGV_event_times', path)
        
    latFault, lonFault, depthFault, mL, rad, ST1, DIP1, RK1, ST2, DIP2, RK2, slip = mostINGV
    depthFault = - depthFault
    
    if minMag is None:
        selectBool = np.ones(latFault.size, dtype = bool)
    else:
        selectBool = mL >= minMag
    if onlyFocal:
        selectBool = selectBool * ~( (ST1==0) * (ST2==0) * (DIP1==0) * (DIP2==0) * (RK1==0) * (RK2==0) )
        
    return( latFault[selectBool], lonFault[selectBool], depthFault[selectBool],
        mL[selectBool], rad[selectBool], ST1[selectBool], DIP1[selectBool],
        RK1[selectBool], ST2[selectBool], DIP2[selectBool], RK2[selectBool],
        time[selectBool], slip[selectBool] )


# In[ ]:


def loadSecondSequence(minMag = None, onlyFocal = False,
                       file = 'focmec-gmt-utm-medi.reloc.dat'):
    table = pd.read_csv(file, delimiter = ' ', skipinitialspace=True)
    latFault = table.LATITUDE.values
    lonFault = table.LONGITUDE.values
    depthFault = table.DEPTH.values * - 1000
    mL =   table.ML.values  
    rad =  np.zeros(latFault.size) * np.nan
    ST1 =  table.STRIKE1.values * np.pi / 180
    ST2 =  table.STRIKE2.values * np.pi / 180
    DIP1 = table.DIP1.values    * np.pi / 180
    DIP2 = table.DIP2.values    * np.pi / 180
    RK1 =  table.RAKE1.values   * np.pi / 180
    RK2 =  table.RAKE2.values   * np.pi / 180
    slip = np.zeros(latFault.size) * np.nan
    
    year = table.YEAR.values
    month = table.MONTH.values
    day = table.DAY.values
    hour = table.HOUR.values
    minute = table.MINUTE.values
    seconds = table.SECONDS.values
    
    datestr = []
    for i in range(year.size):
        if seconds[i] == 60:
            seconds[i] = 0
            minute[i] += 1
        datestr.append('{:04.0f}-'.format(year[i]   ) +
                       '{:02.0f}-'.format(month[i]  ) +
                       '{:02.0f}T'.format(day[i]    ) +
                       '{:02.0f}:'.format(hour[i]   ) +
                       '{:02.0f}:'.format(minute[i] ) +
                       '{:06.3f}' .format(seconds[i])
                      )
                        #aimed format'2009-01-12T20:53:38'
            
    time = np.array(datestr, dtype = 'datetime64')
    
    keep = np.ones(time.shape, dtype = 'bool')
    
    if onlyFocal:
        FA = ~ ((ST1==0) * (ST2==0) * (DIP1==0) * (DIP2==0) * (RK1==0) * (RK2==0) )
        keep = keep * FA
        
    if minMag is not None:
        magPass = mL >= minMag
        keep = keep * magPass

    
    return (latFault[keep], lonFault[keep], depthFault[keep], mL[keep], rad[keep], 
            ST1[keep], ST2[keep], DIP1[keep], DIP2[keep], RK1[keep], RK2[keep],
            time[keep], slip[keep])


# In[2]:


# Loads the GPS displacements and positions that were in Serpelloni's paper.

def loadGPS():
    """return (lon, lat, dx, dy, du, ede, edn, edu) in SI units
    
    Make sure that GPSdata.csv is in the same folder as the jupyter notebook file"""
    
    # There has been inconsistency regarding both encoding and skiprows
    # Some time should be devoted to fixing this in a solid way. 
    try:
        atr = np.loadtxt(fname='GPSdata.csv',delimiter=',',skiprows=1,usecols=(1,2,3,4,7,8,10,11),unpack=True)#, encoding='latin-1') 
    except:
        atr = np.loadtxt(fname='GPSdata.csv',delimiter=',',skiprows=0,usecols=(1,2,3,4,7,8,10,11),unpack=True)#, encoding='latin-1') 
#     except SecondException:
#         atr = np.loadtxt(fname='GPSdata.csv',delimiter=',',skiprows=1,usecols=(1,2,3,4,7,8,10,11),unpack=True) 
#     except ThirdException:
#         atr = np.loadtxt(fname='GPSdata.csv',delimiter=',',skiprows=0,usecols=(1,2,3,4,7,8,10,11),unpack=True) 

        
#     atr=atr.reshape(8,77) 
#     atr = atr[:, (atr[7]!=0)] # remove station with unknown vertical properties
        
    lon = atr[0,]; lat = atr[1,]; dx = atr[2,]; dy = atr[3,]; du = atr[6,] 
    ede = atr[4,]; edn = atr[5,]; edu = atr[7,]
    
    # Error can't be zero. That is null data. I set the error arbitrarily high so that covariance inverse is not affected.
    ede[ede == 0] = 1e6
    edn[edn == 0] = 1e6
    edu[edu == 0] = 1e6
    
        
    #return displacements in meters
    return lon, lat, dx/1000, dy/1000, du/1000, ede/1000, edn/1000, edu/1000


# In[ ]:


# GPS data that I think was used in Cirella et al., from Chelloni.
# I used this to make sure our data was similar. It is. 
def loadGPSChelloni(fileName):
    with open(fileName) as file:
        loaded = pd.read_csv(file, header = None)
        
    values = loaded.values.T
    
    lon = values[1].astype('float64')
    lat = values[2].astype('float64') 
    dx = values[3].astype('float64') / 1000
    edx = values[4].astype('float64') / 1000
    dy = values[5].astype('float64') / 1000
    edy = values[6].astype('float64') / 1000
    dz = values[7].astype('float64') / 1000
    edz = values[8].astype('float64') / 1000
    
    return (lon, lat, dx, dy, dz, edx, edy, edz)


# In[37]:


# def loadGPS():
    
#     data = pd.read_csv( open('GPSdata.csv', encoding='latin-1') ).values
    
#     boo = np.zeros( data[0,:].size, dtype=bool)
    
#     for ind in [1,2,3,4,7,8,10,11]:
#         boo[ind] = True

#     (lon, lat, dx, dy, du, ede, edn, edu) = data.T[boo, :]
    
#     dx=dx/1000
#     dy=dy/1000
#     du=du/1000
#     ede=ede/1000
#     edn=edn/1000
#     edu=edu/1000
    
#     return lon, lat, dx, dy, du, ede, edn, edu


# # Load and save Generated data

# In[4]:


def saveArrays(arrays, fName, path=None):
    """Saves numpy arrays as json files. Pickle is not used.
    
    arrays should be a list or tupple containing all individual arrays. 
    path is a string placed before fName to choose folder."""
    
    if path is not None:
        fName = str(path)+str(fName)
    elif path is None:
        fName = str(fName)
    
    lists = [arr.tolist() for arr in arrays]
    with open(fName, 'w') as outFile:
        json.dump( lists, outFile )

def loadArrays(fName, path=None):
    """Returns what was saved by saveArrays(). 
    
    Returns arrays bundled in a list or tupple. 
    path is a string placed before fName to choose folder"""
    
    if path is not None:
        fName = str(path)+str(fName)
    elif path is None:
        fName = str(fName)  
    
    with open(fName, 'r') as inFile:
        lists = json.load( inFile )
    listArrs = [np.array(arr) for arr in lists]
    
    return listArrs


# In[5]:


# Functions used for storing arrays that are not json saveable, mainly dates.

# Must first convert to string to storage, and later convert back to array.
def toStrings(times):
    """Converts an array or list into a list of strings."""
    stringRep = []
    for data in times:
        stringRep.append(str(data))
    return stringRep

def stringsToArray(strings, dtype = None):
    """Converts a list of strings into an array."""
    Arr=np.zeros( len(strings), dtype=dtype)
    for i in np.arange(len(strings)):
        Arr[i] = strings[i]
    return Arr 

# Functions for saving and loading array of datetime. 
def saveSerialized(array, fName, path=None):
    """Saves array or list of as serialized string."""
    if path is not None:
        fName = str(path)+str(fName)
    elif path is None:
        fName = str(fName)
        
    strings = toStrings(array)
    
    with open(fName, 'w') as outFile:
        json.dump( strings, outFile )
        
def loadSerialized(fName, path=None, retString=False, dtype = None):
    """Loads a string representation of list, returns as array or list."""
    if path is not None:
        fName = str(path)+str(fName)
    elif path is None:
        fName = str(fName)
    
    with open(fName, 'r') as inFile:
        strings = json.loads(inFile.read())

    if retString:
        return strings
    else:
        array = stringsToArray(strings, dtype=dtype)
        return array


# In[ ]:





# In[ ]:




