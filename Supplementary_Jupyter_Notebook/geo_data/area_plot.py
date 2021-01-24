#%%
import matplotlib.pyplot as plt
import numpy as np
import obspy
import matplotlib.patheffects as PathEffects
from mpl_toolkits.basemap import Basemap
import sys
from load_fault_traces import load_fault_traces
import Load_Data_notebooktopy as Load_Data
from netCDF4 import Dataset
from matplotlib.colors import LightSource
import pickle 
class empty:
    pass

color_faults = False

db = pickle.load(open('../stored_variables/'+'dbresults','rb'))
#%%
## Don't delete these commented lines. May need them again later. 
if False: 
    # Dictionary with Boolean arrays to keep track of cluster-fault name relationshp. This is only preliminary! If the cluster labels have changed, the script should be able to keep track of that. 
    # Can either create hypo_faults now if the db results are pretty good, then manually determine the label number to fault name relationship, and save the file, 
    # or load the previously determined hypo_faults variable. 
    hypo_faults = {'Paganica'      :db.labels==0,
                    'Campotosto'    :db.labels==1,
                    'Cittareale'    :db.labels==7,
                    'Mt. San Franco':db.labels==3,
                    'Noise'         :db.labels==-1}
    pickle.dump(hypo_faults, open('../stored_variables/preliminary_clusters','wb'))
hypo_faults = pickle.load(open('../stored_variables/preliminary_clusters','rb'))

# Used the now commented code untill 09/10/2020
# # Keep track of ALL colors while will be used and added to the legend (fault traces and aftershocks)
# fault_color_dict = {
#                     'Paganica':'#32a852',  
#                     'Campotosto':'#4f4fff', 
#                     'Cittareale':'#c74c00',
#                     'Mt. San Franco':'#ff8c00',
#                     # bellow address faults which are only present in fault traces
#                     'Mt. Stabiata':'#a900c7',
#                     'Mt. Castellano-Colle Enzano':'#f2ff00',
#                     'Aterno, San Demetrio systems' : '#03fcf0',

#                     # 'San Demetrio':'#ff0000',
#                     #
#                     'Surface rupture':'#ff0000',
#                     'Noise'     :'#E6E6E6',
#                     }#'#757575'}

# trace_to_label = {'M.Stabbiata':'Mt. Stabiata', # This line is shared with area_plot.py
#                   'M.S.Franco':'Mt. San Franco',
#                   'Aternosys.2':'Paganica',
#                   'Aternosys.6':'Paganica',
#                   'Aternosys.7':'Paganica',
#                   'Aternosys.8':'Paganica',
#                   'Campotosto':'Campotosto',
#                   'Cittareale':'Cittareale',
#                   'Mt.Castellano':'Mt. Castellano-Colle Enzano',
#                   'C.lleEnzano':'Mt. Castellano-Colle Enzano',
#                   'Aternosys.1':'Aterno, San Demetrio systems',
#                   'Aternosys.1':'Aterno, San Demetrio systems',
#                   'Aternosys.3':'Aterno, San Demetrio systems',
#                   'Aternosys.4':'Aterno, San Demetrio systems',
#                   'Aternosys.5':'Aterno, San Demetrio systems',
#                   'Aternosys.9':'Aterno, San Demetrio systems',
#                   'Aternosys.10':'Aterno, San Demetrio systems',
#                   'SanDemetrio':'Aterno, San Demetrio systems',
#                   'D.Demetrio':'Aterno, San Demetrio systems',
#                   }

# Keep track of ALL colors while will be used and added to the legend (fault traces and aftershocks)
psd = 'Paganica - San Demetrio'
fault_color_dict = {
                    'Paganica':'#32a852',  
                    psd: '#03fcf0',
                    'Campotosto':'#4f4fff', 
                    'Cittareale':'#c74c00',
                    'Mt. San Franco':'#ff8c00',
                    # bellow address faults which are only present in fault traces
                    'Mt. Stabiata':'#a900c7',
                    'Mt. Castellano - Colle Enzano':'#f2ff00',

                    # 'San Demetrio':'#ff0000',
                    #
                    'Surface rupture':'#ff0000',
                    'Noise'     :'#E6E6E6',
                    }#'#757575'}
trace_to_label = {'M.Stabbiata':'Mt. Stabiata', # This line is shared with area_plot.py
                  'M.S.Franco':'Mt. San Franco',
                  'Aternosys.2':'Paganica',
                  'Aternosys.6':'Paganica',
                  'Aternosys.7':'Paganica',
                  'Aternosys.8':'Paganica',
                  'Campotosto':'Campotosto',
                  'Cittareale':'Cittareale',
                  'Mt.Castellano':'Mt. Castellano - Colle Enzano',
                  'C.lleEnzano':'Mt. Castellano - Colle Enzano',
                  # 'Aternosys.1':psd,
                  'Aternosys.3':psd,
                  'Aternosys.4':psd,
                  'Aternosys.5':psd,
                  'Aternosys.9':psd,
                  'Aternosys.10':psd,
                  'SanDemetrio':psd,
                  'D.Demetrio':psd,


                  }

with open('traceinfo', 'wb') as file:
    pickle.dump( {'fault_color_dict':fault_color_dict, 'trace_to_label':trace_to_label}, open('traceinfo', 'wb'))
    file.close()

    
def get_cluster_color(labels, label, hypo_faults, fault_color_dict):
    similarity = np.zeros(len(hypo_faults)) # similarity between currently observed cluster and each previously named cluster. 
    name = np.zeros(len(hypo_faults), dtype = object)
    for ikey, key in enumerate(hypo_faults.keys()):
        similarity[ikey] = (hypo_faults[key] * (labels == label)).sum()/(hypo_faults[key].sum()) # Similarity is basically normalized zero lag cross correlation between the boolean aray of this cluster and each previously named cluster
        if similarity[ikey] <= 0:
            name[ikey] = None
        elif similarity[ikey] > 0:
            name[ikey] = key
    if (similarity>0).any():
        fault_name = name[similarity == similarity.max()][0] 
    else:
        fault_name = 'not present'
    color_db    = fault_color_dict.get(fault_name, None)
    return fault_name, color_db
#%%
plt.rcParams.update({'font.size': 16})
if True: 
    width = 80000 # 2000000
    height = 80000 # 2000000
    centlat = 42.5
    centlon = 13.5
    
    llcrnrlat = 42.15 
    llcrnrlon = 13.1
    urcrnrlat = 42.7
    urcrnrlon = 13.65
    
    #%%
    fig, ax = plt.subplots(figsize = (9, 9))
    
    map = Basemap(projection='aeqd',
               lon_0 = centlon,
               lat_0 = centlat,
               # width = width,
               # height = height,
                llcrnrlon = llcrnrlon ,
                 llcrnrlat = llcrnrlat,
                 urcrnrlon = urcrnrlon ,
                 urcrnrlat = urcrnrlat,
              resolution = 'i',
                  ax = ax)
    
    mers = np.arange(12, 15, 0.2)
    map.drawmeridians(mers, labels = np.ones(mers.size), linewidth = 0, fontsize = 10)
    pars = np.arange(40, 44, .2)
    map.drawparallels(pars, labels = np.ones(pars.size), linewidth = 0, fontsize = 10)
    
    if True:
        def get_xyz(fname = 'italy-90m.grd'):
            top = Dataset(fname)
            x_range = top.variables['x_range'][:]
            y_range = top.variables['y_range'][:]
            z_range = top.variables['z_range'][:]
            spacing = top.variables['spacing'][:]
            dimension=top.variables['dimension'][:]
            x = np.arange(x_range[0], x_range[1]+spacing[0], spacing[0])
            y = np.arange(y_range[0], y_range[1]+spacing[1], spacing[1])
            z = top.variables['z'][:]
            # z = z.reshape(y.size, x.size).T
            # y = y[::-1]
            
            z = z.reshape(y.size, x.size).T
            
            z = z[::,::-1]
    
            
            xboo = (x < 14) * (x>13)
            yboo = (y < 43) * (y > 42) 
            xl = x[xboo]
            yl = y[yboo]
            zl = z[xboo,:][:, yboo]
        
            
            d=1
            xl = xl[::d]
            yl = yl[::d]
            zl = zl[::d][:,::d]
            
            
            return xl, yl, zl
        
        lont, latt, zt = get_xyz()
        zt = zt.T
        lont, latt = np.meshgrid(lont, latt)
        
        # n = 1000
        # nx = 1 + int( (map.xmax-map.xmin)/n )
        # ny = 1 + int( (map.ymax-map.ymin)/n )
        nx = 2000; ny = 3000
        topodat = map.transform_scalar(zt, lont[0,:], latt[:,0], nx, ny)
        # im = map.imshow(topodat * 0 + 1, cmap=plt.get_cmap('Greys'), alpha = 1,interpolation = 'bilinear', rasterized=True)
        im = map.imshow(topodat, cmap=plt.get_cmap('terrain'), alpha = 1,interpolation = 'bilinear', rasterized=True)
    
    
        xt, yt = map(lont, latt)
        ls = LightSource(azdeg=315, altdeg=45)
        hs = ls.hillshade(zt, 
                    vert_exag = 2, dx=xt[0,1] - xt[0,0] , dy=yt[0,0] - yt[1,0])
        hs = map.transform_scalar(hs, lont[0,:], latt[:,0], nx, ny)
        map.imshow(hs, cmap = 'gray', alpha = .5,interpolation = 'bilinear', rasterized=True)
#%%
    
    # plt.xlim([0, width]); plt.ylim([0, height])
    # map.drawmapboundary(fill_color=None)
    # map.fillcontinents(color='lightgrey',lake_color=None)
    # map.shadedrelief()
    # map.etopo(scale = 1)
    # map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
    # map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)



    x, y = map(db.lon, db.lat)
    
    # Plot clusters
    # Pain in the ass script to decide which fault is which. Look for the clusters which are most similar to previously named clusters. 
    zorder_dict = {'Noise':.11, 'Campotosto':.12, 'Mt. San Franco':.13, 'Paganica':.14} # Need to know which clusters to plot first, to determine which points are covered by other points. 

    # point_labels = np.zeros(len(db.labels), dtype = object)
    # point_colors = np.zeros(db.labels.shape, dtype = object)
    for ilabel, label in enumerate(np.unique(db.labels)):
        fault_name, color_db = get_cluster_color(db.labels, label, hypo_faults, fault_color_dict)
        this = db.labels == label
        colors = 0 * np.array([1,1,1])# color_db 
        plt.scatter(x[this], y[this], #label = fault_name,
                    s = .3, color = colors, 
                    zorder = zorder_dict.get(fault_name, 100),
                    marker = 'o')#, linewidths = .2, 
#%%
    # plt.figure(4)
    # dat2 = {} 
    # (dat2['lat'], dat2['lon'], dat2['depth'], dat2['MW'], dat2['radius'], 
    #     dat2['ST1'], dat2['DIP1'], dat2['RK1'],
    #      __, __, __, __, __) = Load_Data.loadINGV(use_MW=False) # Load data AGAIN so we can get it with MW instead of ML


    scogdat = np.array([ # lon, lat, depth, strike, dip, rake, Mw # From Scognamiglio et al., 2010 and Chiaraluce et al., 2011
    [13.385371,	    42.349353,	9,	139,	48,	-87,	6.08], # done
    [13.396759,	    42.475676,	9,	154,	57,	-80,	4.98], # LONGITUDE: 13.396759 LATITUDE : 42.475676 DEPTH : 8.781 km
    [13.47844,	    42.309921,	18,	338,	73,	-58,	5.37], # done
    [13.364849,	    42.502472,	9,	322,	46,	-95,	5.21], # MAGNITUDE: 5.2 Mw LONGITUDE: 13.364849 LATITUDE : 42.502472 DEPTH : 8.888 km
    [13.371161,	    42.513332,	6,	137,	48,	-86,	4.96], # LONGITUDE: 13.371161 LATITUDE : 42.513332 DEPTH : 7.602 km 
    [13.206, 	    42.570,     5, 	143, 	52, -76, 	3.48],
    ]) 
    lon, lat, __, st1, dp1, rk1, mags = scogdat.T
         
    # fig = plt.figure(); print('delete')
    from obspy.imaging.beachball import beach
    # great_mag = dat2['MW']>=5.8
    # print(np.arange(great_mag.size)[great_mag])
    # pltfoc = np.logical_or(great_mag, dat2['MW'] == 3.44) 
    focs = np.array([st1, dp1, rk1]).T# * 180 / np.pi 
    # focs = focs[pltfoc,:]
    x, y = map(lon, lat)
    # mags = dat2['MW'][pltfoc]
    scale_beach = 2000 / max(mags)
    for i in range(len(x)):
        iwidth = mags[i] * scale_beach
        beach1 = beach(focs[i], xy = (x[i], y[i]), width = iwidth, linewidth = .1, facecolor = .4 * np.array([1,1,1])  )
        ax = plt.gca()#.add_collection(beach1)
        ax.add_collection(beach1)
        # plt.annotate(str(mags[i]), (x[i], y[i]), size = 30, color = '#ffffff', zorder = 15)
    # ax.autoscale(); print('delete')

    # Legend... Seems unecessary but was requested by reviewer. 
    # xleg, yleg = map(13.6, 42.6)
    xbds = np.array(ax.get_xlim())
    ybds = np.array(ax.get_ylim()) 
    xleg = (xbds[1] - xbds[0]) * .83 + xbds[0] 
    yleg = (ybds[1] - ybds[0]) * .92 + ybds[0]

    beachl1 = beach(np.array([-45, 45, -90]), xy=(xleg, yleg+2000), width = 3.5 * scale_beach, linewidth = .1, facecolor = .4 * np.array([1,1,1]) )
    ax.add_collection(beachl1)
    plt.text(xleg+800, yleg+2000, 'Mw=3.5', horizontalalignment = 'left', verticalalignment = 'center', fontsize = 11)
    
    beachl2 = beach(np.array([-45, 45, -90]), xy=(xleg, yleg-.1), width = 5 * scale_beach, linewidth = .1, facecolor = .4 * np.array([1,1,1]) )
    ax.add_collection(beachl2)
    plt.text(xleg+800, yleg, 'Mw=6.0', horizontalalignment = 'left', verticalalignment = 'center', fontsize = 11)

        
#%%

    color_faults = False
    
    # Plot surface ruptures
    mf, sr = load_fault_traces()  
    for ifault, fault in enumerate(sr):
        x, y = map(fault['lons'], fault['lats'])
        plt.plot(x, y, c = 'r', linewidth = 5, zorder = 10)
        
    # Plot mapped faults. 
    for ifault, fault in enumerate(mf):
        x, y = map(fault['lons'], fault['lats'])
        # print(fault['name'])
        fault_name = trace_to_label.get(fault['name'],'not present') if color_faults else 'do_not_color_this_fault'
        plt.plot(x, y, color = fault_color_dict.get(fault_name, '#ffffff'), linewidth = 1.5, zorder = 10.1)
        
        x, y = map(fault['lonlabel'], fault['latlabel'])
        # plt.annotate(fault['name'], [x,y] , size = 12, color = '#ffffff')


    # Make a legend: only really concerned about portraying colors, so do a trick that makes a dot legend entry for each color used. 
    fault_color_dict.pop('Noise', None)
    markers = [plt.Line2D([-2000000,0],[-200000,0],color=color, marker='o', linestyle='') for color in fault_color_dict.values()]

    if color_faults:
        plt.legend(markers, fault_color_dict.keys(), numpoints=1, prop={'size': 9}, loc = 1)
    
    # xticks = np.arange(13, 14, .2)
    # xtickpos, __ = map(xticks, np.ones(xticks.shape) * centlat)
    # plt.xticks(xtickpos, xticks)
    
    
    # xax = ax.get_xaxis()
    

    map.drawmapscale(13.2,42.185, 13.2, 42.185, 10, fontsize = 11, barstyle = 'simple')
    
    if True:
        # ax.set_rasterized(False)
        plt.savefig('../figures/map.pdf', dpi = 500)

    plt.show()

    