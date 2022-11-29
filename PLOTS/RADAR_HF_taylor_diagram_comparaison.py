#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 13:04:08 2019

@author: tbrivoal
"""
import sys
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 15})
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

from matplotlib import pyplot
from matplotlib import colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.polynomial.polynomial import polyfit
import skill_metrics as sm

# Define radar name, name of the variables etc for each variable (1 radar = 1 column)

list_radar=['BISC2'    ,'GIBR1','GALI2','EBRO' ,'IROI'      , 'TRAD3']
list_var_u=['EWCT'     , 'u'   , 'u'   , 'u'   ,'u'         , 'u']
list_var_v=['NSCT'     , 'v'   , 'v'   ,  'v'  , 'v'        , 'v']
list_time =['TIME'     , 'time', 'time', 'time',  'time'    ,'time']
list_lat  =['LATITUDE' , 'lat' , 'lat' , 'lat' , 'latitude' , 'lat']
list_lon  =['LONGITUDE', 'lon' , 'lon' , 'lon' , 'longitude','lon']

# iterate over radar name  

for radar_num, radar in enumerate(list_radar):

    fig = plt.figure(figsize=(16,16))
    ax1 = fig.add_subplot(111)

    print(radar)

    # Read remapped model files & rename lon / lat with standard lon / lat name
    Ufile_AGRIF = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE/eNEATL36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
    Vfile_AGRIF = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE/eNEATL36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
    Ufile_AGRIF = Ufile_AGRIF.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()
    Vfile_AGRIF = Vfile_AGRIF.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()

    Ufile_TWIN = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_TWIN/eNEATL36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
    Vfile_TWIN = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_TWIN/eNEATL36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
    Ufile_TWIN = Ufile_TWIN.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()
    Vfile_TWIN = Vfile_TWIN.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()

    Ufile_FREE = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_IBI36_FREE/IBI36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
    Vfile_FREE = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_IBI36_FREE/IBI36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
    Ufile_FREE = Ufile_FREE.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()
    Vfile_FREE = Vfile_FREE.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()


    Ufile_ANA = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_IBI36/IBI36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
    Vfile_ANA = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_IBI36/IBI36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
    Ufile_ANA = Ufile_ANA.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()
    Vfile_ANA = Vfile_ANA.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()

    # Read obs & rename lon / lat with standard lon / lat name
    file_OBS = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/201?/*.nc') #, chunks={'time_counter': 50} )#, preprocess = apply_drop_duplicates)
    file_OBS = file_OBS.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}).squeeze()


    # Get variables over the right period & rename time dimension
    U_AGRIF=Ufile_AGRIF.sozocrtx.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'}) 
    V_AGRIF=Vfile_AGRIF.somecrty.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})

    U_TWIN=Ufile_TWIN.sozocrtx.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})
    V_TWIN=Vfile_TWIN.somecrty.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})

    U_FREE=Ufile_FREE.uos.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})
    V_FREE=Vfile_FREE.vos.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})

    U_ANA=Ufile_ANA.uos.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})
    V_ANA=Vfile_ANA.vos.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})

    U_OBS=file_OBS[list_var_u[radar_num]].rename({list_time[radar_num] : 'TIME'}).sel(TIME=slice('2017-06-01', '2018-06-01')) 
    V_OBS=file_OBS[list_var_v[radar_num]].rename({list_time[radar_num] : 'TIME'}).sel(TIME=slice('2017-06-01', '2018-06-01')) 


    # Workaround to make a proper time axis for obs - Round to hour if hourly data, round to minute if higher frequency
    time_counter = pd.to_datetime(U_OBS.TIME.values.astype(str))
    time_counter_round = time_counter.copy().round("H")

    TIME_reindex = pd.date_range(start=time_counter_round[0], end=time_counter_round[-1], freq='H')
    U_OBS['TIME']= time_counter_round
    V_OBS['TIME']= time_counter_round
    try:
        time_counter_round = time_counter.copy().round("H")

        TIME_reindex = pd.date_range(start=time_counter_round[0], end=time_counter_round[-1], freq='H')
        U_OBS['TIME']= time_counter_round
        V_OBS['TIME']= time_counter_round

        U_OBS = U_OBS.reindex(TIME=TIME_reindex, fill_value=np.nan)
        V_OBS = V_OBS.reindex(TIME=TIME_reindex, fill_value=np.nan)
    except:
        time_counter_round = time_counter.copy().round("min")

        TIME_reindex = pd.date_range(start=time_counter_round[0], end=time_counter_round[-1], freq='H')
        U_OBS['TIME']= time_counter_round
        V_OBS['TIME']= time_counter_round

        U_OBS = U_OBS.reindex(TIME=TIME_reindex, fill_value=np.nan)
        V_OBS = V_OBS.reindex(TIME=TIME_reindex, fill_value=np.nan)


    U_OBS_load = U_OBS.values # Just to load data into memory
    V_OBS_load = V_OBS.values

    # Temporal interpolation of the remapped model data + mask data            
    U_AGRIF_intp= U_AGRIF.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS_load),np.nan)
    V_AGRIF_intp= V_AGRIF.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS_load),np.nan)

    U_TWIN_intp= U_TWIN.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS_load),np.nan)
    V_TWIN_intp= V_TWIN.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS_load),np.nan)

    U_FREE_intp= U_FREE.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS_load),np.nan)
    V_FREE_intp= V_FREE.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS_load),np.nan)

    U_ANA_intp= U_ANA.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS_load),np.nan)
    V_ANA_intp= V_ANA.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS_load),np.nan)

    U_AGRIF_intp=U_AGRIF_intp.where(~np.isnan(V_AGRIF_intp),np.nan).where(~np.isnan(V_TWIN_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan)
    V_AGRIF_intp=V_AGRIF_intp.where(~np.isnan(U_AGRIF_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)

    U_TWIN_intp=U_TWIN_intp.where(~np.isnan(V_AGRIF_intp),np.nan).where(~np.isnan(V_TWIN_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan)
    V_TWIN_intp=V_TWIN_intp.where(~np.isnan(U_AGRIF_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)

    U_FREE_intp=U_FREE_intp.where(~np.isnan(V_FREE_intp),np.nan).where(~np.isnan(V_ANA_intp),np.nan).where(~np.isnan(U_ANA_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)
    V_FREE_intp=V_FREE_intp.where(~np.isnan(U_FREE_intp),np.nan).where(~np.isnan(U_ANA_intp),np.nan).where(~np.isnan(V_FREE_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)

    U_ANA_intp=U_ANA_intp.where(~np.isnan(V_FREE_intp),np.nan).where(~np.isnan(V_ANA_intp),np.nan).where(~np.isnan(U_ANA_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)
    V_ANA_intp=V_ANA_intp.where(~np.isnan(U_FREE_intp),np.nan).where(~np.isnan(U_ANA_intp),np.nan).where(~np.isnan(V_FREE_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)


    # mask data where there is less than 30% of the data valid over the time period
    count_U = np.count_nonzero(U_OBS.where(~np.isnan(U_OBS),0).values,axis=0)
    count_V = np.count_nonzero(V_OBS.where(~np.isnan(V_OBS),0).values,axis=0)
    U_AGRIF_intp = U_AGRIF_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_AGRIF_intp = V_AGRIF_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    U_TWIN_intp = U_TWIN_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_TWIN_intp = V_TWIN_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    U_FREE_intp = U_FREE_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_FREE_intp = V_FREE_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    U_ANA_intp = U_ANA_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_ANA_intp = V_ANA_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    U_OBS = U_OBS.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_OBS = V_OBS.where(count_V > int(len(V_AGRIF_intp[:,0,0])*0.3), np.nan)

    # Flatten data
    U_AGRIF_intp_flat = U_AGRIF_intp.values.ravel()
    V_AGRIF_intp_flat = V_AGRIF_intp.values.ravel()
    U_TWIN_intp_flat = U_TWIN_intp.values.ravel()
    V_TWIN_intp_flat = V_TWIN_intp.values.ravel()
    U_FREE_intp_flat = U_FREE_intp.values.ravel()
    V_FREE_intp_flat = V_FREE_intp.values.ravel()
    U_ANA_intp_flat = U_ANA_intp.values.ravel()
    V_ANA_intp_flat = V_ANA_intp.values.ravel()
    U_OBS_flat = U_OBS.where(~np.isnan(U_AGRIF_intp.values),np.nan).values.ravel()
    V_OBS_flat = V_OBS.where(~np.isnan(V_AGRIF_intp.values),np.nan).values.ravel()

    print(U_OBS_flat[~np.isnan(U_OBS_flat)].shape)
    print(U_AGRIF_intp_flat[~np.isnan(U_AGRIF_intp_flat)].shape)
    print(U_TWIN_intp_flat[~np.isnan(U_TWIN_intp_flat)].shape)

    # Compute statistics for taylor diagrams
    taylor_stats1 = sm.taylor_statistics(U_AGRIF_intp_flat[~np.isnan(U_AGRIF_intp_flat)],U_OBS_flat[~np.isnan(U_OBS_flat)])
    taylor_stats2 = sm.taylor_statistics(U_TWIN_intp_flat[~np.isnan(U_TWIN_intp_flat)],U_OBS_flat[~np.isnan(U_OBS_flat)])
#    taylor_stats3 = sm.taylor_statistics(U_FREE_intp_flat[~np.isnan(U_FREE_intp_flat)],U_OBS_flat[~np.isnan(U_OBS_flat)])
    taylor_stats3 = sm.taylor_statistics(U_ANA_intp_flat[~np.isnan(U_ANA_intp_flat)],U_OBS_flat[~np.isnan(U_OBS_flat)])
 
    # Store statistics in arrays
    sdevU = np.around(np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1],
                     taylor_stats2['sdev'][1], taylor_stats3['sdev'][1]]),4)

    crmsdU = np.around(np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1],
                      taylor_stats2['crmsd'][1],taylor_stats3['crmsd'][1]]),4)

    ccoefU = np.around(np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1],
                      taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1]]),4)

    print('sdevU',sdevU)
    print('crmsdU',crmsdU)
    print('ccoefU',ccoefU)

    # Compute statistics for taylor diagrams
    taylor_stats1 = sm.taylor_statistics(V_AGRIF_intp_flat[~np.isnan(V_AGRIF_intp_flat)],V_OBS_flat[~np.isnan(V_OBS_flat)])
    taylor_stats2 = sm.taylor_statistics(V_TWIN_intp_flat[~np.isnan(V_TWIN_intp_flat)],V_OBS_flat[~np.isnan(V_OBS_flat)])
    #taylor_stats3 = sm.taylor_statistics(V_FREE_intp_flat[~np.isnan(V_FREE_intp_flat)],V_OBS_flat[~np.isnan(V_OBS_flat)])
    taylor_stats3 = sm.taylor_statistics(V_ANA_intp_flat[~np.isnan(V_ANA_intp_flat)],V_OBS_flat[~np.isnan(V_OBS_flat)])

    # Store statistics in arrays
    sdevV = np.around(np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1],
                     taylor_stats2['sdev'][1], taylor_stats3['sdev'][1]]),4)

    crmsdV = np.around(np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1],
                      taylor_stats2['crmsd'][1],taylor_stats3['crmsd'][1]]),4)

    ccoefV = np.around(np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1],
                      taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1]]),4)
    print('sdevV',sdevV)
    print('crmsdV',crmsdV)
    print('ccoefV',ccoefV)
    stdmax=max([max(sdevU),max(sdevV)])
    stdmax=round(stdmax + 0.2*stdmax,2)
    rmsmaxU=max(crmsdU)
    rmsmaxV=max(crmsdV)

    print('MAX',stdmax)
    print('MAX rms U',rmsmaxU)
    print('MAX rms V',rmsmaxV)

    # NOW PLOT
    label = ['Non-Dimensional Observation', 'NEST', 'TWIN', 'ANA']
    sm.taylor_diagram(sdevU,crmsdU,ccoefU,
                      markerLabel = label,markerSize=20,markerLabelColor='r', markerColor='r',styleOBS = '-',
                      colOBS = 'r',colCOR = 'k',colRMS = 'r', markerobs = 'o',
                      titleOBS = 'observation U ', axisMax=stdmax, tickRMS=np.around(np.linspace(0,rmsmaxU,4),2))

    sm.taylor_diagram(sdevV,crmsdV,ccoefV,
                      markerLabel = label,markerSize=20,markerLabelColor='b',markerColor='b',styleOBS = '-',
                      colOBS = 'b',colCOR = 'k', colRMS = 'b', markerobs = 'o',
                      titleOBS = 'observation V', axisMax=stdmax, tickRMS=np.around(np.linspace(0,rmsmaxV,4),2))



    # Write plot to file

    plt.savefig('taylor_diagram_'+str(radar)+'.png')

