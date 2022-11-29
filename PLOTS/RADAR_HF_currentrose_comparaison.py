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
#matplotlib.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

from matplotlib import pyplot
from matplotlib import colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import windrose
from windrose import WindroseAxes



# Define radar name, name of the variables etc for each variable (1 radar = 1 column)

list_radar=['BISC2','GIBR1', 'GALI2',  'EBRO', 'IROI', 'TIRLIG', 'TRAD3']
list_var_u=['EWCT', 'u', 'u', 'u','u','EWCT', 'u']
list_var_v=['NSCT', 'v', 'v',  'v' , 'v','NSCT', 'v']
list_time =['TIME', 'time', 'time', 'time',  'time','TIME','time']
list_lat  =['LATITUDE', 'lat', 'lat', 'lat', 'latitude','LATITUDE', 'lat']
list_lon  =['LONGITUDE', 'lon', 'lon', 'lon', 'longitude','LONGITUDE','lon']

# iterate over radar name  

for radar_num, radar in enumerate(list_radar):

    fig = plt.figure(figsize=(16,16))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    #ax1 = plt.subplot(231, projection=proj)
    ax1.set_title('Zonal current (m/s)')
    ax2.set_title('Meridional current (m/s)')

    # Read remapped model files & rename time / lon / lat with standard time / lon / lat name
    print(radar) 
    Ufile_AGRIF = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE/eNEATL36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
    Vfile_AGRIF = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE/eNEATL36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
    Ufile_AGRIF = Ufile_AGRIF.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()
    Vfile_AGRIF = Vfile_AGRIF.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()

    Ufile_TWIN = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_TWIN/eNEATL36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
    Vfile_TWIN = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_TWIN/eNEATL36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
    Ufile_TWIN = Ufile_TWIN.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()
    Vfile_TWIN = Vfile_TWIN.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}) #.squeeze()


    file_OBS = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/201?/*.nc') #, chunks={'time_counter': 50} )#, preprocess = apply_drop_duplicates)
    file_OBS = file_OBS.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}).squeeze()
    U_AGRIF=Ufile_AGRIF.sozocrtx.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'}) 
    V_AGRIF=Vfile_AGRIF.somecrty.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})

    U_TWIN=Ufile_TWIN.sozocrtx.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})
    V_TWIN=Vfile_TWIN.somecrty.squeeze().sel(time_counter=slice('2017-06-01', '2018-06-01')).rename({'time_counter' : 'TIME'})

    # Read obs & rename time / lon / lat with standard time / lon / lat name
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

    # Temporal interpolation of the remapped model data + mask data            
    U_AGRIF_intp= U_AGRIF.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS),np.nan)
    V_AGRIF_intp= V_AGRIF.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS),np.nan)

    U_TWIN_intp= U_TWIN.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS),np.nan)
    V_TWIN_intp= V_TWIN.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS),np.nan)

    U_AGRIF_intp=U_AGRIF_intp.where(~np.isnan(V_AGRIF_intp),np.nan).where(~np.isnan(V_TWIN_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan)
    V_AGRIF_intp=V_AGRIF_intp.where(~np.isnan(U_AGRIF_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)

    U_TWIN_intp=U_TWIN_intp.where(~np.isnan(V_AGRIF_intp),np.nan).where(~np.isnan(V_TWIN_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan)
    V_TWIN_intp=V_TWIN_intp.where(~np.isnan(U_AGRIF_intp),np.nan).where(~np.isnan(U_TWIN_intp),np.nan).where(~np.isnan(V_AGRIF_intp),np.nan)

    # mask data where there is less than 30% of the data valid over the time period
    count_U = np.count_nonzero(U_OBS.where(~np.isnan(U_OBS),0).values,axis=0)
    count_V = np.count_nonzero(V_OBS.where(~np.isnan(V_OBS),0).values,axis=0)
    U_AGRIF_intp = U_AGRIF_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_AGRIF_intp = V_AGRIF_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    U_TWIN_intp = U_TWIN_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_TWIN_intp = V_TWIN_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    U_OBS = U_OBS.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
    V_OBS = V_OBS.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)

    # flatten data 
    U_AGRIF_intp_flat = U_AGRIF_intp.values.ravel()
    V_AGRIF_intp_flat = V_AGRIF_intp.values.ravel()
    U_TWIN_intp_flat = U_TWIN_intp.values.ravel()
    V_TWIN_intp_flat = V_TWIN_intp.values.ravel()
    U_OBS_flat = U_OBS.where(~np.isnan(U_AGRIF_intp),np.nan).values.ravel()
    V_OBS_flat = V_OBS.where(~np.isnan(V_AGRIF_intp),np.nan).values.ravel()

###########################################################################
### PLOTS #################################
###################

    #fig=plt.figure()
    #ax = WindroseAxes.from_ax()
    #print(U_OBS_flat[~np.isnan(U_OBS_flat)].shape)
    #print(V_OBS_flat[~np.isnan(V_OBS_flat)].shape)

    #print(U_AGRIF_intp_flat[~np.isnan(U_AGRIF_intp_flat)].shape)
    #print(V_AGRIF_intp_flat[~np.isnan(V_AGRIF_intp_flat)].shape)

    #print(U_TWIN_intp_flat[~np.isnan(U_TWIN_intp_flat)].shape)
    #print(V_TWIN_intp_flat[~np.isnan(V_TWIN_intp_flat)].shape)
    U = U_OBS_flat[~np.isnan(U_OBS_flat)]
    V = V_OBS_flat[~np.isnan(V_OBS_flat)]
    ws = np.sqrt(U**2 + V**2)
    wd = np.degrees(np.arctan2(U,V)) % 360.0
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(131, projection="windrose")
    ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0, 1.5, 0.3) )
    #ax.set_ylim(0,20)
    ax.set_legend()
    ax.set_title('OBS')
    #plt.tight_layout()

    U = U_AGRIF_intp_flat[~np.isnan(U_AGRIF_intp_flat)]
    V = V_AGRIF_intp_flat[~np.isnan(V_AGRIF_intp_flat)]
    ws = np.sqrt(U**2 + V**2)
    wd = np.degrees(np.arctan2(U,V)) % 360.0
    ax = fig.add_subplot(132, projection="windrose")
    ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0, 1.5, 0.3) )
    #ax.set_ylim(0,20)

    ax.set_legend()
    ax.set_title("AGRIF")
    #plt.tight_layout()

    U = U_TWIN_intp_flat[~np.isnan(U_TWIN_intp_flat)]
    V = V_TWIN_intp_flat[~np.isnan(V_TWIN_intp_flat)]
    ws = np.sqrt(U**2 + V**2)
    wd = np.degrees(np.arctan2(U,V)) % 360.0
    ax = fig.add_subplot(133, projection="windrose")
    ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0, 1.5, 0.3) )
    #ax.set_ylim(0,20)

    ax.set_legend()
    ax.set_title("TWIN")
    #plt.tight_layout()
    plt.savefig('HFRADAR_Windrose_ALL_intp_'+radar+'.png')
    plt.close()


