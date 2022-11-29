#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 13:04:08 2019

@author: tbrivoal
"""
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib import pyplot
from matplotlib import colors as mcolors
import matplotlib.gridspec as gridspec
from scipy.stats import norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["darkblue","blue","white","red","firebrick"])
darkred = plt.cm.Reds(np.linspace(0.99, 1, 2))
darkblue = plt.cm.YlGnBu(np.linspace(0.99, 1, 2))
bwr_cmap1 = cmap(np.linspace(0, 0.5, 252))
bwr_cmap2 = cmap(np.linspace(0.5, 1, 252))
w_cmap = cmap(np.linspace(0.5, 0.5, 78))
# bwr_cmap2 = plt.cm.(np.linspace(0.5, 1, 126))
# combine them and build a new colormap
colors = np.vstack((darkblue,bwr_cmap1,w_cmap,bwr_cmap2,darkred))
mymap =mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

list_radar=['BISC2','GIBR1', 'GALI2',  'EBRO', 'IROI', 'TIRLIG', 'TRAD3']
list_var_u=['EWCT', 'u', 'u', 'u','u','EWCT', 'u']
list_var_v=['NSCT', 'v', 'v',  'v' , 'v','NSCT', 'v']
list_time =['TIME', 'time', 'time', 'time',  'time','TIME','time']
list_lat  =['LATITUDE', 'lat', 'lat', 'lat', 'latitude','LATITUDE', 'lat']
list_lon  =['LONGITUDE', 'lon', 'lon', 'lon', 'longitude','LONGITUDE','lon']


seasons = ["JJA", "SON", "DJF", "MAM"]
season_month_number =[[6, 7, 8], [9, 10, 11], [12, 1, 2], [3, 4, 5]]


#axes=ax1.gridlines( draw_labels=False, linewidth=0)
#axes.ylabels_right = False
#axes.xlabels_top = False
cnt=1
def apply_drop_duplicates(obj):
    data = obj
    for dim in obj.dims: 
        data = data.drop_duplicates(dim=dim)
    return data 

for radar_num, radar in enumerate(list_radar):
    try:
       for iseason in range(len(seasons)): 
            fig = plt.figure(figsize=(8,8))
            proj=ccrs.PlateCarree()
            lon_formatter = LongitudeFormatter(degree_symbol='° ')
            lat_formatter = LatitudeFormatter(degree_symbol='° ')
            
            ax1 = fig.add_subplot(121, projection=proj)
            ax2 = fig.add_subplot(122, projection=proj)
        
            #ax1 = plt.subplot(231, projection=proj)
            ax1.coastlines(resolution='50m')
            ax1.xaxis.set_major_formatter(lon_formatter)
            ax1.yaxis.set_major_formatter(lat_formatter)
            ax1.set_title('AGRIF')
        
            ax2.coastlines(resolution='50m')
            ax2.xaxis.set_major_formatter(lon_formatter)
            ax2.yaxis.set_major_formatter(lat_formatter)
            ax2.set_title('TWIN')
           
 
            print(radar)    
        
            lon_lat = xr.open_dataset('/scratch/work/brivoalt/WORK_ANAMOD/TWIN/mesh_mask_eNEATL36_AGRIF108.nc', drop_variables={"x", "y",})
            lont = lon_lat['glamt'].squeeze()
            latt = lon_lat['gphit'].squeeze()
            lonv = lon_lat['glamv'].squeeze()
            latv = lon_lat['gphiv'].squeeze()    
        
            angle = np.arctan2(((lont-lonv)*np.cos(np.radians(latv))),(latt - latv)) 
            print(np.min(angle))
            Ufile_AGRIF = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE/eNEATL36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
            Vfile_AGRIF = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE/eNEATL36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
            Ufile_TWIN = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_TWIN/eNEATL36_1h_gridU_*') #, preprocess = apply_drop_duplicates)
            Vfile_TWIN = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/EQUIVALENT_MODELE_TWIN/eNEATL36_1h_gridV_*') #, preprocess = apply_drop_duplicates)
        
            Ufile_AGRIF = Ufile_AGRIF.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'})
            Vfile_AGRIF = Vfile_AGRIF.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'})
            Ufile_TWIN = Ufile_TWIN.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'})
            Vfile_TWIN = Vfile_TWIN.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'})
        
            file_OBS = xr.open_mfdataset('../DATA_RADAR/'+str(radar)+'/201?/*.nc') #, chunks={'time_counter': 50} )#, preprocess = apply_drop_duplicates)
            file_OBS = file_OBS.rename({list_lon[radar_num] : 'LON'}).rename({list_lat[radar_num] : 'LAT'}).squeeze()
            U_AGRIF=Ufile_AGRIF.sozocrtx.squeeze().rename({'time_counter' : 'TIME'}) 
            V_AGRIF=Vfile_AGRIF.somecrty.squeeze().rename({'time_counter' : 'TIME'}) 
            U_TWIN=Ufile_TWIN.sozocrtx.squeeze().rename({'time_counter' : 'TIME'}) 
            V_TWIN=Vfile_TWIN.somecrty.squeeze().rename({'time_counter' : 'TIME'}) 
        
            U_OBS=file_OBS[list_var_u[radar_num]].rename({list_time[radar_num] : 'TIME'}) 
            V_OBS=file_OBS[list_var_v[radar_num]].rename({list_time[radar_num] : 'TIME'}) 
            U_AGRIF=U_AGRIF.sel(TIME=U_AGRIF.TIME.dt.month.isin(season_month_number[iseason]))
            V_AGRIF=V_AGRIF.sel(TIME=V_AGRIF.TIME.dt.month.isin(season_month_number[iseason]))
            U_TWIN=U_TWIN.sel(TIME=U_TWIN.TIME.dt.month.isin(season_month_number[iseason]))
            V_TWIN=V_TWIN.sel(TIME=V_TWIN.TIME.dt.month.isin(season_month_number[iseason]))
            U_OBS=U_OBS.sel(TIME=U_OBS.TIME.dt.month.isin(season_month_number[iseason]))
            V_OBS=V_OBS.sel(TIME=V_OBS.TIME.dt.month.isin(season_month_number[iseason]))
  
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
        
            lat_tmp=file_OBS['LAT'].squeeze().values
            lon_tmp=file_OBS['LON'].squeeze().values
            lon, lat = np.meshgrid(lon_tmp,lat_tmp)
        
            latmod_tmp=Ufile_AGRIF['LAT'].squeeze().values
            lonmod_tmp=Ufile_AGRIF['LON'].squeeze().values
            lonmod, latmod = np.meshgrid(lonmod_tmp,latmod_tmp)
        
            
            U_AGRIF_intp = U_AGRIF.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS),np.nan) #.where(LAT == file_OBS[list_lat[radar_num]])
            U_TWIN_intp = U_TWIN.interp(TIME=U_OBS['TIME']).where(~np.isnan(U_OBS),np.nan) #.where(LAT == file_OBS[list_lat[radar_num]])
            V_AGRIF_intp = V_AGRIF.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS),np.nan) #.where(LAT == file_OBS[list_lat[radar_num]])
            V_TWIN_intp = V_TWIN.interp(TIME=V_OBS['TIME']).where(~np.isnan(V_OBS),np.nan) #.where(LAT == file_OBS[list_lat[radar_num]])

            count_U = np.count_nonzero(U_OBS.where(~np.isnan(U_OBS),0).values,axis=0)
            count_V = np.count_nonzero(V_OBS.where(~np.isnan(V_OBS),0).values,axis=0)
            U_AGRIF_intp = U_AGRIF_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
            V_AGRIF_intp = V_AGRIF_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
            U_TWIN_intp = U_TWIN_intp.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
            V_TWIN_intp = V_TWIN_intp.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
            U_OBS = U_OBS.where(count_U > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
            V_OBS = V_OBS.where(count_V > int(len(U_AGRIF_intp[:,0,0])*0.3), np.nan)
        
            matplotlib.rcParams.update({'font.size': 10})
            print(lon.shape, U_AGRIF)
            var1 = np.nanmean(np.sqrt((U_AGRIF_intp - U_OBS)**2),axis=0)
            var2 = np.nanmean(np.sqrt((U_TWIN_intp - U_OBS)**2),axis=0)
            print(var1)
            #maxval = np.nanmax([int(np.nanmax(var1)), int(np.nanmax(var2)), int(np.nanmax(var3)), abs(int(np.nanmin(var1))), abs(int(np.nanmin(var2))) , abs(int(np.nanmin(var3)))])
            maxval =np.nanmax([var1,var2])
            print(maxval)
            im1 = ax1.contourf(lon, lat, var1 ,levels=np.linspace(0,maxval,51), cmap="jet",transform=ccrs.PlateCarree(),extend='both')

        
            cbar1=fig.colorbar(im1, ax=ax1,ticks=np.linspace(0,maxval,11),orientation="vertical",fraction=0.046, pad=0.04) #ticks=[-40,20,0,20,40]) #ticks=[0,1,2,3]) #
            cbar1.set_label('U (m/s)', rotation=270, labelpad=12, fontsize=12,fontweight="bold")
        
            im2 = ax2.contourf(lon, lat, var2 ,levels=np.linspace(0,maxval,51), cmap="jet",transform=ccrs.PlateCarree(),extend='both')
        
            cbar2=fig.colorbar(im2, ax=ax2,ticks=np.linspace(0,maxval,11),orientation="vertical",fraction=0.046, pad=0.04) #ticks=[-40,20,0,20,40]) #ticks=[0,1,2,3]) #
            cbar2.set_label('U (m/s)', rotation=270, labelpad=12, fontsize=12,fontweight="bold")
        #
        
            cnt = cnt + 1
        
            plt.tight_layout()
            plt.savefig('MAP_U_RMS_AGRIF'+radar+'_'+seasons[iseason]+'.png')
            plt.close()
            print('PA DE PROBELEME')
    except: 
        print('PROBELEME')
        #    
        #     
